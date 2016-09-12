module piv_mod
	use kinds_mod
	use settings_mod
	use text_mod
	use array_mod
	use autodiff_mod
	use displacement_mod
	use pair_mod
	use omp_lib
	use cluster_mod
	use netCDF_mod
	implicit none
	
contains
	
	!==========!
	!= Passes =!
	!==========!
	
	subroutine doPass(p,k,N,method,reference)
		class(pair_t),intent(inout)::p
			!! Image pair
		integer,intent(in)::k
			!! Pass index
		integer,dimension(2),intent(in)::N
			!! Window size
		character(*),intent(in)::method
			!! Computation method
		integer,intent(in),optional::reference
		
		integer::tid
		integer::lref
		type(regions_t)::R
		type(ad3_t),dimension(2)::d
		character(:),allocatable::fn
		real(wp)::F
		integer::i,j
		
		lref = -1
		if(present(reference)) lref = reference
		
		!$omp parallel private(i,j,d,tid)
		tid = omp_get_thread_num()
		!$omp barrier
		!$omp do schedule(static,1)
		do j=1,p%Nv(2)
			if(tid==0 .and. j/=p%Nv(2) .and. amRoot()) then
				call showProgress('Correlating '//intToChar(product(p%Nv))//' vectors',real(j-1,wp)/real(p%Nv(2)-1,wp))
			end if
			do i=1,p%Nv(1)
				select case(lref)
				case(0:)
					R = secondPass(i,j,lref)
					F = 0.1_wp
				case default
					R = firstPass(i,j)
					F = correlationFactor
				end select
				
				select case(method)
				case('map')
					d = R%crossCorrelateDirect(F,p%idx,k)
				case('lsq')
					d = R%leastSquares(lsqOrder,p%idx,k)
				end select
				
				p%passes(k)%u(i,j) = d(1)
				p%passes(k)%v(i,j) = d(2)
				
				if(per_pixel) then
					fn = './results/'//prefix//'/vector'
					fn = fn//'-'//intToChar(p%idx)
					fn = fn//'-['//intToChar(i)//','//intToChar(j)//'|'
					fn = fn//''//intToChar(k)//']'
					fn = fn//'.nc'
					call writeVector(fn,d,R)
				end if
			end do
		end do
		!$omp end do
		!$omp barrier
		if(tid==0 .and. amRoot()) call showProgress('Correlating '//intToChar(product(p%Nv))//' vectors',1.0_wp)
		!$omp barrier
		!$omp end parallel
		
	contains
	
		function firstPass(i,j) result(o)
			integer,intent(in)::i,j
				!! Corrdinates of vector
			type(regions_t)::o
				!! Result
			
			integer::il,ih
			integer::jl,jh
			
			il = nint( p%vx(i)-real(N(1),wp)/2.0_wp )
			jl = nint( p%vy(j)-real(N(2),wp)/2.0_wp )
			ih = il+N(1)-1
			jh = jl+N(2)-1
			
			o = regions_t( p%A(il:ih,jl:jh),p%B(il:ih,jl:jh) , [i,j] )
		end function firstPass
	
		function secondPass(i,j,ref) result(o)
			!! @todo
			!! Re-evaluation of algorithm is needed
			integer,intent(in)::i,j
				!! Corrdinates of vector
			integer,intent(in)::ref
				!! Reference pass for window shifting
			type(regions_t)::o
				!! Result
			
			integer,dimension(2)::sp,sm,s
			real(wp),dimension(2)::up
			integer::il,ih
			integer::jl,jh
			
			integer::ilA,ihA,ilB,ihB
			integer::jlA,jhA,jlB,jhB
			logical::ok
			
			up = real([p%passes(ref)%u(i,j),p%passes(ref)%v(i,j)])
			s  = nint(up)
			sm = s/2
			sp = s-sm
			
			il = nint( p%vx(i)-real(N(1),wp)/2.0_wp )
			jl = nint( p%vy(j)-real(N(2),wp)/2.0_wp )
			ih = il+N(1)-1
			jh = jl+N(2)-1
			
			ilA = il-sm(1)
			ihA = ih-sm(1)
			ilB = il+sp(1)
			ihB = ih+sp(1)
			
			jlA = jl-sm(2)
			jhA = jh-sm(2)
			jlB = jl+sp(2)
			jhB = jh+sp(2)
			
			ok = checkBoundsError([ilA,ihA],[jlA,jhA],[ilB,ihB],[jlB,jhB],p%N)
			
			o = regions_t( p%A(ilA:ihA,jlA:jhA) , p%B(ilB:ihB,jlB:jhB) , [i,j] , s )
		end function secondPass
	
		function checkBoundsError(iA,jA,iB,jB,N) result(ok)
			integer,dimension(2),intent(in)::iA,jA
			integer,dimension(2),intent(in)::iB,jB
			integer,dimension(2),intent(in)::N
			logical::ok
			
			ok = .true.
			
			if( any([iA,iB]<1) ) ok = .false.
			if( any([jA,jB]<1) ) ok = .false.
			if( any([iA,iB]>N(1)) ) ok = .false.
			if( any([jA,jB]>N(2)) ) ok = .false.
			
			if(.not.ok) then
				write(*,*) ''
				write(*,*) N(1),N(2)
				write(*,*) 'A('// & 
					& intToChar( iA(1) )//':'//intToChar( iA(2) )//','// &
					& intToChar( jA(1) )//':'//intToChar( jA(2) )//')'
				write(*,*) 'B('// &
					& intToChar( iB(1) )//':'//intToChar( iB(2) )//','// &
					& intToChar( jB(1) )//':'//intToChar( jB(2) )//')'
				write(*,*) ''
				stop 1
			end if
		end function checkBoundsError
	
	end subroutine doPass

	!=============!
	!= Utilities =!
	!=============!

	subroutine filter(p,k,tol)
		class(pair_t),intent(inout)::p
			!! Image pair
		integer,intent(in)::k
			!! Pass index
		real(wp),intent(in)::tol
			!! Tolerance before culling [0,1]
		
		type(ad3_t),dimension(2)::u,m,t
		integer,dimension(2)::N
		integer::i,j
		
		N = [size(p%vx),size(p%vy)]
		
		do j=1,N(2)
			do i=1,N(1)
				u = [p%passes(k)%u(i,j),p%passes(k)%v(i,j)]
				t = [p%passes(0)%u(i,j),p%passes(0)%v(i,j)]
				m = tryMean(i,j)
				if( norm2(real(u-t))/norm2(real(t)) > tol) then
					p%passes(k)%u(i,j) = m(1)
					p%passes(k)%v(i,j) = m(2)
				end if
			end do
		end do
	
	contains
	
		function tryMean(i,j) result(o)
			integer,intent(in)::i,j
			type(ad3_t),dimension(2)::o
			
			integer::c,ii,jj
			
			o = 0.0_wp
			c = 0
			do ii=-1,1,2
				do jj=-1,1,2
					if(i+ii<1   ) cycle
					if(i+ii>N(1)) cycle
					if(j+jj<1   ) cycle
					if(j+jj>N(2)) cycle
					c = c+1
					o = o+[p%passes(k)%u(i+ii,j+jj),p%passes(k)%v(i+ii,j+jj)]
				end do
			end do
			
			o = o/real(c,wp)
		end function tryMean
	
	end subroutine filter

	subroutine writeVector(fn,v,R)
		character(*),intent(in)::fn
		type(ad3_t),dimension(2)::v
		type(regions_t),intent(in)::R
		
		real(wp),dimension(:),allocatable::x,y
		real(wp),dimension(:,:),allocatable::var
		integer,dimension(4)::adN
		integer,dimension(2)::N
		
		adN = get_adN()
		N = adN(2:3)
		
		x = linspace(1.0_wp,real(N(1),wp),N(1))
		y = linspace(1.0_wp,real(N(2),wp),N(2))
		
		call writeGrid(fn, &
			& ['I    ','dudI ','dvdI ', &
			& 'dIdU ','dIdUx','dIdUy', &
			& 'dIdV ','dIdVx','dIdVy', &
			& 'dIdR ','dIdN '], &
			& x,y)
		
		! First Image
		
		call writeStep(fn,0.0_wp,1,'I',real(R%A))
		
		var = v(1)%I(1:N(1),1:N(2),1)
		call writeStep(fn,0.0_wp,1,'dudI',var)
		
		var = v(2)%I(1:N(1),1:N(2),1)
		call writeStep(fn,0.0_wp,1,'dvdI',var)
		
		call writeStep(fn,0.0_wp,1,'dIdU ',der(R%A,ADS_U ))
		call writeStep(fn,0.0_wp,1,'dIdUx',der(R%A,ADS_Ux))
		call writeStep(fn,0.0_wp,1,'dIdUy',der(R%A,ADS_Uy))
		call writeStep(fn,0.0_wp,1,'dIdV ',der(R%A,ADS_V ))
		call writeStep(fn,0.0_wp,1,'dIdVx',der(R%A,ADS_Vx))
		call writeStep(fn,0.0_wp,1,'dIdVy',der(R%A,ADS_Vy))
		call writeStep(fn,0.0_wp,1,'dIdR ',der(R%A,ADS_R ))
		call writeStep(fn,0.0_wp,1,'dIdN ',der(R%A,ADS_N ))
		
		! Second Image
		
		call writeStep(fn,1.0_wp,2,'I',real(R%B))
		
		var = v(1)%I(1:N(1),1:N(2),2)
		call writeStep(fn,1.0_wp,2,'dudI',var)
		
		var = v(2)%I(1:N(1),1:N(2),2)
		call writeStep(fn,1.0_wp,2,'dvdI',var)
		
		call writeStep(fn,1.0_wp,2,'dIdU ',der(R%B,ADS_U ))
		call writeStep(fn,1.0_wp,2,'dIdUx',der(R%B,ADS_Ux))
		call writeStep(fn,1.0_wp,2,'dIdUy',der(R%B,ADS_Uy))
		call writeStep(fn,1.0_wp,2,'dIdV ',der(R%B,ADS_V ))
		call writeStep(fn,1.0_wp,2,'dIdVx',der(R%B,ADS_Vx))
		call writeStep(fn,1.0_wp,2,'dIdVy',der(R%B,ADS_Vy))
		call writeStep(fn,1.0_wp,2,'dIdR ',der(R%B,ADS_R ))
		call writeStep(fn,1.0_wp,2,'dIdN ',der(R%B,ADS_N ))
	end subroutine writeVector

end module piv_mod
