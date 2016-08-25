module piv_mod
	use kinds_mod
	use settings_mod
	use utilities_mod
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
		type(regions_t)::R
		type(ad3_t),dimension(2)::d
		character(:),allocatable::fn
		integer::i,j
		
		!$omp parallel private(i,j,d,tid)
		tid = omp_get_thread_num()
		!$omp barrier
		!$omp do schedule(static,1)
		do j=1,p%Nv(2)
			if(tid==0 .and. j/=p%Nv(2) .and. amRoot()) then
				call showProgress('Correlating '//int2char(product(p%Nv))//' vectors',real(j-1,wp)/real(p%Nv(2)-1,wp))
			end if
			do i=1,p%Nv(1)
				if(present(reference)) then
					if(reference>=0) then
						R = secondPass(i,j,reference)
					else
						R = firstPass(i,j)
					end if
				else
					R = firstPass(i,j)
				end if
				
				select case(method)
				case('map')
					d = R%crossCorrelateDirect(correlationFactor,p%idx)
				case('lsq')
					d = R%leastSquares(lsqOrder)
				end select
				
				p%passes(k)%u(i,j) = d(1)
				p%passes(k)%v(i,j) = d(2)
				
				fn = './results/'//prefix//'/vector'
				fn = fn//'-'//int2char(p%idx)
				fn = fn//'-['//int2char(i)//','//int2char(j)//']'
				fn = fn//'.nc'
				call writeVector(fn,d)
			end do
		end do
		!$omp end do
		!$omp barrier
		if(tid==0 .and. amRoot()) call showProgress('Correlating '//int2char(product(p%Nv))//' vectors',1.0_wp)
		!$omp barrier
		!$omp end parallel
		
	contains
	
		function firstPass(i,j) result(o)
			integer,intent(in)::i,j
				!! Corrdinates of vector
			type(regions_t)::o
				!! Result
			
			real(wp),dimension(2)::ps
			integer::il,ih
			integer::jl,jh
			
			ps = p%L/real(p%N,wp)
			
			il = minloc( abs(p%px-(p%vx(i)-real(N(1),wp)*ps(1)/2.0_wp)) , 1 )
			jl = minloc( abs(p%py-(p%vy(j)-real(N(2),wp)*ps(2)/2.0_wp)) , 1 )
			ih = il+N(1)-1
			jh = jl+N(2)-1
			
			o = regions_t( p%A(il:ih,jl:jh),p%B(il:ih,jl:jh) , [i,j] )
		end function firstPass
	
		function secondPass(i,j,ref) result(o)
			integer,intent(in)::i,j
				!! Corrdinates of vector
			integer,intent(in)::ref
				!! Reference pass for window shifting
			type(regions_t)::o
				!! Result
			
			real(wp),dimension(2)::ps
			integer,dimension(2)::sp,sm,s
			real(wp),dimension(2)::up
			integer::il,ih
			integer::jl,jh
			
			integer::ilA,ihA,ilB,ihB
			integer::jlA,jhA,jlB,jhB
			
			ps = p%L/real(p%N,wp)
			
			up = real([p%passes(ref)%u(i,j),p%passes(ref)%v(i,j)])
			
			s  = nint(up)
			sp = s/2
			sm = s-sp
			
			il = minloc( abs(p%px-(p%vx(i)-real(N(1),wp)*ps(1)/2.0_wp)) , 1 )
			jl = minloc( abs(p%py-(p%vy(j)-real(N(2),wp)*ps(2)/2.0_wp)) , 1 )
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
			
			if( any([ilA,ihA,ilB,ihB]<1) ) return
			if( any([jlA,jhA,jlB,jhB]<1) ) return
			if( any([ilA,ihA,ilB,ihB]>p%N(1)) ) return
			if( any([jlA,jhA,jlB,jhB]>p%N(2)) ) return
			
			o = regions_t( p%A(ilA:ihA,jlA:jhA) , p%B(ilB:ihB,jlB:jhB) , [i,j] , s )
		end function secondPass
	
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

	subroutine writeVector(fn,v)
		character(*),intent(in)::fn
		type(ad3_t),dimension(2)::v
		
		real(wp),dimension(:),allocatable::x,y
		real(wp),dimension(:,:),allocatable::var
		integer,dimension(4)::adN
		integer,dimension(2)::N
		
		adN = get_adN()
		N = adN(2:3)
		
		x = linspace(1.0_wp,real(N(1),wp),N(1))
		y = linspace(1.0_wp,real(N(2),wp),N(2))
		
		call write_grid(fn,['U','V'],x,y)
		
		! First Image
		
		var = v(1)%I(1:N(1),1:N(2),1)
		call write_step(fn,0.0_wp,1,'U',var)
		
		var = v(2)%I(1:N(1),1:N(2),1)
		call write_step(fn,0.0_wp,1,'V',var)
		
		! Second Image
		
		var = v(1)%I(1:N(1),1:N(2),2)
		call write_step(fn,1.0_wp,2,'U',var)
		
		var = v(2)%I(1:N(1),1:N(2),2)
		call write_step(fn,1.0_wp,2,'V',var)
		
	end subroutine writeVector

end module piv_mod
