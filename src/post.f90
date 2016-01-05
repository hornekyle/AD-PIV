program post_prg
	use kinds_mod
	use autodiff_mod
	use utilities_mod
	use plplotlib_mod
	use pair_mod
	use piv_mod
	implicit none
	
	character(128)::cfn
	type(pair_t),dimension(:),allocatable::p
	type(pair_t)::mean_pair
	type(pair_t)::stdev_pair
	integer::k
	
	call get_command_argument(1,cfn)
	call readConfig(cfn)
	
	call setup(device='svgqt',filename='./results/'//prefix//'/output-post-%n.svg',colormap='BlueYellow')
	
	allocate(p(N_pairs))
	do k=1,N_pairs
		call p(k)%readPair(vfn='./results/'//prefix//'/vectors-'//int2char(k)//'.nc')
		call doMask(p(k),0.2_wp)
	end do
	mean_pair  = meanPair(p)
	stdev_pair = stdevPair(p,mean_pair)
	
	call mean_pair%writeVectors('./results/'//prefix//'/mean.nc')
	call stdev_pair%writeVectors('./results/'//prefix//'/stdev.nc')
	
	call plotPair(mean_pair)
	call fullStats(p)
	call plotPair(stdev_pair)
	call show()
	
contains

	subroutine plotPair(self)
		class(pair_t),intent(in)::self
		
		real(wp),dimension(:),allocatable::x,y
		real(wp),dimension(:,:),allocatable::u,v,h
		character(1),dimension(4)::names
		integer,dimension(2)::N,s
		integer::i,k
		
		names = ['U','V','R','N']
		
		do k=lbound(self%passes,1)+1,ubound(self%passes,1)
			
			x = self%vx
			y = self%vy
			u = real(self%passes(k)%u)
			v = real(self%passes(k)%v)
			N = [size(x),size(y)]
			s = N/16+1
			
			call figure()
			call subplot(1,1,1,aspect=1.0_wp)
			call xylim(mixval(x),mixval(y))
			call quiver(x(::s(1)),y(::s(2)),u(::s(1),::s(2)),v(::s(1),::s(2)),lineColor='c')
			call ticks()
			call labels('Position #fix#fn [m]','Position #fiy#fn [m]','Pass '//int2char(k))
			
			call figure()
			
			h = real(self%passes(k)%u)
			call figure()
			call subplot(1,1,1,aspect=1.0_wp)
			call xylim(mixval(self%vx),mixval(self%vy))
			call contourf(self%vx,self%vy,h,20)
			call colorbar2(h,20)
			call ticks()
			call labels('Position #fix#fn [m]','Position #fiy#fn [m]','Pass '//int2char(k)//' u')
			
			h = real(self%passes(k)%v)
			call figure()
			call subplot(1,1,1,aspect=1.0_wp)
			call xylim(mixval(self%vx),mixval(self%vy))
			call contourf(self%vx,self%vy,h,20)
			call colorbar2(h,20)
			call ticks()
			call labels('Position #fix#fn [m]','Position #fiy#fn [m]','Pass '//int2char(k)//' v')
			
			h = real(self%passes(0)%u-self%passes(k)%u)
			call figure()
			call subplot(1,1,1,aspect=1.0_wp)
			call xylim(mixval(self%vx),mixval(self%vy))
			call contourf(self%vx,self%vy,h,20)
			call colorbar2(h,20)
			call ticks()
			call labels('Position #fix#fn [m]','Position #fiy#fn [m]','Pass '//int2char(k)//' #ge#du#u')
			
			h = real(self%passes(0)%v-self%passes(k)%v)
			call figure()
			call subplot(1,1,1,aspect=1.0_wp)
			call xylim(mixval(self%vx),mixval(self%vy))
			call contourf(self%vx,self%vy,h,20)
			call colorbar2(h,20)
			call ticks()
			call labels('Position #fix#fn [m]','Position #fiy#fn [m]','Pass '//int2char(k)//' #ge#dv#u')
		
			do i=1,size(names)
				h = der(self%passes(k)%u,i)
				call figure()
				call subplot(1,1,1,aspect=1.0_wp)
				call xylim(mixval(self%vx),mixval(self%vy))
				call contourf(self%vx,self%vy,h,20)
				call colorbar2(h,20)
				call ticks()
				call labels('Position #fix#fn [m]','Position #fiy#fn [m]','Pass '//int2char(k)//' du/d'//names(i))
				
				h = der(self%passes(k)%v,i)
				call figure()
				call subplot(1,1,1,aspect=1.0_wp)
				call xylim(mixval(self%vx),mixval(self%vy))
				call contourf(self%vx,self%vy,h,20)
				call colorbar2(h,20)
				call ticks()
				call labels('Position #fix#fn [m]','Position #fiy#fn [m]','Pass '//int2char(k)//' dv/d'//names(i))
			end do
			
		end do
	end subroutine plotPair

	subroutine pairStats(self)
		class(pair_t),intent(in)::self
		
		real(wp),dimension(:),allocatable::h,he
		integer::k
		
		do k=1,ubound(self%passes,1)
			h  = flatten(real( self%passes(k)%u ))
			he = flatten(real( self%passes(k)%u-self%passes(0)%u ))
			call doHist( pack(h ,abs(he)<2.0_wp) ,'Displacement #fi#gd#dx#u#fn [px]',k)
			call doHist( pack(he,abs(he)<2.0_wp) ,'Displacement Error #fi#ge#dx#u#fn [px]',k)
			
			h  = flatten(real( self%passes(k)%v ))
			he = flatten(real( self%passes(k)%v-self%passes(0)%v ))
			call doHist( pack(h ,abs(he)<2.0_wp) ,'Displacement #fi#gd#dy#u#fn [px]',k)
			call doHist( pack(he,abs(he)<2.0_wp) ,'Displacement Error #fi#ge#dy#u#fn [px]',k)
			
			h  = flatten(der( self%passes(k)%u , 1 ))
			he = flatten(der( self%passes(k)%u-self%passes(0)%u , 1 ))
			call doHist( h  ,'Displacement Derivative #fid#gd#dx#u/dU#dx#u#fn [px]',k)
			call doHist( he ,'Displacement Error Derivative #fid#ge#dx#u/dU#dx#u#fn [px]',k)
			
			h  = flatten(der( self%passes(k)%u , 2 ))
			he = flatten(der( self%passes(k)%u-self%passes(0)%u , 2 ))
			call doHist( h  ,'Displacement Derivative #fid#gd#dx#u/dU#dy#u#fn [px]',k)
			call doHist( he ,'Displacement Error Derivative #fid#ge#dx#u/dU#dy#u#fn [px]',k)
			
			h  = flatten(der( self%passes(k)%v , 1 ))
			he = flatten(der( self%passes(k)%v-self%passes(0)%v , 1 ))
			call doHist( h  ,'Displacement Derivative #fid#gd#dy#u/dU#dx#u#fn [px]',k)
			call doHist( he ,'Displacement Error Derivative #fid#ge#dy#u/dU#dx#u#fn [px]',k)
			
			h  = flatten(der( self%passes(k)%v , 2 ))
			he = flatten(der( self%passes(k)%v-self%passes(0)%v , 2 ))
			call doHist( h  ,'Displacement Derivative #fid#gd#dy#u/dU#dy#u#fn [px]',k)
			call doHist( he ,'Displacement Error Derivative #fid#ge#dy#u/dU#dy#u#fn [px]',k)
		end do
	
	end subroutine pairStats

	subroutine doHist(h,L,k)
		real(wp),dimension(:),intent(in)::h
		character(*),intent(in)::L
		integer,intent(in)::k
		
		integer::Nb
		
		Nb = nint(0.5_wp*(real(size(h),wp))**(0.5_wp))
		
		call figure()
		call subplot(1,1,1)
		call xylim(mixval(h),[0.0_wp,1.05_wp])
		call hist(h,Nb,relWidth=0.7_wp)
		call xticks(primary=.true.,secondary=.false.)
		call labels(L,'','Pass '//int2char(k)//'; N = '//int2char(size(h)))
	end subroutine doHist

	subroutine fullStats(pairs)
		type(pair_t),dimension(:),intent(in)::pairs
		
		real(wp),dimension(:),allocatable::h,he,hd
		integer::j,k
		
		! u value
		do k=1,size(pairs(1)%passes)-1
			h  = [real(wp)::]
			he = [real(wp)::]
			do j=1,size(pairs)
				h  = [h, real( pack( pairs(j)%passes(k)%u , pairs(j)%passes(k)%mask ) ) ]
				he = [he,real( pack( pairs(j)%passes(0)%u , pairs(j)%passes(k)%mask ) ) ]
			end do
			hd = h-he
			call doHist( h  ,'Displacement #fi#gd#dx#u#fn [px]',k)
			call doHist( hd ,'Displacement Error #fi#ge#dx#u#fn [px]',k)
		end do
		
		! v value
		do k=1,size(pairs(1)%passes)-1
			h  = [real(wp)::]
			he = [real(wp)::]
			do j=1,size(pairs)
				h  = [h, real( pack( pairs(j)%passes(k)%v , pairs(j)%passes(k)%mask ) ) ]
				he = [he,real( pack( pairs(j)%passes(0)%v , pairs(j)%passes(k)%mask ) ) ]
			end do
			hd = h-he
			call doHist( h  ,'Displacement #fi#gd#dy#u#fn [px]',k)
			call doHist( hd ,'Displacement Error #fi#ge#dy#u#fn [px]',k)
		end do
		
		! dudU
		do k=1,size(pairs(1)%passes)-1
			h  = [real(wp)::]
			he = [real(wp)::]
			do j=1,size(pairs)
				h  = [h, der( pack( pairs(j)%passes(k)%u , pairs(j)%passes(k)%mask ) , 1 ) ]
				he = [he,der( pack( pairs(j)%passes(0)%u , pairs(j)%passes(k)%mask ) , 1 ) ]
			end do
			hd = h-he
			call doHist( h  ,'Displacement Derivative #fid#gd#dx#u/dU#fn [px]',k)
			call doHist( hd ,'Displacement Error Derivative #fid#ge#dx#u/dU#fn [px]',k)
		end do
		
		! dvdU
		do k=1,size(pairs(1)%passes)-1
			h  = [real(wp)::]
			he = [real(wp)::]
			do j=1,size(pairs)
				h  = [h, der( pack( pairs(j)%passes(k)%v , pairs(j)%passes(k)%mask ) , 1 ) ]
				he = [he,der( pack( pairs(j)%passes(0)%v , pairs(j)%passes(k)%mask ) , 1 ) ]
			end do
			hd = h-he
			call doHist( h  ,'Displacement Derivative #fid#gd#dy#u/dU#fn [px]',k)
			call doHist( hd ,'Displacement Error Derivative #fid#ge#dy#u/dU#fn [px]',k)
		end do
		
		! dudV
		do k=1,size(pairs(1)%passes)-1
			h  = [real(wp)::]
			he = [real(wp)::]
			do j=1,size(pairs)
				h  = [h, der( pack( pairs(j)%passes(k)%u , pairs(j)%passes(k)%mask ) , 2 ) ]
				he = [he,der( pack( pairs(j)%passes(0)%u , pairs(j)%passes(k)%mask ) , 2 ) ]
			end do
			hd = h-he
			call doHist( h  ,'Displacement Derivative #fid#gd#dx#u/dV#fn [px]',k)
			call doHist( hd ,'Displacement Error Derivative #fid#ge#dx#u/dV#fn [px]',k)
		end do
		
		! dvdV
		do k=1,size(pairs(1)%passes)-1
			h  = [real(wp)::]
			he = [real(wp)::]
			do j=1,size(pairs)
				h  = [h, der( pack( pairs(j)%passes(k)%v , pairs(j)%passes(k)%mask ) , 2 ) ]
				he = [he,der( pack( pairs(j)%passes(0)%v , pairs(j)%passes(k)%mask ) , 2 ) ]
			end do
			hd = h-he
			call doHist( h  ,'Displacement Derivative #fid#gd#dy#u/dV#fn [px]',k)
			call doHist( hd ,'Displacement Error Derivative #fid#ge#dy#u/dV#fn [px]',k)
		end do
		
		! dudR
		do k=1,size(pairs(1)%passes)-1
			h  = [real(wp)::]
			he = [real(wp)::]
			do j=1,size(pairs)
				h  = [h, der( pack( pairs(j)%passes(k)%u , pairs(j)%passes(k)%mask ) , 3 ) ]
				he = [he,der( pack( pairs(j)%passes(0)%u , pairs(j)%passes(k)%mask ) , 3 ) ]
			end do
			hd = h-he
			call doHist( h  ,'Displacement Derivative #fid#gd#dx#u/dR#fn [px]',k)
			call doHist( hd ,'Displacement Error Derivative #fid#ge#dx#u/dR#fn [px]',k)
		end do
		
		! dvdR
		do k=1,size(pairs(1)%passes)-1
			h  = [real(wp)::]
			he = [real(wp)::]
			do j=1,size(pairs)
				h  = [h, der( pack( pairs(j)%passes(k)%v , pairs(j)%passes(k)%mask ) , 3 ) ]
				he = [he,der( pack( pairs(j)%passes(0)%v , pairs(j)%passes(k)%mask ) , 3 ) ]
			end do
			hd = h-he
			call doHist( h  ,'Displacement Derivative #fid#gd#dy#u/dR#fn [px]',k)
			call doHist( hd ,'Displacement Error Derivative #fid#ge#dy#u/dR#fn [px]',k)
		end do
		
		! dudN
		do k=1,size(pairs(1)%passes)-1
			h  = [real(wp)::]
			he = [real(wp)::]
			do j=1,size(pairs)
				h  = [h, der( pack( pairs(j)%passes(k)%u , pairs(j)%passes(k)%mask ) , 4 ) ]
				he = [he,der( pack( pairs(j)%passes(0)%u , pairs(j)%passes(k)%mask ) , 4 ) ]
			end do
			hd = h-he
			call doHist( h  ,'Displacement Derivative #fid#gd#dx#u/dN#fn [px]',k)
			call doHist( hd ,'Displacement Error Derivative #fid#ge#dx#u/dN#fn [px]',k)
		end do
		
		! dvdN
		do k=1,size(pairs(1)%passes)-1
			h  = [real(wp)::]
			he = [real(wp)::]
			do j=1,size(pairs)
				h  = [h, der( pack( pairs(j)%passes(k)%v , pairs(j)%passes(k)%mask ) , 4 ) ]
				he = [he,der( pack( pairs(j)%passes(0)%v , pairs(j)%passes(k)%mask ) , 4 ) ]
			end do
			hd = h-he
			call doHist( h  ,'Displacement Derivative #fid#gd#dy#u/dN#fn [px]',k)
			call doHist( hd ,'Displacement Error Derivative #fid#ge#dy#u/dN#fn [px]',k)
		end do
	end subroutine fullStats

	subroutine doMask(p,tol)
		type(pair_t),intent(inout)::p
		real(wp),intent(in)::tol
		
		type(ad_t),dimension(2)::u,t
		integer::i,j,k
		
		do k=1,size(p%passes)-1
			if(allocated(p%passes(k)%mask)) deallocate(p%passes(k)%mask)
			allocate(p%passes(k)%mask(size(p%vx),size(p%vy)))
			do j=1,size(p%vy)
				do i=1,size(p%vx)
					u = [p%passes(k)%u(i,j),p%passes(k)%v(i,j)]
					t = [p%passes(0)%u(i,j),p%passes(0)%v(i,j)]
					p%passes(k)%mask(i,j) = norm2(real(u-t))/norm2(real(t))<tol
				end do
			end do
		end do
	end subroutine doMask

	function meanPair(p) result(o)
		type(pair_t),dimension(:),intent(in)::p
		type(pair_t)::o
		
		integer,dimension(:,:),allocatable::c
		integer::j,k
		
		o = p(1)
		allocate(c(size(o%vx),size(o%vy)))
		c = 0
		
		do k=1,size(o%passes)-1
			o%passes(k)%u = 0.0_wp
			o%passes(k)%v = 0.0_wp
			o%passes(k)%mask = .true.
		end do
		
		do j=1,size(o%passes)-1
			c = 0
			do k=1,size(p)
				where(p(k)%passes(j)%mask)
					o%passes(j)%u = o%passes(j)%u+p(k)%passes(j)%u
					o%passes(j)%v = o%passes(j)%v+p(k)%passes(j)%v
					c = c+1
				end where
			end do
			o%passes(j)%u = o%passes(j)%u/real(c,wp)
			o%passes(j)%v = o%passes(j)%v/real(c,wp)
		end do
	end function meanPair

	function stdevPair(p,m) result(o)
		type(pair_t),dimension(:),intent(in)::p
		type(pair_t),intent(in)::m
		type(pair_t)::o
		
		integer,dimension(:,:),allocatable::c
		integer::j,k
		
		o = p(1)
		allocate(c(size(o%vx),size(o%vy)))
		c = 0
		
		do k=1,size(o%passes)-1
			o%passes(k)%u = 0.0_wp
			o%passes(k)%v = 0.0_wp
			o%passes(k)%mask = .true.
		end do
		
		do j=1,size(o%passes)-1
			c = 0
			do k=1,size(p)
				where(p(k)%passes(j)%mask)
					o%passes(j)%u = o%passes(j)%u+(p(k)%passes(j)%u-m%passes(j)%u)**2
					o%passes(j)%v = o%passes(j)%v+(p(k)%passes(j)%v-m%passes(j)%v)**2
					c = c+1
				end where
			end do
			o%passes(j)%u = o%passes(j)%u/real(c-1,wp)
			o%passes(j)%v = o%passes(j)%v/real(c-1,wp)
		end do
	end function stdevPair

end program post_prg
 
