program post_prg
	use kinds_mod
	use autodiff_mod
	use utilities_mod
	use plplotlib_mod
	use pair_mod
	use piv_mod
	implicit none
	
	type(pair_t),dimension(:),allocatable::p
	character(:),allocatable::pfn,vfn
	integer::j,k
	
	N_pairs = 4
	allocate(p(N_pairs))
	do k=1,N_pairs
		pfn = 'pair-'//int2char(k)//'.nc'
		vfn = 'vectors-'//int2char(k)//'.nc'
		
		call p(k)%readPair(pfn,vfn)
		deallocate(p(k)%A,p(k)%B)
		
		do j=1,size(p(k)%passes)-1
			call filter(p(k),j,0.1_wp)
		end do
	end do
	
	call setup(device='svgqt',filename='output-post-%n.svg',colormap='BlueYellow')
	call fullStats(p)
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
			
			call figure()
			call subplot(1,1,1,aspect=1.0_wp)
			call xylim(mixval(self%px),mixval(self%py))
			if(allocated(self%vx) .and. allocated(self%vy)) then
				x = self%vx
				y = self%vy
				u = real(self%passes(k)%u)
				v = real(self%passes(k)%v)
				N = [size(x),size(y)]
				s = N/16+1
				call quiver(x(::s(1)),y(::s(2)),u(::s(1),::s(2)),v(::s(1),::s(2)),lineColor='c')
			else
				call contourf(self%px,self%py,real(self%B)-real(self%A),10)
			end if
			call ticks()
			call labels('Position #fix#fn [m]','Position #fiy#fn [m]','Pass '//int2char(k))
			
			call figure()
			
			if(allocated(self%vx) .and. allocated(self%vy) .and. k>0) then
				h = real(self%passes(k)%u)
				call figure()
				call subplot(1,1,1,aspect=1.0_wp)
				call xylim(mixval(self%px),mixval(self%py))
				call contourf(self%vx,self%vy,h,20)
				call colorbar2(h,20)
				call ticks()
				call labels('Position #fix#fn [m]','Position #fiy#fn [m]','Pass '//int2char(k)//' u')
				
				h = real(self%passes(k)%v)
				call figure()
				call subplot(1,1,1,aspect=1.0_wp)
				call xylim(mixval(self%px),mixval(self%py))
				call contourf(self%vx,self%vy,h,20)
				call colorbar2(h,20)
				call ticks()
				call labels('Position #fix#fn [m]','Position #fiy#fn [m]','Pass '//int2char(k)//' v')
				
				h = real(self%passes(0)%u-self%passes(k)%u)
				call figure()
				call subplot(1,1,1,aspect=1.0_wp)
				call xylim(mixval(self%px),mixval(self%py))
				call contourf(self%vx,self%vy,h,20)
				call colorbar2(h,20)
				call ticks()
				call labels('Position #fix#fn [m]','Position #fiy#fn [m]','Pass '//int2char(k)//' #ge#du#u')
				
				h = real(self%passes(0)%v-self%passes(k)%v)
				call figure()
				call subplot(1,1,1,aspect=1.0_wp)
				call xylim(mixval(self%px),mixval(self%py))
				call contourf(self%vx,self%vy,h,20)
				call colorbar2(h,20)
				call ticks()
				call labels('Position #fix#fn [m]','Position #fiy#fn [m]','Pass '//int2char(k)//' #ge#dv#u')
			
				do i=1,size(names)
					h = der(self%passes(k)%u,i)
					call figure()
					call subplot(1,1,1,aspect=1.0_wp)
					call xylim(mixval(self%px),mixval(self%py))
					call contourf(self%vx,self%vy,h,20)
					call colorbar2(h,20)
					call ticks()
					call labels('Position #fix#fn [m]','Position #fiy#fn [m]','Pass '//int2char(k)//' du/d'//names(i))
					
					h = der(self%passes(k)%v,i)
					call figure()
					call subplot(1,1,1,aspect=1.0_wp)
					call xylim(mixval(self%px),mixval(self%py))
					call contourf(self%vx,self%vy,h,20)
					call colorbar2(h,20)
					call ticks()
					call labels('Position #fix#fn [m]','Position #fiy#fn [m]','Pass '//int2char(k)//' dv/d'//names(i))
				end do
			end if
			
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
		
		Nb = nint((real(size(h),wp))**(0.5_wp))
		
		call figure()
		call subplot(1,1,1)
		call xylim(mixval(h),[0.0_wp,1.05_wp])
		call hist(h,Nb)
		call xticks(primary=.true.,secondary=.false.)
		call labels(L,'','Pass '//int2char(k)//'; N = '//int2char(size(h)))
	end subroutine doHist

	subroutine fullStats(pairs)
		type(pair_t),dimension(:),intent(in)::pairs
		
		real(wp),dimension(:),allocatable::h,he,hd
		integer::j,k
		
		h  = [real(wp)::]
		he = [real(wp)::]
		do k=1,size(pairs(1)%passes)-1
			do j=1,size(pairs)
				h  = [h, real( pairs(j)%passes(k)%u )  ]
				he = [he,real( pairs(j)%passes(0)%u ) ]
			end do
		end do
		hd = h-he
		call doHist( h  ,'Displacement #fi#gd#dx#u#fn [px]',0)
		call doHist( hd ,'Displacement Error Derivative #fid#ge#dx#u/dU#dx#u#fn [px]',k)
	end subroutine fullStats
	
end program post_prg
 
