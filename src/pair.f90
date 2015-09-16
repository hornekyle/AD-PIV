module pair_mod
	use kinds_mod
	use utilities_mod
	use autodiff_mod
	use netCDF_mod
	use plplotlib_mod
	implicit none
	
	type::pair_t
		real(wp),dimension(:),allocatable::x,y
		real(wp),dimension(:),allocatable::vx,vy
		type(ad_t),dimension(:,:),allocatable::A,B
		type(ad_t),dimension(:,:),allocatable::u,v
		type(ad_t),dimension(:,:),allocatable::u2,v2
		real(wp),dimension(2)::d
		real(wp)::dt
	contains
		procedure::writeNC
	end type
	
contains

	function newPair(N,L,dt) result(o)
		integer,dimension(2),intent(in)::N
		real(wp),dimension(2),intent(in),optional::L
		real(wp),intent(in),optional::dt
		type(pair_t)::o
		
		allocate(o%x(N(1)))
		allocate(o%y(N(2)))
		allocate(o%A(N(1),N(2)))
		allocate(o%B(N(1),N(2)))
		
		if(present(L)) then
			o%x = linspace(0.0_wp,L(1),N(1))
			o%y = linspace(0.0_wp,L(2),N(2))
			o%d = [o%x(2)-o%x(1),o%y(2)-o%y(1)]
		else
			o%x = 0.0_wp
			o%y = 0.0_wp
			o%d = 0.0_wp
		end if
		
		if(present(dt)) then
			o%dt = dt
		else
			o%dt = 0.0_wp
		end if
	end function newPair

	subroutine writeNC(self,fn)
		class(pair_t),intent(in)::self
		character(*),intent(in)::fn
		
		call write_grid(fn,['I'],self%x,self%y)
		call write_step(fn,0.0_wp,1,'I',real(self%A))
		call write_step(fn,self%dt,2,'I',real(self%B))
	end subroutine writeNC

	subroutine plotPair(p)
		type(pair_t),intent(in)::p
		real(wp),dimension(:),allocatable::h,h1,h2
		
		call figure()
		call subplot(1,1,1,aspect=1.0_wp)
		call xylim(mixval(p%x),mixval(p%y))
		call contourf(p%x,p%y,real(p%A),5)
		if(allocated(p%u)) then
			call quiver(p%vx,p%vy,real(p%u),real(p%v),lineColor='k',lineWidth=2.0_wp)
		end if
		call ticks()
		call labels('Position #fix#fn [m]','Position #fiy#fn [m]','Image A')
	end subroutine plotPair

end module pair_mod
