module settings_mod
	use kinds_mod
	use autodiff_mod
	use config_mod
	implicit none
	
	real(wp)::Ux0 = 5.0_wp
	real(wp)::Uy0 = 5.0_wp
	
	type(ad_t)::Ux
	type(ad_t)::Uy
	
	real(wp)::Lx = 1.0_wp
	real(wp)::Ly = 1.0_wp
	
	integer::velocity_mode = 3
	
	real(wp)::noise_level = 0.1_wp
	
	integer::image_scale = 10
	
	integer::N_Passes = 2
	integer,dimension(2)::buffer_window_size
	integer,dimension(2)::spacing_window_size
	integer,dimension(:,:),allocatable::pass_sizes
	integer,dimension(:),allocatable::pass_guesses
	character(3),dimension(:),allocatable::pass_types
	
contains

	function uf(x) result(o)
		type(ad_t),dimension(2),intent(in)::x
		type(ad_t),dimension(2)::o
		
		type(ad_t)::r,s,c,R0
		
		r = sqrt(sum( [x(1)-Lx,x(2)]**2 ))
		R0 = (Lx+Ly)/2.0_wp
		c = (x(1)-Lx)/r
		s = x(2)/r
		
		select case(velocity_mode)
		case(1)
			o(1) = Ux
			o(2) = Uy
		case(2)
			o(1) =  Ux*real(x(2)/R0)
			o(2) =  Uy
		case(3)
			o(1) =  Ux*real(r/R0*s)
			o(2) = -Uy*real(r/R0*c)
		end select
	end function uf

end module settings_mod
