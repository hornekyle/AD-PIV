module settings_mod
	use kinds_mod
	use autodiff_mod
	use utilities_mod
	use config_mod
	implicit none
	
	character(:),allocatable::prefix
	
	real(wp)::Ux0
	real(wp)::Uy0
	
	type(ad_t)::Ux
	type(ad_t)::Uy
	
	real(wp)::Lx
	real(wp)::Ly
	
	character(16)::velocity_mode
	real(wp)::noise_level
	
	integer::image_scale
	integer::N_pairs
	
	integer::N_Passes
	integer,dimension(2)::buffer_window_size
	integer,dimension(2)::spacing_window_size
	integer,dimension(:),allocatable::pass_guesses
	integer,dimension(:,:),allocatable::pass_sizes
	character(3),dimension(:),allocatable::pass_types
	
contains

	subroutine readConfig(fn)
		character(*),intent(in)::fn
		type(config_t)::cfg
		integer::k
		
		cfg = newConfig(fn)
		
		prefix = trim(cfg%getString('prefix'))
		
		Ux0 = cfg%getReal('Ux0')
		Uy0 = cfg%getReal('Uy0')
		Lx = cfg%getReal('Lx')
		Ly = cfg%getReal('Ly')
		
		velocity_mode = trim(cfg%getString('velocity_mode'))
		noise_level   = cfg%getReal('noise_level')
		
		image_scale = cfg%getInteger('image_scale')
		N_pairs     = cfg%getInteger('N_pairs')
		
		N_passes = cfg%getInteger('N_passes')
		buffer_window_size  = nint(cfg%getVector('buffer_window_size'))
		spacing_window_size = nint(cfg%getVector('spacing_window_size'))
		pass_guesses        = nint(cfg%getVector('pass_guesses'))
		pass_sizes          = nint(cfg%getMatrix('pass_sizes'))
		pass_types = [character(3):: ( cfg%getString('pass_types['//int2char(k)//']') , k=1,N_passes )]
	end subroutine readConfig

	function uf(x) result(o)
		type(ad_t),dimension(2),intent(in)::x
		type(ad_t),dimension(2)::o
		
		type(ad_t)::r,s,c,R0
		
		r = sqrt(sum( [x(1)-Lx,x(2)]**2 ))
		R0 = (Lx+Ly)/2.0_wp
		c = (x(1)-Lx)/r
		s = x(2)/r
		
		select case(velocity_mode)
		case('uniform')
			o(1) = Ux
			o(2) = Uy
		case('shear')
			o(1) =  Ux*real(x(2)/R0)
			o(2) =  Uy
		case('vortex')
			o(1) =  Ux*real(r/R0*s)
			o(2) = -Uy*real(r/R0*c)
		end select
	end function uf

end module settings_mod
