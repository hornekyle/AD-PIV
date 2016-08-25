module settings_mod
	!! Module containting global settings
	use kinds_mod
	use autodiff_mod
	use utilities_mod
	use config_mod
	implicit none
	
	!==============!
	!= Parameters =!
	!==============!
	
	integer,parameter::ADS_U = 1
	integer,parameter::ADS_V = 2
	integer,parameter::ADS_R = 3
	integer,parameter::ADS_N = 4
	
	!=============!
	!= Variables =!
	!=============!
	
	character(:),allocatable::prefix
	
	real(wp)::Ux0
	real(wp)::Uy0
	
	real(wp)::Sx
	real(wp)::Sy
	
	real(wp)::Lx
	real(wp)::Ly
	
	character(16)::velocity_mode
	real(wp)::noise_level
	
	integer::image_scale
	real(wp)::particles_scale
	real(wp)::particle_radius
	integer::N_pairs
	
	integer::N_Passes
	integer,dimension(2)::buffer_window_size
	integer,dimension(2)::spacing_window_size
	integer,dimension(:),allocatable::pass_guesses
	integer,dimension(:,:),allocatable::pass_sizes
	character(3),dimension(:),allocatable::pass_types
	
	real(wp)::correlationFactor = 0.8
	integer::lsqOrder = 1
	
	logical::write_pair = .true.
	
	logical::write_map = .true.
	integer::write_map_k = 1
	
	logical::per_pixel = .true.
	
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
		
		image_scale     = cfg%getInteger('image_scale')
		particles_scale = cfg%getReal('particles_scale')
		particle_radius = cfg%getReal('particle_radius')
		N_pairs         = cfg%getInteger('N_pairs')
		
		N_passes = cfg%getInteger('N_passes')
		buffer_window_size  = nint(cfg%getVector('buffer_window_size'))
		spacing_window_size = nint(cfg%getVector('spacing_window_size'))
		pass_guesses        = nint(cfg%getVector('pass_guesses'))
		pass_sizes          = nint(cfg%getMatrix('pass_sizes'))
		pass_types = [character(3):: ( cfg%getString('pass_types['//int2char(k)//']') , k=1,N_passes )]
		
		write_map = cfg%getLogical('write_map')
		write_map_k = cfg%getInteger('write_map_k')
		per_pixel = cfg%getLogical('per_pixel')
	end subroutine readConfig

	function uf(x) result(o)
		type(ad1_t),dimension(2),intent(in)::x
		type(ad1_t),dimension(2)::o
		
		type(ad1_t)::r,s,c,R0
		
		r = sqrt(sum( [x(1)-Lx,x(2)]**2 ))
		R0 = (Lx+Ly)/2.0_wp
		c = (x(1)-Lx)/r
		s = x(2)/r
		
		select case(velocity_mode)
		case('uniform')
			o(1) = Sx*diff1(Ux0,ADS_U)
			o(2) = Sy*diff1(Uy0,ADS_V)
		case('shear')
			o(1) =  Sx*diff1(Ux0*real(x(2)/R0),ADS_U)
			o(2) =  Sy*diff1(Uy0,ADS_V)
		case('vortex')
			o(1) =  Sx*diff1(real( Ux0*r/R0*s),ADS_U)
			o(2) =  Sy*diff1(real(-Uy0*r/R0*c),ADS_V)
		end select
	end function uf

end module settings_mod
