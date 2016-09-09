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
	
	character(16)::velocity_mode
	real(wp)::noise_level
	
	integer,dimension(2)::image_size
	integer::particle_count
	real(wp)::particle_radius
	integer::N_pairs
	
	integer::N_Passes
	integer,dimension(2)::buffer_window_size
	integer,dimension(2)::spacing_window_size
	integer,dimension(:),allocatable::pass_guesses
	integer,dimension(:,:),allocatable::pass_sizes
	character(3),dimension(:),allocatable::pass_types
	
	real(wp)::correlationFactor = 0.5
	integer::lsqOrder = 1
	
	logical::write_pair = .true.
	logical::write_map  = .true.
	logical::per_pixel  = .true.
	
contains

	subroutine readConfig(fn)
		character(*),intent(in)::fn
		type(config_t)::cfg
		integer::k
		
		cfg = newConfig(fn)
		
		prefix = trim(cfg%getString('prefix'))
		
		Ux0 = cfg%getReal('Ux0')
		Uy0 = cfg%getReal('Uy0')
		
		velocity_mode = trim(cfg%getString('velocity_mode'))
		noise_level   = cfg%getReal('noise_level')
		
		image_size      = nint(cfg%getVector('image_size'))
		particle_count  = cfg%getInteger('particle_count')
		particle_radius = cfg%getReal('particle_radius')
		N_pairs         = cfg%getInteger('N_pairs')
		
		N_passes = cfg%getInteger('N_passes')
		buffer_window_size  = nint(cfg%getVector('buffer_window_size'))
		spacing_window_size = nint(cfg%getVector('spacing_window_size'))
		pass_guesses        = nint(cfg%getVector('pass_guesses'))
		pass_sizes          = nint(cfg%getMatrix('pass_sizes'))
		pass_types = [character(3):: ( cfg%getString('pass_types['//int2char(k)//']') , k=1,N_passes )]
		
		write_map = cfg%getLogical('write_map')
		per_pixel = cfg%getLogical('per_pixel')
	end subroutine readConfig

	function uf(x) result(o)
		type(ad1_t),dimension(2),intent(in)::x
		type(ad1_t),dimension(2)::o
		
		select case(velocity_mode)
		case('uniform')
			o(1) = diff1(Ux0,ADS_U)
			o(2) = diff1(Uy0,ADS_V)
		end select
	end function uf

end module settings_mod
