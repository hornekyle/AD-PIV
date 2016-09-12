module settings_mod
	!! Module containting global settings
	use kinds_mod
	use autodiff_mod
	use config_mod
	use text_mod
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
	
	real(wp),dimension(:),allocatable::Ux0
	real(wp),dimension(:),allocatable::Uy0
	
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
	
	real(wp)::correlationFactor = 0.4
	integer::lsqOrder = 1
	
	logical::write_pair = .true.
	logical::write_map  = .true.
	logical::per_pixel  = .true.
	
contains

	subroutine readConfig(fn)
		character(*),intent(in)::fn
		type(config_t)::cfg
		integer::k
		
		cfg = config_t(fn)
		
		prefix = trim(cfg%getString('prefix'))
		
		Ux0 = cfg%getVector('Ux0')
		Uy0 = cfg%getVector('Uy0')
		
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
		pass_types = [character(3):: ( cfg%getString('pass_types['//intToChar(k)//']') , k=1,N_passes )]
		
		write_map = cfg%getLogical('write_map')
		per_pixel = cfg%getLogical('per_pixel')
	end subroutine readConfig

	function uf(x) result(o)
		type(ad1_t),dimension(2),intent(in)::x
		type(ad1_t),dimension(2)::o
		
		type(ad1_t),dimension(2)::xl
		
		xl = x-real(image_size,wp)/2.0_wp
		
		o(1) = diff1(Ux0(1),ADS_U)+Ux0(2)*xl(1)
		o(2) = diff1(Uy0(1),ADS_V)+Uy0(2)*xl(2)
	end function uf

end module settings_mod
