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
	
	integer,parameter::ADS_U  = 1
	integer,parameter::ADS_Ux = 2
	integer,parameter::ADS_Uy = 3
	
	integer,parameter::ADS_V  = 4
	integer,parameter::ADS_Vx = 5
	integer,parameter::ADS_Vy = 6
	
	integer,parameter::ADS_R  = 7
	integer,parameter::ADS_N  = 8
	
	integer,parameter::ADS_COUNT = 8
	
	!=============!
	!= Variables =!
	!=============!
	
	character(:),allocatable::prefix
	
	real(wp),dimension(:),allocatable::U0
	real(wp),dimension(:,:),allocatable::vU
	
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
	
	real(wp),dimension(:),allocatable::correlationFactors
	
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
		
		U0 = cfg%getVector('U0')
		vU = cfg%getMatrix('vU')
		
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
		correlationFactors = cfg%getVector('correlationFactors')
		
		write_map = cfg%getLogical('write_map')
		per_pixel = cfg%getLogical('per_pixel')
	end subroutine readConfig

	function uf(x) result(o)
		type(ad1_t),dimension(2),intent(in)::x
		type(ad1_t),dimension(2)::o
		
		type(ad1_t),dimension(2)::xl,U0l
		type(ad1_t),dimension(2,2)::vUl
		
		xl  = x-real(image_size,wp)/2.0_wp
		U0l = [ diff1(U0(1),ADS_U) , diff1(U0(2),ADS_V) ]
		vUl(1,1:2) = [ diff1(vU(1,1),ADS_Ux) , diff1(vU(1,2),ADS_Uy) ]
		vUl(2,1:2) = [ diff1(vU(2,1),ADS_Vx) , diff1(vU(2,2),ADS_Vy) ]
		
		o = U0l+matmul(vUl,xl)
	end function uf

end module settings_mod
