module settings_mod
	!! Module containting global settings
	use autodiff_mod
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
	
	
	character(3),dimension(ADS_COUNT),parameter::ADS_CHS = &
		& ['dU ','dUx','dUy','dV ','dVx','dVy','dR ','dN ']
	
	!=============!
	!= Variables =!
	!=============!
	
	real(wp),dimension(:),allocatable::U0
	real(wp),dimension(:,:),allocatable::vU
	
	real(wp)::noise_level
	
	integer,dimension(2)::image_size
	integer::particle_count
	real(wp)::particle_radius
	integer::N_pairs
	integer::pairs_start
	
	integer::N_Passes
	integer,dimension(2)::buffer_window_size
	integer,dimension(2)::spacing_window_size
	integer,dimension(:),allocatable::pass_guesses
	integer,dimension(:,:),allocatable::pass_sizes
	character(3),dimension(:),allocatable::pass_types
	integer,dimension(2)::max_pass_sizes
	
	real(wp),dimension(:),allocatable::correlationFactors
	
	integer::lsqOrder = 1
	
	logical::write_pair = .true.
	logical::write_vectors = .true.
	logical::write_map  = .true.
	logical::per_pixel  = .true.
	
	integer::adN = ADS_COUNT
	
contains

	subroutine readConfig(fn)
		use config_mod
		use text_mod
		character(*),intent(in)::fn
		type(config_t)::cfg
		integer::k
		
		cfg = config_t(fn)
		
		U0 = cfg%getVector('U0')
		vU = cfg%getMatrix('vU')
		
		noise_level   = cfg%getReal('noise_level')
		
		image_size      = nint(cfg%getVector('image_size'))
		particle_count  = cfg%getInteger('particle_count')
		particle_radius = cfg%getReal('particle_radius')
		N_pairs         = cfg%getInteger('N_pairs')
		pairs_start     = cfg%getInteger('pairs_start')
		
		N_passes = cfg%getInteger('N_passes')
		buffer_window_size  = nint(cfg%getVector('buffer_window_size'))
		spacing_window_size = nint(cfg%getVector('spacing_window_size'))
		pass_guesses        = nint(cfg%getVector('pass_guesses'))
		pass_sizes          = nint(cfg%getMatrix('pass_sizes'))
		pass_types = [character(3):: ( cfg%getString('pass_types['//intToChar(k)//']') , k=1,N_passes )]
		correlationFactors = cfg%getVector('correlationFactors')
		
		write_pair = cfg%getLogical('write_pair')
		write_map  = cfg%getLogical('write_map')
		per_pixel  = cfg%getLogical('per_pixel')
		
		max_pass_sizes = [(pass_sizes(1,:)) , maxval(pass_sizes(2,:))]
		if(per_pixel) then
			adN = ADS_COUNT+2*product(max_pass_sizes)
		end if
	end subroutine readConfig

	function uf(x) result(o)
		type(ad_t),dimension(2),intent(in)::x
		type(ad_t),dimension(2)::o
		
		type(ad_t),dimension(2)::xl,U0l
		type(ad_t),dimension(2,2)::vUl
		
		xl  = x-real(image_size,wp)/2.0_wp
		
		U0l(1) = ad_t(U0(1),ADS_COUNT,ADS_U)
		U0l(2) = ad_t(U0(2),ADS_COUNT,ADS_V)
		
		vUl(1,1) = ad_t(vU(1,1),ADS_COUNT,ADS_Ux)
		vUl(1,2) = ad_t(vU(1,2),ADS_COUNT,ADS_Uy)
		vUl(2,1) = ad_t(vU(2,1),ADS_COUNT,ADS_Vx)
		vUl(2,2) = ad_t(vU(2,2),ADS_COUNT,ADS_Vy)
		
		o = U0l+matmul(vUl,xl)
	end function uf

end module settings_mod
