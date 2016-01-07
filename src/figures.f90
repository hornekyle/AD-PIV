program figures_prg
	use kinds_mod
	use autodiff_mod
	use utilities_mod
	use settings_mod
	use pair_mod
	use generator_mod
	use piv_mod
	use plplotlib_mod
	
	character(:),allocatable::fn
	type(pair_t)::p
	type(map_t)::m
	
	call putSettings()
	call execute_command_line('mkdir -p ./results/'//prefix)
	
	p = createFullPair()
	fn = './results/'//prefix//'/pair.nc'
	call p%writePair(fn)
	
	fn = './results/'//prefix//'/output-figures-%n.svg'
	call setup(device='svgqt',filename=fn,colormap='BlueYellow',figSize=[500,400],transparent=.true.)
	
	call plotPair(p)
	
	do k=1,write_map_k-1
		fn = './results/'//prefix//'/map-'//int2char(k)//'.nc'
		call m%readMap(fn)
		call plotMap(m)
	end do
	
	call show()
contains

	subroutine plotPair(p)
		type(pair_t),intent(in)::p
		
		real(wp),dimension(:),allocatable::x,y
		integer::Nl,i,j
		
		Nl = 25
		x = [( real(i,wp) , i=1,size(p%px) )]
		y = [( real(j,wp) , j=1,size(p%py) )]
		
		call figure()
		call subplot(1,1,1,aspect=span(y)/span(x))
		call xylim(mixval(x),mixval(y))
		call contourf(x,y,real(p%A),Nl)
		call colorbar2(real(p%A),Nl)
		call ticks()
		call labels('Position #fix#fn [px]','Position #fiy#fn [px]','Image Intensity #fiI#fn')
		
		call figure()
		call subplot(1,1,1,aspect=span(y)/span(x))
		call xylim(mixval(x),mixval(y))
		call contourf(x,y,der(p%A,1),Nl)
		call colorbar2(der(p%A,1),Nl)
		call ticks()
		call labels('Position #fix#fn [px]','Position #fiy#fn [px]','Image Intensity Derivative #fi #(2265)I/#(2265)U #fn')
		
		call figure()
		call subplot(1,1,1,aspect=span(y)/span(x))
		call xylim(mixval(x),mixval(y))
		call contourf(x,y,der(p%A,2),Nl)
		call colorbar2(der(p%A,2),Nl)
		call ticks()
		call labels('Position #fix#fn [px]','Position #fiy#fn [px]','Image Intensity Derivative #fi #(2265)I/#(2265)V #fn')
		
		call figure()
		call subplot(1,1,1,aspect=span(y)/span(x))
		call xylim(mixval(x),mixval(y))
		call contourf(x,y,der(p%A,3),Nl)
		call colorbar2(der(p%A,3),Nl)
		call ticks()
		call labels('Position #fix#fn [px]','Position #fiy#fn [px]','Image Intensity Derivative #fi #(2265)I/#(2265)R #fn')
		
! 		call figure()
! 		call subplot(1,1,1,aspect=span(y)/span(x))
! 		call xylim(mixval(x),mixval(y))
! 		call contourf(x,y,der(p%A,4),Nl)
! 		call colorbar2(der(p%A,4),Nl)
! 		call ticks()
! 		call labels('Position #fix#fn [px]','Position #fiy#fn [px]','Image Intensity Derivative #fi #(2265)I/#(2265)N #fn')
	end subroutine plotPair

	subroutine plotMap(m)
		type(map_t),intent(in)::m
		
		integer::Nl
		
		Nl = 25
		
		call figure()
		call subplot(1,1,1,aspect=span(m%dy)/span(m%dx))
		call xylim(mixval(m%dx),mixval(m%dy))
		call contourf(m%dx,m%dy,real(m%C),Nl)
		call colorbar2(real(m%C),Nl)
		call ticks()
		call labels('Displacement #fi#gd#dx#u#fn [px]','Displacement #fi#gd#dy#u#fn [px]','Map Intensity #fiC#fn')
		
		call figure()
		call subplot(1,1,1,aspect=span(m%dy)/span(m%dx))
		call xylim(mixval(m%dx),mixval(m%dy))
		call contourf(m%dx,m%dy,der(m%C,1),Nl)
		call colorbar2(der(m%C,1),Nl)
		call ticks()
		call labels('Displacement #fi#gd#dx#u#fn [px]','Displacement #fi#gd#dy#u#fn [px]', &
			& 'Map Intensity Derivative #fi#(2265)C/#(2265)U#fn')
		
		call figure()
		call subplot(1,1,1,aspect=span(m%dy)/span(m%dx))
		call xylim(mixval(m%dx),mixval(m%dy))
		call contourf(m%dx,m%dy,der(m%C,2),Nl)
		call colorbar2(der(m%C,2),Nl)
		call ticks()
		call labels('Displacement #fi#gd#dx#u#fn [px]','Displacement #fi#gd#dy#u#fn [px]', &
			& 'Map Intensity Derivative #fi#(2265)C/#(2265)V#fn')
		
		call figure()
		call subplot(1,1,1,aspect=span(m%dy)/span(m%dx))
		call xylim(mixval(m%dx),mixval(m%dy))
		call contourf(m%dx,m%dy,der(m%C,3),Nl)
		call colorbar2(der(m%C,3),Nl)
		call ticks()
		call labels('Displacement #fi#gd#dx#u#fn [px]','Displacement #fi#gd#dy#u#fn [px]', &
			& 'Map Intensity Derivative #fi#(2265)C/#(2265)R#fn')
		
! 		call figure()
! 		call subplot(1,1,1,aspect=span(m%dy)/span(m%dx))
! 		call xylim(mixval(m%dx),mixval(m%dy))
! 		call contourf(m%dx,m%dy,der(m%C,4),Nl)
! 		call colorbar2(der(m%C,4),Nl)
! 		call ticks()
! 		call labels('Displacement #fi#gd#dx#u#fn [px]','Displacement #fi#gd#dy#u#fn [px]', &
! 			& 'Map Intensity Derivative #fi#(2265)C/#(2265)N#fn')
	end subroutine plotMap

	subroutine putSettings()
		prefix = 'figures'
		
		Ux0 = 5.0_wp
		Uy0 = 5.0_wp
		Lx = 1.0_wp
		Ly = 1.0_wp
		
		velocity_mode = 'uniform'
		noise_level   = 0.0_wp
		
		image_scale = 7
		N_pairs     = 1
		
		N_passes = 1
		buffer_window_size  = nint([32.0_wp,32.0_wp])
		spacing_window_size = nint([32.0_wp,32.0_wp])
		pass_guesses        = nint([-1.0_wp])
		pass_sizes          = nint(reshape([32.0_wp,32.0_wp],[2,1]))
		
		write_map = .true.
		pass_types = ['map']
	end subroutine putSettings

	function createFullPair() result(p)
		type(pair_t)::p
			!! Result
		
		integer,dimension(2)::N
		real(wp),dimension(2)::L
		type(ad_t),dimension(2)::R
		real(wp)::dt
		integer::Np,k
		
		N = 2**image_scale
		L = 0.0512_wp
		Np = nint(0.05_wp*4.0_wp**(image_scale-2))
		dt = 1.0_wp
		R = L/real(N,wp)*[5.0_wp,0.0_wp]*diff(1.0_wp,3)
		
		Sx = L(1)/real(N(1),wp)
		Sy = L(2)/real(N(2),wp)
		Lx = L(1)
		Ly = L(2)
		
		p = generatePair(N,L,Np,dt,R)
		call p%setupPasses(N_passes,buffer_window_size,spacing_window_size)
		
		call doTrue(p)
		
		do k=1,N_passes
			call doPass(p,k,pass_sizes(:,k),pass_types(k),pass_guesses(k))
		end do
	end function createFullPair

end program figures_prg
