program main
	use kinds_mod
	use plplotlib_mod
	use test_mod
	implicit none
	
	call setup(device='pngqt',fileName='plot-%n.png',figSize=[800,600],colormap='BlueYellow',transparent=.false.)
!~ 	call setup(figSize=[500,400],colormap='BlueYellow')
	call test()
	call show()
	
end program main 
