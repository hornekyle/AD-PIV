program main
	use kinds_mod
	use plplotlib_mod
	use test_mod
	implicit none
	
!~ 	call setup(device='svgqt',fileName='plot-%n.svg',figSize=[500,400],colormap='BlueYellow',transparent=.true.)
	call setup(figSize=[500,400],colormap='BlueYellow')
!~ 	call testSingleCorrelation
	call testPairs
	call show()
	
end program main 
