module kinds_mod
	!! Module to manage kinds for the program
	implicit none
	private
	
	!==============!
	!= Real Kinds =!
	!==============!
	
	integer,parameter::sp = selected_real_kind(6)
		!! Single precision
	integer,parameter::dp = selected_real_kind(15)
		!! Double precision
	integer,parameter::ep = selected_real_kind(18)
		!! Extended precision
	integer,parameter::qp = selected_real_kind(32)
		!! Quad precision
	integer,parameter::wp = dp
		!! Set working precision to double
	
	!==================!
	!= Math Constants =!
	!==================!
	
	real(wp),parameter::PI = 4.0_wp*atan(1.0_wp)
		!! Archimedes' constant
	real(wp),parameter::E  = exp(1.0_wp)
		!! Euler's constant
	
	!==================!
	!= Buffer Lengths =!
	!==================!
	
	integer,parameter::str_long = 128
	character(:),parameter::fmt_long = '(1A128)'
	integer,parameter::str_short = 32
	character(:),parameter::fmt_short = '(1A32)'
	
	!===========!
	!= Exports =!
	!===========!
	
	public::sp,dp,ep,qp,wp
	public::PI,E
	
	public::str_long,str_short
	public::fmt_long,fmt_short
	
	public::printTypes
	
contains

	subroutine printTypes
		!! Print the integer kinds for each real type
		write(*,*) 'sp: ',sp
		write(*,*) 'dp: ',dp
		write(*,*) 'ep: ',ep
		write(*,*) 'qp: ',qp
		write(*,*) 'wp: ',wp
	end subroutine printTypes

end module kinds_mod
 
