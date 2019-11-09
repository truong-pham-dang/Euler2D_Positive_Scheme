	subroutine configuration
	common /conf/config
	integer :: config
	write(*,*) 'Enter configuration, config:'
	read(*,*) config 
	end