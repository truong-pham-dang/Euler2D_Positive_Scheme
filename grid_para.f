c**********************************************************************
c
c Subroutine: grid_para
c
c Object: l) To setup grid for the computing domain [0,l]X[0,l] in R^2
c           and time interval [0,T]
c         2) To setup parameters \alpha and \beta
c
c Variables: lx number of grid points in x direction
c ly number of grid points in y direction
c lsteps number of time steps
c T stopping time T
c dt time step
c dx stepsize in x direction
c dy stepsize in y direction
c
c**********************************************************************
c
	subroutine grid_para
	implicit real*8 (a-h,m-z)
	common /gridI/lx,ly,lsteps
	common /gridR/T,dt,dx,dy,xlambda,ylambda
	common /para/alpha,beta
	integer :: lx,ly,lsteps
	real*8  :: T,dt,dx,dy,xlambda,ylambda
	real*8  :: alpha,beta
	write(*,*) ' Enter the Number of Grid Points in x direction '
	read(*,*) lx
	write(*,*) ' Enter the Number of Grid Points in y direction '
	read(*,*) ly
	dx = 1.0d00 / dble(lx)
	dy = 1.0d00 / dble(ly)
	write(*,*) 'Enter Stopping Time, T:'
	read(*,*) T
	write(*,*) 'Enter the minmum of dt/dx and dt/dy:'
	read(*,*) xlambda
	dt=xlambda*dy
	if(xlambda*dx.lt.xlambda*dy)dt=xlambda*dx
	if(T/dt.gt.dble(dint(T/dt)))then
		lsteps=dint(T/dt)+1
	else
		lsteps=dint(T/dt)
	endif
	if(T.eq.0.0d00)then
		dt=0.0d00
	else
		dt=T/dble(lsteps)
	end if
	xlambda=dt/dx
	ylambda=dt/dy
	write(*,*) 'Enter the parameters alpha= and beta='
	read(*,*) alpha,beta
	write(*,*) 'the number of timesteps : ',lsteps
	write(*,*) 'Running...... '
	return
	end