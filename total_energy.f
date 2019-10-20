	real*8 function total_energy(rho,u,v,p)
	real*8, intent(in) :: rho,u,v,p
	real*8, parameter  :: gamma=1.4d0
	total_energy = p/(gamma-1.0d0) + rho*(u*u+v*v)/2.0d0
	end function