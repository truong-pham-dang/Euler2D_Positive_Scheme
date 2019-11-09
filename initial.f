	subroutine initial(U)
	implicit real*8 (a-h,m-z)
	parameter (lmx=400,lmy=400)
	common /gridI/lx,ly,lsteps
      common /gridR/T,dt,dx,dy,xlambda,ylambda
	common /conf/config
      dimension U(4,-1:lmx+3,-1:lmy+3)
	real*8, dimension(4) :: vec1,vec2,vec3,vec4
	integer :: config
	lx_half=lx/2
	ly_half=ly/2
	
c-----Configuration 3--------------------------------------------------
c
	if (config .eq. 3) then
c-----Initial data-----------------------------------------------------
c
c-----First quadrant---------------------------------------------------
	rho1 = 1.5d0
	p1   = 1.5d0
	u1   = 0.0d0
	v1   = 0.0d0 
	e1   = total_energy(rho1,u1,v1,p1)
	vec1 = (/ rho1, rho1*u1, rho1*v1, e1/)
c-----Second quadrant--------------------------------------------------
	rho2 = 0.5323d0
	p2   = 0.3d0
	u2   = 1.206d0
	v2   = 0.0d0
	e2   = total_energy(rho2,u2,v2,p2)
	vec2 = (/ rho2, rho2*u2, rho2*v2, e2/)
c-----Third quadrant---------------------------------------------------
	rho3 = 0.138d0
	p3   = 0.029d0
	u3   = 1.206d0
	v3   = 1.206d0
	e3   = total_energy(rho3,u3,v3,p3)
	vec3 = (/ rho3, rho3*u3, rho3*v3, e3/)
c-----Fourth quadrant--------------------------------------------------
	rho4 = 0.5323d0
	p4   = 0.3d0
	u4   = 0.0d0
	v4   = 1.206d0
	e4   = total_energy(rho4,u4,v4,p4)
	vec4 = (/ rho4, rho4*u4, rho4*v4, e4/)
	endif
c
	if (config .eq. 5) then
c-----Initial data-----------------------------------------------------
c
c-----First quadrant---------------------------------------------------
	rho1 = 1.0d0
	p1   = 1.0d0
	u1   = -0.75d0
	v1   = -0.50d0 
	e1   = total_energy(rho1,u1,v1,p1)
	vec1 = (/ rho1, rho1*u1, rho1*v1, e1/)
c-----Second quadrant--------------------------------------------------
	rho2 = 2.0d0
	p2   = 1.0d0
	u2   = -0.75d0
	v2   = 0.50d0
	e2   = total_energy(rho2,u2,v2,p2)
	vec2 = (/ rho2, rho2*u2, rho2*v2, e2/)
c-----Third quadrant---------------------------------------------------
	rho3 = 1.0d0
	p3   = 1.0d0
	u3   = 0.75d0
	v3   = 0.50d0
	e3   = total_energy(rho3,u3,v3,p3)
	vec3 = (/ rho3, rho3*u3, rho3*v3, e3/)
c-----Fourth quadrant--------------------------------------------------
	rho4 = 3.0d0
	p4   = 1.0d0
	u4   = 0.75d0
	v4   = -0.50d0
	e4   = total_energy(rho4,u4,v4,p4)
	vec4 = (/ rho4, rho4*u4, rho4*v4, e4/)
	endif
c
c-----Initial solution-------------------------------------------------
c
c-----First quadrant---------------------------------------------------
	do j = ly_half+1, ly+3
		do i = lx_half+1, lx+3
			U(1:4,i,j) = vec1(1:4)
		enddo
	enddo
c-----Second quadrant--------------------------------------------------
	do j = ly_half+1, ly+3
		do i = -1, lx_half
			U(1:4,i,j) = vec2(1:4)
		enddo
	enddo
c-----Third quadrant---------------------------------------------------
	do j = -1, ly_half
		do i = -1, lx_half
			U(1:4,i,j) = vec3(1:4)
		enddo
	enddo
c-----Fourth quadrant--------------------------------------------------
	do j = -1, ly_half
		do i = lx_half+1, lx+3
			U(1:4,i,j) = vec4(1:4)
		enddo
      enddo
	end
