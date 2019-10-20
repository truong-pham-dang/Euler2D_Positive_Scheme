
c*************************************************************************
c
c Program: positive_scheme
c
c Object: Using 2nd order accurate positive schemes to solve
c 2-D gas dynamics problems
c
c Variables: U=(ul,u2,u3,u4)^T Conservative Variables
c ul=rho Density
c u2=m Momentum in x direction
c u3=n Momentum in y direction
c u4=e Total energy
c
c***************************************************************************
c
      program positive_scheme
      implicit real*8 (a-h,m-z)
      parameter (lmx=400,lmy=400)
      common /gridI/lx,ly,lsteps
      common /gridR/T,dt,dx,dy,xlambda,ylambda
      dimension U(4,-1:lmx+3,-1:lmy+3),Ul(4,-1:lmx+3,-1:lmy+3)
	integer      :: lx,ly,lsteps
	real*8       :: T,dt,dx,dy,xlambda,ylambda
	real*8       :: T_simul 
	character*10 :: cstep
	real*8, allocatable :: XC(:),YC(:)
      call grid_para
	ALLOCATE(XC(lx+1),YC(ly+1))
	XC = 0.0d0
	YC = 0.0d0
	DO i = 1, lx+1
		XC(i) = (i-1) * dx
	ENDDO
	DO i = 1, ly+1
		YC(i) = (i-1) * dy
	ENDDO
      U  = 0.0d0
      Ul = 0.0d0
	T_simul = 0.0d0
      call initial(U)
      do l=1,lsteps
		T_simul = T_simul + dt
		call evolve(U,Ul)
		call bdy(Ul)
		call evolve(Ul,Ul)
		do j=1,ly+l
		do i=1,lx+l
		do k=1,4
			U(k,i,j)=0.5d0*(U(k,i,j)+Ul(k,i,j))
		enddo
		enddo
		enddo
		call bdy(U)
c-----Write solution to Tecplot----------------------------------------
		if (mod(l,20)==0) then
			write(*,*) 'T = ',T_simul
			CALL WRITE_SOLUTION_TECPLOT
		endif
      enddo
      stop
      CONTAINS
          SUBROUTINE WRITE_SOLUTION_TECPLOT
          IMPLICIT NONE
          INTEGER :: I1,J1
          WRITE(cstep,'(i5.5)') l
		OPEN(71,FILE='pltflow_'//TRIM(CSTEP)//'.plt') 
		WRITE(71,*) 'TITLE="EULER_FLOW"' 
		WRITE(71,*) 'VARIABLES="X","Y","DENSITY"' 
		WRITE(71,*) 'ZONE T="EULER_FLOW",I=',lx+1,',J=',ly+1,',F=POINT' 
		DO J1=1,ly+1 
		DO I1=1,lx+1  
			WRITE(71,7101) XC(I1),YC(J1),U(1,I1,J1)   
		ENDDO
		ENDDO 
		CLOSE(71) 
7101  FORMAT(20(1X,F18.8))
          END SUBROUTINE
      end

