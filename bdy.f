	subroutine bdy(U)
	implicit real*8 (a-h,m-z)
	parameter (lmx=400,lmy=400)
	common /gridI/lx,ly,lsteps
      common /gridR/T,dt,dx,dy,xlambda,ylambda
      dimension U(4,-1:lmx+3,-1:lmy+3)
c
	U(1:4,0,:) = U(1:4,1,:)	
	U(1:4,lx+2,:) = U(1:4,lx+1,:) 
	U(1:4,:,0) = U(1:4,:,1) 
	U(1:4,:,ly+2) = U(1:4,:,ly+1)
c
	U(1:4,-1,:) = U(1:4,0,:)
	U(1:4,lx+3,:) = U(1:4,lx+2,:)	
	U(1:4,:,-1) = U(1:4,:,0)
	U(1:4,:,ly+3) = U(1:4,:,ly+2)
	end