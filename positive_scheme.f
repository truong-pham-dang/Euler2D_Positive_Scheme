
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
      call grid_para
      call initial(U)
      do l=l,lsteps
        call evolve(U,Ul)
        call bdy(Ul)
        call evolve(Ul,Ul)
        do j=l,ly+l
        do i=l,lx+l
        do k=l,4
            U(k,i,j)=0.5d0*(U(k,i,j)+Ul(k,i,j))
        end do
        end do
        end do
        call bdy(U)
      end do
      open(10,file='rho.mat',status='unknown')
      write(10,2000) T,dt,dt/dx,lx,ly
      write(10,1000) ((U(l,i,j),j=l,ly+l),i=l,lx+l)
1000  format(f9.4,f9.4,f9.4,f9.4,f9.4)
2000  format(f9.4,f9.4,f9.4,i5,i5)
      stop
      end

