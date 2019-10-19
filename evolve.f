c************************************************************************
c
c Subroutine: Evolve
c
c Object: To evaluate the future values: rhol,ml,nl,el
c from the current values: rho, m, n, e
c
c Variables: up U_{j+l}
c um U_{j}
c dup U_{j+2}-U_{j+l}
c dum U_{j}-U_{j-l}
c du U_{j+l}-U_{j}
c dw R^{-l}*(U_{j+l}-U_{j})
c dwf Diffusive flux in char fields
c
c fc Central differencing flux
c df=R*dwf Diffusive flux
c f=fc+df Flux in x direction
c g Flux in y direction
c
c*************************************************************************
c
	subroutine evolve (U,Ul)
	implicit real*8 (a-h,m-z)
	common /gridI/lx,ly,lsteps
	common /gridR/T,dt,dx,dy,xlambda,ylambda
	parameter (lmx=400,lmy=400)
	dimension U(4,-1:lmx+3,-1:lmy+3),Ul(4,-1:lmx+3,-1:lmy+3)
	dimension up(4), um(4),du(4),dup(4),dum(4)
	dimension fc(4),df(4),f(4,0:lmx+1,0:lmy+1),g(4,0:lmx+1,0:lmy+1)
	do j=1,ly+1
		do i=0,lx+1
			do k=1,4
				up(k)=U(k,i+1,j)
				um(k)=U(k,i,j)
				dum(k)=U(k,i,j)-U(k,i-1,j)
				dup(k)=U(k,i+2,j)-U(k,i+1,j)
				du(k)=U(k,i+1,j)-U(k,i,j)
			end do
			call central(up,um,fc)
			call diffusiveflux(up,um,du,dup,dum,df)
			do k=1,4
				f(k,i,j)=fc(k)+df(k)
			end do
		end do
	end do
	do j=0,ly+1
		do i=1,lx+1
			do k=1,4
				kl=k
				if(k.eq.2)kl=3
				if(k.eq.3)kl=2
				up(k)=U(kl,i,j+1)
				um(k)=U(kl,i,j)
				dum(k)=U(kl,i,j)-U(kl,i,j-1)
				dup(k)=U(kl,i,j+2)-U(kl,i,j+1)
				du(k)=U(kl,i,j+1)-U(kl,i,j)
			end do
			call central(up,um,fc)
			call diffusiveflux(up,um,du,dup,dum,df)
			do k=1,4
				kl=k
				if(k.eq.2)kl=3
				if(k.eq.3)kl=2
				g(k,i,j)=fc(kl)+df(kl)
			end do
		end do
	end do
	do j=1,ly+1
		do i=1,lx+1
			do k=1,4
			Ul(k,i,j)=U(k,i,j)-xlambda*(f(k,i,j)-f(k,i-1,j))
     c						  -ylambda*(g(k,i,j)-g(k,i,j-1))
			end do
		end do
	end do
	return
	end