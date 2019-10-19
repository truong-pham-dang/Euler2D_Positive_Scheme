c**********************************************************************
c
c Subroutine: diffusiveflux
c
c Object: To calculate the diffusive flux df
c
c Variables: up U_{j+l}
c um U_{j}
c du U_{j+l}-U{j}
c dup U_{j+2}-U{j+l}
c dum U_{j}-U_{j-l}
c df Diffusive flux
c dw R^{-l}*(U_{j+l}-U{j})
c dwf Diffusive flux in char form
c r Right eigenvector matrix R of the A
c ri Left eigenvector matrix R^{-l} fo the A
c eig Eigenvalues of the A
c
c**********************************************************************
c
	subroutine diffusiveflux(up,um,du,dup,dum,df)
	implicit real*8 (a-h,m-z)
	common /para/alpha,beta
	dimension up(4),um(4),du(4),dup(4),dum(4),df(4)
	dimension dw(4),dwf(4),r(4,4),ri(4,4),eig(4)
	call eigs(up,um,r,ri,eig)
	mu=dmax1(dabs(eig(1)),dabs(eig(4)))
	do k=1,4
		dw(k)=ri(k,1)*du(1)+ri(k,2)*du(2)+ri(k,3)*du(3)+ri(k,4)*du(4)
		dwup=ri(k,1)*dup(1)+ri(k,2)*dup(2)
     &                       +ri(k,3)*dup(3)+ri(k,4)*dup(4)
		if(eig(k).ge.0.0d0)then
			dwup=ri(k,1)*dum(1)+ri(k,2)*dum(2)
     &                           +ri(k,3)*dum(3)+ri(k,4)*dum(4)
		end if
		call limiter(dw(k),dwup,phi0,phi1,k)
		dwf(k)=-0.5d0*
     &    (alpha*(1.0d0-phi0)*dabs(eig(k))+beta*(1.0d0-phi1)*mu)*dw(k)
	end do
	do k=1,4
		df(k)=r(k,1)*dwf(1)+r(k,2)*dwf(2)+r(k,3)*dwf(3)+r(k,4)*dwf(4)
	end do
	return
	end