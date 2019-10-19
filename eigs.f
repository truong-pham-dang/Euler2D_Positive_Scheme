c**************************************************************************
c
c Subroutine: eigs
c
c Object: To calculate the right eigenvector matrx, the left
c eigenvector matrx and eigenvalues of Roe matrix
c A=A(U_{j+l},U_{j}), where F(U_{j+l})-F(U_{j})-
c A(U_{j+l},U_{j})*(U_{j+l}-U_{j})
c
c Variables: up U_{j+1}
c um U_{j}
c r Right eigenvector matrix R of the A
c ri Left eigenvector matrix R^{-l} of the A
c eig Eigenvalues of the A
c
c
c**************************************************************************
c
	subroutine eigs(up,um,r,ri,eig)
	implicit real*8 (a-h,m-z)
	parameter (gamma=1.4d0)
	dimension up(4),um(4),r(4,4),ri(4,4),eig(4)
	ul=um(2)/um(1)
	vl=um(3)/um(1)
	Hl=(um(4)+(gamma-1.0d0)*(um(4)-0.5d0*(ul*ul+vl*vl)*um(1)))/um(1)
	u2=up(2)/up(1)
	v2=up(3)/up(1)
	H2=(up(4)+(gamma-1.0d0)*(up(4)-0.5d0*(u2*u2+v2*v2)*up(1)))/up(1)
	w1=dsqrt(um(1))+dsqrt(up(1))
	u=(dsqrt(um(1))*ul+dsqrt(up(1))*u2)/w1
	v=(dsqrt(um(1))*vl+dsqrt(up(1))*v2)/w1
	H=(dsqrt(um(1))*Hl+dsqrt(up(1))*H2)/w1
	q2=u*u+v*v
	c=dsqrt((gamma-1.0d0)*(H-0.5d0*q2))
	r(1,1)=1.0d0
	r(2,1)=u-c
	r(3,1)=v
	r(4,1)=H-u*c
	r(1,2)=0.0d0
	r(2,2)=0.0d0
	r(3,2)=1.0d0
	r(4,2)=v
	r(1,3)=1.0d0
	r(2,3)=u
	r(3,3)=v
	r(4,3)=0.5d0*q2
	r(1,4)=1.0d0
	r(2,4)=u+c
	r(3,4)=v
	r(4,4)=H+u*c
	bl=1.0d0/(H-0.5d0*q2)
	b2=0.5d0*q2*bl
	ri(1,1)=0.5d0*(b2+u/c)
	ri(1,2)=-0.5d0/c-0.5d0*bl*u
	ri(1,3)=-0.50*bl*v
	ri(1,4)=0.5d0*bl
	ri(2,1)=-v
	ri(2,2)=0.0d0
	ri(2,3)=1.0d0
	ri(2,4)=0.0d0
	ri(3,1)=1.0d0-b2
	ri(3,2)=bl*u
	ri(3,3)=bl*v
	ri(3,4)=-bl
	ri(4,1)=0.5d0*(b2-u/c)
	ri(4,2)=0.5d0/c-0.5d0*bl*u
	ri(4,3)=-0.5d0*bl*v
	ri(4,4)=0.5d0*bl
	eig(1)=u-c
	eig(2)=u
	eig(3)=u
	eig(4)=u+c
	return
	end