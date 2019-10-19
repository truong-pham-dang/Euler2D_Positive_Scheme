c**********************************************************************
c
c Subroutine: limiter
c
c Object: To evaluate two limiters phi^0 and phi^l from
c         theta=dw/dwup
c
c Variables: phi0 Limiter in the least dissipative flux
c phil Limiter in the more dissipative flux
c dw R^{-l}*(U_{j+l}-U_{j})
c dwup R^{-l}*(U_{j+2}-U_{j+l}) if eig<0, or
c R^{-l}*(U_{j}-U_{j-l}) otherwise
c
c**********************************************************************
c
	subroutine limiter(dw,dwup,phi0,phil,k)
	implicit real*8 (a-h,m-z)
c--------Superbee------------------------------------------------------
	phi0=0.0d0
	if(dw.eq.0.0d0.and.dwup.gt.0.0d0)phi0=2.0d0
	if(dw*dwup.gt.0.0d0)then
		theta=dwup/dw
		if(theta.le.0.5d0)then
			phi0=2.0d0*theta
		else if(theta.le.1.0d0.and.theta.gt.0.5d0)then
			phi0=1.0d0
		else if(theta.le.2.0d0.and.theta.gt.1.0d0)then
			phi0=theta
		else
			phi0=2.0d0
		end if
	end if
c-------VanLeer--------------------------------------------------------
	phi=0.0d0
	if(dw.eq.0.0d0.and.dwup.gt.0.0d0)phi=2.0d0
	if(dw*dwup.gt.0.0d0)then
		theta=dwup/dw
		phi=2.0d0*theta/(1.0d0+theta)
	end if
c-------MinMod---------------------------------------------------------
	phil=0.0d0
	if(dw.eq.0.0d0.and.dwup.gt.0.0d0)phil=1.0d0
	if(dw*dwup.gt.0.0d0)then
		phil=1.0d0
		if(dwup/dw.le.1.0d0)then
			phil=dwup/dw
		end if
	end if
c----------------------------------------------------------------------
	if(k.eq.1.or.k.eq.4)then
		phi0=phi
	end if
	return
	end