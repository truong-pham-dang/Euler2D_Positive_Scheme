c**********************************************************************
c
c Subroutine: central
c
c Object: To calculate the central differencing flux
c
c Variables: fc=0.5d0*(F(U_{j+l})+F(U_{j}))
c
c**********************************************************************
c
	subroutine central(up,um,fc)
	implicit real*8 (a-h,m-z)
	parameter (gamma=1.4d0)
	dimension up(4),um(4),fc(4)
      pl=(gamma-1.0d0)*(um(4)-0.5d0*(um(2)**2+um(3)**2)/um(1))
	pr=(gamma-1.0d0)*(up(4)-0.5d0*(up(2)**2+up(3)**2)/up(1))
	fc(1)=0.5d0*(um(2)+up(2))
	fc(2)=0.5d0*(um(2)**2/um(1)+pl+up(2)**2/up(1)+pr)
	fc(3)=0.5d0*(um(2)*um(3)/um(1)+up(2)*up(3)/up(1))
	fc(4)=0.5d0*((um(4)+pl)*um(2)/um(1)+(up(4)+pr)*up(2)/up(1))
	return
	end