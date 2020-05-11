module read_model !!! Homogeneous !!! 


implicit none

public :: get_value, get_value_aniso


contains

! #######################################################
subroutine get_value (r,theta_rad,phi_rad,rho,vp,vs,Qmu,moho)


implicit none

integer, optional, intent(IN) :: moho
doubleprecision, intent(IN) :: r,theta_rad,phi_rad
doubleprecision, intent(OUT) :: rho,vp,vs,Qmu


vs = 5000.d0
vp = 8000.d0
rho = 3000.d0
Qmu = 300.d0


end subroutine get_value

! #######################################################
subroutine get_value_aniso (r,theta_rad,phi_rad,rho,A,C,F,L,M,Gc,Gs,Hc,Hs,Bc,Bs,Ec,Es,Qmu,moho)


implicit none

integer, optional, intent(IN) :: moho
doubleprecision, intent(IN) :: r,theta_rad,phi_rad
doubleprecision, intent(OUT) :: rho,A,C,F,L,M,Gc,Gs,Hc,Hs,Bc,Bs,Ec,Es,Qmu

doubleprecision :: vp,vs, vpv,vph,vsv,vsh, eta_aniso
doubleprecision, parameter :: anom = 0.02


vs = 5000.d0
vp = 8000.d0
rho = 3000.d0
Qmu = 300.d0

vsv = vs*(1.d0+anom);   vsh = vs*(1.d0-anom)
vpv = vp*(1.d0+anom);   vph = vp*(1.d0-anom)
eta_aniso = 1.d0

A = rho*vph**2
C = rho*vpv**2
L = rho*vsv**2
M = rho*vsh**2
F = eta_aniso*(A-2.d0*L)
Gc=0.d0; Gs=0.d0; Hc=0.d0; Hs=0.d0; Bc=0.d0; Bs=0.d0; Ec=0.d0; Es=0.d0


end subroutine get_value_aniso

! #######################################################
end module read_model
