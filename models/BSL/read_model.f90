module read_model !!! BSL !!! 


implicit none

private

logical :: ever_read_moho = .false., ever_read_A3d = .false.

public :: get_value, get_value_aniso

contains

! #######################################################
subroutine get_value (r,theta_rad,phi_rad,rho,vp,vs,Qmu,moho)

use module_A3d

implicit none

integer, intent(IN) :: moho
doubleprecision, intent(IN) :: r,theta_rad,phi_rad
doubleprecision, intent(INOUT) :: rho,vp,vs,Qmu

print *, "BSL MODELS ARE ANISOTROPIC !!!"
print *, "Provide A3d.dat with parameter #2 (xi) = 1 if you want an isotropic medium."
stop

end subroutine get_value

! #######################################################
subroutine get_value_aniso (r,theta_rad,phi_rad,rho,A,C,F,L,M,Gc,Gs,Hc,Hs,Bc,Bs,Ec,Es,Qmu,moho)

use module_A3d

implicit none

integer, intent(IN) :: moho
doubleprecision, intent(IN) :: r,theta_rad,phi_rad
doubleprecision, intent(OUT) :: rho,A,C,F,L,M,Gc,Gs,Hc,Hs,Bc,Bs,Ec,Es,Qmu

integer :: ifanis=1, mh
doubleprecision :: rr, vs,vp, xi,fi,eta_aniso, vpv,vph,vsv,vsh

mh=moho
rr=r/1000.0d0

if (ever_read_A3d==.false.) then
    call init_A3d
    ever_read_A3d = .true.
endif
if (mh==1 .or. mh==-1) then
 if (ever_read_moho==.false.) then
     call init_crust
     ever_read_moho = .true.
 endif
endif

Gc=0.d0; Gs=0.d0; Hc=0.d0; Hs=0.d0; Bc=0.d0; Bs=0.d0; Ec=0.d0; Es=0.d0
if (mh==1) then 
    call get_crust(rr,theta_rad,phi_rad,rho,vpv,vph,vsv,vsh,eta_aniso,Qmu)
else
    if (ever_read_moho .and. rr>6346.5)   rr = 6346.5   ! This is to push the mantle up to the Moho in the case of a thin crust (< 24.5km)
    call A3d_full(rr,theta_rad,phi_rad,ifanis,vs,vp,rho,Qmu,mh,xi,fi,eta_aniso,Gc,Gs,Hc,Hs,Bc,Bs)
    vsv = sqrt(3.d0/(xi+2.d0))*vs
    vsh = sqrt(xi)*vsv
    vph = sqrt(5.d0/(fi+4.d0))*vp
    vpv = sqrt(fi)*vph
    rho = rho
endif
A = rho*vph**2
C = rho*vpv**2
L = rho*vsv**2
M = rho*vsh**2
F = eta_aniso*(A-2.d0*L)
Gc = L*Gc
Gs = L*Gs

end subroutine get_value_aniso


! ########################################################
end module read_model
