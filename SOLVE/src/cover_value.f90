subroutine cover_value (x,y,z,rho,A,C,F,L,M,Gc,Gs,Hc,Hs,Bc,Bs,Ec,Es,Qmu,cov,cover)


use angles
use read_model
use tensor_util

implicit none

integer, intent(IN) :: cov
doubleprecision, intent(IN) :: x,y,z
doubleprecision, intent(INOUT) :: rho,A,C,F,L,M,Gc,Gs,Hc,Hs,Bc,Bs,Ec,Es,Qmu
doubleprecision, dimension(1:6,0:1), intent(IN) :: cover

integer :: i,j
doubleprecision :: x0,y0,z0,rho0,A0,C0,F0,L0,M0,Gc0,Gs0,Hc0,Hs0,Bc0,Bs0,Ec0,Es0,Qmu0, &
                   x1,y1,z1,rho1,A1,C1,F1,L1,M1,Gc1,Gs1,Hc1,Hs1,Bc1,Bs1,Ec1,Es1,Qmu1, &
                   lam, mu, coord, delta, wt
doubleprecision, dimension(1:6,1:6) :: Cij


! Defining coord, delta, the coordinates of the border points and the taper coefficient
x0 = x;   y0 = y;   z0 = z
x1 = x;   y1 = y;   z1 = z
if (cov==1 .or. cov==6) then
   z0 = cover(cov,0)
   z1 = cover(cov,1)
   if (z1-z0>=0.d0) then
      coord = z-z0;   delta = z1-z0
   else
      coord = z0-z;   delta = z0-z1
   endif
else if (cov==2 .or. cov==4) then
   y0 = cover(cov,0)
   y1 = cover(cov,1)
   if (y1-y0>=0.d0) then
      coord = y-y0;   delta = y1-y0
   else
      coord = y0-y;   delta = y0-y1
   endif
else if (cov==3 .or. cov==5) then
   x0 = cover(cov,0)
   x1 = cover(cov,1)
   if (x1-x0>=0.d0) then
      coord = x-x0;   delta = x1-x0
   else
      coord = x0-x;   delta = x0-x1
   endif
endif
call wtcoef (coord,0.d0,delta,2.d0*delta,2.d0*delta,wt)


! Elastic properties at the border points

! Point 0
call get_value_aniso (x0,y0,z0,rho0,A0,C0,F0,L0,M0,Gc0,Gs0,Hc0,Hs0,Bc0,Bs0,Ec0,Es0,Qmu0)
           Cij(:,:) = 0.d0
           Cij(1,1) = A0+Bc0+Ec0
           Cij(2,2) = A0-Bc0+Ec0
           Cij(3,3) = C0
           Cij(1,2) = A0-2.d0*M0-Ec0
           Cij(1,3) = F0+Hc0
           Cij(2,3) = F0-Hc0 
           Cij(1,6) = s2*(Bs0/2.d0+Es0) 
           Cij(2,6) = s2*(Bs0/2.d0-Es0)
           Cij(3,6) = s2*Hs0
           Cij(4,4) = 2.d0*(L0-Gc0)
           Cij(5,5) = 2.d0*(L0+Gc0)
           Cij(6,6) = 2.d0*(M0-Ec0)
           Cij(4,5) = 2.d0*Gs0
           do i = 2,6
              do j = 1,i-1
                 Cij(i,j) = Cij(j,i)
              enddo
           enddo
           lam = lambda_from_Cij(Cij)
           mu = mu_from_Cij(Cij)
A0 = lam+2.d0*mu;   C0 = A0;   F0 = lam;   L0 = mu;   M0 = mu
Gc0=0.d0; Gs0=0.d0; Hc0=0.d0; Hs0=0.d0; Bc0=0.d0; Bs0=0.d0; Ec0=0.d0; Es0=0.d0

! Point 1
call get_value_aniso (x1,y1,z1,rho1,A1,C1,F1,L1,M1,Gc1,Gs1,Hc1,Hs1,Bc1,Bs1,Ec1,Es1,Qmu1)


! Computing the interpolated elastic properties
rho = rho0 + (rho1-rho0)*wt
A = A0 + (A1-A0)*wt
C = C0 + (C1-C0)*wt
F = F0 + (F1-F0)*wt
L = L0 + (L1-L0)*wt
M = M0 + (M1-M0)*wt
Gc = Gc0 + (Gc1-Gc0)*wt
Gs = Gs0 + (Gs1-Gs0)*wt
Hc = Hc0 + (Hc1-Hc0)*wt
Hs = Hs0 + (Hs1-Hs0)*wt
Bc = Bc0 + (Bc1-Bc0)*wt
Bs = Bs0 + (Bs1-Bs0)*wt
Ec = Ec0 + (Ec1-Ec0)*wt
Es = Es0 + (Es1-Es0)*wt
Qmu = Qmu0 + (Qmu1-Qmu0)*wt


return
end subroutine cover_value

