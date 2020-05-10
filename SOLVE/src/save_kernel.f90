subroutine save_kernel (Tdomain,rg,delta_chi)


use sdomains
use angles

implicit none

include 'mpif.h'

type (Domain), intent (INOUT) :: Tdomain
integer, intent(IN) :: rg
doubleprecision, intent(INOUT) :: delta_chi

integer :: n, x,y,z, ngll1,ngll2,ngll3, i, code
doubleprecision :: u,v,w, xa,ya,za, r,theta,phi, lat_deg,phi_deg, delta, &
                   rho,lambda,mu, A,L,M, kern_beta,kern_xi
doubleprecision, dimension(0:2) :: tmp
doubleprecision, dimension (0:2,0:2) :: Rot
character*60 :: fname1, fname2


write (fname1,"(a,I3.3)") "kernel_beta_",rg
open (31,file=trim(fname1))
if (Tdomain%aniso) then
    write (fname2,"(a,I3.3)") "kernel_xi_",rg
    open (32,file=trim(fname2))
endif


delta_chi = 0.d0

do n = 0,Tdomain%n_elem-1

    if (.not. Tdomain%specel(n)%PML) then

        ngll1 = Tdomain%specel(n)%ngllx
        ngll2 = Tdomain%specel(n)%nglly
        ngll3 = Tdomain%specel(n)%ngllz
        do z = 0,ngll3-1 
         do y = 0,ngll2-1 
          do x = 0,ngll1-1

              ! Taking the cartesian coordinates of the GLL
              i = Tdomain%specel(n)%Iglobnum(x,y,z)
              u = Tdomain%Globcoord(0,i)
              v = Tdomain%Globcoord(1,i)
              w = Tdomain%Globcoord(2,i)
              if (Tdomain%curve) then
                  ! Coordinates in the real chunk
                  xa = u;   ya = v;   za = w
                  Rot = Tdomain%rot
                  u = Rot(0,0)*xa + Rot(0,1)*ya + Rot(0,2)*za
                  v = Rot(1,0)*xa + Rot(1,1)*ya + Rot(1,2)*za
                  w = Rot(2,0)*xa + Rot(2,1)*ya + Rot(2,2)*za
                  ! Convert the cartesian coordinates into spherical coordinates
                  call cart2sph (u,v,w,r,theta,phi)
                  lat_deg = 90.d0 - 180.d0*theta/pi
                  phi_deg = 180.d0*phi/pi
                  if (phi_deg>180.d0)   phi_deg = phi_deg - 360.d0
                  tmp(0) = lat_deg;   tmp(1) = phi_deg;   tmp(2) = r/1000.d0
              else
                  tmp(0) = u/1000.d0;   tmp(1) = v/1000.d0;   tmp(2) = w/1000.d0
              endif

              !!! Computing beta and xi kernels !!!
              if (Tdomain%aniso) then
                  A = Tdomain%specel(n)%save_TIparam(1,x,y,z)
                  L = Tdomain%specel(n)%save_TIparam(4,x,y,z)
                  M = Tdomain%specel(n)%save_TIparam(5,x,y,z)
                  kern_beta = 2.d0 * (Tdomain%specel(n)%Kern_aniso(4,x,y,z) + &
                                      Tdomain%specel(n)%Kern_aniso(5,x,y,z) - &
                                      2.d0*L*Tdomain%specel(n)%Kern_aniso(3,x,y,z)/(A-2.d0*L))
                  kern_xi = (2.d0*L*Tdomain%specel(n)%Kern_aniso(5,x,y,z) - &
                             M*Tdomain%specel(n)%Kern_aniso(4,x,y,z) + &
                             2.d0*L*M*Tdomain%specel(n)%Kern_aniso(3,x,y,z)/(A-2.d0*L)) / (2.d0*L+M)
              else
                  lambda = Tdomain%specel(n)%Lambda(x,y,z)
                  mu = Tdomain%specel(n)%Mu(x,y,z)
                  kern_beta = 2.d0 * (Tdomain%specel(n)%Kern_mu(x,y,z) - &
                                      2.d0*mu*Tdomain%specel(n)%Kern_lambda(x,y,z)/lambda)
              endif

              !!! check de la derivee partielle !!!
!              if (Tdomain%aniso) then
!                  delta = Tdomain%specel(n)%anomaly(x,y,z)/(M/L) - 1.d0
!                  delta_chi = delta_chi + &
!                              delta * kern_xi * &
!                              Tdomain%specel(n)%Jacob(x,y,z) * &
!                              Tdomain%specel(n)%wgtx(x) * Tdomain%specel(n)%wgty(y) * Tdomain%specel(n)%wgtz(z)
!              else
!                  rho = Tdomain%specel(n)%Density(x,y,z) 
!                  delta = Tdomain%specel(n)%anomaly(x,y,z)/dsqrt(mu/rho) - 1.d0
!                  delta_chi = delta_chi + &
!                              delta * kern_beta * &
!                              Tdomain%specel(n)%Jacob(x,y,z) * &
!                              Tdomain%specel(n)%wgtx(x) * Tdomain%specel(n)%wgty(y) * Tdomain%specel(n)%wgtz(z)
!              endif

              !!! Storing the value of the kernel at the interior GLLs !!!
              if (x/=0 .and. y/=0 .and. z/=0 .and. x/=ngll1-1 .and. y/=ngll2-1 .and. z/=ngll3-1) then
                  write (31,*) tmp(0:2), kern_beta
                  if (Tdomain%aniso)   write (32,*) tmp(0:2), kern_xi
              endif

          enddo
         enddo
        enddo

    endif

enddo


close (31)
if (Tdomain%aniso)   close (32)
call MPI_BARRIER(MPI_COMM_WORLD, code)
if (rg==0) then
    call system("cat kernel_beta_* > kern_beta")
    call system("rm kernel_beta_*")
    if (Tdomain%aniso) then
        call system("cat kernel_xi_* > kern_xi")
        call system("rm kernel_xi_*")
    endif
endif


return
end subroutine save_kernel
