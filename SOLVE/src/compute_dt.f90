subroutine Compute_dt (Tdomain,dt)


use sdomains
use tensor_util

implicit none

type (Domain), intent (IN) :: Tdomain
doubleprecision, intent (INOUT) :: dt

integer :: i,j,k, n, ngllx,nglly,ngllz, idef0,idef1, u,v,w
doubleprecision :: dx,dxmin, vp, lam,mu,rho, ratio,ratio_max
doubleprecision, dimension(1:6,1:6) :: C

doubleprecision, parameter :: courant_max = 0.35


ratio_max = 0.d0

do n = 0, Tdomain%n_elem -1 

    ngllx = Tdomain%specel(n)%ngllx
    nglly = Tdomain%specel(n)%nglly
    ngllz = Tdomain%specel(n)%ngllz

    do k = 0,ngllz-2
       do j = 0,nglly-2
          do i = 0,ngllx-2

              dxmin = 1e10

              idef0 = Tdomain%specel(n)%Iglobnum (i,j,k)
              idef1 = Tdomain%specel(n)%Iglobnum (i+1,j,k)
              dx = (Tdomain%Globcoord(0,idef0) - Tdomain%Globcoord(0,idef1))**2 + &
                   (Tdomain%Globcoord(1,idef0) - Tdomain%Globcoord(1,idef1))**2 + & 
                   (Tdomain%Globcoord(2,idef0) - Tdomain%Globcoord(2,idef1))**2
              if (dx < dxmin ) dxmin = dx
              idef1 = Tdomain%specel(n)%Iglobnum (i,j+1,k)
              dx = (Tdomain%Globcoord(0,idef0) - Tdomain%Globcoord(0,idef1))**2 + &
                   (Tdomain%Globcoord(1,idef0) - Tdomain%Globcoord(1,idef1))**2 + & 
                   (Tdomain%Globcoord(2,idef0) - Tdomain%Globcoord(2,idef1))**2
              if (dx < dxmin ) dxmin = dx
              idef1 = Tdomain%specel(n)%Iglobnum (i,j,k+1)
              dx = (Tdomain%Globcoord(0,idef0) - Tdomain%Globcoord(0,idef1))**2 + &
                   (Tdomain%Globcoord(1,idef0) - Tdomain%Globcoord(1,idef1))**2 + & 
                   (Tdomain%Globcoord(2,idef0) - Tdomain%Globcoord(2,idef1))**2
              if (dx < dxmin ) dxmin = dx

              idef1 = Tdomain%specel(n)%Iglobnum (i+1,j+1,k)
              dx = (Tdomain%Globcoord(0,idef0) - Tdomain%Globcoord(0,idef1))**2 + &
                   (Tdomain%Globcoord(1,idef0) - Tdomain%Globcoord(1,idef1))**2 + & 
                   (Tdomain%Globcoord(2,idef0) - Tdomain%Globcoord(2,idef1))**2
              if (dx < dxmin ) dxmin = dx
              idef1 = Tdomain%specel(n)%Iglobnum (i+1,j,k+1)
              dx = (Tdomain%Globcoord(0,idef0) - Tdomain%Globcoord(0,idef1))**2 + &
                   (Tdomain%Globcoord(1,idef0) - Tdomain%Globcoord(1,idef1))**2 + & 
                   (Tdomain%Globcoord(2,idef0) - Tdomain%Globcoord(2,idef1))**2
              if (dx < dxmin ) dxmin = dx
              idef1 = Tdomain%specel(n)%Iglobnum (i,j+1,k+1)
              dx = (Tdomain%Globcoord(0,idef0) - Tdomain%Globcoord(0,idef1))**2 + &
                   (Tdomain%Globcoord(1,idef0) - Tdomain%Globcoord(1,idef1))**2 + & 
                   (Tdomain%Globcoord(2,idef0) - Tdomain%Globcoord(2,idef1))**2
              if (dx < dxmin ) dxmin = dx

              idef1 = Tdomain%specel(n)%Iglobnum (i+1,j+1,k+1)
              dx = (Tdomain%Globcoord(0,idef0) - Tdomain%Globcoord(0,idef1))**2 + &
                   (Tdomain%Globcoord(1,idef0) - Tdomain%Globcoord(1,idef1))**2 + & 
                   (Tdomain%Globcoord(2,idef0) - Tdomain%Globcoord(2,idef1))**2
              if (dx < dxmin ) dxmin = dx

              if (Tdomain%specel(n)%PML) then
                  lam = Tdomain%specel(n)%Lambda(i,j,k)
                  mu = Tdomain%specel(n)%Mu(i,j,k)
              else
                  if ((Tdomain%aniso) .and. (Tdomain%n_sls==0)) then
                      w = 0
                      do u = 1,6
                          do v = u,6
                              C(u,v) = Tdomain%specel(n)%Cij(w,i,j,k)
                              w = w + 1
                          enddo
                      enddo
                      do u = 2,6
                          do v = 1,u-1
                              C(u,v) = C(v,u)
                          enddo
                      enddo
                      lam = lambda_from_Cij(C)
                      mu = mu_from_Cij(C)
                  else
                      lam = Tdomain%specel(n)%Lambda(i,j,k)
                      mu = Tdomain%specel(n)%Mu(i,j,k)
                  endif
              endif
              rho = Tdomain%specel(n)%Density(i,j,k)
              vp = (lam + 2.d0*mu) / rho

              ratio = dsqrt(vp/dxmin)
              if (ratio>ratio_max)   ratio_max = ratio
          enddo
       enddo
    enddo

enddo

dt = courant_max/ratio_max


return
end subroutine Compute_dt
