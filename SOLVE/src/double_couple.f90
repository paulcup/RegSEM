subroutine double_couple(Tdomain,rg)


use sdomains
use angles

implicit none

type (domain), intent(inout) :: Tdomain
integer, intent(in) :: rg

integer :: n, elem, x,y,z, ngllx,nglly,ngllz, mat, a,b, i,j,k
doubleprecision :: ct,st,cp,sp, xi,eta,zeta, coord, xa,ya,za, xb,yb,zb, num,denom,prod, &
                   x0,y0,z0, x1,y1,z1, x2,y2,z2, x3,y3,z3, x4,y4,z4, x5,y5,z5, x6,y6,z6, x7,y7,z7
doubleprecision, dimension(:), allocatable :: xpol,ypol,zpol, xdh,ydh,zdh
doubleprecision, dimension(0:2,0:2) :: Pcs, Psc, TMP, M, Rot, tRot


do n = 0,Tdomain%n_source-1
    if (Tdomain%sSource(n)%i_type_source == 2 .and. rg == Tdomain%sSource(n)%proc) then

        M = Tdomain%sSource(n)%Moment
        if (Tdomain%curve) then
            Tdomain%sSource(n)%realcolat = pi*Tdomain%sSource(n)%realcolat/180.d0
            Tdomain%sSource(n)%reallong = pi*Tdomain%sSource(n)%reallong/180.d0
            ct=dcos(Tdomain%sSource(n)%realcolat)
            st=dsin(Tdomain%sSource(n)%realcolat)
            cp=dcos(Tdomain%sSource(n)%reallong)
            sp=dsin(Tdomain%sSource(n)%reallong)
            ! matrice de passage du systeme cartesien au systeme spherique
            Pcs(0,0) = st*cp; Pcs(0,1) = ct*cp; Pcs(0,2) = -sp
            Pcs(1,0) = st*sp; Pcs(1,1) = ct*sp; Pcs(1,2) = cp
            Pcs(2,0) = ct   ; Pcs(2,1) = -st  ; Pcs(2,2) = 0.d0
            ! matrice de passage du systeme spherique au systeme cartesien (inverse de Pcs)
            Psc(0,0) = st*cp; Psc(0,1) = st*sp; Psc(0,2) = ct
            Psc(1,0) = ct*cp; Psc(1,1) = ct*sp; Psc(1,2) = -st
            Psc(2,0) = -sp  ; Psc(2,1) = cp   ; Psc(2,2) = 0.0d0
            ! Calcul du tenseur moment sismique dans le systeme cartesien
            do a = 0,2
             do b = 0,2
                 TMP(a,b) = 0.d0
                 do k = 0,2
                     TMP(a,b) = TMP(a,b) + M(a,k)*Psc(k,b)
                 enddo
             enddo
            enddo
            do a = 0,2
             do b = 0,2
                 M(a,b) = 0.d0
                 do k = 0,2
                     M(a,b) = M(a,b) + Pcs(a,k)*TMP(k,b)
                 enddo
             enddo
            enddo
            ! Rotation (du chunk reel au chunk de reference) du tenseur moment sismique
            Rot = Tdomain%rot
            tRot = transpose(Rot)
            do a = 0,2
             do b = 0,2
                 TMP(a,b) = 0.d0
                 do k = 0,2
                     TMP(a,b) = TMP(a,b) + M(a,k)*Rot(k,b)
                 enddo
             enddo
            enddo
            do a = 0,2
             do b = 0,2
                 M(a,b) = 0.d0
                 do k = 0,2
                     M(a,b) = M(a,b) + tRot(a,k)*TMP(k,b)
                 enddo
             enddo
            enddo
        endif

        do a = 0,2
         do b = 0,2
             TMP(a,b) = 0.d0
             do k = 0,2
                 TMP(a,b) = TMP(a,b) + Tdomain%sSource(n)%InvGrad(k,a)*M(b,k)
             enddo
         enddo
        enddo

        elem = Tdomain%Ssource(n)%elem
        ngllx = Tdomain%specel(elem)%ngllx
        nglly = Tdomain%specel(elem)%nglly
        ngllz = Tdomain%specel(elem)%ngllz
        allocate (xpol(0:ngllx-1));   xpol = 1.d0
        allocate (ypol(0:nglly-1));   ypol = 1.d0
        allocate (zpol(0:ngllz-1));   zpol = 1.d0
        allocate (xdh(0:ngllx-1))
        allocate (ydh(0:nglly-1))
        allocate (zdh(0:ngllz-1))
        xi = Tdomain%sSource(n)%refcoord(0)
        eta = Tdomain%sSource(n)%refcoord(1)
        zeta = Tdomain%sSource(n)%refcoord(2)
        mat = Tdomain%specel(elem)%mat_index
        do x = 0,ngllx-1
            coord = Tdomain%sSubdomain(mat)%GLLcx(x)
            num = 0.d0;   denom = 1.d0
            do a = 0,ngllx-1
                if (a/=x) then
                   xa = Tdomain%sSubdomain(mat)%GLLcx(a)
                   xpol(x) = xpol(x) * (xi-xa)/(coord-xa)
                   denom = denom * (coord-xa)
                   prod = 1.d0
                   do b = 0,ngllx-1
                       if ((b/=x) .and. (b/=a)) then
                           xb = Tdomain%sSubdomain(mat)%GLLcx(b)
                           prod = prod * (xi-xb)
                       endif
                   enddo
                   num = num + prod
                endif
            enddo
            xdh(x) = num/denom
        enddo
        do y = 0,nglly-1
            coord = Tdomain%sSubdomain(mat)%GLLcy(y)
            num = 0.d0;   denom = 1.d0
            do a = 0,nglly-1
                if (a/=y) then
                   ya = Tdomain%sSubdomain(mat)%GLLcy(a)
                   ypol(y) = ypol(y) * (eta-ya)/(coord-ya)
                   denom = denom * (coord-ya)
                   prod = 1.d0
                   do b = 0,nglly-1
                       if ((b/=y) .and. (b/=a)) then
                           yb = Tdomain%sSubdomain(mat)%GLLcy(b)
                           prod = prod * (eta-yb)
                       endif
                   enddo
                   num = num + prod
                endif
            enddo
            ydh(y) = num/denom
        enddo
        do z = 0,ngllz-1
            coord = Tdomain%sSubdomain(mat)%GLLcz(z)
            num = 0.d0;   denom = 1.d0
            do a = 0,ngllz-1
                if (a/=z) then
                   za = Tdomain%sSubdomain(mat)%GLLcz(a)
                   zpol(z) = zpol(z) * (zeta-za)/(coord-za)
                   denom = denom * (coord-za)
                   prod = 1.d0
                   do b = 0,ngllz-1
                       if ((b/=z) .and. (b/=a)) then
                           zb = Tdomain%sSubdomain(mat)%GLLcz(b)
                           prod = prod * (zeta-zb)
                       endif
                   enddo
                   num = num + prod
                endif
            enddo
            zdh(z) = num/denom
        enddo

        allocate (Tdomain%sSource(n)%coeff(0:ngllx-1, 0:nglly-1, 0:ngllz-1, 0:2))
        do x = 0,ngllx-1
         do y = 0,nglly-1
          do z = 0,ngllz-1
              Tdomain%sSource(n)%coeff(x,y,z,:) = xdh(x)*ypol(y)*zpol(z)*TMP(0,:) + &
                                                  xpol(x)*ydh(y)*zpol(z)*TMP(1,:) + &
                                                  xpol(x)*ypol(y)*zdh(z)*TMP(2,:)
          enddo
         enddo
        enddo

        deallocate (xpol,ypol,zpol, xdh,ydh,zdh)

    endif
enddo

return
end subroutine double_couple
