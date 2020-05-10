subroutine forceAdjoint(Tdomain,rg)


use sdomains
use angles

implicit none

type (domain), intent(inout) :: Tdomain
integer, intent(in) :: rg

integer :: n, elem, x,y,z, ngllx,nglly,ngllz, mat, a, i,j,k
doubleprecision :: ct,st,cp,sp, xi,eta,zeta, coord, xa,ya,za
doubleprecision, dimension(:), allocatable :: xpol,ypol,zpol
doubleprecision, dimension(0:2) :: TMP
doubleprecision, dimension(0:2,0:2) :: Pcs, Rot, tRot


do n = 0,Tdomain%n_source_adj-1
    if (rg == Tdomain%sAdjoint(n)%proc) then

        if (Tdomain%curve) then
            Tdomain%sAdjoint(n)%realcolat = pi*Tdomain%sAdjoint(n)%realcolat/180.d0
            Tdomain%sAdjoint(n)%reallong = pi*Tdomain%sAdjoint(n)%reallong/180.d0
            ct=dcos(Tdomain%sAdjoint(n)%realcolat)
            st=dsin(Tdomain%sAdjoint(n)%realcolat)
            cp=dcos(Tdomain%sAdjoint(n)%reallong)
            sp=dsin(Tdomain%sAdjoint(n)%reallong)
            ! Passage de spherique a cartesien
            Pcs(0,0) = st*cp; Pcs(0,1) = -ct*cp; Pcs(0,2) = -sp
            Pcs(1,0) = st*sp; Pcs(1,1) = -ct*sp; Pcs(1,2) = cp
            Pcs(2,0) = ct   ; Pcs(2,1) = st  ; Pcs(2,2) = 0.d0
            do i = 0,Tdomain%sTimeParam%ntime-1
                do a = 0,2
                    TMP(a) = 0.d0
                    do k = 0,2
                        TMP(a) = TMP(a) + Pcs(a,k)*Tdomain%sAdjoint(n)%timefunc(k,i)
                    enddo
                enddo
                ! Rotation (du chunk reel au chunk de reference) du vecteur force
                Rot = Tdomain%rot
                tRot = transpose(Rot)
                Tdomain%sAdjoint(n)%timefunc(0,i) = tRot(0,0)*TMP(0) + tRot(0,1)*TMP(1) + tRot(0,2)*TMP(2)
                Tdomain%sAdjoint(n)%timefunc(1,i) = tRot(1,0)*TMP(0) + tRot(1,1)*TMP(1) + tRot(1,2)*TMP(2)
                Tdomain%sAdjoint(n)%timefunc(2,i) = tRot(2,0)*TMP(0) + tRot(2,1)*TMP(1) + tRot(2,2)*TMP(2)
            enddo
        endif

        elem = Tdomain%sAdjoint(n)%elem
        ngllx = Tdomain%specel(elem)%ngllx
        nglly = Tdomain%specel(elem)%nglly
        ngllz = Tdomain%specel(elem)%ngllz
        allocate (xpol(0:ngllx-1));   xpol = 1
        allocate (ypol(0:nglly-1));   ypol = 1
        allocate (zpol(0:ngllz-1));   zpol = 1
        xi = Tdomain%sAdjoint(n)%refcoord(0)
        eta = Tdomain%sAdjoint(n)%refcoord(1)
        zeta = Tdomain%sAdjoint(n)%refcoord(2)
        mat = Tdomain%specel(elem)%mat_index
        do x = 0,ngllx-1
            coord = Tdomain%sSubdomain(mat)%GLLcx(x)
            do a = 0,ngllx-1
                if (a/=x) then
                   xa = Tdomain%sSubdomain(mat)%GLLcx(a)
                   xpol(x) = xpol(x) * (xi-xa)/(coord-xa)
                endif
            enddo
        enddo
        do y = 0,nglly-1
            coord = Tdomain%sSubdomain(mat)%GLLcy(y)
            do a = 0,nglly-1
                if (a/=y) then
                   ya = Tdomain%sSubdomain(mat)%GLLcy(a)
                   ypol(y) = ypol(y) * (eta-ya)/(coord-ya)
                endif
            enddo
        enddo
        do z = 0,ngllz-1
            coord = Tdomain%sSubdomain(mat)%GLLcz(z)
            do a = 0,ngllz-1
                if (a/=z) then
                   za = Tdomain%sSubdomain(mat)%GLLcz(a)
                   zpol(z) = zpol(z) * (zeta-za)/(coord-za)
                endif
            enddo
        enddo

        allocate (Tdomain%sAdjoint(n)%coeff(0:ngllx-1, 0:nglly-1, 0:ngllz-1))
        do x = 0,ngllx-1
         do y = 0,nglly-1
          do z = 0,ngllz-1
              Tdomain%sAdjoint(n)%coeff(x,y,z) = xpol(x)*ypol(y)*zpol(z)
          enddo
         enddo
        enddo

        deallocate (xpol,ypol,zpol)

    endif
enddo

return
end subroutine forceAdjoint
