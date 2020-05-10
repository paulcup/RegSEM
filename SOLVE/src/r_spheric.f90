subroutine r_spheric(ngllx,nglly,ngllz,Tdomain,mat,n,Rsph)


use sdomains

implicit none

type(domain), target, intent (INOUT) :: Tdomain
integer, intent(IN) :: n, ngllx,nglly,ngllz, mat
doubleprecision, dimension(0:ngllx-1,0:nglly-1,0:ngllz-1), intent(OUT) :: Rsph

integer :: i,j,k, i_count
doubleprecision :: xi,eta,zeta, zp, f
doubleprecision, dimension(:), allocatable :: rayon


allocate (rayon(0:Tdomain%n_nodes-1))
do i = 0,Tdomain%n_nodes-1
   j = Tdomain%specel(n)%Control_Nodes(i)
   rayon(i) = Tdomain%Rsph_Nodes(j)
enddo

do k = 0,ngllz - 1
    zeta =  Tdomain%sSubdomain(mat)%GLLcz (k)
    do j = 0,nglly - 1
        eta =  Tdomain%sSubdomain(mat)%GLLcy (j)
        do i = 0,ngllx - 1
            xi = Tdomain%sSubdomain(mat)%GLLcx (i)

            if (Tdomain%n_nodes==27) then
                zp = 0
                do i_count = 0,Tdomain%n_nodes-1
                   f = Comp_shapefunc(i_count,xi,eta,zeta)
                   zp = zp + rayon(i_count)*f
                enddo
            else if (Tdomain%n_nodes==8) then
                zp = 0.125 * (rayon(0)*(1-xi)*(1-eta)*(1-zeta) + rayon(1)*(1+xi)*(1-eta)*(1-zeta) + &
                              rayon(2)*(1+xi)*(1+eta)*(1-zeta) + rayon(3)*(1-xi)*(1+eta)*(1-zeta) + &
                              rayon(4)*(1-xi)*(1-eta)*(1+zeta) + rayon(5)*(1+xi)*(1-eta)*(1+zeta) + &
                              rayon(6)*(1+xi)*(1+eta)*(1+zeta) + rayon(7)*(1-xi)*(1+eta)*(1+zeta))
            else
                print *,"YOU CAN'T HAVE ELLIPTICITY WITH 20 CTRL_PTS FOR NOW"
                stop
            endif
            Rsph(i,j,k) = zp

        enddo
    enddo
enddo

deallocate (rayon)

end subroutine r_spheric

