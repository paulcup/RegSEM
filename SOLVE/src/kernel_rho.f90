subroutine kernel_rho (Tdomain,dt,ntime)


use sdomains

implicit none

type (Domain), intent (INOUT) :: Tdomain
doubleprecision, intent(IN) :: dt
integer, intent (IN) :: ntime

integer :: n, x,y,z, ngll1,ngll2,ngll3
doubleprecision :: rho
doubleprecision, dimension(:,:,:), allocatable :: kern_rho_loc


do n = 0,Tdomain%n_elem-1
    if (.not. Tdomain%specel(n)%PML) then
        ngll1 = Tdomain%specel(n)%ngllx
        ngll2 = Tdomain%specel(n)%nglly
        ngll3 = Tdomain%specel(n)%ngllz
        allocate (kern_rho_loc(0:ngll1-1,0:ngll2-1,0:ngll3-1))
        do z = 1,ngll3-2
         do y = 1,ngll2-2
          do x = 1,ngll1-2
              rho = Tdomain%specel(n)%Density(x,y,z)
              kern_rho_loc(x,y,z) = -rho * &
                                    ( Tdomain%specel(n)%sSimu(1)%Displ(x,y,z,0)*Tdomain%specel(n)%sSimu(0)%Accel(x,y,z,0) &
                                    + Tdomain%specel(n)%sSimu(1)%Displ(x,y,z,1)*Tdomain%specel(n)%sSimu(0)%Accel(x,y,z,1) &
                                    + Tdomain%specel(n)%sSimu(1)%Displ(x,y,z,2)*Tdomain%specel(n)%sSimu(0)%Accel(x,y,z,2) )
              Tdomain%specel(n)%Kern_rho(x,y,z) = Tdomain%specel(n)%Kern_rho(x,y,z) + dt*kern_rho_loc(x,y,z)
          enddo
         enddo
        enddo
        !!! Reste a chopper les valeurs des differents champs sur les faces, les edges et les vertices !!!
        deallocate (kern_rho_loc)
    endif
enddo


return
end subroutine kernel_rho
