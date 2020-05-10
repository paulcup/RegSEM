subroutine kernel_iso (Tdomain,dt,ntime)


use sdomains

implicit none

type (Domain), intent (INOUT) :: Tdomain
doubleprecision, intent(IN) :: dt
integer, intent (IN) :: ntime

integer :: n, x,y,z, ngll1,ngll2,ngll3, i
doubleprecision :: mu,lambda
doubleprecision, dimension(0:1) :: div
doubleprecision, dimension(:,:,:), allocatable :: kern_lambda_loc, kern_mu_loc


do n = 0,Tdomain%n_elem-1
    if (.not. Tdomain%specel(n)%PML) then
        ngll1 = Tdomain%specel(n)%ngllx
        ngll2 = Tdomain%specel(n)%nglly
        ngll3 = Tdomain%specel(n)%ngllz
        allocate (kern_lambda_loc(0:ngll1-1,0:ngll2-1,0:ngll3-1))
        allocate (kern_mu_loc(0:ngll1-1,0:ngll2-1,0:ngll3-1))
        do z = 0,ngll3-1 
         do y = 0,ngll2-1 
          do x = 0,ngll1-1
!!! CAUTION: kern_lambda_loc is the value of the kernel at t-dt/2 !!!
!!! This is because save_strain comes from internal_forces.f90 where the displacement from the prediction phase is used !!!
              lambda = Tdomain%specel(n)%Lambda(x,y,z)
              do i = 0,1
                  div(i) = Tdomain%specel(n)%sSimu(i)%save_strain(x,y,z,0) + &
                           Tdomain%specel(n)%sSimu(i)%save_strain(x,y,z,1) + &
                           Tdomain%specel(n)%sSimu(i)%save_strain(x,y,z,2)
              enddo
              kern_lambda_loc(x,y,z) = -lambda * div(1) * div(0)
              Tdomain%specel(n)%Kern_lambda(x,y,z) = Tdomain%specel(n)%Kern_lambda(x,y,z) + dt*kern_lambda_loc(x,y,z)
!!! CAUTION: kern_mu_loc is the value of the kernel at t-dt/2 !!!
!!! This is because save_strain come from internal_forces.f90 where the displacement from the prediction phase is used !!!
              mu = Tdomain%specel(n)%Mu(x,y,z)
              kern_mu_loc(x,y,z) = -2.d0*mu * &
                                   ( Tdomain%specel(n)%sSimu(1)%save_strain(x,y,z,0)*Tdomain%specel(n)%sSimu(0)%save_strain(x,y,z,0) + &
                                     Tdomain%specel(n)%sSimu(1)%save_strain(x,y,z,1)*Tdomain%specel(n)%sSimu(0)%save_strain(x,y,z,1) + &
                                     Tdomain%specel(n)%sSimu(1)%save_strain(x,y,z,2)*Tdomain%specel(n)%sSimu(0)%save_strain(x,y,z,2) + &
                                     2*Tdomain%specel(n)%sSimu(1)%save_strain(x,y,z,3)*Tdomain%specel(n)%sSimu(0)%save_strain(x,y,z,3) + &
                                     2*Tdomain%specel(n)%sSimu(1)%save_strain(x,y,z,4)*Tdomain%specel(n)%sSimu(0)%save_strain(x,y,z,4) + &
                                     2*Tdomain%specel(n)%sSimu(1)%save_strain(x,y,z,5)*Tdomain%specel(n)%sSimu(0)%save_strain(x,y,z,5) )
              Tdomain%specel(n)%Kern_mu(x,y,z) = Tdomain%specel(n)%Kern_mu(x,y,z) + dt*kern_mu_loc(x,y,z)
          enddo
         enddo
        enddo
        deallocate (kern_lambda_loc,kern_mu_loc)
    endif
enddo


return
end subroutine kernel_iso
