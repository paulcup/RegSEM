subroutine kernel_aniso (Tdomain,dt,ntime)


use sdomains
use angles

implicit none

type (Domain), intent (INOUT) :: Tdomain
doubleprecision, intent(IN) :: dt
integer, intent (IN) :: ntime

integer :: n, x,y,z, ngll1,ngll2,ngll3, i,k, idef, a,b
doubleprecision :: r,theta_ref,phi_ref, u,v,w, ct,st,cp,sp
doubleprecision, dimension(0:2,0:2) :: Pcs, Psc, TMP, M
doubleprecision, dimension(:,:,:), allocatable :: kern_A_loc, kern_C_loc, kern_F_loc, kern_L_loc, kern_N_loc


do n = 0,Tdomain%n_elem-1
    if (.not. Tdomain%specel(n)%PML) then
        ngll1 = Tdomain%specel(n)%ngllx
        ngll2 = Tdomain%specel(n)%nglly
        ngll3 = Tdomain%specel(n)%ngllz
        allocate (kern_A_loc(0:ngll1-1,0:ngll2-1,0:ngll3-1))
        allocate (kern_C_loc(0:ngll1-1,0:ngll2-1,0:ngll3-1))
        allocate (kern_F_loc(0:ngll1-1,0:ngll2-1,0:ngll3-1))
        allocate (kern_L_loc(0:ngll1-1,0:ngll2-1,0:ngll3-1))
        allocate (kern_N_loc(0:ngll1-1,0:ngll2-1,0:ngll3-1))
        do z = 0,ngll3-1 
         do y = 0,ngll2-1 
          do x = 0,ngll1-1

              !!! Expressing the adjoint and regular strain tensors in spherical coordinates !!!
              idef = Tdomain%specel(n)%Iglobnum(x,y,z)
              u = Tdomain%Globcoord(0,idef)
              v = Tdomain%Globcoord(1,idef)
              w = Tdomain%Globcoord(2,idef)
              call cart2sph(u,v,w,r,theta_ref,phi_ref)
              ct=dcos(theta_ref)
              st=dsin(theta_ref)
              cp=dcos(phi_ref)
              sp=dsin(phi_ref)
              ! matrice de passage du systeme cartesien au systeme spherique
              Pcs(0,0) = st*cp; Pcs(0,1) = ct*cp; Pcs(0,2) = -sp
              Pcs(1,0) = st*sp; Pcs(1,1) = ct*sp; Pcs(1,2) = cp
              Pcs(2,0) = ct   ; Pcs(2,1) = -st  ; Pcs(2,2) = 0.0d0
              ! matrice de passage du systeme spherique au systeme cartesien (inverse de Pcs)
              Psc(0,0) = st*cp; Psc(0,1) = st*sp; Psc(0,2) = ct
              Psc(1,0) = ct*cp; Psc(1,1) = ct*sp; Psc(1,2) = -st
              Psc(2,0) = -sp  ; Psc(2,1) = cp   ; Psc(2,2) = 0.0d0
              do i = 0,1
                  M(0,0) = Tdomain%specel(n)%sSimu(i)%save_strain(x,y,z,0)
                  M(1,1) = Tdomain%specel(n)%sSimu(i)%save_strain(x,y,z,1)
                  M(2,2) = Tdomain%specel(n)%sSimu(i)%save_strain(x,y,z,2)
                  M(0,1) = Tdomain%specel(n)%sSimu(i)%save_strain(x,y,z,3);   M(1,0) = M(0,1)
                  M(0,2) = Tdomain%specel(n)%sSimu(i)%save_strain(x,y,z,4);   M(2,0) = M(0,2)
                  M(1,2) = Tdomain%specel(n)%sSimu(i)%save_strain(x,y,z,5);   M(2,1) = M(1,2)
                  do a = 0,2
                   do b = 0,2
                       TMP(a,b) = 0.0d0
                       do k = 0,2 
                           TMP(a,b) = TMP(a,b) + M(a,k)*Pcs(k,b)
                       enddo
                   enddo
                  enddo
                  do a = 0,2
                   do b = 0,2
                       M(a,b) = 0.0d0
                       do k = 0,2
                           M(a,b) = M(a,b) + Psc(a,k)*TMP(k,b)
                       enddo
                   enddo
                  enddo
                  Tdomain%specel(n)%sSimu(i)%save_strain(x,y,z,0) = M(0,0)
                  Tdomain%specel(n)%sSimu(i)%save_strain(x,y,z,1) = M(1,1)
                  Tdomain%specel(n)%sSimu(i)%save_strain(x,y,z,2) = M(2,2)
                  Tdomain%specel(n)%sSimu(i)%save_strain(x,y,z,3) = M(0,1)
                  Tdomain%specel(n)%sSimu(i)%save_strain(x,y,z,4) = M(0,2)
                  Tdomain%specel(n)%sSimu(i)%save_strain(x,y,z,5) = M(1,2)
              enddo

!!! CAUTION: kern_A_loc is the value of the kernel at t-dt/2 !!!
!!! This is because save_strain comes from internal_forces.f90 where the displacement from the prediction phase is used !!!
              kern_A_loc(x,y,z) = -Tdomain%specel(n)%save_TIparam(1,x,y,z) * &
                                  ( Tdomain%specel(n)%sSimu(1)%save_strain(x,y,z,1)*Tdomain%specel(n)%sSimu(0)%save_strain(x,y,z,1) + &
                                    Tdomain%specel(n)%sSimu(1)%save_strain(x,y,z,2)*Tdomain%specel(n)%sSimu(0)%save_strain(x,y,z,2) + &
                                    Tdomain%specel(n)%sSimu(1)%save_strain(x,y,z,1)*Tdomain%specel(n)%sSimu(0)%save_strain(x,y,z,2) + &
                                    Tdomain%specel(n)%sSimu(1)%save_strain(x,y,z,2)*Tdomain%specel(n)%sSimu(0)%save_strain(x,y,z,1) )
              Tdomain%specel(n)%Kern_aniso(1,x,y,z) = Tdomain%specel(n)%Kern_aniso(1,x,y,z) + dt*kern_A_loc(x,y,z)
!!! CAUTION: kern_C_loc is the value of the kernel at t-dt/2 !!!
!!! This is because save_strain come from internal_forces.f90 where the displacement from the prediction phase is used !!!
              kern_C_loc(x,y,z) = -Tdomain%specel(n)%save_TIparam(2,x,y,z) * &
                                  Tdomain%specel(n)%sSimu(1)%save_strain(x,y,z,0)*Tdomain%specel(n)%sSimu(0)%save_strain(x,y,z,0)
              Tdomain%specel(n)%Kern_aniso(2,x,y,z) = Tdomain%specel(n)%Kern_aniso(2,x,y,z) + dt*kern_C_loc(x,y,z)
!!! CAUTION: kern_F_loc is the value of the kernel at t-dt/2 !!!
!!! This is because save_strain comes from internal_forces.f90 where the displacement from the prediction phase is used !!!
              kern_F_loc(x,y,z) = -Tdomain%specel(n)%save_TIparam(3,x,y,z) * &
                                  ( Tdomain%specel(n)%sSimu(1)%save_strain(x,y,z,0)*Tdomain%specel(n)%sSimu(0)%save_strain(x,y,z,1) + &
                                    Tdomain%specel(n)%sSimu(1)%save_strain(x,y,z,0)*Tdomain%specel(n)%sSimu(0)%save_strain(x,y,z,2) + &
                                    Tdomain%specel(n)%sSimu(1)%save_strain(x,y,z,1)*Tdomain%specel(n)%sSimu(0)%save_strain(x,y,z,0) + &
                                    Tdomain%specel(n)%sSimu(1)%save_strain(x,y,z,2)*Tdomain%specel(n)%sSimu(0)%save_strain(x,y,z,0) )
              Tdomain%specel(n)%Kern_aniso(3,x,y,z) = Tdomain%specel(n)%Kern_aniso(3,x,y,z) + dt*kern_F_loc(x,y,z)
!!! CAUTION: kern_L_loc is the value of the kernel at t-dt/2 !!!
!!! This is because save_strain comes from internal_forces.f90 where the displacement from the prediction phase is used !!!
              kern_L_loc(x,y,z) = -4.d0*Tdomain%specel(n)%save_TIparam(4,x,y,z) * &
                                  ( Tdomain%specel(n)%sSimu(1)%save_strain(x,y,z,4)*Tdomain%specel(n)%sSimu(0)%save_strain(x,y,z,4) + &
                                    Tdomain%specel(n)%sSimu(1)%save_strain(x,y,z,3)*Tdomain%specel(n)%sSimu(0)%save_strain(x,y,z,3) )
              Tdomain%specel(n)%Kern_aniso(4,x,y,z) = Tdomain%specel(n)%Kern_aniso(4,x,y,z) + dt*kern_L_loc(x,y,z)
!!! CAUTION: kern_N_loc is the value of the kernel at t-dt/2 !!!
!!! This is because save_strain comes from internal_forces.f90 where the displacement from the prediction phase is used !!!
              kern_N_loc(x,y,z) = -2.d0*Tdomain%specel(n)%save_TIparam(5,x,y,z) * &
                                  ( 2.d0*Tdomain%specel(n)%sSimu(1)%save_strain(x,y,z,5)*Tdomain%specel(n)%sSimu(0)%save_strain(x,y,z,5) - &
                                    Tdomain%specel(n)%sSimu(1)%save_strain(x,y,z,1)*Tdomain%specel(n)%sSimu(0)%save_strain(x,y,z,2) - &
                                    Tdomain%specel(n)%sSimu(1)%save_strain(x,y,z,2)*Tdomain%specel(n)%sSimu(0)%save_strain(x,y,z,1) )
              Tdomain%specel(n)%Kern_aniso(5,x,y,z) = Tdomain%specel(n)%Kern_aniso(5,x,y,z) + dt*kern_N_loc(x,y,z)
          enddo
         enddo
        enddo
        deallocate (kern_A_loc,kern_C_loc,kern_F_loc,kern_L_loc,kern_N_loc)
    endif
enddo


return
end subroutine kernel_aniso
