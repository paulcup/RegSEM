module sedges


use ssimus

type :: edge

type(simuE), dimension (:), pointer :: sSimu

logical :: PML, Abs, ocean, mirror

integer :: ngll
integer, dimension (:), pointer :: Iglobnum_Edge

doubleprecision, dimension (:), pointer  :: MassMat, Mocean
doubleprecision, dimension (:,:), pointer :: DumpMass, DumpVx, DumpVy, DumpVz
doubleprecision, dimension (:,:,:), pointer :: verticale, M33ocean, random_coeff

end type

contains

! ############################################################
subroutine Prediction_Edge_Veloc (E, alpha, bega, gam1, dt, i_simu)
implicit none

type (Edge), intent (INOUT) :: E
integer, intent (IN) :: i_simu
doubleprecision, intent (IN) :: alpha, bega, gam1, dt

integer :: ngll

ngll = E%ngll

E%sSimu(i_simu)%Forces(:,:) = E%sSimu(i_simu)%Displ(:,:) + &
                              dt * E%sSimu(i_simu)%Veloc(:,:) + &
                              dt**2 * (0.5-bega) * E%sSimu(i_simu)%Accel(:,:)
E%sSimu(i_simu)%V0(:,:) = E%sSimu(i_simu)%Veloc(:,:)
E%sSimu(i_simu)%Forces(:,:) = alpha * E%sSimu(i_simu)%Forces(:,:) + &
                              (1-alpha) * E%sSimu(i_simu)%Displ(:,:)

return
end subroutine Prediction_Edge_Veloc

! ###########################################################
subroutine Correction_Edge_Veloc (E, i_simu, bega, gam1, dt)
implicit none

type (Edge), intent (INOUT) :: E
doubleprecision, intent (IN) :: bega, gam1, dt
integer, intent (IN) :: i_simu

integer :: i,k, ngll
doubleprecision, dimension(0:2) :: tmp

ngll = E%ngll

if (E%ocean) then
    do i = 1,ngll-2
        tmp(:) = E%sSimu(i_simu)%Forces(i,:)
        do k = 0,2
            E%sSimu(i_simu)%Forces(i,k) = E%M33ocean(i,k,0) * tmp(0) + &
                                          E%M33ocean(i,k,1) * tmp(1) + &
                                          E%M33ocean(i,k,2) * tmp(2)
        enddo
    enddo
else
    do k = 0,2
        E%sSimu(i_simu)%Forces(:,k) = E%MassMat(:) * E%sSimu(i_simu)%Forces(:,k)
    enddo
endif
E%sSimu(i_simu)%Veloc(:,:) = E%sSimu(i_simu)%V0(:,:) + &
                             dt * E%sSimu(i_simu)%Forces(:,:)
E%sSimu(i_simu)%Accel(:,:) = E%sSimu(i_simu)%Accel(:,:) + &
                             gam1 / dt * (E%sSimu(i_simu)%Veloc(:,:)-E%sSimu(i_simu)%V0(:,:))
E%sSimu(i_simu)%Displ(:,:) = E%sSimu(i_simu)%Displ(:,:) + &
                             bega * dt * (E%sSimu(i_simu)%Veloc(:,:)+E%sSimu(i_simu)%V0(:,:))

return
end subroutine Correction_Edge_Veloc

! ############################################################
subroutine Correction_Edge_PML_Veloc (E, i_simu, dt)

implicit none

type (Edge), intent (INOUT) :: E
doubleprecision, intent (IN) :: dt
integer, intent (IN) :: i_simu

integer :: i

do i = 0,2
    E%sSimu(i_simu)%Veloc1(:,i) = E%DumpVx(:,0) * E%sSimu(i_simu)%Veloc1(:,i) + &
                                  dt * E%DumpVx(:,1) * E%sSimu(i_simu)%Forces1(:,i)
    E%sSimu(i_simu)%Veloc2(:,i) = E%DumpVy(:,0) * E%sSimu(i_simu)%Veloc2(:,i) + &
                                  dt * E%DumpVy(:,1) * E%sSimu(i_simu)%Forces2(:,i)
    E%sSimu(i_simu)%Veloc3(:,i) = E%DumpVz(:,0) * E%sSimu(i_simu)%Veloc3(:,i) + &
                                  dt * E%DumpVz(:,1) * E%sSimu(i_simu)%Forces3(:,i)
enddo

E%sSimu(i_simu)%Veloc = E%sSimu(i_simu)%Veloc1 + E%sSimu(i_simu)%Veloc2 + E%sSimu(i_simu)%Veloc3

if (E%Abs) then
    E%sSimu(i_simu)%Veloc = 0.d0
endif

return
end subroutine Correction_Edge_PML_Veloc

! ############################################
end module sedges
