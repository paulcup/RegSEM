module sfaces


use ssimus

type :: face

type(simuF), dimension (:), pointer :: sSimu

logical :: PML, Abs, ocean, mirror

integer :: ngll1, ngll2
integer, dimension (:,:), pointer :: Iglobnum_Face

doubleprecision, dimension (:,:), pointer  :: MassMat, Mocean
doubleprecision, dimension (:,:,:), pointer :: DumpVx, DumpVy, DumpVz, DumpMass
doubleprecision, dimension (:,:,:,:), pointer :: verticale, M33ocean

end type
 
contains

! ############################################################
subroutine Prediction_Face_Veloc (F, alpha, bega, gam1, dt, i_simu)

implicit none

type (Face), intent (INOUT) :: F
integer, intent (IN) :: i_simu
doubleprecision, intent (IN) :: alpha, bega, gam1, dt
  
integer :: ngll1, ngll2


ngll1 = F%ngll1 ; ngll2 =F%ngll2
 
F%sSimu(i_simu)%Forces(:,:,:) = F%sSimu(i_simu)%Displ(:,:,:) + &
                                dt * F%sSimu(i_simu)%Veloc(:,:,:) + &
                                dt**2 * (0.5 - bega) * F%sSimu(i_simu)%Accel(:,:,:)
F%sSimu(i_simu)%V0(:,:,:) = F%sSimu(i_simu)%Veloc(:,:,:)
F%sSimu(i_simu)%Forces(:,:,:) = alpha * F%sSimu(i_simu)%Forces(:,:,:) + &
                                (1-alpha) * F%sSimu(i_simu)%Displ(:,:,:) 

return
end subroutine Prediction_Face_Veloc

! ###########################################################
subroutine Correction_Face_Veloc (F, i_simu, bega, gam1, dt)

implicit none

type (Face), intent (INOUT) :: F
doubleprecision, intent (IN) :: bega, gam1, dt
integer, intent (IN) :: i_simu

integer :: i,j,k, ngll1,ngll2
doubleprecision, dimension(0:2) :: tmp


ngll1 = F%ngll1 ; ngll2 = F%ngll2

if (F%ocean) then
    do j = 1,ngll2-2
     do i = 1,ngll1-2
         tmp(:) = F%sSimu(i_simu)%Forces(i,j,:)
         do k = 0,2
             F%sSimu(i_simu)%Forces(i,j,k) = F%M33ocean(i,j,k,0)*tmp(0) + &
                                             F%M33ocean(i,j,k,1)*tmp(1) + &
                                             F%M33ocean(i,j,k,2)*tmp(2)
         enddo
     enddo
    enddo
else
    do k = 0,2
        F%sSimu(i_simu)%Forces(:,:,k) = F%MassMat(:,:) * F%sSimu(i_simu)%Forces(:,:,k)
    enddo
endif
F%sSimu(i_simu)%Veloc(:,:,:) = F%sSimu(i_simu)%V0(:,:,:) + &
                               dt * F%sSimu(i_simu)%Forces(:,:,:)
F%sSimu(i_simu)%Accel(:,:,:) = F%sSimu(i_simu)%Accel(:,:,:) + &
                               gam1 / dt * (F%sSimu(i_simu)%Veloc(:,:,:)-F%sSimu(i_simu)%V0(:,:,:))
F%sSimu(i_simu)%Displ(:,:,:) = F%sSimu(i_simu)%Displ(:,:,:) + &
                               bega * dt * (F%sSimu(i_simu)%Veloc(:,:,:)+F%sSimu(i_simu)%V0(:,:,:))

return
end subroutine Correction_Face_Veloc

! ##########################################################
subroutine Correction_Face_PML_Veloc (F, i_simu, dt)

implicit none

type (Face), intent (INOUT) :: F
doubleprecision, intent (IN) :: dt
integer, intent (IN) :: i_simu

integer :: i


do i = 0,2
    F%sSimu(i_simu)%Veloc1(:,:,i) = F%DumpVx(:,:,0) * F%sSimu(i_simu)%Veloc1(:,:,i) + &
                                    dt * F%DumpVx(:,:,1) * F%sSimu(i_simu)%Forces1(:,:,i)
    F%sSimu(i_simu)%Veloc2(:,:,i) = F%DumpVy(:,:,0) * F%sSimu(i_simu)%Veloc2(:,:,i) + &
                                    dt * F%DumpVy(:,:,1) * F%sSimu(i_simu)%Forces2(:,:,i)
    F%sSimu(i_simu)%Veloc3(:,:,i) = F%DumpVz(:,:,0) * F%sSimu(i_simu)%Veloc3(:,:,i) + &
                                    dt * F%DumpVz(:,:,1) * F%sSimu(i_simu)%Forces3(:,:,i)
enddo

F%sSimu(i_simu)%Veloc = F%sSimu(i_simu)%Veloc1 + F%sSimu(i_simu)%Veloc2 + F%sSimu(i_simu)%Veloc3

if (F%Abs) then
    F%sSimu(i_simu)%Veloc = 0.d0
endif

return
end subroutine Correction_Face_PML_Veloc

! ############################################
end module sfaces
