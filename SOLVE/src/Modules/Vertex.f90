module svertices


use ssimus

type :: vertex

type(simuV), dimension (:), pointer :: sSimu

logical :: PML, Abs, ocean, mirror

integer :: Iglobnum_Vertex

doubleprecision :: MassMat, Mocean
doubleprecision, dimension(:), pointer :: DumpVx, DumpVy, DumpVz, DumpMass
doubleprecision, dimension(:,:), pointer :: verticale, M33ocean, random_coeff

end type

contains

! ############################################################
subroutine Prediction_Vertex_Veloc (V, alpha, bega, gam1, dt, i_simu)

implicit none

type (Vertex), intent (INOUT) :: V
integer, intent (IN) :: i_simu
doubleprecision, intent (IN) :: alpha, bega, gam1, dt


V%sSimu(i_simu)%Forces = V%sSimu(i_simu)%Displ + &
                         dt * V%sSimu(i_simu)%Veloc + &
                         dt**2 * (0.5-bega) * V%sSimu(i_simu)%Accel
V%sSimu(i_simu)%V0 = V%sSimu(i_simu)%Veloc
V%sSimu(i_simu)%Forces = alpha * V%sSimu(i_simu)%Forces + &
                         (1-alpha) * V%sSimu(i_simu)%Displ

return
end subroutine Prediction_Vertex_Veloc

! ###########################################################
subroutine Correction_Vertex_Veloc (V, i_simu, bega, gam1, dt)

implicit none

type (Vertex), intent (INOUT) :: V
doubleprecision, intent (IN) :: bega, gam1, dt
integer, intent (IN) :: i_simu

integer :: k
doubleprecision, dimension(0:2) :: tmp

if (V%ocean) then
    tmp(:) = V%sSimu(i_simu)%Forces(:)
    do k = 0,2
        V%sSimu(i_simu)%Forces(k) = V%M33ocean(k,0) * tmp(0) + &
                                    V%M33ocean(k,1) * tmp(1) + &
                                    V%M33ocean(k,2) * tmp(2)
    enddo
else
    do k = 0,2
        V%sSimu(i_simu)%Forces(k) = V%MassMat * V%sSimu(i_simu)%Forces(k)
    enddo
endif
V%sSimu(i_simu)%Veloc = V%sSimu(i_simu)%V0 + &
                        dt * V%sSimu(i_simu)%Forces
V%sSimu(i_simu)%Accel = V%sSimu(i_simu)%Accel + &
                        gam1 / dt * (V%sSimu(i_simu)%Veloc-V%sSimu(i_simu)%V0)
V%sSimu(i_simu)%Displ = V%sSimu(i_simu)%Displ + &
                        bega * dt * (V%sSimu(i_simu)%Veloc+V%sSimu(i_simu)%V0)

return
end subroutine Correction_Vertex_Veloc

! ############################################################
subroutine Correction_Vertex_PML_Veloc (V, i_simu, dt)

implicit none

type (Vertex), intent (INOUT) :: V
doubleprecision, intent (IN) :: dt
integer, intent (IN) :: i_simu

integer :: i


do i = 0,2
    V%sSimu(i_simu)%Veloc1(i) = V%DumpVx(0) * V%sSimu(i_simu)%Veloc1(i) + &
                                dt * V%DumpVx(1) * V%sSimu(i_simu)%Forces1(i)
    V%sSimu(i_simu)%Veloc2(i) = V%DumpVy(0) * V%sSimu(i_simu)%Veloc2(i) + &
                                dt * V%DumpVy(1) * V%sSimu(i_simu)%Forces2(i)
    V%sSimu(i_simu)%Veloc3(i) = V%DumpVz(0) * V%sSimu(i_simu)%Veloc3(i) + &
                                dt * V%DumpVz(1) * V%sSimu(i_simu)%Forces3(i)
enddo

V%sSimu(i_simu)%Veloc = V%sSimu(i_simu)%Veloc1 + V%sSimu(i_simu)%Veloc2 + V%sSimu(i_simu)%Veloc3

if (V%Abs) then
    V%sSimu(i_simu)%Veloc = 0.d0
endif

return
end subroutine Correction_Vertex_PML_Veloc

! ##########################################################
end module svertices
