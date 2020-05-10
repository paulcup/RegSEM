module ssimus


type :: simu

doubleprecision, dimension (:,:,:,:), pointer :: Forces, Veloc,Displ,Accel, V0

! PML
doubleprecision, dimension (:,:,:,:), pointer :: Diagonal_Stress, Diagonal_Stress1, Diagonal_Stress2, Diagonal_Stress3
doubleprecision, dimension (:,:,:,:), pointer :: Residual_Stress, Residual_Stress1, Residual_Stress2, Residual_Stress3
doubleprecision, dimension (:,:,:,:), pointer :: Forces1,Forces2,Forces3, Veloc1,Veloc2,Veloc3

! Attenuation
doubleprecision, dimension (:,:,:), pointer :: epsilondev_xx_,epsilondev_yy_,epsilondev_xy_,epsilondev_xz_,epsilondev_yz_
doubleprecision, dimension (:,:,:,:), pointer :: R_xx_,R_yy_,R_xy_,R_xz_,R_yz_

! Adjoint
doubleprecision, dimension (:,:,:,:), pointer :: save_strain

end type


type :: simuF

doubleprecision, dimension (:,:,:), pointer :: Forces, Veloc,Displ,Accel, V0

! PML
doubleprecision, dimension (:,:,:), pointer :: Forces1,Forces2,Forces3, Veloc1,Veloc2,Veloc3

end type


type :: simuE

doubleprecision, dimension (:,:), pointer :: Forces, Veloc,Displ,Accel, V0

! PML
doubleprecision, dimension (:,:), pointer :: Forces1,Forces2,Forces3, Veloc1,Veloc2,Veloc3

end type


type :: simuV

doubleprecision, dimension (:), pointer :: Forces, Veloc,Displ,Accel, V0

! PML
doubleprecision, dimension (:), pointer :: Forces1,Forces2,Forces3, Veloc1,Veloc2,Veloc3

end type


end module ssimus
