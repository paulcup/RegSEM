module sadjoints


type :: Adjoint

integer :: elem, proc
integer, dimension(0:2) :: gll
doubleprecision :: depth,realcolat,reallong,refcolat,reflong
doubleprecision, dimension(0:2) :: refcoord
doubleprecision, dimension(:,:), pointer :: timefunc
doubleprecision, dimension(0:2,0:2) :: InvGrad
doubleprecision, dimension (:,:,:), pointer :: coeff

end type 


end module sadjoints
