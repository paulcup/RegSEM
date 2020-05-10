module sreceivers


type :: receiver

integer :: elem, proc
doubleprecision :: refcolat,reflong,elevation,realcolat,reallong, cosgamma,singamma
doubleprecision, dimension(0:2) :: gll
doubleprecision, dimension(0:2,0:2) :: Passage, Pass
doubleprecision, dimension(:,:), pointer :: StoreTrace
doubleprecision, dimension(:,:,:), pointer :: pol
doubleprecision, dimension(:,:,:,:), pointer :: coeff
character (len=5) :: sta_name

end type
 
end module sreceivers
