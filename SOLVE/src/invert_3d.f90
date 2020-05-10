subroutine invert_3d (a,Det)


implicit none

doubleprecision, dimension (0:2,0:2), intent (INOUT) :: a
doubleprecision, intent (OUT) :: Det 

doubleprecision, dimension(0:2,0:2) :: Inverse_A


Det = a(0,0)* (a(1,1)*a(2,2) - a(2,1)*a(1,2) ) + a(1,0) * (a(0,2) * a(2,1) - a(0,1) * a(2,2) ) + a(2,0)  * (a(0,1)*a(1,2) - a(0,2)*a(1,1) ) 

if (Det==0.d0) then
    write (*,*) "Non inversible matrix !!!"
    stop
endif

Inverse_A (0,0) = a(1,1)*a(2,2) - a(2,1)*a(1,2)
Inverse_A (1,0) = a(0,2) * a(2,1) - a(0,1) * a(2,2) 
Inverse_A (2,0) = a(0,1)*a(1,2) - a(0,2)*a(1,1)

Inverse_A (0,1) = a(1,2)*a(2,0) - a(1,0)*a(2,2)
Inverse_A (1,1) = a(0,0) * a(2,2) - a(0,2) * a(2,0) 
Inverse_A (2,1) = a(0,2)*a(1,0) - a(0,0)*a(1,2)

Inverse_A (0,2) = a(1,0)*a(2,1) - a(1,1)*a(2,0)
Inverse_A (1,2) = a(0,1) * a(2,0) - a(2,1) * a(0,0) 
Inverse_A (2,2) = a(1,1)*a(0,0) - a(0,1)*a(1,0)

A = TRANSPOSE (Inverse_A)
A = A/Det 

return
end subroutine invert_3d
