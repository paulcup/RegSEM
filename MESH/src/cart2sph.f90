subroutine cart2sph (x,y,z,r,theta,phi)


implicit none

doubleprecision, intent(IN) :: x,y,z
doubleprecision, intent(OUT) :: r,theta,phi
doubleprecision :: dx
doubleprecision, parameter :: pi = 3.141592653589793238462643383279502884197


r = dsqrt(x**2 + y**2 + z**2)
if (r==0) then
    theta = 0.d0;   phi = 0.d0
else
    dx = z/r
    if (dx >= 1.d0) then
        theta = 0.d0
    else if (dx <= -1.d0) then
        theta = pi
    else
        theta = dacos(dx)
    endif
    if ((theta==0.d0) .or. (theta==pi)) then
        phi = 0.d0
    else
        dx = x/(r*dsin(theta))
        if (dx > 1.d0) then
            phi = 0.d0
        else if (dx < -1.d0) then
            phi = pi
        else
            phi = dacos(dx)
            if (y < 0.d0)   phi = 2.d0*pi - phi
        endif
    endif
endif


end subroutine cart2sph
