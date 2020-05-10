doubleprecision function pow (x,vp,npm,dx,A,np)


implicit none

doubleprecision :: x,vp,dx,A
integer :: npm,np

doubleprecision :: rnpm,pp1


rnpm = dble(npm)

pp1 = dx
pp1 = 1.d0/pp1
pow = A * vp * pp1 * (x/rnpm)**np

return
end function
