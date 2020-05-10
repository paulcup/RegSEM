module read_mirror !!! An example of subroutine get_mirror !!! 


implicit none

public :: get_mirror


contains

! #######################################################
subroutine get_mirror(x,y,z,mirror)


implicit none

double precision :: x,y,z
integer :: mirror
double precision :: xc,yc,zc,r,dist


r = 112000*5
xc = 0.
yc = 0.
zc = 6371000.

dist = dsqrt((x-xc)**2.d0+(y-yc)**2.d0+(z-zc)**2.d0)

if(dist<=r)then
	mirror=1
else
	mirror=0
endif


end subroutine get_mirror


end module read_mirror
