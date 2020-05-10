module angles


doubleprecision, parameter :: pi = 3.141592653589793238462643383279502884197, &
                              Rterre = 6371000.d0


contains

! #######################################################
subroutine cart2sph (x,y,z,r,theta,phi)


implicit none

doubleprecision, intent(IN) :: x,y,z
doubleprecision, intent(OUT) :: r,theta,phi
doubleprecision :: dx


r = dsqrt(x**2 + y**2 + z**2)
if (r==0.d0) then
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


return
end subroutine

! #######################################################
subroutine vect_sph2cart (position_pt,n,input_vect,output_vect)


implicit none

integer, intent(IN) :: n
doubleprecision, dimension(0:2), intent(IN) :: position_pt
doubleprecision, dimension(0:n-1,0:2), intent(INOUT) :: input_vect, output_vect

integer :: i,j
doubleprecision :: r,theta,phi, ct,st,cp,sp
doubleprecision, dimension(0:2,0:2) :: Pcs


call cart2sph(position_pt(0),position_pt(1),position_pt(2),r,theta,phi)
ct=dcos(theta)
st=dsin(theta)
cp=dcos(phi)
sp=dsin(phi)
!   matrice de passage
Pcs(0,0) = st*cp; Pcs(0,1) = ct*cp; Pcs(0,2) = -sp
Pcs(1,0) = st*sp; Pcs(1,1) = ct*sp; Pcs(1,2) = cp
Pcs(2,0) = ct   ; Pcs(2,1) = -st  ; Pcs(2,2) = 0.d0
!   projection
do i = 0,2
   output_vect(:,i) = 0.d0
   do j = 0,2
      output_vect(:,i) = output_vect(:,i) + Pcs(i,j)*input_vect(:,j)
   enddo
enddo


return
end subroutine

! #######################################################
subroutine arc_gd_cercle (longA,longB,colatA,colatB,cosdelta,sindelta)


implicit none

doubleprecision, intent(IN) :: longA,longB,colatA,colatB
doubleprecision, intent(OUT) :: cosdelta,sindelta


cosdelta = dcos(dabs(longA-longB)) * dsin(colatA) * dsin(colatB) &
           + dcos(colatA) * dcos(colatB)

sindelta = dsqrt(1.d0-cosdelta**2)


return
end subroutine

! #######################################################
subroutine barycentre3d (nb_pts,pts,res)


implicit none

integer, intent(IN) :: nb_pts
doubleprecision, dimension(1:nb_pts,1:3), intent(IN) :: pts
doubleprecision, dimension(1:3), intent(OUT) :: res

integer :: i,j


do i = 1,3
   res(i) = 0.d0
   do j = 1,nb_pts
      res(i) = res(i) + pts(j,i)
   enddo
enddo
res = res/dble(nb_pts)


return
end subroutine barycentre3d

! #######################################################
subroutine wtcoef (f,f1,f2,f3,f4,wt)


   implicit none

   doubleprecision, intent(in) ::  f,f1,f2,f3,f4
   doubleprecision, intent(out)::  wt


   if (f3.gt.f4) stop 'wtcoef: f3>f4 '
   if (f1.gt.f2) stop 'wtcoef: f1>f2 '
   if (f.le.f3.and.f.ge.f2) then
      wt = 1.
   else if (f.gt.f4.or.f.lt.f1 ) then
      wt = 0.
   else if (f.gt.f3.and.f.le.f4) then
      wt = 0.5*(1.0+cos(pi*(f-f3)/(f4-f3)))
   else if (f.ge.f1.and.f.lt.f2) then
      wt = 0.5*(1.0+cos(pi*(f-f2)/(f2-f1)))
   endif


return
end subroutine wtcoef

! #######################################################
doubleprecision function deg2rad(val)

 doubleprecision :: val

 deg2rad = pi*val/180.d0

end function

! ########################################################
doubleprecision function rad2deg(val)

 doubleprecision :: val

 rad2deg = 180.d0*val/pi

end function

! ########################################################
doubleprecision function distance3d(u,v)

 implicit none

 doubleprecision, dimension(0:2), intent(IN) :: u,v

 distance3d = dsqrt ( (u(0)-v(0))**2 + (u(1)-v(1))**2 + (u(2)-v(2))**2 )

end function distance3d

! ########################################################
end module angles
