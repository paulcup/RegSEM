subroutine find_normales (Tdomain, n)


use sdomains

implicit none

type (Domain), intent (INOUT) :: Tdomain
integer, intent (IN) :: n

integer :: mat, ngllx,nglly,ngllz, i,k
doubleprecision :: a,b,c
doubleprecision, dimension(0:2) :: u,v, corner0,corner1,corner2,corner3
doubleprecision, dimension(0:2,0:2) :: TMP


mat = Tdomain%specel(n)%mat_index
ngllx = Tdomain%specel(n)%ngllx
nglly = Tdomain%specel(n)%nglly
ngllz = Tdomain%specel(n)%ngllz

do k = 0,1
  if (k==0) then
    if (Tdomain%sSubDomain(mat)%Left) then
        i = Tdomain%specel(n)%Iglobnum(ngllx-1,0,0)
        corner0(:) = Tdomain%Globcoord(:,i)
        i = Tdomain%specel(n)%Iglobnum(ngllx-1,nglly-1,0)
        corner1(:) = Tdomain%Globcoord(:,i)
        i = Tdomain%specel(n)%Iglobnum(ngllx-1,0,ngllz-1)
        corner2(:) = Tdomain%Globcoord(:,i)
        i = Tdomain%specel(n)%Iglobnum(ngllx-1,nglly-1,ngllz-1)
        corner3(:) = Tdomain%Globcoord(:,i)
    else
        i = Tdomain%specel(n)%Iglobnum(0,0,0)
        corner0(:) = Tdomain%Globcoord(:,i)
        i = Tdomain%specel(n)%Iglobnum(0,nglly-1,0)
        corner1(:) = Tdomain%Globcoord(:,i)
        i = Tdomain%specel(n)%Iglobnum(0,0,ngllz-1)
        corner2(:) = Tdomain%Globcoord(:,i)
        i = Tdomain%specel(n)%Iglobnum(0,nglly-1,ngllz-1)
        corner3(:) = Tdomain%Globcoord(:,i)
    endif
  else
    if (Tdomain%sSubDomain(mat)%Forward) then
        i = Tdomain%specel(n)%Iglobnum(0,nglly-1,0)
        corner0(:) = Tdomain%Globcoord(:,i)
        i = Tdomain%specel(n)%Iglobnum(ngllx-1,nglly-1,0)
        corner1(:) = Tdomain%Globcoord(:,i)
        i = Tdomain%specel(n)%Iglobnum(0,nglly-1,ngllz-1)
        corner2(:) = Tdomain%Globcoord(:,i)
        i = Tdomain%specel(n)%Iglobnum(ngllx-1,nglly-1,ngllz-1)
        corner3(:) = Tdomain%Globcoord(:,i)
    else
        i = Tdomain%specel(n)%Iglobnum(0,0,0)
        corner0(:) = Tdomain%Globcoord(:,i)
        i = Tdomain%specel(n)%Iglobnum(ngllx-1,0,0)
        corner1(:) = Tdomain%Globcoord(:,i)
        i = Tdomain%specel(n)%Iglobnum(0,0,ngllz-1)
        corner2(:) = Tdomain%Globcoord(:,i)
        i = Tdomain%specel(n)%Iglobnum(ngllx-1,0,ngllz-1)
        corner3(:) = Tdomain%Globcoord(:,i)
    endif
  endif
  u = corner0 - corner3
  if (u(0)==0.d0) then
      v = u
      u = corner1 - corner2
  else
      v = corner1 - corner2
  endif
  if (k==0 .and. u(0)==0.d0) then
      Tdomain%specel(n)%Normales(0,k) = 1.d0
      Tdomain%specel(n)%Normales(1,k) = 0.d0
      Tdomain%specel(n)%Normales(2,k) = 0.d0
      !!! (u(1)==0 .and. v(1)==0) indicates that face 2 or 4 is in the (x,z) plan.
      !!! It is impossible in the cases we consider.
      !!! (v(0)==0 .and. v(1)==0) indicates that v==z.
      !!! It is impossible if Xi <= 90
  else if (k==1 .and. u(1)==0.d0 .and. v(1)==0.d0) then
      Tdomain%specel(n)%Normales(0,k) = 0.d0
      Tdomain%specel(n)%Normales(1,k) = 1.d0
      Tdomain%specel(n)%Normales(2,k) = 0.d0
      !!! u(0)==0 implies that v(0)==0 and indicates that face 1 or 3 is in the (y,z) plan.
      !!! It is impossible in the cases we consider.
      !!! (v(0)==0 .and. v(1)==0) indicates that v==z.
      !!! It is impossible if Eta <= 90
  else
      a = u(1)/u(0)
      b = u(2)/u(0)
      c = (v(0)*u(2)-v(2)*u(0)) / (v(1)*u(0)-v(0)*u(1))
      Tdomain%specel(n)%Normales(2,k) = -1. / dsqrt((a**2+1.d0)*c**2 + 2.d0*a*b*c + b**2 + 1.d0)
      Tdomain%specel(n)%Normales(1,k) = c*Tdomain%specel(n)%Normales(2,k)
      Tdomain%specel(n)%Normales(0,k) = -a*Tdomain%specel(n)%Normales(1,k) &
                                        -b*Tdomain%specel(n)%Normales(2,k)
  endif
enddo

Tdomain%specel(n)%Normales(0,2) = Tdomain%specel(n)%Normales(1,0) * Tdomain%specel(n)%Normales(2,1) &
                                  - Tdomain%specel(n)%Normales(2,0) * Tdomain%specel(n)%Normales(1,1)
Tdomain%specel(n)%Normales(1,2) = Tdomain%specel(n)%Normales(2,0) * Tdomain%specel(n)%Normales(0,1) &
                                  - Tdomain%specel(n)%Normales(0,0) * Tdomain%specel(n)%Normales(2,1)
Tdomain%specel(n)%Normales(2,2) = Tdomain%specel(n)%Normales(0,0) * Tdomain%specel(n)%Normales(1,1) &
                                  - Tdomain%specel(n)%Normales(1,0) * Tdomain%specel(n)%Normales(0,1)
if (Tdomain%sSubDomain(mat)%Down) then
    if (Tdomain%specel(n)%Normales(2,2)>0) Tdomain%specel(n)%Normales(:,2) = -Tdomain%specel(n)%Normales(:,2)
else
    if (Tdomain%specel(n)%Normales(2,2)<0) Tdomain%specel(n)%Normales(:,2) = -Tdomain%specel(n)%Normales(:,2)
endif

TMP = Tdomain%specel(n)%Normales
call invert_3d(TMP,a)
Tdomain%specel(n)%Inv_Normales = TMP


end subroutine find_normales
