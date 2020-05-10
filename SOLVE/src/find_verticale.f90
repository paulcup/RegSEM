subroutine find_verticale (Tdomain,n)


use sdomains

implicit none

type (Domain), intent (INOUT) :: Tdomain
integer, intent (IN) :: n

integer :: num, ngllx,nglly,ngllz, i,j,k, nf,ne,nv
doubleprecision, dimension(0:2) :: u,v, ptA,ptB,ptC,ptD, vect


ngllx = Tdomain%specel(n)%ngllx
nglly = Tdomain%specel(n)%nglly
ngllz = Tdomain%specel(n)%ngllz

nf = Tdomain%specel(n)%Near_Faces(5)
k = ngllz-1
do j = 1,nglly-2
 do i = 1,ngllx-2
     num = Tdomain%specel(n)%Iglobnum(i-1,j,k)
     ptA(:) = Tdomain%Globcoord(:,num)
     num = Tdomain%specel(n)%Iglobnum(i+1,j,k)
     ptB(:) = Tdomain%Globcoord(:,num)
     num = Tdomain%specel(n)%Iglobnum(i,j-1,k)
     ptC(:) = Tdomain%Globcoord(:,num)
     num = Tdomain%specel(n)%Iglobnum(i,j+1,k)
     ptD(:) = Tdomain%Globcoord(:,num)
     u(:) = ptB(:) - ptA(:)
     v(:) = ptD(:) - ptC(:)
     vect(0) = u(1)*v(2) - u(2)*v(1)
     vect(1) = u(2)*v(0) - u(0)*v(2)
     vect(2) = u(0)*v(1) - u(1)*v(0)
     call normalisation(vect)
     call vect2matrix(vect,Tdomain%sFace(nf)%verticale(i,j,0:2,0:2))
 enddo
enddo

ne = Tdomain%specel(n)%Near_Edges(5)
k = ngllz-1
j = 0
do i = 1,ngllx-2
    num = Tdomain%specel(n)%Iglobnum(i-1,j,k)
    ptA(:) = Tdomain%Globcoord(:,num)
    num = Tdomain%specel(n)%Iglobnum(i+1,j,k)
    ptB(:) = Tdomain%Globcoord(:,num)
    num = Tdomain%specel(n)%Iglobnum(i,j,k)
    ptC(:) = Tdomain%Globcoord(:,num)
    num = Tdomain%specel(n)%Iglobnum(i,j+1,k)
    ptD(:) = Tdomain%Globcoord(:,num)
    u(:) = ptB(:) - ptA(:)
    v(:) = ptD(:) - ptC(:)
    vect(0) = u(1)*v(2) - u(2)*v(1)
    vect(1) = u(2)*v(0) - u(0)*v(2)
    vect(2) = u(0)*v(1) - u(1)*v(0)
    call normalisation(vect)
    call vect2matrix(vect,Tdomain%sEdge(ne)%verticale(i,0:2,0:2))
enddo

ne = Tdomain%specel(n)%Near_Edges(9)
k = ngllz-1
j = nglly-1
do i = 1,ngllx-2
    num = Tdomain%specel(n)%Iglobnum(i-1,j,k)
    ptA(:) = Tdomain%Globcoord(:,num)
    num = Tdomain%specel(n)%Iglobnum(i+1,j,k)
    ptB(:) = Tdomain%Globcoord(:,num)
    num = Tdomain%specel(n)%Iglobnum(i,j-1,k)
    ptC(:) = Tdomain%Globcoord(:,num)
    num = Tdomain%specel(n)%Iglobnum(i,j,k)
    ptD(:) = Tdomain%Globcoord(:,num)
    u(:) = ptB(:) - ptA(:)
    v(:) = ptD(:) - ptC(:)
    vect(0) = u(1)*v(2) - u(2)*v(1)
    vect(1) = u(2)*v(0) - u(0)*v(2)
    vect(2) = u(0)*v(1) - u(1)*v(0)
    call normalisation(vect)
    call vect2matrix(vect,Tdomain%sEdge(ne)%verticale(i,0:2,0:2))
enddo

ne = Tdomain%specel(n)%Near_Edges(8)
k = ngllz-1
i = ngllx-1
do j = 1,nglly-2
    num = Tdomain%specel(n)%Iglobnum(i-1,j,k)
    ptA(:) = Tdomain%Globcoord(:,num)
    num = Tdomain%specel(n)%Iglobnum(i,j,k)
    ptB(:) = Tdomain%Globcoord(:,num)
    num = Tdomain%specel(n)%Iglobnum(i,j-1,k)
    ptC(:) = Tdomain%Globcoord(:,num)
    num = Tdomain%specel(n)%Iglobnum(i,j+1,k)
    ptD(:) = Tdomain%Globcoord(:,num)
    u(:) = ptB(:) - ptA(:)
    v(:) = ptD(:) - ptC(:)
    vect(0) = u(1)*v(2) - u(2)*v(1)
    vect(1) = u(2)*v(0) - u(0)*v(2)
    vect(2) = u(0)*v(1) - u(1)*v(0)
    call normalisation(vect)
    call vect2matrix(vect,Tdomain%sEdge(ne)%verticale(j,0:2,0:2))
enddo

ne = Tdomain%specel(n)%Near_Edges(11)
k = ngllz-1
i = 0
do j = 1,nglly-2
    num = Tdomain%specel(n)%Iglobnum(i,j,k)
    ptA(:) = Tdomain%Globcoord(:,num)
    num = Tdomain%specel(n)%Iglobnum(i+1,j,k)
    ptB(:) = Tdomain%Globcoord(:,num)
    num = Tdomain%specel(n)%Iglobnum(i,j-1,k)
    ptC(:) = Tdomain%Globcoord(:,num)
    num = Tdomain%specel(n)%Iglobnum(i,j+1,k)
    ptD(:) = Tdomain%Globcoord(:,num)
    u(:) = ptB(:) - ptA(:)
    v(:) = ptD(:) - ptC(:)
    vect(0) = u(1)*v(2) - u(2)*v(1)
    vect(1) = u(2)*v(0) - u(0)*v(2)
    vect(2) = u(0)*v(1) - u(1)*v(0)
    call normalisation(vect)
    call vect2matrix(vect,Tdomain%sEdge(ne)%verticale(j,0:2,0:2))
enddo

nv = Tdomain%specel(n)%Near_Vertices(4)
k = ngllz-1
j = 0
i = 0
num = Tdomain%specel(n)%Iglobnum(i,j,k)
ptA(:) = Tdomain%Globcoord(:,num)
num = Tdomain%specel(n)%Iglobnum(i+1,j,k)
ptB(:) = Tdomain%Globcoord(:,num)
num = Tdomain%specel(n)%Iglobnum(i,j,k)
ptC(:) = Tdomain%Globcoord(:,num)
num = Tdomain%specel(n)%Iglobnum(i,j+1,k)
ptD(:) = Tdomain%Globcoord(:,num)
u(:) = ptB(:) - ptA(:)
v(:) = ptD(:) - ptC(:)
vect(0) = u(1)*v(2) - u(2)*v(1)
vect(1) = u(2)*v(0) - u(0)*v(2)
vect(2) = u(0)*v(1) - u(1)*v(0)
call normalisation(vect)
call vect2matrix(vect,Tdomain%sVertex(nv)%verticale(0:2,0:2))

nv = Tdomain%specel(n)%Near_Vertices(5)
k = ngllz-1
j = 0
i = ngllx-1
num = Tdomain%specel(n)%Iglobnum(i-1,j,k)
ptA(:) = Tdomain%Globcoord(:,num)
num = Tdomain%specel(n)%Iglobnum(i,j,k)
ptB(:) = Tdomain%Globcoord(:,num)
num = Tdomain%specel(n)%Iglobnum(i,j,k)
ptC(:) = Tdomain%Globcoord(:,num)
num = Tdomain%specel(n)%Iglobnum(i,j+1,k)
ptD(:) = Tdomain%Globcoord(:,num)
u(:) = ptB(:) - ptA(:)
v(:) = ptD(:) - ptC(:)
vect(0) = u(1)*v(2) - u(2)*v(1)
vect(1) = u(2)*v(0) - u(0)*v(2)
vect(2) = u(0)*v(1) - u(1)*v(0)
call normalisation(vect)
call vect2matrix(vect,Tdomain%sVertex(nv)%verticale(0:2,0:2))

nv = Tdomain%specel(n)%Near_Vertices(6)
k = ngllz-1
j = nglly-1
i = ngllx-1
num = Tdomain%specel(n)%Iglobnum(i-1,j,k)
ptA(:) = Tdomain%Globcoord(:,num)
num = Tdomain%specel(n)%Iglobnum(i,j,k)
ptB(:) = Tdomain%Globcoord(:,num)
num = Tdomain%specel(n)%Iglobnum(i,j-1,k)
ptC(:) = Tdomain%Globcoord(:,num)
num = Tdomain%specel(n)%Iglobnum(i,j,k)
ptD(:) = Tdomain%Globcoord(:,num)
u(:) = ptB(:) - ptA(:)
v(:) = ptD(:) - ptC(:)
vect(0) = u(1)*v(2) - u(2)*v(1)
vect(1) = u(2)*v(0) - u(0)*v(2)
vect(2) = u(0)*v(1) - u(1)*v(0)
call normalisation(vect)
call vect2matrix(vect,Tdomain%sVertex(nv)%verticale(0:2,0:2))

nv = Tdomain%specel(n)%Near_Vertices(7)
k = ngllz-1
j = nglly-1
i = 0
num = Tdomain%specel(n)%Iglobnum(i,j,k)
ptA(:) = Tdomain%Globcoord(:,num)
num = Tdomain%specel(n)%Iglobnum(i+1,j,k)
ptB(:) = Tdomain%Globcoord(:,num)
num = Tdomain%specel(n)%Iglobnum(i,j-1,k)
ptC(:) = Tdomain%Globcoord(:,num)
num = Tdomain%specel(n)%Iglobnum(i,j,k)
ptD(:) = Tdomain%Globcoord(:,num)
u(:) = ptB(:) - ptA(:)
v(:) = ptD(:) - ptC(:)
vect(0) = u(1)*v(2) - u(2)*v(1)
vect(1) = u(2)*v(0) - u(0)*v(2)
vect(2) = u(0)*v(1) - u(1)*v(0)
call normalisation(vect)
call vect2matrix(vect,Tdomain%sVertex(nv)%verticale(0:2,0:2))


contains


 subroutine normalisation(vect)

 implicit none

 doubleprecision, dimension(0:2), intent(INOUT) :: vect
 doubleprecision :: norm

 norm = dsqrt (vect(0)**2 + vect(1)**2 + vect(2)**2)
 vect(:) = vect(:) / norm

 end subroutine normalisation


 subroutine vect2matrix(vect,mat)

 implicit none

 doubleprecision, dimension(0:2), intent(IN) :: vect
 doubleprecision, dimension(0:2,0:2), intent(OUT) :: mat
 integer :: i,j

 do i = 0,2
  do j = 0,2
      mat(i,j) = vect(i)*vect(j)
  enddo
 enddo

 end subroutine vect2matrix


end subroutine find_verticale
