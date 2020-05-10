!-----------------------------------------------------------------
module module_polinterp   ! from the Numerical Recipes
!-----------------------------------------------------------------
  implicit none
  public :: polint, polin2, polint3D
  private :: locate
contains
!-----------------------------------------------------------------
SUBROUTINE polint3D(x1a,x2a,x3a,ya,n,x1,x2,x3,y)
!-----------------------------------------------------------------
USE nrtype
IMPLICIT NONE
REAL(DP), DIMENSION(:), INTENT(IN) :: x1a,x2a,x3a
REAL(DP), DIMENSION(:,:,:), INTENT(IN) :: ya
REAL(DP), INTENT(IN) :: x1,x2,x3
REAL(DP), INTENT(OUT) :: y
INTEGER(I4B), INTENT(IN) :: n
INTEGER(I4B) :: i_count, k, klo,khi, jlo,jhi, ilo,ihi
REAL(DP) :: dy
REAL(DP), DIMENSION(1:n) :: ytmp
ilo = locate(x1a,x1)
ilo = min(max(ilo-(n-1)/2,1),size(x1a)+1-n)
ihi = ilo+n-1
jlo = locate(x2a,x2)
jlo = min(max(jlo-(n-1)/2,1),size(x2a)+1-n)
jhi = jlo+n-1
klo = locate(x3a,x3)
klo = min(max(klo-(n-1)/2,1),size(x3a)+1-n)
khi = klo+n-1
i_count = 1
do k = klo,khi
   call polin2 (x1a(ilo:ihi), x2a(jlo:jhi), ya(ilo:ihi,jlo:jhi,k), &
                x1, x2, ytmp(i_count), dy)
   i_count = i_count+1
enddo
call polint (x3a(klo:khi), ytmp, x3, y, dy)
!-----------------------------------------------------------------
END SUBROUTINE polint3D
!-----------------------------------------------------------------
!-----------------------------------------------------------------
SUBROUTINE polin2(x1a,x2a,ya,x1,x2,y,dy)
!-----------------------------------------------------------------
USE nrtype; USE nrutil, ONLY : assert_eq
IMPLICIT NONE
REAL(DP), DIMENSION(:), INTENT(IN) :: x1a,x2a
REAL(DP), DIMENSION(:,:), INTENT(IN) :: ya
REAL(DP), INTENT(IN) :: x1,x2
REAL(DP), INTENT(OUT) :: y,dy
INTEGER(I4B) :: j,m,ndum
REAL(DP), DIMENSION(size(x1a)) :: ymtmp
REAL(DP), DIMENSION(size(x2a)) :: yntmp
m=assert_eq(size(x1a),size(ya,1),'polin2: m')
ndum=assert_eq(size(x2a),size(ya,2),'polin2: ndum')
do j=1,m
   yntmp=ya(j,:)
   call polint(x2a,yntmp,x2,ymtmp(j),dy)
end do
call polint(x1a,ymtmp,x1,y,dy)
!-----------------------------------------------------------------
END SUBROUTINE polin2
!-----------------------------------------------------------------
!-----------------------------------------------------------------
SUBROUTINE polint(xa,ya,x,y,dy)
!-----------------------------------------------------------------
USE nrtype; USE nrutil, ONLY : assert_eq,iminloc,nrerror
IMPLICIT NONE
REAL(DP), DIMENSION(:), INTENT(IN) :: xa,ya
REAL(DP), INTENT(IN) :: x
REAL(DP), INTENT(OUT) :: y,dy
INTEGER(I4B) :: m,n,ns
REAL(DP), DIMENSION(size(xa)) :: c,d,den,ho
n=assert_eq(size(xa),size(ya),'polint')
c=ya
d=ya
ho=xa-x
ns=iminloc(dabs(x-xa))
y=ya(ns)
ns=ns-1
do m=1,n-1 
   den(1:n-m)=ho(1:n-m)-ho(1+m:n)
   if (any(den(1:n-m) == 0.0))   call nrerror('polint: calculation failure')
   den(1:n-m)=(c(2:n-m+1)-d(1:n-m))/den(1:n-m)
   d(1:n-m)=ho(1+m:n)*den(1:n-m)
   c(1:n-m)=ho(1:n-m)*den(1:n-m)
   if (2*ns < n-m) then
      dy=c(ns+1)
   else
      dy=d(ns)
      ns=ns-1
   end if
   y=y+dy
end do
!-----------------------------------------------------------------
END SUBROUTINE polint
!-----------------------------------------------------------------
!-----------------------------------------------------------------
FUNCTION locate(xx,x)
!-----------------------------------------------------------------
  USE nrtype
  IMPLICIT NONE
  REAL(DP), DIMENSION(:), INTENT(IN) :: xx
  REAL(DP), INTENT(IN) :: x
  INTEGER(I4B) :: locate
  INTEGER(I4B) :: n,jl,jm,ju
  LOGICAL :: ascnd
  n=size(xx)
  ascnd = (xx(n) >= xx(1))
  jl=0
  ju=n+1
  do
     if (ju-jl <= 1) exit
     jm=(ju+jl)/2
     if (ascnd .eqv. (x >= xx(jm))) then
        jl=jm
     else
        ju=jm
     end if
  end do
  if (x == xx(1)) then
     locate=1
  else if (x == xx(n)) then
     locate=n-1
  else
     locate=jl
  end if
!-----------------------------------------------------------------
END FUNCTION locate
!-----------------------------------------------------------------
!-----------------------------------------------------------------
end module module_polinterp
!-----------------------------------------------------------------
