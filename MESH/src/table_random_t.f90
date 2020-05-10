subroutine table_random_t (rg, n, mat, coord_node25, moho, topo, rot, Rterre, pi, ellipt, Rsph)


implicit none

logical, intent(IN) :: ellipt
integer, intent(IN) :: rg, n, mat, moho, topo
doubleprecision, intent(IN) :: Rterre, pi
doubleprecision, intent(IN), optional :: Rsph
doubleprecision, dimension(0:2), intent(IN) :: coord_node25
doubleprecision, dimension(0:2,0:2), intent(IN) :: rot
integer :: i,j,k
doubleprecision :: xa,ya,za, xs,ys,zs, R,theta,phi


! Coordonnees cartesiennes du noeud 25 dans le chunk de ref
xa = coord_node25(0);   ya = coord_node25(1);   za = coord_node25(2)
! Rotation du chunk de ref vers le chunk reel
xs = rot(0,0)*xa + rot(0,1)*ya + rot(0,2)*za
ys = rot(1,0)*xa + rot(1,1)*ya + rot(1,2)*za
zs = rot(2,0)*xa + rot(2,1)*ya + rot(2,2)*za
! Passage en coordonnees spheriques
call cart2sph(xs,ys,zs,R,theta,phi)
theta = 90.d0 - 180.d0*theta/pi
phi = 180.d0*phi/pi

if (ellipt)   R = Rsph
if (topo==1 .or. moho==1 .or. (abs(R-Rterre)<1)) then   ! C'est un element en surface
   write (25,*) rg, n
   if (mat==0) then
      write (25,*) theta, phi, "T"
   else
      write (25,*) theta, phi, "F"
   endif
endif


end subroutine table_random_t
