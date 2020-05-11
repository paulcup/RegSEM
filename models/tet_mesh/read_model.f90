module read_model !!! from a tetgen mesh with regions !!! 

use read_input
use kdtree2_module

implicit none

public :: get_value, get_value_aniso
private :: mesh1, tree, results, nb_regions, density, Pspeed, Sspeed, &
           inv_bary_mat, first_time

type(mesh) :: mesh1
type(kdtree2), pointer :: tree
type(kdtree2_result), allocatable :: results(:)
integer :: nb_regions
doubleprecision, dimension(:), allocatable :: density, Pspeed, Sspeed
doubleprecision, dimension(:,:,:), allocatable :: inv_bary_mat
logical :: first_time = .true.

contains

! #######################################################
subroutine get_value (r,theta_rad,phi_rad,rho,vp,vs,Qmu,moho)


implicit none

integer, optional, intent(IN) :: moho
doubleprecision, intent(IN) :: r,theta_rad,phi_rad
doubleprecision, intent(OUT) :: rho,vp,vs,Qmu

integer :: i,j,k, ii,jj,kk, region
doubleprecision :: x, lambda, mu
doubleprecision, dimension(1:3) :: shift, pts, vect, bary_coord
doubleprecision, dimension(1:3,1:3) :: bary_mat, tmp_mat
doubleprecision, dimension(1:4,1:3) :: vertices

integer, parameter :: n = 6   ! Nb of nearest nodes for each voxet point

shift(1) = 0.d0;   shift(2) = 0.d0;   shift(3) = -6000.d0


if (first_time) then

   first_time = .false.

   write (mesh1%meshtype,'(a)') "tetgen"
   write (mesh1%meshname,'(a)') "../model_brep_90_simplify_meshed_new_options"
   call load_mesh (mesh1)

   open (21,file="elastic_regions")
      read (21,*) nb_regions
      allocate (density(1:nb_regions), Pspeed(1:nb_regions), Sspeed(1:nb_regions))
      do i = 1,nb_regions
         read (21,*) density(i), lambda, mu
         Pspeed(i) = dsqrt((lambda+2.d0*mu)/density(i))
         Sspeed(i) = dsqrt(mu/density(i))
      enddo
   close (21)

   ! Defining the barycentric basis associated with each element
   allocate (inv_bary_mat(1:mesh1%nb_ele,1:3,1:3))
   do i = 1,mesh1%nb_ele
      do j = 1,4
         k = mesh1%ele_node(i,j)
         vertices(j,1) = mesh1%node_coord(1,k)
         vertices(j,2) = mesh1%node_coord(2,k)
         vertices(j,3) = mesh1%node_coord(3,k)
      enddo
      do j = 1,3
       do k = 1,3
          bary_mat(k,j) = vertices(j+1,k) - vertices(1,k)
       enddo
      enddo
      call invert_3x3mat (bary_mat, x, tmp_mat)
      inv_bary_mat(i,:,:) = tmp_mat(:,:)
   enddo

   tree => kdtree2_create(mesh1%node_coord)
   allocate (results(1:n))

endif

! Finding the n nearest nodes of the current voxet point
pts(1) = r+shift(1);   pts(2) = theta_rad+shift(2);   pts(3) = phi_rad+shift(3)
call kdtree2_n_nearest (tree, pts, n, results)
region = -1
scan_nearest_nodes : do jj = 1,n
   i = results(jj)%idx
   ! Finding the element of the current voxet point
   scan_elem : do j = 1,mesh1%node_valence(i)
      k = mesh1%node_ele(i,j)
      ii = mesh1%ele_node(k,1)
      vect(1) = pts(1) - mesh1%node_coord(1,ii)
      vect(2) = pts(2) - mesh1%node_coord(2,ii)
      vect(3) = pts(3) - mesh1%node_coord(3,ii)
      tmp_mat(:,:) = inv_bary_mat(k,:,:)
      bary_coord = matmul(tmp_mat,vect)
      if (all(bary_coord>=0.d0) .and. (sum(bary_coord)<=1.d0)) then
         region = mesh1%rmark(k)
         exit scan_nearest_nodes
      endif
   enddo scan_elem
enddo scan_nearest_nodes
! Si on est en dehors du maillage...
if (region==-1)   region = 3 ! alors pour l'overthrust on est tout
                             ! en-dessous, i.e. dans la region 3

vs = Sspeed(region)
vp = Pspeed(region)
rho = density(region)
Qmu = 300.d0


end subroutine get_value

! #######################################################
subroutine get_value_aniso (r,theta_rad,phi_rad,rho,A,C,F,L,M,Gc,Gs,Hc,Hs,Bc,Bs,Ec,Es,Qmu,moho)


implicit none

integer, optional, intent(IN) :: moho
doubleprecision, intent(IN) :: r,theta_rad,phi_rad
doubleprecision, intent(INOUT) :: rho,A,C,F,L,M,Gc,Gs,Hc,Hs,Bc,Bs,Ec,Es,Qmu


stop 'get_value_aniso IS NOT IMPLEMENTED IN THIS CASEi YET !!!'


end subroutine get_value_aniso

!------------------------------------------------------------------------------------
subroutine invert_3x3mat(A,det,invA)


doubleprecision, dimension(1:3,1:3), intent(IN) :: A
doubleprecision, intent(OUT) :: det
doubleprecision, dimension(1:3,1:3), intent(OUT) :: invA

doubleprecision, dimension(1:3,1:3) :: tmp


det = A(1,1) * (A(2,2)*A(3,3) - A(3,2)*A(2,3)) + &
      A(2,1) * (A(1,3)*A(3,2) - A(1,2)*A(3,3)) + &
      A(3,1) * (A(1,2)*A(2,3) - A(1,3)*A(2,2))
if (det==0.d0) then
   stop 'NON INVERTIBLE 3x3 MATRIX'
else
   tmp(1,1) = A(2,2)*A(3,3) - A(3,2)*A(2,3)
   tmp(2,1) = A(1,3)*A(3,2) - A(1,2)*A(3,3)
   tmp(3,1) = A(1,2)*A(2,3) - A(1,3)*A(2,2)
   tmp(1,2) = A(2,3)*A(3,1) - A(2,1)*A(3,3)
   tmp(2,2) = A(1,1)*A(3,3) - A(1,3)*A(3,1)
   tmp(3,2) = A(1,3)*A(2,1) - A(1,1)*A(2,3)
   tmp(1,3) = A(2,1)*A(3,2) - A(2,2)*A(3,1)
   tmp(2,3) = A(1,2)*A(3,1) - A(3,2)*A(1,1)
   tmp(3,3) = A(2,2)*A(1,1) - A(1,2)*A(2,1)
endif
invA = transpose(tmp)
invA = invA/det


end subroutine invert_3x3mat

! #######################################################
end module read_model
