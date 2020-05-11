module read_model !!!  1x1x1km random cubes embedded in a homogeneous medium !!!
! Routine get_region_param is used here instead of get_value_aniso.

implicit none

public :: get_value, get_value_aniso, get_region_param
private :: nb_region, region, region_bary, tab_rho,tab_lambda,tab_mu, first_time

integer :: nb_region, region
doubleprecision, dimension(:), allocatable :: tab_rho,tab_lambda,tab_mu
doubleprecision, dimension(:,:), allocatable :: region_bary
logical :: first_time = .true.

contains

! #######################################################
subroutine get_value (r,theta_rad,phi_rad,rho,vp,vs,Qmu,moho)

implicit none

integer, optional, intent(IN) :: moho
doubleprecision, intent(IN) :: r,theta_rad,phi_rad
doubleprecision, intent(INOUT) :: rho,vp,vs,Qmu


stop 'get_value IS NOT IMPLEMENTED IN THIS MODULE !!!'


end subroutine get_value

! #######################################################
subroutine get_value_aniso (r,theta_rad,phi_rad,rho,A,C,F,L,M,Gc,Gs,Hc,Hs,Bc,Bs,Ec,Es,Qmu,moho)

implicit none

integer, optional, intent(INOUT) :: moho
doubleprecision, intent(IN) :: r,theta_rad,phi_rad
doubleprecision, intent(INOUT) :: rho,A,C,F,L,M,Gc,Gs,Hc,Hs,Bc,Bs,Ec,Es,Qmu


moho = 10   ! then get_region_param is used in define_arrays.f90


end subroutine get_value_aniso

! ########################################################
subroutine get_region_param (bary,gll,rho,C)

use angles, only : distance3d

implicit none

integer, intent(IN) :: gll
doubleprecision, dimension(1:3), intent(IN) :: bary
doubleprecision, intent(OUT) :: rho
doubleprecision, dimension(1:6,1:6), intent(OUT) :: C

integer :: i,j, ios, n_reg
doubleprecision :: lambda,mu


if (first_time) then
   first_time = .false.
   open (11,file="elastic_regions",form="formatted",status="old",iostat=ios)
   if (ios>0)   stop 'TEXT FILE elastic_regions IS NEEDED TO GET THE ELASTIC PARAMETERS OF THE MEDIUM'
   read (11,*) nb_region
   allocate (tab_rho(0:nb_region))
   allocate (tab_lambda(0:nb_region))
   allocate (tab_mu(0:nb_region))
   do i = 0,nb_region
      read (11,*) tab_rho(i), tab_lambda(i), tab_mu(i)
   enddo 
   close (11)
   open (12,file="bary_region",form="formatted",status="old",iostat=ios)
   if (ios>0)   stop 'TEXT FILE bary_region IS NEEDED TO GET THE ELASTIC PARAMETERS OF THE MEDIUM'
   read (12,*) n_reg
   if (n_reg/=nb_region)   stop 'bary_region AND elastic_regions ARE INCOMPATIBLE: DIFFERENT NUMBER OF REGIONS'
   allocate (region_bary(1:nb_region,1:3))
   do i = 1,nb_region
      read (12,*) j, region_bary(i,1:3)
   enddo
   close (12)
   region_bary = region_bary/1000.d0
endif

if (gll==0) then
   region = 0
   if (bary(1)>13.d0 .and. bary(1)<113.d0 .and. bary(2)>13.d0.and. bary(2)<113.d0 .and. bary(3)>13.d0 .and. bary(3)<113.d0) then
      search : do i = 1,nb_region
         if (distance3d(bary(1:3),region_bary(i,1:3))<0.1) then
            region = i
            exit search
         endif
      enddo search
   endif
endif

rho = tab_rho(region)
lambda = tab_lambda(region)
mu = tab_mu(region)
C(:,:) = 0.d0
C(1,1) = 2.d0*mu + lambda;   C(2,2) = C(1,1);   C(3,3) = C(1,1)
C(1,2) = lambda;   C(1,3) = lambda;   C(2,3) = lambda
C(4,4) = 2.d0*mu;   C(5,5) = C(4,4);   C(6,6) = C(4,4)   ! Mandel's notation
do i = 2,6
   do j = 1,i-1
      C(i,j) = C(j,i)
   enddo
enddo


end subroutine get_region_param

! #######################################################
end module read_model
