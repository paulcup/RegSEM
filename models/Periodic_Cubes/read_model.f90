module read_model !!!  A periodic medium made of eight 1x1x1km random cubes !!!
! Routine get_region_param is used here instead of get_value_aniso.

implicit none

public :: get_value, get_value_aniso, get_region_param, count_region
private :: nb_region, region_vect, tab_rho,tab_lambda,tab_mu, first_time

integer, parameter :: nb_region = 8
integer, dimension(1:nb_region) :: count_region = 0
integer, dimension(1:nb_region,1:3) :: region_vect
doubleprecision, dimension(1:nb_region) :: tab_rho,tab_lambda,tab_mu
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
subroutine get_region_param (bary,rho,C)

implicit none

doubleprecision, dimension(1:3), intent(IN) :: bary
doubleprecision, intent(OUT) :: rho
doubleprecision, dimension(1:6,1:6), intent(OUT) :: C

integer :: i,j, ios, n_reg, integer_part, region
integer, dimension(1:3) :: vect
doubleprecision :: lambda,mu


if (first_time) then
   first_time = .false.
   region_vect(1,1:3) = (/ 0, 0, 0 /)
   region_vect(2,1:3) = (/ 1, 0, 0 /)
   region_vect(3,1:3) = (/ 0, 1, 0 /)
   region_vect(4,1:3) = (/ 1, 1, 0 /)
   region_vect(5,1:3) = (/ 0, 0, 1 /)
   region_vect(6,1:3) = (/ 1, 0, 1 /)
   region_vect(7,1:3) = (/ 0, 1, 1 /)
   region_vect(8,1:3) = (/ 1, 1, 1 /)
   open (11,file="elastic_regions",form="formatted",status="old",iostat=ios)
   if (ios>0)   stop 'TEXT FILE elastic_regions IS NEEDED TO GET THE ELASTIC PARAMETERS OF THE MEDIUM'
   read (11,*) n_reg
   if (n_reg/=nb_region)   stop 'get_region_param AND elastic_param ARE INCOMPATIBLE: DIFFERENT NUMBER OF REGIONS'
   do i = 1,nb_region
      read (11,*) tab_rho(i), tab_lambda(i), tab_mu(i)
   enddo 
   close (11)
endif

do i = 1,3
   integer_part = int(bary(i))
   vect(i) = mod(integer_part,2)
enddo

region = 0
do i = 1,8
   if (vect(1)==region_vect(i,1) .and. vect(2)==region_vect(i,2) .and. vect(3)==region_vect(i,3)) then
      region = i
      count_region(i) = count_region(i) + 1
   endif
enddo
if (region==0)   stop 'PROBLEM IN get_region_param'

rho = tab_rho(region)
lambda = tab_lambda(region)
mu = tab_mu(region)
C(:,:) = 0.d0
C(1,1) = 2.d0*mu + lambda;   C(2,2) = 2.d0*mu + lambda;   C(3,3) = 2.d0*mu + lambda
C(1,2) = lambda;   C(1,3) = lambda;   C(2,3) = lambda
C(4,4) = 2.d0*mu;   C(5,5) = 2.d0*mu;   C(6,6) = 2.d0*mu   ! Mandel's notation
do i = 2,6
   do j = 1,i-1
      C(i,j) = C(j,i)
   enddo
enddo


end subroutine get_region_param

! #######################################################
end module read_model
