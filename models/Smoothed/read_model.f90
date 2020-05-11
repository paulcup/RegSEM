module read_model !!! Polynomial interpolation of a smoothed medium !!!
! The routine to be used depends on the first number in Csmooth_rho.res:
! 1->get_value, 2->get_value_aniso, other->get_elastic_param

implicit none

public :: get_value, get_value_aniso, get_elastic_param, deload
private :: nb_pt_dir, model_bounds, rho_smooth, C_smooth, xseries,yseries,zseries, &
           first_time, order

integer, parameter :: order = 3   ! 2=linear, 3=quadratic, 4=cubic, etc...

integer, dimension(1:3) :: nb_pt_dir
doubleprecision, dimension(:), allocatable :: xseries,yseries,zseries
doubleprecision, dimension(1:3,1:2) :: model_bounds
doubleprecision, dimension(:,:,:), allocatable :: rho_smooth
doubleprecision, dimension(:,:,:,:), allocatable :: C_smooth
logical :: first_time = .true.

contains

! #######################################################
subroutine get_value (x,y,z,rho,vp,vs,Qmu,moho)

use module_polinterp

implicit none

integer, optional, intent(IN) :: moho
doubleprecision, intent(INOUT) :: x,y,z
doubleprecision, intent(OUT) :: rho,vp,vs,Qmu

integer :: sortie, i,j,k, ios, i_count
doubleprecision :: lambda, mu
doubleprecision, dimension(1:3) :: delta_dir


! Shifting the input point
x = x + 400.d0
y = y + 400.d0
z = z - 6200.d0

if (first_time) then
   first_time = .false.
   open (11,file="Csmooth_rho.res",form="unformatted",access="sequential",status="old",iostat=ios)
   if (ios>0)   stop 'FILE Csmooth_rho.res IS NEEDED TO GET THE ELASTIC PARAMETERS OF THE MEDIUM'
   read (11) sortie
   if (sortie==2) then
      print *, 'THE INPUT FILE IS FOR A TI MEDIUM'
      stop 'PLEASE USE get_value_aniso INSTEAD OF get_value !!!'
   else if (sortie/=1) then
      print *, 'THE INPUT FILE IS FOR A FULLY ANISOTROPIC MEDIUM'
      stop 'PLEASE USE get_elastic_param INSTEAD OF get_value !!!'
   endif
   do i = 1,3
      read (11) model_bounds(i,1), model_bounds(i,2), nb_pt_dir(i)
      delta_dir(i) = (model_bounds(i,2)-model_bounds(i,1)) / (nb_pt_dir(i)-1)
   enddo
   allocate (C_smooth(1:nb_pt_dir(1),1:nb_pt_dir(2),1:nb_pt_dir(3),1:2))
   allocate (rho_smooth(1:nb_pt_dir(1),1:nb_pt_dir(2),1:nb_pt_dir(3)))
   do k = 1,nb_pt_dir(3)
    do j = 1,nb_pt_dir(2)
     do i = 1,nb_pt_dir(1)
        read (11) rho_smooth(i,j,k)
        do i_count = 1,2
           read (11) C_smooth(i,j,k,i_count)   ! NB: the 1st parameter is Vp, the 2nd is Vs (in m/s)
        enddo
     enddo
    enddo
   enddo
   close (11)
   allocate (xseries(1:nb_pt_dir(1)))
   allocate (yseries(1:nb_pt_dir(2)))
   allocate (zseries(1:nb_pt_dir(3)))
   do i = 1,nb_pt_dir(1)
      xseries(i) = model_bounds(1,1) + (i-1)*delta_dir(1)
   enddo
   do i = 1,nb_pt_dir(2)
      yseries(i) = model_bounds(2,1) + (i-1)*delta_dir(2)
   enddo
   do i = 1,nb_pt_dir(3)
      zseries(i) = model_bounds(3,1) + (i-1)*delta_dir(3)
   enddo
endif

! Different solutions depending on what happens beyond the model bounds:
! 1) Homogeneous material:
!if (x<model_bounds(1,1) .or. x>model_bounds(1,2) .or. &
!    y<model_bounds(2,1) .or. y>model_bounds(2,2) .or. &
!    z<model_bounds(3,1) .or. z>model_bounds(3,2)) then
!   rho = 2999.76550425812;   lambda = 59746012739.2350;   mu = 36418757245.0581
!   vp = dsqrt((lambda+2.d0*mu)/rho);   vs = dsqrt(mu/rho)
!else
!   call polint3D (xseries, yseries, zseries, C_smooth(:,:,:,1), order, x, y, z, vp)
!   call polint3D (xseries, yseries, zseries, C_smooth(:,:,:,2), order, x, y, z, vs)
!   call polint3D (xseries, yseries, zseries, rho_smooth(:,:,:), order, x, y, z, rho)
!endif
! 2) Taking values at the closest bound
if (x<model_bounds(1,1))   x = model_bounds(1,1)
if (x>model_bounds(1,2))   x = model_bounds(1,2)
if (y<model_bounds(2,1))   y = model_bounds(2,1)
if (y>model_bounds(2,2))   y = model_bounds(2,2)
if (z<model_bounds(3,1))   z = model_bounds(3,1)
if (z>model_bounds(3,2))   z = model_bounds(3,2)
call polint3D (xseries, yseries, zseries, C_smooth(:,:,:,1), order, x, y, z, vp)
call polint3D (xseries, yseries, zseries, C_smooth(:,:,:,2), order, x, y, z, vs)
call polint3D (xseries, yseries, zseries, rho_smooth(:,:,:), order, x, y, z, rho)

Qmu = 300.d0


end subroutine get_value

! #######################################################
subroutine get_value_aniso (x,y,z,rho,A,C,F,L,M,Gc,Gs,Hc,Hs,Bc,Bs,Ec,Es,Qmu,moho)

implicit none

integer, optional, intent(INOUT) :: moho
doubleprecision, intent(IN) :: x,y,z
doubleprecision, intent(INOUT) :: rho,A,C,F,L,M,Gc,Gs,Hc,Hs,Bc,Bs,Ec,Es,Qmu


if (first_time) then
   open (11,file="Csmooth_rho.res",form="unformatted",access="sequential",status="old",iostat=ios)
   if (ios>0)   stop 'FILE Csmooth_rho.res IS NEEDED TO GET THE ELASTIC PARAMETERS OF THE MEDIUM'
   read (11) sortie
   if (sortie==1) then
      print *, 'THE INPUT FILE IS FOR AN ISOTROPIC MEDIUM'
      stop 'PLEASE USE get_value !!!'
   else if (sortie==2) then
      print *, 'THE INPUT FILE IS FOR A TI MEDIUM'
      stop 'get_value_aniso IS NOT IMPLEMENTED FOR THAT CASE YET !!!'
   endif
   close (11)
endif

moho = 10   ! then get_elastic_param is used in define_arrays.f90


end subroutine get_value_aniso

! ########################################################
subroutine get_elastic_param (x,y,z,rho,C,dom_bounds)

use module_polinterp

implicit none

doubleprecision, intent(INOUT) :: x,y,z
doubleprecision, intent(OUT) :: rho
doubleprecision, dimension(1:3,1:2), intent(IN) :: dom_bounds   ! not used yet
doubleprecision, dimension(1:6,1:6), intent(OUT) :: C

integer :: sortie, i,j,k, ios, i_count
doubleprecision :: lambda, mu
doubleprecision, dimension(1:3) :: delta_dir
doubleprecision, dimension(1:21) :: tmp_C


if (first_time) then
   first_time = .false.
   open (11,file="Csmooth_rho.res",form="unformatted",access="sequential",status="old",iostat=ios)
   if (ios>0)   stop 'FILE Csmooth_rho.res IS NEEDED TO GET THE ELASTIC PARAMETERS OF THE MEDIUM'
   read (11) sortie
   if (sortie==1) then
      print *, 'THE INPUT FILE IS FOR AN ISOTROPIC MEDIUM'
      stop 'PLEASE USE get_value INSTEAD OF get_elastic_param !!!'
   else if (sortie==2) then
      print *, 'THE INPUT FILE IS FOR A TI MEDIUM'
      stop 'PLEASE USE get_value_aniso INSTEAD OF get_elastic_param !!!'
   endif
   do i = 1,3
      read (11) model_bounds(i,1), model_bounds(i,2), nb_pt_dir(i)
      delta_dir(i) = (model_bounds(i,2)-model_bounds(i,1)) / (nb_pt_dir(i)-1)
   enddo
   allocate (C_smooth(1:nb_pt_dir(1),1:nb_pt_dir(2),1:nb_pt_dir(3),1:21))
   allocate (rho_smooth(1:nb_pt_dir(1),1:nb_pt_dir(2),1:nb_pt_dir(3)))
   do k = 1,nb_pt_dir(3)
    do j = 1,nb_pt_dir(2)
     do i = 1,nb_pt_dir(1)
        read (11) rho_smooth(i,j,k)
        do i_count = 1,21
           read (11) C_smooth(i,j,k,i_count)   ! NB: Cij values in Csmooth_rho.res are in Mandel's notation
        enddo
     enddo
    enddo
   enddo
   close (11)
   allocate (xseries(1:nb_pt_dir(1)))
   allocate (yseries(1:nb_pt_dir(2)))
   allocate (zseries(1:nb_pt_dir(3)))
   do i = 1,nb_pt_dir(1)
      xseries(i) = model_bounds(1,1) + (i-1)*delta_dir(1)
   enddo
   do i = 1,nb_pt_dir(2)
      yseries(i) = model_bounds(2,1) + (i-1)*delta_dir(2)
   enddo
   do i = 1,nb_pt_dir(3)
      zseries(i) = model_bounds(3,1) + (i-1)*delta_dir(3)
   enddo
endif

! Different solutions depending on what happens beyond the model bounds:
! 1) Homogeneous material:
!if (x<model_bounds(1,1) .or. x>model_bounds(1,2) .or. &
!    y<model_bounds(2,1) .or. y>model_bounds(2,2) .or. &
!    z<model_bounds(3,1) .or. z>model_bounds(3,2)) then
!   rho = 3003.74107037824;   lambda = 55260199143.1541;   mu = 42994835872.9560
!   tmp_C(:) = 0.d0
!   tmp_C(1) = 2.d0*mu + lambda;   tmp_C(7) = 2.d0*mu + lambda;   tmp_C(12) = 2.d0*mu + lambda
!   tmp_C(16) = 2.d0*mu;   tmp_C(19) = 2.d0*mu;   tmp_C(21) = 2.d0*mu
!   tmp_C(2) = lambda;   tmp_C(3) = lambda;   tmp_C(8) = lambda
!else
!   do i_count = 1,21
!      call polint3D (xseries, yseries, zseries, C_smooth(:,:,:,i_count), order, x, y, z, tmp_C(i_count))
!   enddo
!   call polint3D (xseries, yseries, zseries, rho_smooth(:,:,:), order, x, y, z, rho)
!endif
! 2) Taking values at the closest bound
if (x<model_bounds(1,1))   x = model_bounds(1,1)
if (x>model_bounds(1,2))   x = model_bounds(1,2)
if (y<model_bounds(2,1))   y = model_bounds(2,1)
if (y>model_bounds(2,2))   y = model_bounds(2,2)
if (z<model_bounds(3,1))   z = model_bounds(3,1)
if (z>model_bounds(3,2))   z = model_bounds(3,2)
do i_count = 1,21
   call polint3D (xseries, yseries, zseries, C_smooth(:,:,:,i_count), order, x, y, z, tmp_C(i_count))
enddo
call polint3D (xseries, yseries, zseries, rho_smooth(:,:,:), order, x, y, z, rho)

C(1,1:6) = tmp_C(1:6)
C(2,2:6) = tmp_C(7:11)
C(3,3:6) = tmp_C(12:15)
C(4,4:6) = tmp_C(16:18)
C(5,5:6) = tmp_C(19:20)
C(6,6) = tmp_C(21)
do i = 2,6
   do j = 1,i-1
      C(i,j) = C(j,i)
   enddo
enddo


end subroutine get_elastic_param

! #######################################################
subroutine deload

implicit none


deallocate (rho_smooth, C_smooth, xseries, yseries, zseries)


end subroutine deload


end module read_model
