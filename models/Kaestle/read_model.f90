module read_model !!! Polynomial interpolation of E. Kaestle's model (JGR 2018) !!!

implicit none

public :: get_value, get_value_aniso, deload
private :: load_model, sw_vel, latitude,longitude,depth, &
           first_time, order, n_lat,n_long,n_depth

integer, parameter :: order = 3   ! 2=linear, 3=quadratic, 4=cubic, etc...

integer, parameter :: n_lat = 110, n_long = 109, n_depth = 79
doubleprecision, dimension(1:n_lat) :: latitude = &
(/ 40.05, 40.15, 40.25, 40.35, 40.45, 40.55, 40.65, 40.75, 40.85, 40.95, &
   41.05, 41.15, 41.25, 41.35, 41.45, 41.55, 41.65, 41.75, 41.85, 41.95, &
   42.05, 42.15, 42.25, 42.35, 42.45, 42.55, 42.65, 42.75, 42.85, 42.95, &
   43.05, 43.15, 43.25, 43.35, 43.45, 43.55, 43.65, 43.75, 43.85, 43.95, &
   44.05, 44.15, 44.25, 44.35, 44.45, 44.55, 44.65, 44.75, 44.85, 44.95, &
   45.05, 45.15, 45.25, 45.35, 45.45, 45.55, 45.65, 45.75, 45.85, 45.95, &
   46.05, 46.15, 46.25, 46.35, 46.45, 46.55, 46.65, 46.75, 46.85, 46.95, &
   47.05, 47.15, 47.25, 47.35, 47.45, 47.55, 47.65, 47.75, 47.85, 47.95, &
   48.05, 48.15, 48.25, 48.35, 48.45, 48.55, 48.65, 48.75, 48.85, 48.95, &
   49.05, 49.15, 49.25, 49.35, 49.45, 49.55, 49.65, 49.75, 49.85, 49.95, &
   50.05, 50.15, 50.25, 50.35, 50.45, 50.55, 50.65, 50.75, 50.85, 50.95 /)
doubleprecision, dimension(1:n_long) :: longitude = &
(/ 4.057,  4.205,  4.352,  4.5,    4.648,  4.795,  4.943,  5.09,   5.238, &
   5.385,  5.533,  5.68,   5.828,  5.975,  6.123,  6.27,   6.418,  6.566, &
   6.713,  6.861,  7.008,  7.156,  7.303,  7.451,  7.598,  7.746,  7.893, &
   8.041,  8.189,  8.336,  8.484,  8.631,  8.779,  8.926,  9.074,  9.221, &
   9.369,  9.516,  9.664,  9.811,  9.959, 10.107, 10.254, 10.402, 10.549, &
  10.697, 10.844, 10.992, 11.139, 11.287, 11.434, 11.582, 11.73,  11.877, &
  12.025, 12.172, 12.32,  12.467, 12.615, 12.762, 12.91,  13.057, 13.205, &
  13.352, 13.5,   13.648, 13.795, 13.943, 14.09,  14.238, 14.385, 14.533, &
  14.68,  14.828, 14.975, 15.123, 15.27,  15.418, 15.566, 15.713, 15.861, &
  16.008, 16.156, 16.303, 16.451, 16.598, 16.746, 16.893, 17.041, 17.189, &
  17.336, 17.484, 17.631, 17.779, 17.926, 18.074, 18.221, 18.369, 18.516, &
  18.664, 18.811, 18.959, 19.107, 19.254, 19.402, 19.549, 19.697, 19.844, &
  19.992 /)
doubleprecision, dimension(1:n_depth) :: depth = &
(/ 0.00000000e+00,  1.00000000e-01,  2.00000000e-01,  3.00000000e-01, &
   4.00000000e-01,  5.00000000e-01,  6.00000000e-01,  1.00000000e+00, &
   1.50000000e+00,  2.00000000e+00,  2.50000000e+00,  3.00000000e+00, &
   4.00000000e+00,  5.00000000e+00,  6.00000000e+00,  7.00000000e+00, &
   8.00000000e+00,  9.00000000e+00,  1.00000000e+01,  1.20000000e+01, &
   1.40000000e+01,  1.60000000e+01,  1.80000000e+01,  2.00000000e+01, &
   2.20000000e+01,  2.40000000e+01,  2.60000000e+01,  2.80000000e+01, &
   3.00000000e+01,  3.20000000e+01,  3.40000000e+01,  3.60000000e+01, &
   3.80000000e+01,  4.00000000e+01,  4.20000000e+01,  4.40000000e+01, &
   4.60000000e+01,  4.80000000e+01,  5.00000000e+01,  5.20000000e+01, &
   5.40000000e+01,  5.60000000e+01,  5.80000000e+01,  6.00000000e+01, &
   6.20000000e+01,  6.40000000e+01,  6.60000000e+01,  6.80000000e+01, &
   7.00000000e+01,  7.20000000e+01,  7.40000000e+01,  7.60000000e+01, &
   7.80000000e+01,  8.00000000e+01,  8.20000000e+01,  8.50000000e+01, &
   9.00000000e+01,  9.50000000e+01,  1.00000000e+02,  1.05000000e+02, &
   1.10000000e+02,  1.15000000e+02,  1.20000000e+02,  1.25000000e+02, &
   1.30000000e+02,  1.35000000e+02,  1.40000000e+02,  1.45000000e+02, &
   1.50000000e+02,  1.55000000e+02,  1.60000000e+02,  1.65000000e+02, &
   1.70000000e+02,  1.75000000e+02,  1.80000000e+02,  1.85000000e+02, &
   1.90000000e+02,  1.95000000e+02,  2.00000000e+02 /)

doubleprecision, dimension(:,:,:), allocatable :: sw_vel
logical :: first_time = .true.

contains

! #######################################################
subroutine get_value (r,theta_rad,phi_rad,rho,vp,vs,Qmu,moho)

use angles
use module_polinterp

implicit none

integer, optional, intent(IN) :: moho
doubleprecision, intent(IN) :: r,theta_rad,phi_rad
doubleprecision, intent(INOUT) :: rho,vp,vs,Qmu

doubleprecision :: x,y,z

doubleprecision, parameter :: PS_ratio = 1.8, vs_min = 1800.d0


if (first_time) then
   first_time = .false.
   call load_model
endif

x = 90.d0 - rad2deg(theta_rad)
y = rad2deg(phi_rad)
z = (Rterre-r)/1000.d0

if (x<latitude(1))       x = latitude(1)
if (x>latitude(n_lat))   x = latitude(n_lat)
if (y<longitude(1))      y = longitude(1)
if (y>longitude(n_long)) y = longitude(n_long)
if (z<depth(1))          z = depth(1)
if (z>depth(n_depth))    z = depth(n_depth)
call polint3D (latitude, longitude, depth, sw_vel(:,:,:), order, x, y, z, vs)

vs = vs*1000.d0
vs = max(vs,vs_min)
vp = PS_ratio*vs
rho = 2670.d0
Qmu = 300.d0


end subroutine get_value

! #######################################################
subroutine get_value_aniso (r,theta_rad,phi_rad,rho,A,C,F,L,M,Gc,Gs,Hc,Hs,Bc,Bs,Ec,Es,Qmu,moho)

implicit none

integer, optional, intent(IN) :: moho
doubleprecision, intent(IN) :: r,theta_rad,phi_rad
doubleprecision, intent(INOUT) :: rho,A,C,F,L,M,Gc,Gs,Hc,Hs,Bc,Bs,Ec,Es,Qmu


stop "KAESTLE'S MODEL IS ISOTROPIC !!!"


end subroutine get_value_aniso

! #######################################################
subroutine load_model

implicit none

integer :: ios, i,j,k
real :: x,y,z


open (11, file="~/RegSEM/models/Kaestle/swmod_alps.txt", form="formatted", status="old", iostat=ios)
if (ios/=0)   stop "Cannot open swmod_alps.txt"

allocate (sw_vel(1:n_lat,1:n_long,1:n_depth))
do i = 1,n_lat
    do j = 1,n_long
        do k = 1,n_depth
            read(11,*) x, y, z, sw_vel(i,j,k)
        enddo
    enddo
enddo


end subroutine load_model

! #######################################################
subroutine deload

implicit none


deallocate (sw_vel)


end subroutine deload


end module read_model
