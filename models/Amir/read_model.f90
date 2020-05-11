module read_model


implicit none

public :: get_value, get_value_aniso

private

logical :: first_time = .true.
integer, parameter :: N_DEPTH_MAX = 10000, &
                      n_lat = 6, &   ! Number of different latitudes in the .dat file
                      n_param = 7   ! Number of elastic parameters in the .dat file
character(len=30) :: inputfile = "aus_model_smooth.dat"

integer :: n_profile, n_depth
integer, dimension(1:n_lat) :: n_long
doubleprecision, dimension(1:n_lat) :: lat
doubleprecision, dimension(:), allocatable :: depth
doubleprecision, dimension(:,:), allocatable :: long
doubleprecision, dimension(:,:,:), allocatable :: profile


contains

!------------------------------------------------------------------------------------
subroutine get_value (r,theta_rad,phi_rad,rho,vp,vs,Qmu,moho)


implicit none

integer, optional, intent(IN) :: moho
doubleprecision, intent(IN) :: r,theta_rad,phi_rad
doubleprecision, intent(INOUT) :: rho,vp,vs,Qmu

print *,"AMIR'S MODELS ARE ANISOTROPIC"
stop


end subroutine get_value

!------------------------------------------------------------------------------------
subroutine get_value_aniso (r,theta_rad,phi_rad,rho,A,C,F,L,M,Gc,Gs,Hc,Hs,Bc,Bs,Ec,Es,Qmu,moho)


use angles

implicit none

integer, optional, intent(IN) :: moho
doubleprecision, intent(IN) :: r,theta_rad,phi_rad
doubleprecision, intent(OUT) :: rho,A,C,F,L,M,Gc,Gs,Hc,Hs,Bc,Bs,Ec,Es,Qmu

integer :: i,k, ii,jj,kk, i_profile
integer, dimension(0:1) :: j
doubleprecision :: d, wi,wk
doubleprecision, dimension(0:1) :: wj
doubleprecision, dimension(1:n_param) :: val
doubleprecision, dimension(0:1,0:1,0:1) :: w


if (first_time) then
   call read_init_model
   do i = 1,n_lat   ! On travaille avec des colatitudes en radian ordonnees de maniere croissante
      lat(i) = 90.d0 - lat(i)
      lat(i) = deg2rad(lat(i))
      do jj = 1,n_long(i)   ! On travaille avec des longitudes positives en radian ordonnees de maniere croissante
         if (long(i,jj)<0.d0)   long(i,jj) = 360.d0 + long(i,jj)
         long(i,jj) = deg2rad(long(i,jj))
      enddo
   enddo
   first_time = .false.
endif

d = Rterre - r
if (d<=depth(1)) then
   k = 1;   wk = 1.d0
else if (d>=depth(n_depth)) then
   k = n_depth-1;   wk = 0.d0
else
   call dichotomy (n_depth, depth, d, 1, n_depth, k, wk)
endif

if (theta_rad<=lat(1)) then
   i = 1;   wi = 1.d0
else if (theta_rad>=lat(n_lat)) then
   i = n_lat-1;   wi = 0.d0
else
   call dichotomy (n_lat, lat, theta_rad, 1, n_lat, i, wi)
endif

do ii = 0,1   ! Et si on passe Greenwich ?!?
 if (phi_rad<=long(i+ii,1)) then
    j(ii) = 1;   wj(ii) = 1.d0
 else if (phi_rad>=long(i+ii,n_long(i+ii))) then
    j(ii) = n_long(i+ii)-1;   wj(ii) = 0.d0
 else
    call dichotomy (n_long(i+ii), long(i+ii,1:n_long(i+ii)), phi_rad, 1, n_long(i+ii), j(ii), wj(ii))
 endif
enddo

! Linear interpolation
w(0,0,0) = wi*wj(0)*wk
w(1,0,0) = (1-wi)*wj(1)*wk
w(0,1,0) = wi*(1-wj(0))*wk
w(1,1,0) = (1-wi)*(1-wj(1))*wk
w(0,0,1) = wi*wj(0)*(1-wk)
w(1,0,1) = (1-wi)*wj(1)*(1-wk)
w(0,1,1) = wi*(1-wj(0))*(1-wk)
w(1,1,1) = (1-wi)*(1-wj(1))*(1-wk)
val(:) = 0.d0
do kk = 0,1
 do ii = 0,1
  do jj = 0,1
     i_profile = j(ii) + jj
     if (i+ii>1)   i_profile = i_profile + sum(n_long(1:i+ii-1))
     val(:) = val(:) + w(ii,jj,kk)*profile(i_profile,k+kk,:)
  enddo
 enddo
enddo
rho=val(1); A=val(2); C=val(3); F=val(4); L=val(5); M=val(6); Qmu=val(7)
Gc=0.d0; Gs=0.d0; Hc=0.d0; Hs=0.d0; Bc=0.d0; Bs=0.d0; Ec=0.d0; Es=0.d0


end subroutine get_value_aniso

!------------------------------------------------------------------------------------
subroutine read_init_model
! format of .dat file:
!  lat  long  depth  rho    Vsv   Vsh   Vpv   Vph    eta      Qmu
! (deg) (deg) (km) (kg/m3) (m/s) (m/s) (m/s) (m/s) (no dim) (no dim)
! Longitude changes first, then latitude and then depth.
! Latitude and depth must be regular.
! Longitude can be irregular (i.e different for each latitude).
! Depth must increase.


implicit none

integer :: ios, i,j
doubleprecision :: current, trash
doubleprecision, parameter :: epsil = 0.000001
doubleprecision, dimension(1:n_param) :: tmp
doubleprecision, dimension(1:N_DEPTH_MAX) :: tmp2


open (11,file=trim(inputfile),form="formatted",status="old",iostat=ios)
if (ios>0)   stop 'UNABLE TO OPEN the inputfile !!!'

read (11,*) current
lat(1) = current
j = 0
do i = 2,n_lat+1
   do while (dabs(current-lat(i-1))<epsil)
      read (11,*) current
      j = j + 1
   enddo
   n_long(i-1) = j
   if (i<=n_lat)   lat(i) = current
   j = 0
enddo

rewind(11)
allocate (long(1:n_lat,maxval(n_long)))
do i = 1,n_lat
 do j = 1,n_long(i)
    read (11,*) trash, long(i,j)
 enddo
enddo
n_profile = sum(n_long)

rewind(11)
j = 1
do while (ios==0)
   do i = 1,n_profile
      read (11,fmt=*,iostat=ios) trash, trash, current
   enddo
   tmp2(j) = current
   j = j + 1
enddo
n_depth = j - 2
allocate (depth(1:n_depth))
depth(1:n_depth) = tmp2(1:n_depth)*1000.d0

rewind(11)
allocate (profile(1:n_profile,1:n_depth,1:n_param))
do i = 1,n_depth
 do j = 1,n_profile
    read (11,*) trash, trash, trash, tmp(1:n_param)
    profile(j,i,1) = tmp(1)   ! rho = rho
    profile(j,i,2) = tmp(1)*tmp(5)**2   ! A = rho*Vph**2
    profile(j,i,3) = tmp(1)*tmp(4)**2   ! C = rho*Vpv**2
    profile(j,i,5) = tmp(1)*tmp(2)**2   ! L = rho*Vsv**2
    profile(j,i,6) = tmp(1)*tmp(3)**2   ! N = rho*Vsh**2
    profile(j,i,4) = tmp(6)*(profile(j,i,2)-2.d0*profile(j,i,5))   ! F = eta*(A-2*L)
    profile(j,i,7) = tmp(7)   ! Qmu = Qmu
 enddo
enddo

close (11)


end subroutine

!------------------------------------------------------------------------------------
recursive subroutine dichotomy(n,array,val,debut,fin,iii,w)


implicit none

integer, intent(IN) :: n, debut,fin
doubleprecision, dimension(1:n), intent(IN) :: array
doubleprecision, intent(IN) :: val
integer, intent(OUT) :: iii
doubleprecision, intent(OUT) :: w

integer :: i, a,b


a = debut;   b = fin;   i = int(real(a+b)/2.)
if (val>array(i)) then
   a = i
else
   b = i
endif
if (a==b-1) then
   iii = a
   w = (array(b)-val)/(array(b)-array(a))
else
   call dichotomy (n, array, val, a, b, iii, w)
endif


end subroutine


end module read_model
