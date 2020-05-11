module read_model !!! CUB !!!

use angles

private

type :: nodes_modele
doubleprecision, dimension(:,:), pointer :: param
end type

type :: long_modele
type(nodes_modele), dimension(:), pointer :: longitude
end type

type(long_modele), dimension(:), allocatable :: latitude
doubleprecision, dimension(:), pointer :: mem_lat,mem_long
integer, parameter :: n_lat = 91, n_long = 180
logical :: first_time = .true.

public :: get_value, get_value_aniso, deload

contains

! #######################################################
subroutine get_value (r,theta_rad,phi_rad,rho,vp,vs,Qmu,moho)

implicit none

integer, optional, intent(IN) :: moho
doubleprecision, intent(IN) :: r,theta_rad,phi_rad
doubleprecision, intent(INOUT) :: rho,vp,vs,Qmu

stop 'CUB IS AN ANISOTROPIC MODEL !!!'

end subroutine get_value

! #######################################################
subroutine get_value_aniso (r,theta_rad,phi_rad,rho,A,C,F,L,M,Gc,Gs,Hc,Hs,Bc,Bs,Ec,Es,Qmu,moho)

implicit none

integer, optional, intent(IN) :: moho
doubleprecision, intent(IN) :: r,theta_rad,phi_rad
doubleprecision, intent(OUT) :: rho,A,C,F,L,M,Gc,Gs,Hc,Hs,Bc,Bs,Ec,Es,Qmu

integer :: i,j, ok
doubleprecision :: depth, Wcolat,Wlong,Wrad, theta,phi, vpv,vph,vsv,vsh,eta_aniso
doubleprecision, dimension(1:9) :: param
doubleprecision, dimension(1:2,1:9) :: pt_intermediaire
doubleprecision, dimension(1:8,1:9) :: pt
doubleprecision, parameter :: d = 6071.d0, mille = 1000.d0

if (first_time) then
    call load_model
    first_time = .false.
endif

if (r>d*1000.d0) then

    ok = 1

    theta = 180.d0*theta_rad/pi
    if (theta>180.d0)   theta = 179.999999
    if (theta<0.d0)   theta = 0.000001
    phi = 180.d0*phi_rad/pi

    depth = (Rterre-r)/1000.d0
    if (depth<=0)   depth = 0.001
    if (depth>=d)   depth = d - 0.001

    i = 0
    do while (mem_lat(i) - theta > 2.d0)
        i = i + 1
    enddo

    j = 0
    do while (phi - mem_long(j) > 2.d0)
        j = j + 1
    enddo

    if (moho==1) then
        call lookformohoup(i,j,pt(1,1:9),pt(2,1:9))
        if (j==n_long-1) then
            call lookformohoup(i,0,pt(3,1:9),pt(4,1:9))
        else
            call lookformohoup(i,j+1,pt(3,1:9),pt(4,1:9))
        endif
        call lookformohoup(i+1,j,pt(5,1:9),pt(6,1:9))
        if (j==n_long-1) then
            call lookformohoup(i+1,0,pt(7,1:9),pt(8,1:9))
        else
            call lookformohoup(i+1,j+1,pt(7,1:9),pt(8,1:9))
        endif

        Wcolat = abs(theta-mem_lat(i))/2.d0
        Wlong = abs(phi-mem_long(j))/2.d0
        pt_intermediaire(2,1:9) = Wcolat * (Wlong*pt(8,1:9) + (1-Wlong)*pt(6,1:9)) + &
                                  (1-Wcolat) * (Wlong*pt(4,1:9) + (1-Wlong)*pt(2,1:9))
        if (depth>=pt_intermediaire(2,1)) then
            param(1:9) = pt_intermediaire(2,1:9)
            ok = 0
        else if (depth>=pt(2,1) .or. depth>=pt(4,1) .or. depth>=pt(6,1) .or. depth>=pt(8,1)) then
            pt_intermediaire(1,1:9) = Wcolat * (Wlong*pt(7,1:9) + (1-Wlong)*pt(5,1:9)) + &
                                      (1-Wcolat) * (Wlong*pt(3,1:9) + (1-Wlong)*pt(1,1:9))
            Wrad = abs(depth-pt_intermediaire(1,1)) / abs(pt_intermediaire(2,1)-pt_intermediaire(1,1))
            param(1:9) = Wrad*pt_intermediaire(2,1:9) + (1-Wrad)*pt_intermediaire(1,1:9)
            ok = 0
        endif
    endif

    if (moho==-1) then
        call lookformohodown(i,j,pt(1,1:9),pt(2,1:9))
        if (j==n_long-1) then
            call lookformohodown(i,0,pt(3,1:9),pt(4,1:9))
        else
            call lookformohodown(i,j+1,pt(3,1:9),pt(4,1:9))
        endif
        call lookformohodown(i+1,j,pt(5,1:9),pt(6,1:9))
        if (j==n_long-1) then
            call lookformohodown(i+1,0,pt(7,1:9),pt(8,1:9))
        else
            call lookformohodown(i+1,j+1,pt(7,1:9),pt(8,1:9))
        endif

        Wcolat = abs(theta-mem_lat(i))/2.d0
        Wlong = abs(phi-mem_long(j))/2.d0
        pt_intermediaire(1,1:9) = Wcolat * (Wlong*pt(7,1:9) + (1-Wlong)*pt(5,1:9)) + &
                                  (1-Wcolat) * (Wlong*pt(3,1:9) + (1-Wlong)*pt(1,1:9))
        if (depth<=pt_intermediaire(1,1)) then
            param(1:9) = pt_intermediaire(1,1:9)
            ok = 0
        else if (depth<=pt(1,1) .or. depth<=pt(3,1) .or. depth<=pt(5,1) .or. depth<=pt(7,1)) then
            pt_intermediaire(2,1:9) = Wcolat * (Wlong*pt(8,1:9) + (1-Wlong)*pt(6,1:9)) + &
                                      (1-Wcolat) * (Wlong*pt(4,1:9) + (1-Wlong)*pt(2,1:9))
            Wrad = abs(depth-pt_intermediaire(1,1)) / abs(pt_intermediaire(2,1)-pt_intermediaire(1,1))
            param(1:9) = Wrad*pt_intermediaire(2,1:9) + (1-Wrad)*pt_intermediaire(1,1:9)
            ok = 0
        endif
    endif

    if (ok==1) then
        call lookfordepths(i,j,pt(1,1:9),pt(2,1:9),depth)
        if (j==n_long-1) then
            call lookfordepths(i,0,pt(3,1:9),pt(4,1:9),depth)
        else
            call lookfordepths(i,j+1,pt(3,1:9),pt(4,1:9),depth)
        endif
        call lookfordepths(i+1,j,pt(5,1:9),pt(6,1:9),depth)
        if (j==n_long-1) then
            call lookfordepths(i+1,0,pt(7,1:9),pt(8,1:9),depth)
        else
            call lookfordepths(i+1,j+1,pt(7,1:9),pt(8,1:9),depth)
        endif

        Wcolat = abs(theta-mem_lat(i))/2.d0
        Wlong = abs(phi-mem_long(j))/2.d0
        pt_intermediaire(1,1:9) = Wcolat * (Wlong*pt(7,1:9) + (1-Wlong)*pt(5,1:9)) + &
                                  (1-Wcolat) * (Wlong*pt(3,1:9) + (1-Wlong)*pt(1,1:9))
        pt_intermediaire(2,1:9) = Wcolat * (Wlong*pt(8,1:9) + (1-Wlong)*pt(6,1:9)) + &
                                  (1-Wcolat) * (Wlong*pt(4,1:9) + (1-Wlong)*pt(2,1:9))
        Wrad = abs(depth-pt_intermediaire(1,1)) / abs(pt_intermediaire(2,1)-pt_intermediaire(1,1))
        param(1:9) = Wrad*pt_intermediaire(2,1:9) + (1-Wrad)*pt_intermediaire(1,1:9)
    endif

    rho = param(6)
    vpv = param(2)
    vph = param(3)
    vsv = param(4)
    vsh = param(5)
    eta_aniso = param(7)
    Qmu = param(8)

else

    call prem_aniso(r/Rterre,rho,vpv,vph,vsv,vsh,eta_aniso,Qmu)

endif

rho = rho*mille;   vpv = vpv*mille;   vsv = vsv*mille;   vph = vph*mille;   vsh = vsh*mille
A = rho*vph**2
C = rho*vpv**2
L = rho*vsv**2
M = rho*vsh**2
F = eta_aniso*(A-2.d0*L)
Gc=0.d0; Gs=0.d0; Hc=0.d0; Hs=0.d0; Bc=0.d0; Bs=0.d0; Ec=0.d0; Es=0.d0

return
end subroutine get_value_aniso

! #######################################################
subroutine load_model

implicit none

integer :: ios, i,j,k, n, nb_nodes

allocate (latitude(0:n_lat-1))

open (11, file="/data2/paulcup/RegSEM_training/DATA/CUB/goodCUB.asc", form="formatted", status="old", iostat=ios)
if (ios/=0) then
    print *,"Cannot open goodCUB.asc"
    stop
endif

allocate (mem_lat(0:n_lat-1))
allocate (mem_long(0:n_long-1))

do i = 0,n_lat-1
    allocate (latitude(i)%longitude(0:n_long-1))
    do j = 0,n_long-1
        read(11,*) mem_lat(i), mem_long(j), nb_nodes
        mem_lat(i) = 90.d0 - mem_lat(i)
        if (mem_long(j)<0.d0)   mem_long(j) = 360.d0 + mem_long(j)
        allocate (latitude(i)%longitude(j)%param(1:nb_nodes,1:9))
        do k = 1,nb_nodes
            read(11,*) (latitude(i)%longitude(j)%param(k,n), n=1,9)
        enddo
    enddo
enddo

end subroutine load_model

! #######################################################
subroutine deload

implicit none

integer :: i,j

do i = 0,n_lat-1
    do j = 0,n_long-1
        deallocate (latitude(i)%longitude(j)%param)
    enddo
    deallocate (latitude(i)%longitude)
enddo
deallocate (latitude)

end subroutine deload

! #######################################################
subroutine lookfordepths (latnum,longnum,ptA,ptB,depth)

implicit none

integer, intent(IN) :: latnum,longnum
doubleprecision, intent(IN) :: depth
doubleprecision, dimension(1:9), intent(OUT) :: ptA, ptB

integer :: k

k = 1
do while (latitude(latnum)%longitude(longnum)%param(k,1) < depth)
    k = k + 1
enddo
ptA(:) = latitude(latnum)%longitude(longnum)%param(k-1,:)
ptB(:) = latitude(latnum)%longitude(longnum)%param(k,:)

return
end subroutine

! ########################################################
subroutine lookformohoup (latnum,longnum,ptA,ptB)

implicit none

integer, intent(IN) :: latnum,longnum
doubleprecision, dimension(1:9), intent(OUT) :: ptA, ptB

integer :: k
doubleprecision :: previous

previous = latitude(latnum)%longitude(longnum)%param(1,1)
k = 2
do while (abs(latitude(latnum)%longitude(longnum)%param(k,1)-previous)>0.001)
    previous = latitude(latnum)%longitude(longnum)%param(k,1)
    k = k + 1
enddo
ptA(:) = latitude(latnum)%longitude(longnum)%param(k-2,:)
ptB(:) = latitude(latnum)%longitude(longnum)%param(k-1,:)

return
end subroutine

! ########################################################
subroutine lookformohodown (latnum,longnum,ptA,ptB)

implicit none

integer, intent(IN) :: latnum,longnum
doubleprecision, dimension(1:9), intent(OUT) :: ptA, ptB

integer :: k
doubleprecision :: previous

previous = latitude(latnum)%longitude(longnum)%param(1,1)
k = 2
do while (abs(latitude(latnum)%longitude(longnum)%param(k,1)-previous)>0.001)
    previous = latitude(latnum)%longitude(longnum)%param(k,1)
    k = k + 1
enddo
ptA(:) = latitude(latnum)%longitude(longnum)%param(k,:)
ptB(:) = latitude(latnum)%longitude(longnum)%param(k+1,:)

return
end subroutine

! #######################################################
subroutine prem_aniso(x0,ro,vpv,vph,vsv,vsh,eta,Qmu)

doubleprecision                   :: x0,ro,vpv,vph,vsv,vsh,eta,Qmu

integer                           :: i
doubleprecision                   :: x
doubleprecision, parameter        :: Rt = 6371000.d0
doubleprecision, dimension(14)    :: r,q
doubleprecision, dimension(13,4)  :: d,pv,ph,sv,sh,et
logical                           :: pastrouve

r(1) = 0.d0;   r(2) = 1221.5d0;   r(3) = 3480.d0;   r(4) = 3630.d0
r(5) = 5600.d0;   r(6) = 5701.d0;   r(7) = 5771.d0;   r(8) = 5971.d0
r(9) = 6151.d0;   r(10) = 6291.d0;   r(11) = 6346.6d0
r(12) = 6356.d0;   r(13) = 6368.d0;   r(14) = 6371.d0

q(1) = 84.6d0;   q(2) = 1.d5;   q(3) = 355.d0;   q(4) = 355.d0;   q(5) = 355.d0;
q(6) = 165.d0;   q(7) = 165.d0;   q(8) = 165.d0;   q(9) = 70.d0;   q(10) = 191.d0;
q(11) = 300.d0;   q(12) = 300.d0;   q(13) = 300.d0

d(:,:) = 0.d0
d(1,1) = 13.0885d0;                     d(1,3) = -8.8381d0
d(2,1) = 12.5815d0; d(2,2) = -1.2638d0; d(2,3) = -3.6426d0; d(2,4) = -5.5281d0
d(3,1) = 7.9565d0 ; d(3,2) = -6.4761;   d(3,3) = 5.5283d0;  d(3,4) = -3.0807d0
d(4,1) = 7.9565d0 ; d(4,2) = -6.4761;   d(4,3) = 5.5283d0;  d(4,4) = -3.0807d0
d(5,1) = 7.9565d0 ; d(5,2) = -6.4761;   d(5,3) = 5.5283d0;  d(5,4) = -3.0807d0
d(6,1) = 5.3197d0 ; d(6,2) = -1.4836d0
d(7,1) = 11.2494d0; d(7,2) = -8.0298d0
d(8,1) = 7.1089d0 ; d(8,2) = -3.8045d0
d(9,1) = 2.6910d0 ; d(9,2) = 0.6924d0
d(10,1) = 2.6910d0; d(10,2) = 0.6924d0
d(11,1) = 2.9d0
d(12,1) = 2.6d0
d(13,1) = 2.6d0

pv(:,:) = 0.d0
pv(1,1) = 11.2622d0; pv(1,3) = -6.3640d0
pv(2,1) = 11.0487d0; pv(2,2) = -4.0362d0;  pv(2,3) = 4.8023d0;  pv(2,4) = -13.5732d0
pv(3,1) = 15.3891d0; pv(3,2) = -5.3181d0;  pv(3,3) = 5.5242d0;  pv(3,4) = -2.5514d0
pv(4,1) = 24.952d0 ; pv(4,2) = -40.4673d0; pv(4,3) = 51.4832d0; pv(4,4) = -26.6419d0
pv(5,1) = 29.2766d0; pv(5,2) = -23.6027d0; pv(5,3) = 5.5242d0;  pv(5,4) = -2.5514d0
pv(6,1) = 19.0957d0; pv(6,2) = -9.8672d0
pv(7,1) = 39.7027d0; pv(7,2) = -32.6166d0
pv(8,1) = 20.3926d0; pv(8,2) = -12.2569d0
pv(9,1) = 0.8317d0 ; pv(9,2) = 7.2180d0
pv(10,1) = 0.8317d0; pv(10,2) = 7.2180d0
pv(11,1) = 6.8d0
pv(12,1) = 5.8d0
pv(13,1) = pv(12,1)

ph(:,:) = pv(:,:)
ph(9,1) = 3.5908d0 ; ph(9,2) = 4.6172d0
ph(10,1) = 3.5908d0; ph(10,2) = 4.6172d0

sv(:,:) = 0.d0
sv(1,1) = 3.6678d0; sv(1,3) = -4.4475d0

sv(3,1) = 6.9254d0;  sv(3,2) = 1.4672d0;   sv(3,3) = -2.0834d0; sv(3,4) = 0.9783d0
sv(4,1) = 11.1671d0; sv(4,2) = -13.7818d0; sv(4,3) = 17.4575d0; sv(4,4) = -9.2777d0
sv(5,1) = 22.3459d0; sv(5,2) = -17.2473d0; sv(5,3) = -2.0834d0; sv(5,4) = 0.9783d0
sv(6,1) = 9.9839;    sv(6,2) = -4.9324
sv(7,1) = 22.3512d0; sv(7,2) = -18.5856d0
sv(8,1) = 8.9496d0;  sv(8,2) = -4.4597
sv(9,1) = 5.8582d0;  sv(9,2) = -1.4678d0
sv(10,1) = 5.8582d0; sv(10,2) = -1.4678d0
sv(11,1) = 3.9d0
sv(12,1) = 3.2d0
sv(13,1) = sv(12,1)

sh(:,:) = sv(:,:)
sh(9,1) = -1.0839d0; sh(9,2)  = 5.7176d0 
sh(10,1) = -1.0839d0; sh(10,2)  = 5.7176d0

et(:,:) = 0.0d0
et(:,1) = 1.d0
et(9,1) = 3.3687d0 ; et(9,2) = -2.4778d0
et(10,1) = 3.3687d0; et(10,2) = -2.4778d0

r(:) = r(:) * 1000.d0
x    = x0   * Rt

if (x >= r(14)) then
   x0 = 1.d0
   ro = d(13,1) + x0*( d(13,2) + x0*( d(13,3) + x0*( d(13,4) )))
   vpv = pv(13,1) + x0*( pv(13,2) + x0*( pv(13,3) + x0*( pv(13,4) )))
   vsv = sv(13,1) + x0*( sv(13,2) + x0*( sv(13,3) + x0*( sv(13,4) )))
   vph = ph(13,1) + x0*( ph(13,2) + x0*( ph(13,3) + x0*( ph(13,4) )))
   vsh = sh(13,1) + x0*( sh(13,2) + x0*( sh(13,3) + x0*( sh(13,4) )))
   eta = et(13,1) + x0*( et(13,2) + x0*( et(13,3) + x0*( et(13,4) )))
   Qmu = q(13)
   pastrouve = .false.
else
   pastrouve = .true.
endif

do i = 1,13
 if (pastrouve .and. (x >= r(i)) .and. (x < r(i+1))) then
  ro  = d (i,1) + x0*( d (i,2) + x0*( d (i,3) + x0*( d (i,4) )))
  vpv = pv(i,1) + x0*( pv(i,2) + x0*( pv(i,3) + x0*( pv(i,4) )))
  vsv = sv(i,1) + x0*( sv(i,2) + x0*( sv(i,3) + x0*( sv(i,4) )))
  vph = ph(i,1) + x0*( ph(i,2) + x0*( ph(i,3) + x0*( ph(i,4) )))
  vsh = sh(i,1) + x0*( sh(i,2) + x0*( sh(i,3) + x0*( sh(i,4) )))
  eta = et(i,1) + x0*( et(i,2) + x0*( et(i,3) + x0*( et(i,4) )))
  Qmu = q(i)
  pastrouve = .false.
 endif
enddo

end subroutine prem_aniso

! ########################################################
end module read_model
