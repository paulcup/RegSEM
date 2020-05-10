subroutine define_arrays (Tdomain, rg)


use sdomains
use angles
use tensor_util
use read_model
use module_A3d

implicit none

include 'mpif.h'

type (domain), intent (INOUT), target :: Tdomain
integer, intent(IN) :: rg

integer :: n, mat, cov, ngllx,nglly,ngllz, ngll1,ngll2, ngll, ngllPML, ngllocean, i,j,k, n_elem, nf,ne,nv, &
           idef, code, shift, I_give_to, I_take_from, n_rings, ii,jj, meshtype, nb_interfaces
integer, parameter :: etiquette=100
integer, dimension(mpi_status_size) :: statut
doubleprecision :: vp,vs,rho,Qmu, dx, x,y,z, r,theta,phi, A,C,L,M,F, Gc,Gs,Hc,Hs,Bc,Bs,Ec,Es, &
                   lambda,mu, xa,ya,za, xi,eta, Jac_surf, det, dt, lat,long, epsil, &
                   rho_anom,vp_anom,vs_anom,Qmu_anom, A_anom,C_anom,F_anom,L_anom,M_anom, &
                   Gc_anom,Gs_anom,Hc_anom,Hs_anom,Bc_anom,Bs_anom,Ec_anom,Es_anom
doubleprecision, parameter :: rhowater=1000.d0
doubleprecision, dimension (1:3) :: barycentre, tmp
doubleprecision, dimension (0:1,0:1) :: LocInvGrad_surf
doubleprecision, dimension (1:3,1:2) :: boundary
doubleprecision, dimension (0:2,0:2) :: Rot
doubleprecision, dimension (1:3,1:3) :: Rot2
doubleprecision, dimension (1:6,1:6) :: Cij
doubleprecision, dimension (1:8,1:3) :: vertices
doubleprecision, dimension (:), allocatable :: timestep, rad
doubleprecision, dimension (:,:), allocatable :: coord
doubleprecision, dimension (:,:,:), allocatable :: xix,xiy,xiz, etax,etay,etaz, zetax,zetay,zetaz, Jac, Rsph, &
                                                   Rlam,Rmu,RKmod, Whei, LocMassMat, wx,wy,wz, temp_PMLx,temp_PMLy
character*60 :: modelfile
character*20, parameter :: myfmt = "(4f10.3)"


Tdomain%save_model = .true.


!!! Initialisation !!!

meshtype = Tdomain%mesh
if (meshtype==1 .or. meshtype==2 .or. meshtype==3 .or. meshtype==6) then
   allocate (rad(0:10))
   rad(0)=1221500; rad(1)=3480000; rad(2)=3630000; rad(3)=5600000
   rad(4)=5701000; rad(5)=5771000; rad(6)=5971000; rad(7)=6151000
   rad(8)=6291000; rad(9)=6346600; rad(10)=6356000
   if (meshtype==1)   nb_interfaces = 11
   if (meshtype==2)   nb_interfaces = 8
   if (meshtype==3)   nb_interfaces = 7
   if (meshtype==6)   nb_interfaces = 10
else if (meshtype==5) then
   allocate (rad(0:5))
   rad(0)=1221500; rad(1)=3480000; rad(2)=3630000; rad(3)=5701000
   rad(4)=5971000; rad(5)=6311000
   nb_interfaces = 6
else if (meshtype==4) then   ! To be set for any user-defined 1D model (model=4)
   allocate (rad(0:1))
!   rad(0)=3480000; rad(1)=6341000 ! the CMB + a 30km-deep interface
   rad(0)=3480000; rad(1)=3490000
   nb_interfaces = 2
else if (meshtype==0) then   ! To be set for any cuboid chunk (curve=F)
   allocate (rad(0:1))
!   rad(0)=-1000; rad(1)=17000 ! the bottom + an interface at z = 17km
   rad(0)=-1000; rad(1)=-500
   nb_interfaces = 2
endif
epsil = 1.d0   ! Precision (in meter) when matching the elastic model with the mesh at the interfaces

if (Tdomain%curve) then
   Rot = transpose(Tdomain%rot)
   Rot2(1,1:3) = Rot(0,0:2);   Rot2(2,1:3) = Rot(1,0:2);   Rot2(3,1:3) = Rot(2,0:2)
   Rot = Tdomain%rot
endif

! Determination des coordonnees extremes du sous-domaine
boundary(:,1) = huge(x);   boundary(:,2) = -huge(x)
do n = 0,Tdomain%n_elem-1
   ngllx = Tdomain%specel(n)%ngllx
   nglly = Tdomain%specel(n)%nglly
   ngllz = Tdomain%specel(n)%ngllz
   vertices(1,1:3) = Tdomain%Globcoord(0:2,Tdomain%specel(n)%Iglobnum(0,0,0))
   vertices(2,1:3) = Tdomain%Globcoord(0:2,Tdomain%specel(n)%Iglobnum(ngllx-1,0,0))
   vertices(3,1:3) = Tdomain%Globcoord(0:2,Tdomain%specel(n)%Iglobnum(0,nglly-1,0))
   vertices(4,1:3) = Tdomain%Globcoord(0:2,Tdomain%specel(n)%Iglobnum(ngllx-1,nglly-1,0))
   vertices(5,1:3) = Tdomain%Globcoord(0:2,Tdomain%specel(n)%Iglobnum(0,0,ngllz-1))
   vertices(6,1:3) = Tdomain%Globcoord(0:2,Tdomain%specel(n)%Iglobnum(ngllx-1,0,ngllz-1))
   vertices(7,1:3) = Tdomain%Globcoord(0:2,Tdomain%specel(n)%Iglobnum(0,nglly-1,ngllz-1))
   vertices(8,1:3) = Tdomain%Globcoord(0:2,Tdomain%specel(n)%Iglobnum(ngllx-1,nglly-1,ngllz-1))
   tmp(:) = minval(vertices,dim=1)
   do i = 1,3
      boundary(i,1) = min(boundary(i,1),tmp(i))
   enddo
   tmp(:) = maxval(vertices,dim=1)
   do i = 1,3
      boundary(i,2) = max(boundary(i,2),tmp(i))
   enddo
enddo


!!! Introducing elastic properties !!!

do n = 0,Tdomain%n_elem-1
   
   mat = Tdomain%specel(n)%mat_index
   cov = Tdomain%specel(n)%cov_index
   ngllx = Tdomain%specel(n)%ngllx
   nglly = Tdomain%specel(n)%nglly
   ngllz = Tdomain%specel(n)%ngllz

   if (Tdomain%ellipticity) then
       allocate (Rsph(0:ngllx-1,0:nglly-1,0:ngllz-1))
       call r_spheric (ngllx,nglly,ngllz,Tdomain,mat,n,Rsph)
   endif
   if (Tdomain%specel(n)%ocean) then
       allocate (Tdomain%specel(n)%hocean(0:ngllx-1,0:nglly-1))
       call find_verticale (Tdomain,n)
   endif

   ! Calcul du barycentre
!   vertices(1,1:3) = Tdomain%Globcoord(0:2,Tdomain%specel(n)%Iglobnum(0,0,0))
!   vertices(2,1:3) = Tdomain%Globcoord(0:2,Tdomain%specel(n)%Iglobnum(ngllx-1,0,0))
!   vertices(3,1:3) = Tdomain%Globcoord(0:2,Tdomain%specel(n)%Iglobnum(0,nglly-1,0))
!   vertices(4,1:3) = Tdomain%Globcoord(0:2,Tdomain%specel(n)%Iglobnum(ngllx-1,nglly-1,0))
!   vertices(5,1:3) = Tdomain%Globcoord(0:2,Tdomain%specel(n)%Iglobnum(0,0,ngllz-1))
!   vertices(6,1:3) = Tdomain%Globcoord(0:2,Tdomain%specel(n)%Iglobnum(ngllx-1,0,ngllz-1))
!   vertices(7,1:3) = Tdomain%Globcoord(0:2,Tdomain%specel(n)%Iglobnum(0,nglly-1,ngllz-1))
!   vertices(8,1:3) = Tdomain%Globcoord(0:2,Tdomain%specel(n)%Iglobnum(ngllx-1,nglly-1,ngllz-1))
!   call barycentre3d (8,vertices,barycentre)

   do k = 0,ngllz-1
    do j = 0,nglly-1
     do i = 0,ngllx-1

        ! Taking the cartesian coordinates of the GLL
        idef = Tdomain%specel(n)%Iglobnum(i,j,k)
        x = Tdomain%Globcoord(0,idef)
        y = Tdomain%Globcoord(1,idef)
        z = Tdomain%Globcoord(2,idef)
        if (Tdomain%curve) then
           ! Coordinates in the real chunk
           xa = x;   ya = y;   za = z
           x = Rot(0,0)*xa + Rot(0,1)*ya + Rot(0,2)*za
           y = Rot(1,0)*xa + Rot(1,1)*ya + Rot(1,2)*za
           z = Rot(2,0)*xa + Rot(2,1)*ya + Rot(2,2)*za
           ! Convert the cartesian coordinates into spherical coordinates
           call cart2sph (x,y,z,r,theta,phi)
           if (Tdomain%ellipticity)   r = Rsph(i,j,k)
           if (Tdomain%specel(n)%ocean .and. k==ngllz-1) then
              Tdomain%specel(n)%hocean(i,j) = Rterre - r
              if (Tdomain%specel(n)%hocean(i,j)<0.d0)   Tdomain%specel(n)%hocean(i,j) = 0.d0
           endif
           ! Matching the elastic model with the mesh at the spherical interfaces
           if (k==0 .or. k==ngllz-1) then
              interface : do idef = 0,nb_interfaces-1
                 dx = dabs(r-rad(idef))
                 if (dx<epsil) then
                    r = rad(idef)
                    if (k==0)   r = r + epsil
                    if (k==ngllz-1)   r = r - epsil
                    exit interface
                 endif
              enddo interface
           endif
        else
           ! Matching the elastic model with the mesh at the interfaces
           if (k==0 .or. k==ngllz-1) then
              interface_plane : do idef = 0,nb_interfaces-1
                 dx = dabs(z-rad(idef))
                 if (dx<epsil) then
                    z = rad(idef)
                    if (k==0)   z = z + epsil
                    if (k==ngllz-1)   z = z - epsil
                    exit interface_plane
                 endif
              enddo interface_plane
           endif
        endif

        ! Getting the elastic properties
        if (Tdomain%aniso) then
           if (Tdomain%curve) then
              call get_value_aniso (r,theta,phi,rho,A,C,F,L,M,Gc,Gs,Hc,Hs,Bc,Bs,Ec,Es,Qmu, &
                                    Tdomain%specel(n)%moho_position)
!              if (n==Tdomain%n_elem-1 .and. i==ngllx-1 .and. j==nglly-1 .and. k==ngllz-1) then
!                 call deload   ! Desallocation (pour CUB)
!              endif
           else
              call get_value_aniso (x,y,z,rho,A,C,F,L,M,Gc,Gs,Hc,Hs,Bc,Bs,Ec,Es,Qmu, &
                                    Tdomain%specel(n)%moho_position)
           endif
           if (Tdomain%specel(n)%moho_position == 10) then
!              call get_region_param (barycentre/1000.d0,i+j+k,rho,Cij)   ! Uncomment L125-133
!              call get_homog_param(rho,Cij)
!              call get_elastic_param (x,y,z,rho,Cij,boundary)
!              if (n==Tdomain%n_elem-1 .and. i==ngllx-1 .and. j==nglly-1 .and. k==ngllz-1) then
!                 call deload   ! Desallocation (pour Upscaled/polynomials et Smoothed)
!              endif
           else
              if (Tdomain%curve) then
                 Cij(:,:) = 0.d0
                 Cij(1,1) = C
                 Cij(2,2) = A+Bc+Ec
                 Cij(3,3) = A-Bc+Ec
                 Cij(1,2) = F+Hc
                 Cij(1,3) = F-Hc
                 Cij(2,3) = A-2.d0*M-Ec
                 Cij(1,4) = s2*Hs
                 Cij(2,4) = s2*(Bs/2.d0+Es)
                 Cij(3,4) = s2*(Bs/2.d0-Es)
                 Cij(4,4) = 2.d0*(M-Ec)
                 Cij(5,5) = 2.d0*(L-Gc)
                 Cij(6,6) = 2.d0*(L+Gc)
                 Cij(5,6) = 2.d0*Gs
              else
                 Cij(:,:) = 0.d0
                 Cij(1,1) = A+Bc+Ec
                 Cij(2,2) = A-Bc+Ec
                 Cij(3,3) = C
                 Cij(1,2) = A-2.d0*M-Ec
                 Cij(1,3) = F+Hc
                 Cij(2,3) = F-Hc
                 Cij(1,6) = s2*(Bs/2.d0+Es)
                 Cij(2,6) = s2*(Bs/2.d0-Es)
                 Cij(3,6) = s2*Hs
                 Cij(4,4) = 2.d0*(L-Gc)
                 Cij(5,5) = 2.d0*(L+Gc)
                 Cij(6,6) = 2.d0*(M-Ec)
                 Cij(4,5) = 2.d0*Gs
              endif
              do ii = 2,6
                 do jj = 1,ii-1
                    Cij(ii,jj) = Cij(jj,ii)
                 enddo
              enddo
              if (Tdomain%adjoint .and. (.not. Tdomain%specel(n)%PML)) then
                 Tdomain%specel(n)%save_TIparam(1,i,j,k) = A
                 Tdomain%specel(n)%save_TIparam(2,i,j,k) = C
                 Tdomain%specel(n)%save_TIparam(3,i,j,k) = F
                 Tdomain%specel(n)%save_TIparam(4,i,j,k) = L
                 Tdomain%specel(n)%save_TIparam(5,i,j,k) = M
!                 call get_value_aniso_anom (r,theta,phi,rho_anom,A_anom,C_anom,F_anom,L_anom,M_anom, &
!                                            Gc_anom,Gs_anom,Hc_anom,Hs_anom,Bc_anom,Bs_anom,Ec_anom,Es_anom, &
!                                            Qmu_anom, Tdomain%specel(n)%moho_position)
!                 Tdomain%specel(n)%anomaly(i,j,k) = M_anom/L_anom
              endif
           endif
           ! Si l'elem est PML on considere un milieu isotrope moyen
           if (Tdomain%specel(n)%PML) then
              Tdomain%specel(n)%Lambda(i,j,k) = lambda_from_Cij(Cij)
              Tdomain%specel(n)%Mu(i,j,k) = mu_from_Cij(Cij)
           else ! Et sinon:
              ! Expression de Cij en cartesien
              if (Tdomain%curve)   call c_4tensor(Cij,theta,phi)
              ! Expression de Cij dans le chunk de reference
              if (Tdomain%curve)   call rot_4tensor(Cij,Rot2)
              ! Sauvegarde des 21 coeffs
              idef = 0
              do ii = 1,6
                 do jj = ii,6
                    Tdomain%specel(n)%Cij(idef,i,j,k) = Cij(ii,jj)
                    idef = idef + 1
                 enddo
              enddo
              if (Tdomain%n_sls>0) then
                 lambda = lambda_from_Cij(Cij)
                 mu = mu_from_Cij(Cij)
                 Tdomain%specel(n)%Cij( 0,i,j,k) = Tdomain%specel(n)%Cij( 0,i,j,k) - lambda-2.d0*mu
                 Tdomain%specel(n)%Cij( 6,i,j,k) = Tdomain%specel(n)%Cij( 6,i,j,k) - lambda-2.d0*mu
                 Tdomain%specel(n)%Cij(11,i,j,k) = Tdomain%specel(n)%Cij(11,i,j,k) - lambda-2.d0*mu
                 Tdomain%specel(n)%Cij( 1,i,j,k) = Tdomain%specel(n)%Cij( 1,i,j,k) - lambda
                 Tdomain%specel(n)%Cij( 2,i,j,k) = Tdomain%specel(n)%Cij( 2,i,j,k) - lambda
                 Tdomain%specel(n)%Cij( 7,i,j,k) = Tdomain%specel(n)%Cij( 7,i,j,k) - lambda
                 Tdomain%specel(n)%Cij(15,i,j,k) = Tdomain%specel(n)%Cij(15,i,j,k) - 2.d0*mu
                 Tdomain%specel(n)%Cij(18,i,j,k) = Tdomain%specel(n)%Cij(18,i,j,k) - 2.d0*mu
                 Tdomain%specel(n)%Cij(20,i,j,k) = Tdomain%specel(n)%Cij(20,i,j,k) - 2.d0*mu
                 Tdomain%specel(n)%Lambda(i,j,k) = lambda
                 Tdomain%specel(n)%Mu(i,j,k) = mu
              endif
           endif
        else
           if (Tdomain%curve) then
              call get_value (r,theta,phi,rho,vp,vs,Qmu,Tdomain%specel(n)%moho_position)
!              if (n==Tdomain%n_elem-1 .and. i==ngllx-1 .and. j==nglly-1 .and. k==ngllz-1) then
!                 call deload   ! Desallocation (pour CUB)
!              endif
           else
              call get_value (x,y,z,rho,vp,vs,Qmu,Tdomain%specel(n)%moho_position)
!              if (n==Tdomain%n_elem-1 .and. i==ngllx-1 .and. j==nglly-1 .and. k==ngllz-1) then
!                 call deload   ! Desallocation (pour Smoothed)
!              endif
           endif
           if (Tdomain%adjoint .and. (.not. Tdomain%specel(n)%PML)) then
!              call get_value_anom (r,theta,phi,rho_anom,vp_anom,vs_anom,Qmu_anom, &
!                                   Tdomain%specel(n)%moho_position)
!              Tdomain%specel(n)%anomaly(i,j,k) = vs_anom
           endif
           Tdomain%specel(n)%Lambda(i,j,k) = (vp**2 - 2.d0*vs**2) * rho
           Tdomain%specel(n)%Mu(i,j,k) = vs**2 * rho
        endif
        Tdomain%specel(n)%Density(i,j,k) = rho
        if ((.not.Tdomain%specel(n)%PML) .and. Tdomain%n_sls>0)   Tdomain%specel(n)%Q(i,j,k) = Qmu
     enddo
    enddo
   enddo

   if (Tdomain%ellipticity)   deallocate (Rsph)

enddo

if (Tdomain%ellipticity)   deallocate(Tdomain%Rsph_Nodes)

if (Tdomain%save_model) then
   call save_model (Tdomain,rg)
   call mpi_barrier(mpi_comm_world,code)
   if (rg==0) then
      call merge_model_files (Tdomain%n_proc, 'Vp')
      call merge_model_files (Tdomain%n_proc, 'Vs')
   endif
endif


!!! Calculating the time step following the Courant criteria !!!

call compute_dt (Tdomain,dt)
allocate (timestep(0:Tdomain%n_proc-1))
call mpi_allgather(dt,1,mpi_double_precision,timestep,1,mpi_double_precision,mpi_comm_world,code)
do i = 0,Tdomain%n_proc-1
   if (timestep(i) <= dt)   dt = timestep(i)
enddo
deallocate (timestep)
Tdomain%sTimeParam%ntime = int(Tdomain%sTimeParam%Duration/dt) + 1
Tdomain%sTimeParam%dt = Tdomain%sTimeParam%Duration/(Tdomain%sTimeParam%ntime-1)
print *, "THE NUMBER OF ITERATION IS", Tdomain%sTimeParam%ntime
print *, "THE TIME STEP IS", Tdomain%sTimeParam%dt


!!! Computing the source signal and allocate the traces !!!

call def_timefunc (Tdomain,rg)

if (Tdomain%save_trace) then
   do n = 0,Tdomain%n_receivers-1
      if (rg == Tdomain%sReceiver(n)%proc) then
         allocate (Tdomain%sReceiver(n)%StoreTrace (0:Tdomain%sTimeParam%ntime-1, 0:2))
         Tdomain%sReceiver(n)%StoreTrace = 0.d0
      endif
   enddo
endif


!!! Computing the Mass Matrix !!!

do n = 0,Tdomain%n_elem-1

   mat = Tdomain%specel(n)%mat_index
   ngllx = Tdomain%specel(n)%ngllx
   nglly = Tdomain%specel(n)%nglly
   ngllz = Tdomain%specel(n)%ngllz

   allocate (Jac (0:ngllx-1,0:nglly-1,0:ngllz-1))
   allocate (xix (0:ngllx-1,0:nglly-1,0:ngllz-1)) 
   allocate (xiy (0:ngllx-1,0:nglly-1,0:ngllz-1))
   allocate (xiz (0:ngllx-1,0:nglly-1,0:ngllz-1))    
  
   allocate (etax (0:ngllx-1,0:nglly-1,0:ngllz-1))
   allocate (etay (0:ngllx-1,0:nglly-1,0:ngllz-1))
   allocate (etaz (0:ngllx-1,0:nglly-1,0:ngllz-1))

   allocate (zetax (0:ngllx-1,0:nglly-1,0:ngllz-1))
   allocate (zetay (0:ngllx-1,0:nglly-1,0:ngllz-1))
   allocate (zetaz (0:ngllx-1,0:nglly-1,0:ngllz-1))

   allocate (Whei (0:ngllx-1,0:nglly-1,0:ngllz-1))

   allocate (RKmod (0:ngllx-1,0:nglly-1,0:ngllz-1))
   allocate (Rlam (0:ngllx-1,0:nglly-1,0:ngllz-1))
   allocate (Rmu (0:ngllx-1,0:nglly-1,0:ngllz-1))

   do k = 0,ngllz -1 
       do j = 0,nglly-1 
           do i = 0,ngllx-1
               Whei (i,j,k) = Tdomain%sSubdomain(mat)%GLLwx(i) * Tdomain%sSubdomain(mat)%GLLwy(j) &
                              * Tdomain%sSubdomain(mat)%GLLwz(k)
               if (Tdomain%specel(n)%PML .eqv. .false.) then
                  Tdomain%specel(n)%wgtx(i) = Tdomain%sSubdomain(mat)%GLLwx(i)
                  Tdomain%specel(n)%wgty(j) = Tdomain%sSubdomain(mat)%GLLwy(j)
                  Tdomain%specel(n)%wgtz(k) = Tdomain%sSubdomain(mat)%GLLwz(k)
               endif
           enddo
       enddo
   enddo
      
   xix = Tdomain%specel(n)%InvGrad(:,:,:,0,0)
   xiy = Tdomain%specel(n)%InvGrad(:,:,:,1,0)
   xiz = Tdomain%specel(n)%InvGrad(:,:,:,2,0)
      
   etax = Tdomain%specel(n)%InvGrad(:,:,:,0,1)
   etay = Tdomain%specel(n)%InvGrad(:,:,:,1,1)
   etaz = Tdomain%specel(n)%InvGrad(:,:,:,2,1)

   zetax = Tdomain%specel(n)%InvGrad(:,:,:,0,2)
   zetay = Tdomain%specel(n)%InvGrad(:,:,:,1,2)
   zetaz = Tdomain%specel(n)%InvGrad(:,:,:,2,2)

   Jac  = Tdomain%specel(n)%Jacob

   Tdomain%specel(n)%MassMat = Whei*Tdomain%specel(n)%Density*Jac

   if (Tdomain%specel(n)%ocean) then
       if (Tdomain%n_nodes/=27) then
           print *,"THE OCEAN OPTION REQUIRES 27 CTRL_PTS FOR NOW."
           stop
       endif
       allocate (Tdomain%specel(n)%Mocean(0:ngllx-1,0:nglly-1))
       allocate (coord(0:Tdomain%n_nodes-1,0:2))
       do i = 0,Tdomain%n_nodes-1
           j = Tdomain%specel(n)%Control_Nodes(i)
           coord(i,0:2) = Tdomain%Coord_Nodes(0:2,j)
       enddo
       do j = 0,nglly-1
        eta = Tdomain%sSubdomain(mat)%GLLcy(j)
        do i = 0,ngllx-1
            xi = Tdomain%sSubdomain(mat)%GLLcx(i)
            LocInvGrad_surf(:,:) = 0.d0
            do k = 4,7
                do ii = 0,1
                 do jj = 0,1
                     LocInvGrad_surf(ii,jj) = LocInvGrad_surf(ii,jj) + &
                                              Comp_derivshapefunc(k,xi,eta,1.d0,ii) * coord(k,jj)
                 enddo
                enddo
            enddo
            do k = 16,19
                do ii = 0,1
                 do jj = 0,1
                     LocInvGrad_surf(ii,jj) = LocInvGrad_surf(ii,jj) + &
                                              Comp_derivshapefunc(k,xi,eta,1.d0,ii) * coord(k,jj)
                 enddo
                enddo
            enddo
            do ii = 0,1
             do jj = 0,1
                 LocInvGrad_surf(ii,jj) = LocInvGrad_surf(ii,jj) + &
                                          Comp_derivshapefunc(25,xi,eta,1.d0,ii) * coord(25,jj)
             enddo
            enddo
            Jac_surf = LocInvgrad_surf(0,0) * LocInvgrad_surf(1,1) &
                       - LocInvgrad_surf(1,0) * LocInvgrad_surf(0,1)
            Tdomain%specel(n)%Mocean(i,j) = rhowater * Tdomain%specel(n)%hocean(i,j) * &
                                            Tdomain%sSubdomain(mat)%GLLwx(i) * &
                                            Tdomain%sSubdomain(mat)%GLLwy(j) * &
                                            Jac_surf
        enddo
       enddo
       deallocate (coord,Tdomain%specel(n)%hocean)
   endif
      
   if (.not.Tdomain%specel(n)%PML) then

!     Rlam = Tdomain%specel(n)%Lambda
!     Rmu  = Tdomain%specel(n)%Mu
!     RKmod = Rlam + 2.* Rmu
!
!     Tdomain%specel(n)%Acoeff(:,:,:,0) = -Whei*(RKmod*xix**2+Rmu*(xiy**2+xiz**2))*Jac
!     Tdomain%specel(n)%Acoeff(:,:,:,1) = -Whei*(RKmod*xix*etax+Rmu*(xiy*etay+xiz*etaz))*Jac
!     Tdomain%specel(n)%Acoeff(:,:,:,2) = -Whei*(rKmod*xix*zetax+Rmu*(xiy*zetay+xiz*zetaz))*Jac 
!     Tdomain%specel(n)%Acoeff(:,:,:,3) = -Whei*(Rlam+Rmu)*xix*xiy*Jac
!     Tdomain%specel(n)%Acoeff(:,:,:,4) = -Whei*(Rlam*xix*etay+Rmu*xiy*etax)*Jac  
!     Tdomain%specel(n)%Acoeff(:,:,:,5) = -Whei*(Rlam*xix*zetay+Rmu*xiy*zetax)*Jac
!     Tdomain%specel(n)%Acoeff(:,:,:,6) = -Whei*(Rlam+Rmu)*xix*xiz*Jac
!     Tdomain%specel(n)%Acoeff(:,:,:,7) = -Whei*(Rlam*xix*etaz+Rmu*xiz*etax)*Jac
!     Tdomain%specel(n)%Acoeff(:,:,:,8) = -Whei*(Rlam*xix*zetaz+rmu*xiz*zetax)*Jac
!
!     Tdomain%specel(n)%Acoeff(:,:,:,9) = -Whei*(RKmod*etax**2+Rmu* (etay**2+etaz**2))*Jac
!     Tdomain%specel(n)%Acoeff(:,:,:,10) = -Whei*(RKmod*etax*zetax+Rmu* (etay*zetay+etaz*zetaz))*Jac
!     Tdomain%specel(n)%Acoeff(:,:,:,11) = -Whei*(Rlam*etax*xiy+Rmu*etay*xix)*Jac
!     Tdomain%specel(n)%Acoeff(:,:,:,12) = -Whei*(Rlam+Rmu)*etay*etax*Jac
!     Tdomain%specel(n)%Acoeff(:,:,:,13) = -Whei*(Rlam*etax*zetay+Rmu*etay*zetax)*Jac
!     Tdomain%specel(n)%Acoeff(:,:,:,14) = -Whei*(Rlam*etax*xiz+Rmu*etaz*xix)*Jac
!     Tdomain%specel(n)%Acoeff(:,:,:,15) = -Whei*(Rlam+Rmu)*etaz*etax*Jac
!     Tdomain%specel(n)%Acoeff(:,:,:,16) = -Whei*(Rlam*etax*zetaz+Rmu*etaz*zetax)*Jac
!
!     Tdomain%specel(n)%Acoeff(:,:,:,17) = -Whei*(RKmod*zetax**2+Rmu*  (zetay**2+zetaz**2))*Jac
!     Tdomain%specel(n)%Acoeff(:,:,:,18) = -Whei*(Rlam*zetax*xiy+Rmu*zetay*xix)*Jac
!     Tdomain%specel(n)%Acoeff(:,:,:,19) = -Whei*(Rlam*zetax*etay+Rmu*zetay*etax)*Jac
!     Tdomain%specel(n)%Acoeff(:,:,:,20) = -Whei*(Rlam+Rmu)*zetax*zetay*Jac
!     Tdomain%specel(n)%Acoeff(:,:,:,21) = -Whei*(Rlam*zetax*xiz+Rmu*zetaz*xix)*Jac
!     Tdomain%specel(n)%Acoeff(:,:,:,22) = -Whei*(Rlam*zetax*etaz+Rmu*zetaz*etax)*Jac
!     Tdomain%specel(n)%Acoeff(:,:,:,23) = -Whei*(Rlam+Rmu)*zetax*zetaz*Jac
!
!     Tdomain%specel(n)%Acoeff(:,:,:,24) = -Whei*(RKmod*xiy**2+Rmu* (xix**2+xiz**2))*Jac
!     Tdomain%specel(n)%Acoeff(:,:,:,25) = -Whei*(RKmod*xiy*etay+Rmu* (xix*etax+xiz*etaz))*Jac
!     Tdomain%specel(n)%Acoeff(:,:,:,26) = -Whei*(RKmod*xiy*zetay+Rmu*  (xix*zetax+xiz*zetaz))*Jac
!     Tdomain%specel(n)%Acoeff(:,:,:,27) = -Whei*(Rlam+Rmu)*xiy*xiz*Jac
!     Tdomain%specel(n)%Acoeff(:,:,:,28) = -Whei*(Rlam*etaz*xiy+Rmu*etay*xiz)*Jac
!     Tdomain%specel(n)%Acoeff(:,:,:,29) = -Whei*(Rlam*zetaz*xiy+Rmu*zetay*xiz)*Jac
!
!     Tdomain%specel(n)%Acoeff(:,:,:,30) = -Whei*(RKmod*etay**2+Rmu* (etax**2+etaz**2))*Jac
!     Tdomain%specel(n)%Acoeff(:,:,:,31) = -Whei*(RKmod*zetay*etay+Rmu* (zetax*etax+zetaz*etaz))*Jac
!     Tdomain%specel(n)%Acoeff(:,:,:,32) = -Whei*(Rlam*etay*xiz+Rmu*etaz*xiy)*Jac
!     Tdomain%specel(n)%Acoeff(:,:,:,33) = -Whei*(Rlam+Rmu)*etay*etaz*Jac
!     Tdomain%specel(n)%Acoeff(:,:,:,34) = -Whei*(Rlam*zetaz*etay+Rmu*zetay*etaz)*Jac
!
!     Tdomain%specel(n)%Acoeff(:,:,:,35) = -Whei*(RKmod*zetay**2+Rmu* (zetax**2+zetaz**2))*Jac
!     Tdomain%specel(n)%Acoeff(:,:,:,36) = -Whei*(Rlam*xiz*zetay+Rmu*xiy*zetaz)*Jac
!     Tdomain%specel(n)%Acoeff(:,:,:,37) = -Whei*(Rlam*zetay*etaz+Rmu*zetaz*etay)*Jac
!     Tdomain%specel(n)%Acoeff(:,:,:,38) = -Whei*(Rlam+Rmu)*zetay*zetaz*Jac
!
!     Tdomain%specel(n)%Acoeff(:,:,:,39) = -Whei*(RKmod*xiz**2+Rmu*  (xix**2+xiy**2))*Jac 
!     Tdomain%specel(n)%Acoeff(:,:,:,40) = -Whei*(RKmod*xiz*etaz+Rmu* (xix*etax+xiy*etay))*Jac
!     Tdomain%specel(n)%Acoeff(:,:,:,41) = -Whei*(RKmod*xiz*zetaz+Rmu* (xix*zetax+xiy*zetay))*Jac
!                                       
!     Tdomain%specel(n)%Acoeff(:,:,:,42) = -Whei*(RKmod*etaz**2+Rmu* (etax**2+etay**2))*Jac
!     Tdomain%specel(n)%Acoeff(:,:,:,43) = -Whei*(RKmod*zetaz*etaz+Rmu* (zetax*etax+zetay*etay))*Jac
!                                       
!     Tdomain%specel(n)%Acoeff(:,:,:,44) = -Whei*(RKmod*zetaz**2+Rmu* (zetax**2+zetay**2))*Jac
!
!     deallocate (Tdomain%specel(n)%InvGrad)
                                       
   else

     RLam = Tdomain%specel(n)%Lambda
     RMu = Tdomain%specel(n)%Mu
     RKmod = RLam + 2*RMu

     call define_Acoeff_PML_iso(ngllx,nglly,ngllz,Rkmod,Rmu,Rlam,xix,xiy,xiz, &
            etax,etay,etaz,zetax,zetay,zetaz,Whei,Jac,Tdomain%specel(n)%Acoeff)

     allocate (wx (0:ngllx-1,0:nglly-1,0:ngllz-1))
     allocate (wy (0:ngllx-1,0:nglly-1,0:ngllz-1))
     allocate (wz (0:ngllx-1,0:nglly-1,0:ngllz-1))

     call define_alpha_PML(Tdomain%sSubDomain(mat)%Px,0,Tdomain%sSubDomain(mat)%Left, &
               ngllx,nglly,ngllz,ngllx,Tdomain%n_glob_points,Tdomain%GlobCoord,       &
               Tdomain%sSubDomain(mat)%GLLcx,RKmod(:,0,0),                            &
               Tdomain%specel(n)%Density(:,0,0),Tdomain%specel(n)%Iglobnum(0,int((nglly-1)/2),int((ngllz-1)/2)),   &
               Tdomain%specel(n)%Iglobnum(ngllx-1,int((nglly-1)/2),int((ngllz-1)/2)),Tdomain%sSubdomain(mat)%Apow, &
               Tdomain%sSubdomain(mat)%npow,wx)

     call define_alpha_PML(Tdomain%sSubDomain(mat)%Py,1,Tdomain%sSubDomain(mat)%Forward, &
               ngllx,nglly,ngllz,nglly,Tdomain%n_glob_points,Tdomain%GlobCoord,          &
               Tdomain%sSubDomain(mat)%GLLcy,RKmod(0,:,0),                               &
               Tdomain%specel(n)%Density(0,:,0),Tdomain%specel(n)%Iglobnum(int((ngllx-1)/2),0,int((ngllz-1)/2)),   &
               Tdomain%specel(n)%Iglobnum(int((ngllx-1)/2),nglly-1,int((ngllz-1)/2)),Tdomain%sSubdomain(mat)%Apow, &
               Tdomain%sSubdomain(mat)%npow,wy)

     call define_alpha_PML(Tdomain%sSubDomain(mat)%Pz,2,Tdomain%sSubDomain(mat)%Down, &
               ngllx,nglly,ngllz,ngllz,Tdomain%n_glob_points,Tdomain%GlobCoord,       &
               Tdomain%sSubDomain(mat)%GLLcz,RKmod(0,0,:),                            &
               Tdomain%specel(n)%Density(0,0,:),Tdomain%specel(n)%Iglobnum(int((ngllx-1)/2),int((nglly-1)/2),0),   &
               Tdomain%specel(n)%Iglobnum(int((ngllx-1)/2),int((nglly-1)/2),ngllz-1),Tdomain%sSubdomain(mat)%Apow, &
               Tdomain%sSubdomain(mat)%npow,wz)

     if (Tdomain%MPML) then
         allocate(temp_PMLx(0:ngllx-1,0:nglly-1,0:ngllz-1))
         allocate(temp_PMLy(0:ngllx-1,0:nglly-1,0:ngllz-1))
         temp_PMLx(:,:,:) = wx(:,:,:)
         temp_PMLy(:,:,:) = wy(:,:,:)
         wx(:,:,:) = wx(:,:,:)+Tdomain%MPML_coeff*(wy(:,:,:)+wz(:,:,:))
         wy(:,:,:) = wy(:,:,:)+Tdomain%MPML_coeff*(temp_PMLx(:,:,:)+wz(:,:,:))
         wz(:,:,:) = wz(:,:,:)+Tdomain%MPML_coeff*(temp_PMLx(:,:,:)+temp_PMLy(:,:,:))
         deallocate(temp_PMLx,temp_PMLy)
     endif
    
     dt = Tdomain%sTimeParam%dt 
     call define_PML_DumpInit(ngllx,nglly,ngllz,dt,wx,Tdomain%specel(n)%Density,whei,Jac, &
                              Tdomain%specel(n)%DumpSx,Tdomain%specel(n)%DumpMass(:,:,:,0))
     call define_PML_DumpInit(ngllx,nglly,ngllz,dt,wy,Tdomain%specel(n)%Density,whei,Jac, &
                              Tdomain%specel(n)%DumpSy,Tdomain%specel(n)%DumpMass(:,:,:,1))
     call define_PML_DumpInit(ngllx,nglly,ngllz,dt,wz,Tdomain%specel(n)%Density,whei,Jac, &
                              Tdomain%specel(n)%DumpSz,Tdomain%specel(n)%DumpMass(:,:,:,2))

     if (Tdomain%curve) then
         call find_normales(Tdomain,n)
     endif

     deallocate (wx,wy,wz)
     deallocate (Tdomain%specel(n)%InvGrad)
   endif

   deallocate (Jac, xix, xiy, xiz, etax, etay, etaz, zetax, zetay, zetaz, Whei, RKmod, Rmu, Rlam)  

enddo


!!! Communications in the processor !!!

do n = 0,Tdomain%n_elem-1
    call get_Mass_Elem2Obj(Tdomain,n)
enddo


!!! Invert Mass Matrix expression !!!

do n = 0,Tdomain%n_elem-1
   ngllx = Tdomain%specel(n)%ngllx
   nglly = Tdomain%specel(n)%nglly
   ngllz = Tdomain%specel(n)%ngllz

   if (Tdomain%specel(n)%PML) then
      call define_PML_DumpEnd(n,rg,ngllx,nglly,ngllz,Tdomain%specel(n)%MassMat,   &
                              Tdomain%specel(n)%DumpMass(:,:,:,0),Tdomain%specel(n)%DumpVx)
      call define_PML_DumpEnd(n,rg,ngllx,nglly,ngllz,Tdomain%specel(n)%MassMat,   &
                              Tdomain%specel(n)%DumpMass(:,:,:,1),Tdomain%specel(n)%DumpVy)
      call define_PML_DumpEnd(n,rg,ngllx,nglly,ngllz,Tdomain%specel(n)%MassMat,   &
                              Tdomain%specel(n)%DumpMass(:,:,:,2),Tdomain%specel(n)%DumpVz)
      deallocate (Tdomain%specel(n)%DumpMass)
   endif

   allocate (LocMassMat(1:ngllx-2,1:nglly-2,1:ngllz-2))
   LocMassMat(:,:,:) = Tdomain%specel(n)%MassMat(1:ngllx-2,1:nglly-2,1:ngllz-2)
   LocMassmat = 1.d0/LocMassMat
   deallocate (Tdomain%specel(n)%MassMat) 
   allocate (Tdomain%specel(n)%MassMat(1:ngllx-2,1:nglly-2,1:ngllz-2) )
   Tdomain%specel(n)%MassMat = LocMassMat
   deallocate (LocMassMat)
enddo


!!! MPI communications !!!

! if (Tdomain%n_proc>1) then
do n = 0,Tdomain%n_proc-1
    ngll = 0
    ngllPML = 0
    ngllocean = 0
    do i = 0,Tdomain%sComm(n)%nb_faces-1
        nf = Tdomain%sComm(n)%faces(i)
        do j = 1,Tdomain%sFace(nf)%ngll2-2
            do k = 1,Tdomain%sFace(nf)%ngll1-2
                Tdomain%sComm(n)%Give(ngll) = Tdomain%sFace(nf)%MassMat(k,j)
                ngll = ngll + 1
            enddo
        enddo
        if (Tdomain%sFace(nf)%PML) then
            do j = 1,Tdomain%sFace(nf)%ngll2-2
                do k = 1,Tdomain%sFace(nf)%ngll1-2
                    Tdomain%sComm(n)%GivePML(ngllPML,0:2) = Tdomain%sFace(nf)%DumpMass(k,j,0:2)
                    ngllPML = ngllPML + 1
                enddo
            enddo
        endif
    enddo
    do i = 0,Tdomain%sComm(n)%nb_edges-1
        ne = Tdomain%sComm(n)%edges(i)
        do j = 1,Tdomain%sEdge(ne)%ngll-2
            Tdomain%sComm(n)%Give(ngll) = Tdomain%sEdge(ne)%MassMat(j)
            ngll = ngll + 1
        enddo
        if (Tdomain%sEdge(ne)%ocean) then
            do j = 1,Tdomain%sEdge(ne)%ngll-2
                Tdomain%sComm(n)%Giveocean(ngllocean) = Tdomain%sEdge(ne)%Mocean(j)
                ngllocean = ngllocean + 1
            enddo
        endif
        if (Tdomain%sEdge(ne)%PML) then
            do j = 1,Tdomain%sEdge(ne)%ngll-2
                Tdomain%sComm(n)%GivePML(ngllPML,0:2) = Tdomain%sEdge(ne)%DumpMass(j,0:2)
                ngllPML = ngllPML + 1
            enddo
        endif
    enddo
    do i = 0,Tdomain%sComm(n)%nb_vertices-1
        nv =  Tdomain%sComm(n)%vertices(i)
        Tdomain%sComm(n)%Give(ngll) = Tdomain%svertex(nv)%MassMat
        ngll = ngll + 1
        if (Tdomain%sVertex(nv)%ocean) then
            Tdomain%sComm(n)%Giveocean(ngllocean) = Tdomain%sVertex(nv)%Mocean
            ngllocean = ngllocean + 1
        endif
        if (Tdomain%sVertex(nv)%PML) then
            Tdomain%sComm(n)%GivePML(ngllPML,0:2) = Tdomain%sVertex(nv)%DumpMass(0:2)
            ngllPML = ngllPML + 1
        endif
    enddo
enddo

n = Tdomain%n_proc
do shift = 1,n-1
    I_give_to = rg + shift
    if (I_give_to > n-1)   I_give_to = I_give_to - n
    I_take_from = rg - shift
    if (I_take_from < 0)   I_take_from = I_take_from + n
    if (mod(n,shift)==0 .and. shift/=1) then
        n_rings = shift
    else if (mod(n,n-shift)==0 .and. shift/=n-1) then
        n_rings = n-shift
    else if (mod(n,2)==0 .and. mod(shift,2)==0) then
        n_rings = 2
    else
        n_rings = 1
    endif
    do i = 0,n_rings-1
        if (rg==i) then
            if (Tdomain%sComm(I_give_to)%ngll>0) then
             call MPI_SEND (Tdomain%sComm(I_give_to)%Give, Tdomain%sComm(I_give_to)%ngll, &
                            MPI_DOUBLE_PRECISION, I_give_to, etiquette, MPI_COMM_WORLD, code)
            endif
            if (Tdomain%sComm(I_take_from)%ngll>0) then
             call MPI_RECV (Tdomain%sComm(I_take_from)%Take, Tdomain%sComm(I_take_from)%ngll, &
                            MPI_DOUBLE_PRECISION, I_take_from, etiquette, MPI_COMM_WORLD, statut, code)
            endif
        else
            do j = 0,n/n_rings-1
                if (rg == i + j*n_rings) then
                    if (Tdomain%sComm(I_take_from)%ngll>0) then
                     call MPI_RECV (Tdomain%sComm(I_take_from)%Take, Tdomain%sComm(I_take_from)%ngll, &
                                    MPI_DOUBLE_PRECISION, I_take_from, etiquette, MPI_COMM_WORLD, statut, code)
                    endif
                    if (Tdomain%sComm(I_give_to)%ngll>0) then
                     call MPI_SEND (Tdomain%sComm(I_give_to)%Give, Tdomain%sComm(I_give_to)%ngll, &
                                    MPI_DOUBLE_PRECISION, I_give_to, etiquette, MPI_COMM_WORLD, code)
                    endif
                endif
            enddo
        endif
    enddo
    do i = 0,n_rings-1
        if (rg==i) then
            if (Tdomain%sComm(I_give_to)%ngllocean>0) then
             call MPI_SEND (Tdomain%sComm(I_give_to)%Giveocean, Tdomain%sComm(I_give_to)%ngllocean, &
                            MPI_DOUBLE_PRECISION, I_give_to, etiquette, MPI_COMM_WORLD, code)
            endif
            if (Tdomain%sComm(I_take_from)%ngllocean>0) then
             call MPI_RECV (Tdomain%sComm(I_take_from)%Takeocean, Tdomain%sComm(I_take_from)%ngllocean, &
                            MPI_DOUBLE_PRECISION, I_take_from, etiquette, MPI_COMM_WORLD, statut, code)
            endif
        else
            do j = 0,n/n_rings-1
                if (rg == i + j*n_rings) then
                    if (Tdomain%sComm(I_take_from)%ngllocean>0) then
                     call MPI_RECV (Tdomain%sComm(I_take_from)%Takeocean, Tdomain%sComm(I_take_from)%ngllocean, &
                                    MPI_DOUBLE_PRECISION, I_take_from, etiquette, MPI_COMM_WORLD, statut, code)
                    endif
                    if (Tdomain%sComm(I_give_to)%ngllocean>0) then
                     call MPI_SEND (Tdomain%sComm(I_give_to)%Giveocean, Tdomain%sComm(I_give_to)%ngllocean, &
                                    MPI_DOUBLE_PRECISION, I_give_to, etiquette, MPI_COMM_WORLD, code)
                    endif
                endif
            enddo
        endif
    enddo
    do i = 0,n_rings-1
        if (rg==i) then
            if (Tdomain%sComm(I_give_to)%ngllPML>0) then
             call MPI_SEND (Tdomain%sComm(I_give_to)%GivePML, 3*Tdomain%sComm(I_give_to)%ngllPML, &
                            MPI_DOUBLE_PRECISION, I_give_to, etiquette, MPI_COMM_WORLD, code)
            endif
            if (Tdomain%sComm(I_take_from)%ngllPML>0) then
             call MPI_RECV (Tdomain%sComm(I_take_from)%TakePML, 3*Tdomain%sComm(I_take_from)%ngllPML, &
                            MPI_DOUBLE_PRECISION, I_take_from, etiquette, MPI_COMM_WORLD, statut, code)
            endif
        else
            do j = 0,n/n_rings-1
                if (rg == i + j*n_rings) then
                    if (Tdomain%sComm(I_take_from)%ngllPML>0) then
                     call MPI_RECV (Tdomain%sComm(I_take_from)%TakePML, 3*Tdomain%sComm(I_take_from)%ngllPML, &
                                    MPI_DOUBLE_PRECISION, I_take_from, etiquette, MPI_COMM_WORLD, statut, code)
                    endif
                    if (Tdomain%sComm(I_give_to)%ngllPML>0) then
                     call MPI_SEND (Tdomain%sComm(I_give_to)%GivePML, 3*Tdomain%sComm(I_give_to)%ngllPML, &
                                    MPI_DOUBLE_PRECISION, I_give_to, etiquette, MPI_COMM_WORLD, code)
                    endif
                endif
            enddo
        endif
    enddo
enddo

do n = 0,Tdomain%n_proc-1
    call Comm_Mass_Obj(Tdomain,n,ngll,ngllPML,ngllocean)

    ! Compilation options  "-O0 -C" sometimes doesn't like these following deallocate!?!
    if (Tdomain%sComm(n)%ngll>0) then
        deallocate (Tdomain%sComm(n)%Give)
        deallocate (Tdomain%sComm(n)%Take)
    endif
    if (Tdomain%sComm(n)%ngllocean>0) then
        deallocate (Tdomain%sComm(n)%Giveocean)
        deallocate (Tdomain%sComm(n)%Takeocean)
    endif
    if (Tdomain%sComm(n)%ngllPML>0) then
        deallocate (Tdomain%sComm(n)%GivePML)
        deallocate (Tdomain%sComm(n)%TakePML)
    endif
enddo
!endif


!!! Dumping factor for PML !!!

do nf = 0,Tdomain%n_face-1
    if (Tdomain%sFace(nf)%PML) then
        ngll1 = Tdomain%sFace(nf)%ngll1;   ngll2 = Tdomain%sFace(nf)%ngll2
        call define_PML_Face_DumpEnd(ngll1,ngll2,Tdomain%sFace(nf)%Massmat,  &
             Tdomain%sFace(nf)%DumpMass(:,:,0),Tdomain%sFace(nf)%DumpVx)
        call define_PML_Face_DumpEnd(ngll1,ngll2,Tdomain%sFace(nf)%Massmat,  &
             Tdomain%sFace(nf)%DumpMass(:,:,1),Tdomain%sFace(nf)%DumpVy)
        call define_PML_Face_DumpEnd(ngll1,ngll2,Tdomain%sFace(nf)%Massmat,  &
             Tdomain%sFace(nf)%DumpMass(:,:,2),Tdomain%sFace(nf)%DumpVz)
        deallocate (Tdomain%sFace(nf)%DumpMass)
    endif
    if (Tdomain%sFace(nf)%ocean) then
        ngll1 = Tdomain%sFace(nf)%ngll1
        ngll2 = Tdomain%sFace(nf)%ngll2
        do i = 1,ngll1-2
         do j = 1,ngll2-2
             do ii = 0,2
              do jj = 0,2
                  Tdomain%sFace(nf)%M33ocean(i,j,ii,jj) = Tdomain%sFace(nf)%Mocean(i,j) * &
                                                          Tdomain%sFace(nf)%verticale(i,j,ii,jj)
                  if (ii==jj)   Tdomain%sFace(nf)%M33ocean(i,j,ii,jj) = Tdomain%sFace(nf)%M33ocean(i,j,ii,jj) + &
                                                                        Tdomain%sFace(nf)%MassMat(i,j)
              enddo
             enddo
             call invert_3d (Tdomain%sFace(nf)%M33ocean(i,j,0:2,0:2),det)
         enddo
        enddo
        deallocate (Tdomain%sFace(nf)%Mocean, Tdomain%sFace(nf)%verticale)
    endif
    Tdomain%sFace(nf)%MassMat = 1./Tdomain%sFace(nf)%MassMat
enddo

do ne = 0,Tdomain%n_edge-1
    if (Tdomain%sEdge(ne)%PML) then
        ngll = Tdomain%sEdge(ne)%ngll
        call define_PML_Edge_DumpEnd(ngll,Tdomain%sEdge(ne)%Massmat,    &
             Tdomain%sEdge(ne)%DumpMass(:,0),Tdomain%sEdge(ne)%DumpVx)
        call define_PML_Edge_DumpEnd(ngll,Tdomain%sEdge(ne)%Massmat,    &
             Tdomain%sEdge(ne)%DumpMass(:,1),Tdomain%sEdge(ne)%DumpVy)
        call define_PML_Edge_DumpEnd(ngll,Tdomain%sEdge(ne)%Massmat,    &
             Tdomain%sEdge(ne)%DumpMass(:,2),Tdomain%sEdge(ne)%DumpVz)
        deallocate (Tdomain%sEdge(ne)%DumpMass)
    endif
    if (Tdomain%sEdge(ne)%ocean) then
        ngll = Tdomain%sEdge(ne)%ngll
        do i = 1,ngll-2
            do ii = 0,2
             do jj = 0,2
                 Tdomain%sEdge(ne)%M33ocean(i,ii,jj) = Tdomain%sEdge(ne)%Mocean(i) * &
                                                       Tdomain%sEdge(ne)%verticale(i,ii,jj)
                 if (ii==jj)   Tdomain%sEdge(ne)%M33ocean(i,ii,jj) = Tdomain%sEdge(ne)%M33ocean(i,ii,jj) + &
                                                                     Tdomain%sEdge(ne)%MassMat(i)
             enddo
            enddo
            call invert_3d (Tdomain%sEdge(ne)%M33ocean(i,0:2,0:2),det)
        enddo
        deallocate (Tdomain%sEdge(ne)%Mocean, Tdomain%sEdge(ne)%verticale)
    endif
    Tdomain%sEdge(ne)%MassMat = 1./ Tdomain%sEdge(ne)%MassMat
enddo

do nv = 0,Tdomain%n_vertex-1
    if (Tdomain%sVertex(nv)%PML) then
        call define_PML_Vertex_DumpEnd(Tdomain%sVertex(nv)%Massmat,    &
             Tdomain%sVertex(nv)%DumpMass(0),Tdomain%sVertex(nv)%DumpVx)
        call define_PML_Vertex_DumpEnd(Tdomain%sVertex(nv)%Massmat,    &
             Tdomain%sVertex(nv)%DumpMass(1),Tdomain%sVertex(nv)%DumpVy)
        call define_PML_Vertex_DumpEnd(Tdomain%sVertex(nv)%Massmat,    &
             Tdomain%sVertex(nv)%DumpMass(2),Tdomain%sVertex(nv)%DumpVz)
        deallocate (Tdomain%sVertex(nv)%DumpMass)
    endif
    if (Tdomain%sVertex(nv)%ocean) then
        do ii = 0,2
         do jj = 0,2
             Tdomain%sVertex(nv)%M33ocean(ii,jj) = Tdomain%sVertex(nv)%Mocean * &
                                                   Tdomain%sVertex(nv)%verticale(ii,jj)
             if (ii==jj)   Tdomain%sVertex(nv)%M33ocean(ii,jj) = Tdomain%sVertex(nv)%M33ocean(ii,jj) + &
                                                                 Tdomain%sVertex(nv)%MassMat
         enddo
        enddo
        call invert_3d (Tdomain%sVertex(nv)%M33ocean(0:2,0:2),det)
        deallocate (Tdomain%sVertex(nv)%verticale)
    endif
    Tdomain%sVertex(nv)%MassMat = 1.d0/Tdomain%sVertex(nv)%MassMat
enddo

return


end subroutine define_arrays
!------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------
subroutine define_Acoeff_PML_iso(ngllx,nglly,ngllz,Rkmod,Rmu,Rlam,xix,xiy,xiz, &
                                 etax,etay,etaz,zetax,zetay,zetaz,Whei,Jac,Acoeff)

    implicit none

    integer, intent(in) :: ngllx,nglly,ngllz
    doubleprecision, dimension(0:ngllx-1,0:nglly-1,0:ngllz-1), intent(in) :: Rkmod,Rmu,Rlam, &
                                         xix,xiy,xiz,etax,etay,etaz,zetax,zetay,zetaz,Whei,Jac
    doubleprecision, dimension(0:ngllx-1,0:nglly-1,0:ngllz-1,0:35), intent(out) :: Acoeff


    Acoeff(:,:,:,0) = RKmod *xix
    Acoeff(:,:,:,1) = RKmod *etax
    Acoeff(:,:,:,2) = RKmod *zetax
    Acoeff(:,:,:,3) = RLam *xiy
    Acoeff(:,:,:,4) = RLam *etay
    Acoeff(:,:,:,5) = RLam *zetay
    Acoeff(:,:,:,6) = RLam *xiz
    Acoeff(:,:,:,7) = RLam *etaz
    Acoeff(:,:,:,8) = RLam *zetaz
    Acoeff(:,:,:,9) = RLam *xix
    Acoeff(:,:,:,10) = RLam *etax
    Acoeff(:,:,:,11) = RLam *zetax
    Acoeff(:,:,:,12) = RKmod *xiy
    Acoeff(:,:,:,13) = RKmod *etay
    Acoeff(:,:,:,14) = RKmod *zetay
    Acoeff(:,:,:,15) = RKmod *xiz
    Acoeff(:,:,:,16) = RKmod *etaz
    Acoeff(:,:,:,17) = RKmod *zetaz
    Acoeff(:,:,:,18) = RMu *xix
    Acoeff(:,:,:,19) = RMu *etax
    Acoeff(:,:,:,20) = RMu *zetax
    Acoeff(:,:,:,21) = RMu *xiy
    Acoeff(:,:,:,22) = RMu *etay
    Acoeff(:,:,:,23) = RMu *zetay
    Acoeff(:,:,:,24) = RMu *xiz
    Acoeff(:,:,:,25) = RMu *etaz
    Acoeff(:,:,:,26) = RMu *zetaz
    Acoeff(:,:,:,27) = -Whei * xix * Jac
    Acoeff(:,:,:,28) = -Whei * xiy * Jac
    Acoeff(:,:,:,29) = -Whei * xiz * Jac
    Acoeff(:,:,:,30) = -Whei * etax * Jac
    Acoeff(:,:,:,31) = -Whei * etay * Jac
    Acoeff(:,:,:,32) = -Whei * etaz * Jac
    Acoeff(:,:,:,33) = -Whei * zetax * Jac
    Acoeff(:,:,:,34) = -Whei * zetay * Jac
    Acoeff(:,:,:,35) = -Whei * zetaz * Jac

    return

end subroutine define_Acoeff_PML_iso
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
subroutine define_alpha_PML(lattenu,dir,ldir_attenu,ngllx,nglly,ngllz,ngll,n_pts,   &
                            Coord,GLLc,Rkmod,Density,ind_min,ind_max,Apow,npow,alpha)
! determining attenuation profile in a PML layer (see Festa et al, GJI 2005).
! dir = attenuation's direction, ldir_attenu = the logical giving the orientation.

    use angles

    implicit none

    logical, intent(in) :: lattenu,ldir_attenu
    integer, intent(in) :: dir,ngllx,nglly,ngllz,ngll,n_pts,ind_min,ind_max,npow
    doubleprecision, dimension(0:2,0:n_pts-1), intent(in) :: Coord
    doubleprecision, dimension(0:ngll-1), intent(in) :: GLLc,RKmod,Density
    doubleprecision, intent(in) :: Apow
    doubleprecision, dimension(0:ngllx-1,0:nglly-1,0:ngllz-1), intent(out) :: alpha

    integer :: i
    doubleprecision :: dh
    doubleprecision, dimension(0:2) :: u,v
    doubleprecision, dimension(0:ngll-1) :: ri,vp
    doubleprecision, external :: pow

    if(.not.lattenu) then   ! no attenuation in the dir-direction
        alpha(:,:,:) = 0d0
    else  ! yes, attenuation in this dir-direction
!        dh = Coord(dir,ind_min)
!        dh = dabs(Coord(dir,ind_max)-dh)
        dh = distance3d(Coord(0:2,ind_min),Coord(0:2,ind_max))
        if(ldir_attenu) then  ! Left in x, Forward in y, Down in z 
            ri(:) = 0.5d0*(1d0+GLLc(ngll-1:0:-1))*dble(ngll-1)
        else  ! Right in x, Backward in y, Up in z 
            ri(:) = 0.5d0*(1d0+GLLc(0:ngll-1))*dble(ngll-1)
        endif
        vp(:) = dsqrt(Rkmod(:)/Density(:))
        select case(dir)
            case(0)  ! dir = x
                do i = 0,ngll-1
                    alpha(i,0:,0:) = pow(ri(i),vp(i),ngll-1,dh,Apow,npow)
                enddo
            case(1)  ! dir = y
                do i = 0,ngll-1
                    alpha(0:,i,0:) = pow(ri(i),vp(i),ngll-1,dh,Apow,npow)
                enddo
            case(2)  ! dir = z
                do i = 0,ngll-1
                    alpha(0:,0:,i) = pow(ri(i),vp(i),ngll-1,dh,Apow,npow)
                enddo
        end select
    endif

    return

end subroutine define_alpha_PML
!----------------------------------------------------------------------------------

!----------------------------------------------------------------------------------
subroutine define_PML_DumpInit(ngllx,nglly,ngllz,dt,alpha,density,whei,jac,DumpS,DumpMass)
! defining parameters related to stresses and mass matrix elements, in the case of
! a PML, along a given splitted direction.

    implicit none

    integer, intent(in) :: ngllx,nglly,ngllz
    doubleprecision, intent(in) :: dt
    doubleprecision, dimension(0:ngllx-1,0:nglly-1,0:ngllz-1), intent(in) :: alpha
    doubleprecision, dimension(0:ngllx-1,0:nglly-1,0:ngllz-1), intent(in) :: density,whei,Jac
    doubleprecision, dimension(0:ngllx-1,0:nglly-1,0:ngllz-1,0:1), intent(out) :: DumpS
    doubleprecision, dimension(0:ngllx-1,0:nglly-1,0:ngllz-1), intent(out) :: DumpMass
    doubleprecision, dimension(0:ngllx-1,0:nglly-1,0:ngllz-1) :: Id

    Id = 1d0

    DumpS(:,:,:,1) = Id + 0.5d0*dt*alpha
    DumpS(:,:,:,1) = 1d0/DumpS(:,:,:,1)
    DumpS(:,:,:,0) = (Id - 0.5d0*dt*alpha)*DumpS(:,:,:,1)

    DumpMass(:,:,:) = 0.5d0*Density(:,:,:)*Whei(:,:,:)*Jac(:,:,:)*alpha(:,:,:)*dt

    return

end subroutine define_PML_DumpInit
!----------------------------------------------------------------------------------

!----------------------------------------------------------------------------------
subroutine define_PML_DumpEnd(n,rg,ngllx,nglly,ngllz,Massmat,DumpMass,DumpV)

    implicit none

    integer, intent(in) :: ngllx,nglly,ngllz,n,rg
    doubleprecision, dimension(0:ngllx-1,0:nglly-1,0:ngllz-1), intent(in) :: MassMat
    doubleprecision, dimension(1:ngllx-2,1:nglly-2,1:ngllz-2,0:1), intent(out) :: DumpV
    doubleprecision, dimension(0:ngllx-1,0:nglly-1,0:ngllz-1), intent(in) :: DumpMass
    doubleprecision, dimension(1:ngllx-2,1:nglly-2,1:ngllz-2) :: LocMassMat

    LocMassMat(:,:,:) = MassMat(1:ngllx-2,1:nglly-2,1:ngllz-2)
    DumpV(:,:,:,1) = LocMassMat + DumpMass(1:ngllx-2,1:nglly-2,1:ngllz-2)
    DumpV(:,:,:,1) = 1d0/DumpV(:,:,:,1)
    DumpV(:,:,:,0) = LocMassMat - DumpMass(1:ngllx-2,1:nglly-2,1:ngllz-2)
    DumpV(:,:,:,0) = DumpV(:,:,:,0) * DumpV(:,:,:,1)

    return

end subroutine define_PML_DumpEnd
!----------------------------------------------------------------------------------

!----------------------------------------------------------------------------------
subroutine define_PML_Face_DumpEnd(ngll1,ngll2,Massmat,DumpMass,DumpV)

    implicit none

    integer, intent(in) :: ngll1,ngll2
    doubleprecision, dimension(1:ngll1-2,1:ngll2-2), intent(in) :: MassMat
    doubleprecision, dimension(1:ngll1-2,1:ngll2-2,0:1), intent(out) :: DumpV
    doubleprecision, dimension(1:ngll1-2,1:ngll2-2), intent(in) :: DumpMass

    DumpV(:,:,1) = MassMat + DumpMass
    DumpV(:,:,1) = 1d0/DumpV(:,:,1)
    DumpV(:,:,0) = MassMat - DumpMass
    DumpV(:,:,0) = DumpV(:,:,0) * DumpV(:,:,1)

    return

end subroutine define_PML_Face_DumpEnd
!----------------------------------------------------------------------------------

!----------------------------------------------------------------------------------

subroutine define_PML_Edge_DumpEnd(ngll,Massmat,DumpMass,DumpV)

    implicit none

    integer, intent(in) :: ngll
    doubleprecision, dimension(1:ngll-2), intent(in) :: MassMat
    doubleprecision, dimension(1:ngll-2,0:1), intent(out) :: DumpV
    doubleprecision, dimension(1:ngll-2), intent(in) :: DumpMass

    DumpV(:,1) = MassMat + DumpMass
    DumpV(:,1) = 1d0/DumpV(:,1)
    DumpV(:,0) = MassMat - DumpMass
    DumpV(:,0) = DumpV(:,0) * DumpV(:,1)

    return

end subroutine define_PML_Edge_DumpEnd
!----------------------------------------------------------------------------------

!----------------------------------------------------------------------------------
subroutine define_PML_Vertex_DumpEnd(Massmat,DumpMass,DumpV)

    implicit none

    doubleprecision, intent(in) :: MassMat
    doubleprecision, dimension(0:1), intent(out) :: DumpV
    doubleprecision, intent(in) :: DumpMass

    DumpV(1) = MassMat + DumpMass
    DumpV(1) = 1d0/DumpV(1)
    DumpV(0) = MassMat - DumpMass
    DumpV(0) = DumpV(0) * DumpV(1)

    return

end subroutine define_PML_Vertex_DumpEnd
!----------------------------------------------------------------------------------
