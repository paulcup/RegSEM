program mesher


use module_ellipticity


implicit none

type :: processor
integer, dimension(:), pointer :: Obj
end type

type :: souvenir
type(processor), dimension(:), pointer :: rank 
end type

type(souvenir), dimension(:), pointer :: memory
doubleprecision :: x_len,y_len,z_len, dx,dy,dz, rayon, R, X, Y, D, ratio, Clat,Clong, alpha, dRmax, &
                   xa,ya,za, xs,ys,zs, theta,phi, limit_inf,limit_med, dist, Ar,Atheta,Aphi, &
                   epsil_xi,epsil_eta,epsil_r, latmin,latmax,dlat, lonmin,lonmax,dlong, Rterre, P2
doubleprecision, parameter :: pi = 3.1415926535897932384626433832795028841974, &
                              seuil = 0.2, &
                              vs_surface = 3200.d0
doubleprecision, dimension(0:3) :: weight, rad_loc
doubleprecision, dimension(0:11) :: rad, vs
doubleprecision, dimension(:), allocatable :: radius, xco,yco,zco, dR, Rsph
doubleprecision, dimension(0:2,0:2) :: rot
doubleprecision, dimension(1:6,0:1) :: cover_def
doubleprecision, dimension(0:7,0:2) :: coord_vertices
doubleprecision, dimension(:,:), allocatable :: coord_nodes, dzmin, dzmax, &
                                                moho,depth_moho, topo,altitude_topo
real, dimension(:,:), allocatable :: coord_grid
integer :: i,j,k, n,m, nel, i_count, num, n_elem, n_elem_xy, n_points, n_pts_xy, nx,ny,nz_tot, &
           edgecut, nparts, proc, ok, n_vertices, n_faces, n_edges, n_nodes, nf,ne,nv, &
           neighbor, neighbor_face, neighbor_edge, nods_per_elem, n_layers, model, nb_couches_model, &
           i1,i2,j1,j2,k1,k2, n_out, n_other_proc, nxtmp,nytmp,nztmp
integer, parameter :: n_dim = 3, etype = 3, numflag = 0, wgtflag = 2, options = 0
integer, dimension(0:2) :: neigh_local
integer, dimension(0:25) :: neigh
integer, dimension(:), allocatable :: nz, elmnts, counter, corner, neighbor_corner, tmp1D, &
                                      Material, Cover, moho_position, topo_position, &
                                      dxadj, dxadjncy, part, vwgt, adjwgt, &
                                      nelem_in_proc, which_vertices, nf_shared, ne_shared, nv_shared
integer, dimension(0:2,0:2) :: other_proc
integer, dimension(:,:), allocatable :: elmnts_local, which_elem_in_proc, nodes, tmp2D, elmnts_indexes, &
                                        faces, faces_shared, edges, edges_shared, vertices_shared
logical :: curve, any_PML, free_surface, free_bottom, topo_log, ellipticity, random_trac, mirror, flag
logical, dimension(:), allocatable :: L_Proc, L_Nodes
character(len=1) :: yes_or_no
character(len=13) :: meshfilename

integer, parameter :: n_cover = 4
logical, parameter :: output_chunk = .false., &
                      write_vtk = .true.


!!! Inputs !!!
write (*,*) "Do you want to take into account the Earth sphericity? [y/n]"
read (*,*) yes_or_no
if (yes_or_no == "y") then
   curve = .true.
   write (*,*) "Introduce the coordinates (lat,long) of the center"
   read (*,*) Clat, Clong
   write (*,*) "Introduce the lengths of the horizontal edges: Xi and Eta (max: 90 degrees)"
   read (*,*) x_len, y_len
   write (*,*) "Introduce the horizontal space step"
   read (*,*) dx;   dy = dx
   write (*,*) "Do you want PREM interfaces? (1)"
   write (*,*) "            PREM interfaces + 3D Moho? (2)" ! Les 3 premieres interfaces de PREM sont alors virees
   write (*,*) "            PREM interfaces + 3D Moho - 220km interface? (3)" ! Les 4 premieres interfaces de PREM sont alors virees
   write (*,*) "            another 1D model? (4)"
   write (*,*) "            SVEMum_SAW24B16LM? (5)"
   write (*,*) "            PREM interfaces - 15km interface? (6)"
   read (*,*) model
   if (model==4) then
      ! Ici on peut construire une boule a couches de rayon quelconque
      ! Une topo de surface pourra etre prise en compte
      ! Les oceans pourront aussi etre consideres. Il faudra alors ajuster le parametre Rterre dans angles.f90
      ! Les vitesses du modele qui sera associe a ce maillage devront croitre avec la profondeur
      write (*,*) "Introduce the number of layers"
      read (*,*) n_layers
      write (*,"(a11,i2.2,a18)") " Introduce ", n_layers+1, " radii (in meter)"
      allocate (radius(0:n_layers))
      read (*,*) (radius(i),i=0,n_layers)
      Rterre = radius(n_layers)
      ellipticity = .false.
   else
      if ((model==2 .or. model==3) .and. (dx<0.199999))then
         print *,"The horizontal step has to be larger than 0.2 degree."
         print *,"Otherwise you'll probably have an under-sampled crust."
         stop ! Sinon on aurait plusieurs couches d'elements dans la croute
      endif
      call define_model(model,nb_couches_model,rad,vs)
      Rterre = 6371000.d0
      rad(nb_couches_model-1) = Rterre
      vs(nb_couches_model-1) = vs_surface
      write (*,*) "Introduce the depth (in meter) of the bottom"
      read (*,*) rayon
      rayon = Rterre - rayon
      bottom : do i = 0,nb_couches_model-1
         if (rayon<rad(i)) then
            n_layers = nb_couches_model - i
            exit bottom
         endif
      enddo bottom
      call bottom_remarks(rayon,rad,i,dx,Rterre,model)
      allocate (radius(0:n_layers))
      radius(0) = rayon
      do j = 1,n_layers
         radius(j) = rad(i)
         i = i + 1
      enddo
      write (*,*) "Do you want to take into account the Earth ellipticity? [y/n]"
      read (*,*) yes_or_no
      ellipticity = .false.
      if (yes_or_no == "y") then
         ellipticity = .true.
         call init_ellipticity()
      endif
   endif
   write (*,*) "Do you want a topography at the surface? [y/n]"
   read (*,*) yes_or_no
   topo_log = .false.
   if (yes_or_no == "y") topo_log = .true.
   write (*,*) "How many nodes do you want for the geometry of each element: 8, 20 or 27?"
   read (*,*) nods_per_elem
else
   curve = .false.
   nods_per_elem = 8
   model = 0
   ellipticity = .false.
   topo_log = .false.
   write (*,*) "Introduce the lengths of the edges: x_len, y_len, z_len (in meter)"
   read (*,*) x_len, y_len, z_len
   write (*,*) "Introduce the space steps: dx, dy, dz"
   read (*,*) dx, dy, dz
endif

write (*,*) "Do you want PML? [y/n]"
read (*,*) yes_or_no
any_PML = .false.
if (yes_or_no == "y") then
   any_PML = .true.
   write (*,*) "Introduce the polynomial order for a PML element, and then for a normal element"
   read (*,*) i, j;   ratio = (i+1.)/(j+1.)
   write (*,*) "Do you want a free surface at the top? [y/n]"
   read (*,*) yes_or_no
   free_surface = .false.;   free_bottom = .false.
   if (yes_or_no == "y") free_surface = .true.
   if (iargc()==1) free_bottom = .true.
endif

write (*,*) "Introduce the number of parts to partition the mesh"
read (*,*) nparts

random_trac = .false.
if (curve .and. nods_per_elem==27) then
   write (*,*) "Do you want a file for random tractions? [y/n]"
   read (*,*) yes_or_no
   if (yes_or_no == "y") then
      random_trac = .true.
      write (*,*) "Introduce the direction of the tractions: A_r, A_theta, A_phi"
      read (*,*) Ar, Atheta, Aphi
      open (25,file="random_trac.table")
      write (25,*) Ar, Atheta, Aphi
   endif
endif

mirror = .false.
write (*,*) "Do you want a time-reversal mirror? [y/n]"
read (*,*) yes_or_no
if (yes_or_no == "y") then
   mirror = .true.
   write (*,*) "Define the mirror with the number of elements that you"
   write (*,*) "want to skip in each direction (-x,x,-y,y,-z,z)?"
   write (*,*) "NB: enter -2 for the directions with no mirror"
   read (*,*) i1,i2,j1,j2,k1,k2
else
   i1=-2; i2=-2; j1=-2; j2=-2; k1=-2; k2=-2
endif


!!! Defining the number of elements !!!
nx = int(x_len/dx)
ny = int(y_len/dy)
if (nx*dx < x_len-dx/2.d0)   nx = nx + 1
if (ny*dy < y_len-dy/2.d0)   ny = ny + 1
dx = x_len/nx
dy = y_len/ny
if (curve) then
   dRmax = deg2rad(dx*radius(n_layers))
   allocate (nz(0:n_layers))
   allocate (dR(0:n_layers))
   nz_tot = 0
   do i = 0,n_layers-1
      z_len = radius(i+1) - radius(i)
      if (model==4) then
         dz = dRmax
      else
         dz = dRmax * vs(nb_couches_model-(n_layers-i))/vs(nb_couches_model-1)
         if (egal(dz,0.d0,1.d0) .eqv. .true.)   dz = dRmax ! Quand on se trouve dans le noyau
      endif
      X = z_len/dz
      Y = X - aint(X)
      if (Y/X > seuil) then
         nz(i) = int(X) + 1
      else
         nz(i) = int(X)
      endif
      dR(i) = z_len/nz(i)
      nz_tot = nz_tot + nz(i)
   enddo
   nz(n_layers) = 1
   dR(n_layers) = 0
else
   nz_tot = int(z_len/dz)
   if (nz_tot*dz < z_len-dz/2.d0)   nz_tot = nz_tot + 1
   dz = z_len/nz_tot
endif
n_points = (nx+1) * (ny+1) * (nz_tot+1)
n_pts_xy = (nx+1) * (ny+1)
n_elem = nx * ny * nz_tot
n_elem_xy = nx * ny


!!! Associating to each element 8 vertices by using a global numbering !!!
allocate (elmnts(0:8*n_elem-1))
allocate (elmnts_indexes(0:n_elem-1,0:2))
do nel = 0, n_elem-1 
   k = nel / n_elem_xy
   j = nel - k*n_elem_xy; j = j/nx
   i = nel - k*n_elem_xy - j*nx 
   i_count = i + j*(nx+1) + k*n_pts_xy
   elmnts_indexes(nel,0) = i
   elmnts_indexes(nel,1) = j
   elmnts_indexes(nel,2) = k
   elmnts(nel*8:nel*8+7) = (/ i_count, i_count+1, i_count+nx+2, i_count+nx+1, i_count+n_pts_xy, &
   i_count+n_pts_xy+1, i_count+n_pts_xy+nx+2, i_count+n_pts_xy+nx+1 /)
enddo


!!! Defining the coordinates of each vertex !!!
allocate (xco(0:n_points-1))
allocate (yco(0:n_points-1))
allocate (zco(0:n_points-1))
if (curve) then
   y_len = -y_len/2.d0
   x_len = -x_len/2.d0
   i_count = 0
   do n = 0,n_layers
      rayon = radius(n)
      dz = dR(n)
      do k = 0,nz(n)-1
         R = rayon + k*dz
         do j = 0,ny
            Y = y_len + j*dy
            do i = 0,nx
               X = x_len + i*dx
               xco(i_count) = X
               yco(i_count) = Y
               zco(i_count) = R
               i_count = i_count + 1
            enddo
         enddo
      enddo
   enddo
else
   do k = 0,nz_tot
      do j = 0,ny
         do i = 0,nx
            i_count = i + j*(nx+1) + k*n_pts_xy
            xco(i_count) = i * dx
            yco(i_count) = j * dy
            zco(i_count) = k * dz
         enddo
      enddo
   enddo
endif


!!! Defining the cover !!!
nxtmp = nx-1-n_cover
nytmp = ny-1-n_cover
nztmp = nz_tot-1-n_cover

j = 0;   i = 0
k = 1;   i_count = i + j*(nx+1) + k*n_pts_xy
cover_def(1,0) = zco(i_count)
k = n_cover+1;    i_count = i + j*(nx+1) + k*n_pts_xy
cover_def(1,1) = zco(i_count)

k = 0;   i = 0
j = 1;   i_count = i + j*(nx+1) + k*n_pts_xy
cover_def(2,0) = yco(i_count)
j = n_cover+1;   i_count = i + j*(nx+1) + k*n_pts_xy
cover_def(2,1) = yco(i_count)

k = 0;   j = 0
i = nx-1;   i_count = i + j*(nx+1) + k*n_pts_xy
cover_def(3,0) = xco(i_count)
i = nxtmp;   i_count = i + j*(nx+1) + k*n_pts_xy
cover_def(3,1) = xco(i_count)

k = 0;   i = 0
j = ny-1;   i_count = i + j*(nx+1) + k*n_pts_xy
cover_def(4,0) = yco(i_count)
j = nytmp;   i_count = i + j*(nx+1) + k*n_pts_xy
cover_def(4,1) = yco(i_count)

k = 0;   j = 0
i = 1;   i_count = i + j*(nx+1) + k*n_pts_xy
cover_def(5,0) = xco(i_count)
i = n_cover+1;   i_count = i + j*(nx+1) + k*n_pts_xy
cover_def(5,1) = xco(i_count)

j = 0;   i = 0
k = nz_tot-1;   i_count = i + j*(nx+1) + k*n_pts_xy
cover_def(6,0) = zco(i_count)
k = nztmp;    i_count = i + j*(nx+1) + k*n_pts_xy
cover_def(6,1) = zco(i_count)


!!! Defining the rotation matrix !!!
if (curve) then
   Clat = 90.d0 - Clat
   if (Clong<0.d0) Clong = 360.d0 + Clong
   alpha = 0.d0
   Clat = deg2rad(Clat); Clong = deg2rad(Clong); alpha = deg2rad(alpha);
   rot(0,0) = dcos(alpha)*dcos(Clat)*dcos(Clong)-dsin(alpha)*dsin(Clong)
   rot(0,1) = -dsin(alpha)*dcos(Clat)*dcos(Clong)-dcos(alpha)*dsin(Clong)
   rot(0,2) = dsin(Clat)*dcos(Clong)
   rot(1,0) = dcos(alpha)*dcos(Clat)*dsin(Clong)+dsin(alpha)*dcos(Clong)
   rot(1,1) = -dsin(alpha)*dcos(Clat)*dsin(Clong)+dcos(alpha)*dcos(Clong)
   rot(1,2) = dsin(Clat)*dsin(Clong)
   rot(2,0) = -dcos(alpha)*dsin(Clat)
   rot(2,1) = dsin(alpha)*dsin(Clat)
   rot(2,2) = dcos(Clat)
endif


!!! Defining the topography of the Moho !!!
!!! On agit comme s'il y avait 27 ctrl_pts car c'est le cas general (8 et 20 sont des sous-cas de 27) !!!
if (model==2 .or. model==3) then
    open (11, file="Moho.asc", status="old")
    read (11,*) latmin, latmax, dlat
    if (mod(latmax-latmin,dlat)/=0) then
        print *,"In the file Moho.asc the latitude step doesn't fit the latitude limits"
        stop
    endif
    read (11,*) lonmin, lonmax, dlong
    if (lonmin<0.d0)   lonmin = 360.d0 + lonmin
    if (lonmax<0.d0)   lonmax = 360.d0 + lonmax
    if (mod(lonmax-lonmin,dlong)/=0) then
        print *,"In the file Moho.asc the longitude step doesn't fit the longitude limits"
        stop
    endif
    m = int((latmax-latmin)/dlat) + 1
    n = int((lonmax-lonmin)/dlong) + 1
    allocate (moho(0:m-1,0:n-1))
    do i = 0,m-1
     do j = 0,n-1
         read (11,*) moho(i,j)
     enddo
    enddo
    moho(:,:) = moho(:,:)*1000.d0
    allocate (depth_moho(0:2*nx,0:2*ny))
    open (17,file="moho.out")
    do j = 0,2*ny
        Y = y_len + j*dy/2.d0
        do i = 0,2*nx
            X = x_len + i*dx/2.d0
            ! Passage en coordonnees cartesiennes
            xa = dtan(deg2rad(X))
            ya = dtan(deg2rad(Y))
            D = dsqrt(1.d0 + xa**2 + ya**2)
            xa = xa/D;   ya = ya/D;   za = 1.d0/D
            ! Rotation du chunk de ref vers le chunk reel
            xs = rot(0,0)*xa + rot(0,1)*ya + rot(0,2)*za
            ys = rot(1,0)*xa + rot(1,1)*ya + rot(1,2)*za
            zs = rot(2,0)*xa + rot(2,1)*ya + rot(2,2)*za
            ! Passage en spherique
            call cart2sph(xs,ys,zs,alpha,theta,phi)
            ! Profondeur du Moho en (theta,phi)
            theta = theta*180.d0/pi;   phi = phi*180.d0/pi
            call read_moho(theta,phi,moho,m,n,latmin,latmax,dlat,lonmin,lonmax,dlong,depth_moho(i,j))
            write (17,*) 90.d0-theta, phi, depth_moho(i,j)
        enddo
    enddo
    close(17)
    deallocate (moho)
    allocate (dzmin(0:1,0:nparts-1));   dzmin = 1000000.d0
    allocate (dzmax(0:1,0:nparts-1));   dzmax = 0.d0
    limit_inf = radius(n_layers-2) + dR(n_layers-2)/10.d0
    limit_med = radius(n_layers-1) - dR(n_layers-2)/10.d0
endif


!!! Defining the topography of the surface!!!
!!! On procede comme pour le Moho !!!
if (topo_log) then
    open (11, file="Topo.asc", status="old")
    read (11,*) latmin, latmax, dlat
    if (mod(latmax-latmin,dlat)/=0) then
        print *,"In the file Topo.asc the latitude step doesn't fit the latitude limits"
        stop
    endif
    read (11,*) lonmin, lonmax, dlong
    if (lonmin<0.d0)   lonmin = 360.d0 - lonmin
    if (lonmax<0.d0)   lonmax = 360.d0 - lonmax
    if (mod(lonmax-lonmin,dlong)/=0) then
        print *,"In the file Topo.asc the longitude step doesn't fit the longitude limits"
        stop
    endif
    m = int((latmax-latmin)/dlat) + 1
    n = int((lonmax-lonmin)/dlong) + 1
    allocate (topo(0:m-1,0:n-1))
    do i = 0,m-1
     do j = 0,n-1
         read (11,*) topo(i,j)
     enddo
    enddo
    allocate (altitude_topo(0:2*nx,0:2*ny))
    open (17,file="topo.out")
    do j = 0,2*ny
        Y = y_len + j*dy/2.d0
        do i = 0,2*nx
            X = x_len + i*dx/2.d0
            ! Passage en coordonnees cartesiennes
            xa = dtan(deg2rad(X))
            ya = dtan(deg2rad(Y))
            D = dsqrt(1.d0 + xa**2 + ya**2)
            xa = xa/D;   ya = ya/D;   za = 1/D
            ! Rotation du chunk de ref vers le chunk reel
            xs = rot(0,0)*xa + rot(0,1)*ya + rot(0,2)*za
            ys = rot(1,0)*xa + rot(1,1)*ya + rot(1,2)*za
            zs = rot(2,0)*xa + rot(2,1)*ya + rot(2,2)*za
            ! Passage en spherique
            call cart2sph(xs,ys,zs,alpha,theta,phi)
            ! Topo en (theta,phi). On utilise la routine ecrite pour le Moho.
            theta = theta*180.d0/pi;   phi = phi*180.d0/pi
            call read_moho(theta,phi,topo,m,n,latmin,latmax,dlat,lonmin,lonmax,dlong,altitude_topo(i,j))
            write(17,*) 90.d0-theta, phi, altitude_topo(i,j)
        enddo
    enddo
    close(17)
    deallocate (topo)
    limit_med = radius(n_layers-1) - dR(n_layers-2)/10.d0
endif


!!! Signing the elements above and below the Moho !!!
!!! On suppose ici qu'il n'y a qu'une seule couche d'elements entre le Moho et la surface !!!
allocate (moho_position(0:n_elem-1));   moho_position = 0
if (model==2 .or. model==3) then
    do k = 0,nz_tot-1
      do j = 0,ny-1
        do i = 0,nx-1
            i_count = i + j*nx + k*n_elem_xy
            if (k==nz_tot-2) then
                moho_position(i_count) = -1
            else if (k==nz_tot-1) then
                moho_position(i_count) = 1
            endif
        enddo
      enddo
    enddo
endif


!!! Signing the elements below the 3D surface !!!
allocate (topo_position(0:n_elem-1));   topo_position = 0
if (topo_log) then
    do k = 0,nz_tot-1
      do j = 0,ny-1
        do i = 0,nx-1
            i_count = i + j*nx + k*n_elem_xy
            if (k==nz_tot-1)   topo_position(i_count) = 1
        enddo
      enddo
    enddo
endif


!!! To each element we associate a number which refers to a material !!!
allocate (Material(0:n_elem-1));   Material = 0
if (any_PML) then
    do k = 0,nz_tot-1
      do j = 0,ny-1
        do i = 0,nx-1
            i_count = i + j*nx + k*n_elem_xy

            if (k==0 .and. (.not.free_bottom)) then
                Material(i_count) = 1
            else if (k==nz_tot-1 .and. (.not.free_surface)) then
                Material(i_count) = 6
            else if (j==0) then
                Material(i_count) = 2
            else if (j==ny-1) then
                Material(i_count) = 4
            else if (i==0) then
                Material(i_count) = 5
            else if (i==nx-1) then
                Material(i_count) = 3
            endif

            if (k==0 .and. j==0 .and. (.not.free_bottom)) then
                Material(i_count) = 7
            else if (k==0 .and. i==nx-1 .and. (.not.free_bottom)) then
                Material(i_count) = 8
            else if (k==0 .and. j==ny-1 .and. (.not.free_bottom)) then
                Material(i_count) = 9
            else if (k==0 .and. i==0 .and. (.not.free_bottom)) then
                Material(i_count) = 10
            else if (i==0 .and. j==0) then
                Material(i_count) = 11
            else if (i==nx-1 .and. j==0) then
                Material(i_count) = 12
            else if (i==nx-1 .and. j==ny-1) then
                Material(i_count) = 13
            else if (i==0 .and. j==ny-1) then
                Material(i_count) = 14
            else if (k==nz_tot-1 .and. j==0 .and. (.not.free_surface)) then
                Material(i_count) = 15
            else if (k==nz_tot-1 .and. i==nx-1 .and. (.not.free_surface)) then
                Material(i_count) = 16
            else if (k==nz_tot-1 .and. j==ny-1 .and. (.not.free_surface)) then
                Material(i_count) = 17
            else if (k==nz_tot-1 .and. i==0 .and. (.not.free_surface)) then
                Material(i_count) = 18
            endif

            if (k==0 .and. j==0 .and. i==0 .and. (.not.free_bottom)) then
                Material(i_count) = 19
            else if (k==0 .and. i==nx-1 .and. j==0 .and. (.not.free_bottom)) then
                Material(i_count) = 20
            else if (k==0 .and. j==ny-1 .and. i==nx-1 .and. (.not.free_bottom)) then
                Material(i_count) = 21
            else if (k==0 .and. i==0 .and. j==ny-1 .and. (.not.free_bottom)) then
                Material(i_count) = 22
            else if (k==nz_tot-1 .and. j==0 .and. i==0 .and. (.not.free_surface)) then
                Material(i_count) = 23
            else if (k==nz_tot-1 .and. i==nx-1 .and. j==0 .and. (.not.free_surface)) then
                Material(i_count) = 24
            else if (k==nz_tot-1 .and. j==ny-1 .and. i==nx-1 .and. (.not.free_surface)) then
                Material(i_count) = 25
            else if (k==nz_tot-1 .and. i==0 .and. j==ny-1 .and. (.not.free_surface)) then
                Material(i_count) = 26
            endif

        enddo
      enddo
    enddo
endif
allocate (vwgt(0:n_elem-1))
do i = 0,3
    weight(i) = anint(ratio**i)
enddo
do nel = 0, n_elem-1
    if (Material(nel)<1) then
        vwgt(nel) = weight(0)
    else if (Material(nel)<7) then
        vwgt(nel) = weight(1)
    else if (Material(nel)<19) then
        vwgt(nel) = weight(2)
    else
        vwgt(nel) = weight(3)
    endif
enddo


!!! To each element we associate a number which refers to its position in the cover !!!
allocate (Cover(0:n_elem-1));   Cover = 0
if (any_PML) then
    do k = 1,nz_tot-2
      do j = 1,ny-2
        do i = 1,nx-2
            i_count = i + j*nx + k*n_elem_xy 

            if (k<=n_cover .and. (.not.free_bottom)) then
                Cover(i_count) = 1
            else if (k>=nztmp .and. (.not.free_surface)) then
                Cover(i_count) = 6
            else if (j<=n_cover) then
                Cover(i_count) = 2
            else if (j>=nytmp) then
                Cover(i_count) = 4
            else if (i<=n_cover) then
                Cover(i_count) = 5
            else if (i>=nxtmp) then
                Cover(i_count) = 3
            endif

!            if (k<=n_cover .and. j<=n_cover .and. (.not.free_bottom)) then
!                Cover(i_count) = 7
!            else if (k<=n_cover .and. i>=nxtmp .and. (.not.free_bottom)) then
!                Cover(i_count) = 8
!            else if (k<=n_cover .and. j>=nytmp .and. (.not.free_bottom)) then
!                Cover(i_count) = 9
!            else if (k<=n_cover .and. i<=n_cover .and. (.not.free_bottom)) then
!                Cover(i_count) = 10
!            else if (i<=n_cover .and. j<=n_cover) then
!                Cover(i_count) = 11
!            else if (i>=nxtmp .and. j<=n_cover) then
!                Cover(i_count) = 12
!            else if (i>=nxtmp .and. j>=nytmp) then
!                Cover(i_count) = 13
!            else if (i<=n_cover .and. j>=nytmp) then
!                Cover(i_count) = 14
!            else if (k>=nztmp .and. j<=n_cover .and. (.not.free_surface)) then
!                Cover(i_count) = 15
!            else if (k>=nztmp .and. i>=nxtmp .and. (.not.free_surface)) then
!                Cover(i_count) = 16
!            else if (k>=nztmp .and. j>=nytmp .and. (.not.free_surface)) then
!                Cover(i_count) = 17
!            else if (k>=nztmp .and. i<=n_cover .and. (.not.free_surface)) then
!                Cover(i_count) = 18
!            endif
!
!            if (k<=n_cover .and. j<=n_cover .and. i<=n_cover .and. (.not.free_bottom)) then
!                Cover(i_count) = 19
!            else if (k<=n_cover .and. i>=nxtmp .and. j<=n_cover .and. (.not.free_bottom)) then
!                Cover(i_count) = 20
!            else if (k<=n_cover .and. j>=nytmp .and. i>=nxtmp .and. (.not.free_bottom)) then
!                Cover(i_count) = 21
!            else if (k<=n_cover .and. i<=n_cover .and. j>=nytmp .and. (.not.free_bottom)) then
!                Cover(i_count) = 22
!            else if (k>=nztmp .and. j<=n_cover .and. i<=n_cover .and. (.not.free_surface)) then
!                Cover(i_count) = 23
!            else if (k>=nztmp .and. i>=nxtmp .and. j<=n_cover .and. (.not.free_surface)) then
!                Cover(i_count) = 24
!            else if (k>=nztmp .and. j>=nytmp .and. i>=nxtmp .and. (.not.free_surface)) then
!                Cover(i_count) = 25
!            else if (k>=nztmp .and. i<=n_cover .and. j>=nytmp .and. (.not.free_surface)) then
!                Cover(i_count) = 26
!            endif

        enddo
      enddo
    enddo
endif


!!! Partitioning the domain by using METIS !!!
allocate (dxadj(0:n_elem))
allocate (dxadjncy(0:8*n_elem-1))
call METIS_MeshToDual (n_elem, n_points, elmnts, etype, numflag, dxadj, dxadjncy)
allocate (adjwgt(0:dxadj(n_elem)-1));   adjwgt = 0
allocate (part(0:n_elem-1))
if (nparts==1) then
    part(:) = 0
else
    call METIS_PartGraphKway (n_elem, dxadj, dxadjncy(0:dxadj(n_elem)-1), vwgt, adjwgt, &
                              wgtflag, numflag, nparts, options, edgecut, part)
    deallocate (vwgt, adjwgt)
endif
!!! Now there will be two different numberings: !!!
!!! a local one (i.e. relative to a processor) and a global one (i.e. which concerns the whole domain) !!!


!!! Defining the number of elements per processor !!!
allocate (nelem_in_proc(0:nparts-1))
nelem_in_proc = 0
do nel = 0,n_elem-1
    nelem_in_proc(part(nel)) = nelem_in_proc(part(nel)) + 1
enddo


!!! Defining for each processor which elements are inside !!!
!!! "Counter" refers to the local numberings and "nel" to the global numbering !!!
!!! Note that the elements of a processor are naturally sorted according to the global numbering !!!
allocate (counter(0:nparts-1));   counter = 0
allocate (which_elem_in_proc(0:nparts-1,0:maxval(nelem_in_proc)-1));   which_elem_in_proc(:,:) = -1
do nel = 0,n_elem-1
    num = part(nel)
    which_elem_in_proc(num,counter(num)) = nel
    counter(num) = counter(num) + 1
enddo
deallocate (counter)


!!! Allocating "memory" !!!
!!! Allow to establish the correspondence between objects shared by different processors !!!
allocate (memory(0:n_elem-1))
do nel = 0,n_elem-1
    if (part(nel) /= nparts-1) then
        allocate (memory(nel)%rank(part(nel)+1:nparts-1))
        do proc = part(nel)+1,nparts-1
            allocate (memory(nel)%rank(proc)%Obj(0:17))
            memory(nel)%rank(proc)%Obj(:) = -1
        enddo
    endif
enddo


!!! LET'S CONSIDER A PROCESSOR AND CREATE THE FIELDS REQUIRED TO BUILD ITS MESHFILE !!!
call system ("rm -f mesh4spec.???")
meshfilename(1:10) = "mesh4spec."
if(write_vtk)   allocate(coord_grid(0:(nx+1)*(ny+1)*(nz_tot+1)-1,0:2))
do proc = 0,nparts-1

 !!! Defining the vertices which belong to the processor !!!
 !!! These vertices will be sorted according to the global numbering !!!
 allocate (which_vertices(0:8*nelem_in_proc(proc)-1))
 n_vertices = 0
 do n = 0,nelem_in_proc(proc)-1
     nel = which_elem_in_proc(proc,n)
     do i = 0,7
         num = elmnts(8*nel+i)
         ok = 1
         check : do j = 0,n_vertices-1
             if (which_vertices(j)==num) then
                 ok = 0
                 exit check
             endif
         enddo check
         if (ok==1) then
             which_vertices(n_vertices) = num
             n_vertices = n_vertices+1
         endif
     enddo
 enddo
 allocate (tmp1D(0:n_vertices-1))
 tmp1D(0:n_vertices-1) = which_vertices(0:n_vertices-1)
 deallocate (which_vertices)
 allocate (which_vertices(0:n_vertices-1))
 which_vertices(0:n_vertices-1) = tmp1D(0:n_vertices-1)
 deallocate (tmp1D)
 call sort(which_vertices,n_vertices)

 !!! Associating to each element 8 vertices by using a local numbering !!!
 allocate (elmnts_local(0:nelem_in_proc(proc)-1,0:7))
 do n = 0,nelem_in_proc(proc)-1
     nel = which_elem_in_proc(proc,n)
     do i = 0,7
         num = elmnts(8*nel+i)
         build : do j = 0,n_vertices-1
             if (which_vertices(j)==num) then
                 elmnts_local(n,i) = j
                 exit build
             endif
         enddo build
     enddo
 enddo

 !!! Associating to each element 6 faces by using a local numbering !!!
 !!! Defining the faces shared with another processor !!!
 allocate (faces(0:nelem_in_proc(proc)-1,0:5))
 allocate (faces_shared(0:nparts-1,0:6*nelem_in_proc(proc)-1))
 allocate (nf_shared(0:nparts-1));   nf_shared(:) = 0
 allocate (corner(0:3))
 allocate (neighbor_corner(0:3))
 n_faces = 0
 do n = 0,nelem_in_proc(proc)-1   ! n is the number of the considered element in the local numbering
     nel = which_elem_in_proc(proc,n)   ! nel is the number of the considered element in the global numbering
     do nf = 0,5   ! nf indicates which face of the element we're considering
         select case (nf)   ! Here we pick up the vertices (in the global numbering) which define the face
         case (0)
          corner(0) = elmnts(8*nel)
          corner(1) = elmnts(8*nel+1)
          corner(2) = elmnts(8*nel+2)
          corner(3) = elmnts(8*nel+3)
         case (1)
          corner(0) = elmnts(8*nel)
          corner(1) = elmnts(8*nel+1)
          corner(2) = elmnts(8*nel+4)
          corner(3) = elmnts(8*nel+5)
         case (2)
          corner(0) = elmnts(8*nel+1)
          corner(1) = elmnts(8*nel+2)
          corner(2) = elmnts(8*nel+5)
          corner(3) = elmnts(8*nel+6)
         case (3)
          corner(0) = elmnts(8*nel+2)
          corner(1) = elmnts(8*nel+3)
          corner(2) = elmnts(8*nel+6)
          corner(3) = elmnts(8*nel+7)
         case (4)
          corner(0) = elmnts(8*nel)
          corner(1) = elmnts(8*nel+3)
          corner(2) = elmnts(8*nel+4)
          corner(3) = elmnts(8*nel+7)
         case (5)
          corner(0) = elmnts(8*nel+4)
          corner(1) = elmnts(8*nel+5)
          corner(2) = elmnts(8*nel+6)
          corner(3) = elmnts(8*nel+7)
         end select
         find0 : do i = dxadj(nel), dxadj(nel+1)-1   ! We look at a neighbor of the element
             neighbor = dxadjncy(i)
             find1 : do j = 0,3   ! Does the face belong to this neighbor ?
                 num = corner(j)
                 ok = 0
                 find2 : do k = 0,7
                     if (elmnts(8*neighbor+k)==num) then
                         neighbor_corner(j) = k
                         ok = 1
                         exit find2
                     endif
                 enddo find2
                 if (ok==0) exit find1   ! NO, so let's see another neighbor
                 if (j==3) then   ! YES
                     call sort(neighbor_corner,4)
                     if (neighbor_corner(0)==0) then   ! So which face of the neighbor is it ?
                         if (neighbor_corner(3)==3)   neighbor_face = 0
                         if (neighbor_corner(3)==5)   neighbor_face = 1
                         if (neighbor_corner(3)==7)   neighbor_face = 4
                     else if (neighbor_corner(0)==1) then
                         neighbor_face = 2
                     else if (neighbor_corner(0)==2) then
                         neighbor_face = 3
                     else if (neighbor_corner(0)==4) then
                         neighbor_face = 5
                     else
                         print *,"Coherency Pb between faces and vertices of an element"
                     endif
                     if (part(neighbor)==proc) then   ! The neighbor is on the same processor than the element
                         if (neighbor>nel) then   ! The neighbor is an element we've never seen
                             faces(n,nf) = n_faces
                             n_faces = n_faces + 1
                         else   ! The neighbor is an element we've ever seen
                             g2l : do i_count = 0,n-1
                                 if (which_elem_in_proc(proc,i_count)==neighbor) then
                                     faces(n,nf) = faces(i_count,neighbor_face)
                                     exit g2l
                                 endif
                             enddo g2l
                         endif   
                     else   ! The neighbor is not on the same processor than the element
                         faces(n,nf) = n_faces
                         n_faces = n_faces + 1
                         num = part(neighbor)
                         !!! Ensuring the correspondence between the shared faces!!!
                         if (num<proc) then   ! We've ever seen the processor of the neighbor
                             faces_shared(num,memory(neighbor)%rank(proc)%Obj(neighbor_face)) = faces(n,nf)
                         else   ! We've never seen the processor of the neighbor
                             faces_shared(num,nf_shared(num)) = faces(n,nf)
                             memory(nel)%rank(num)%Obj(nf) = nf_shared(num)
                         endif
                         nf_shared(num) = nf_shared(num) + 1
                     endif
                     exit find0
                 endif
             enddo find1
         enddo find0
         if (ok==0) then   ! The face is not shared by a neighbor
             faces(n,nf) = n_faces
             n_faces = n_faces + 1
         endif
     enddo
 enddo
 deallocate (corner, neighbor_corner)
 allocate (tmp2D(0:nparts-1,0:maxval(nf_shared)-1))
 tmp2D(0:nparts-1,0:maxval(nf_shared)-1) = faces_shared(0:nparts-1,0:maxval(nf_shared)-1)
 deallocate (faces_shared)
 allocate (faces_shared(0:nparts-1,0:maxval(nf_shared)-1))
 faces_shared(0:nparts-1,0:maxval(nf_shared)-1) = tmp2D(0:nparts-1,0:maxval(nf_shared)-1)
 deallocate (tmp2D)

 !!! Associating to each element 12 edges by using a local numbering !!!
 !!! Defining the edges shared with others processors !!!
 allocate (edges(0:nelem_in_proc(proc)-1,0:11))
 allocate (edges_shared(0:nparts-1,0:12*nelem_in_proc(proc)-1))
 allocate (ne_shared(0:nparts-1));   ne_shared(:) = 0
 allocate (corner(0:1))
 allocate (neighbor_corner(0:1))
 n_edges = 0
 do n = 0,nelem_in_proc(proc)-1   ! n is the number of the considered element in the local numbering
     nel = which_elem_in_proc(proc,n)   ! nel is the number of the considered element in the global numbering
     !!! Sorting the 26 neighbors of the element (in the global numbering) !!!
     neigh(0) = nel-n_elem_xy-nx-1;   neigh(1) = neigh(0)+1;   neigh(2) = neigh(0)+2
     neigh(3) = nel-n_elem_xy-1;   neigh(4) = neigh(3)+1;   neigh(5) = neigh(3)+2
     neigh(6) = nel-n_elem_xy+nx-1;   neigh(7) = neigh(6)+1;   neigh(8) = neigh(6)+2
     neigh(9) = nel-nx-1;   neigh(10) = neigh(9)+1;   neigh(11) = neigh(9)+2
     neigh(12) = nel-1;   neigh(13) = nel+1
     neigh(14) = nel+nx-1;   neigh(15) = neigh(14)+1;   neigh(16) = neigh(14)+2
     neigh(17) = nel+n_elem_xy-nx-1;   neigh(18) = neigh(17)+1;   neigh(19) = neigh(17)+2
     neigh(20) = nel+n_elem_xy-1;   neigh(21) = neigh(20)+1;   neigh(22) = neigh(20)+2
     neigh(23) = nel+n_elem_xy+nx-1;   neigh(24) = neigh(23)+1;   neigh(25) = neigh(23)+2
     do ne = 0,11   ! ne indicates which edge of the element we're considering
         select case (ne)   ! Here we pick up the vertices (in the global numbering) which define the edge
         case (0)
          corner(0) = elmnts(8*nel)
          corner(1) = elmnts(8*nel+1)
          neigh_local(0:2) = (/ 1, 4, 10 /)
         case (1)
          corner(0) = elmnts(8*nel+1)
          corner(1) = elmnts(8*nel+2)
          neigh_local(0:2) = (/ 4, 5, 13 /)
         case (2)
          corner(0) = elmnts(8*nel+2)
          corner(1) = elmnts(8*nel+3)
          neigh_local(0:2) = (/ 4, 7, 15 /)
         case (3)
          corner(0) = elmnts(8*nel)
          corner(1) = elmnts(8*nel+3)
          neigh_local(0:2) = (/ 3, 4, 12 /)
         case (4)
          corner(0) = elmnts(8*nel+1)
          corner(1) = elmnts(8*nel+5)
          neigh_local(0:2) = (/ 10, 11, 13 /)
         case (5)
          corner(0) = elmnts(8*nel+4)
          corner(1) = elmnts(8*nel+5)
          neigh_local(0:2) = (/ 10, 18, 21 /)
         case (6)
          corner(0) = elmnts(8*nel)
          corner(1) = elmnts(8*nel+4)
          neigh_local(0:2) = (/ 9, 10, 12 /)
         case (7)
          corner(0) = elmnts(8*nel+2)
          corner(1) = elmnts(8*nel+6)
          neigh_local(0:2) = (/ 13, 15, 16 /)
         case (8)
          corner(0) = elmnts(8*nel+5)
          corner(1) = elmnts(8*nel+6)
          neigh_local(0:2) = (/ 13, 21, 22 /)
         case (9)
          corner(0) = elmnts(8*nel+6)
          corner(1) = elmnts(8*nel+7)
          neigh_local(0:2) = (/ 15, 21, 24 /)
         case (10)
          corner(0) = elmnts(8*nel+3)
          corner(1) = elmnts(8*nel+7)
          neigh_local(0:2) = (/ 12, 14, 15 /)
         case (11)
          corner(0) = elmnts(8*nel+4)
          corner(1) = elmnts(8*nel+7)
          neigh_local(0:2) = (/ 12, 20, 21 /)
         end select
         flag = .true.;   n_other_proc = 0
         findbis0 : do i = 0,2   ! We look at a neighbor of the element
             neighbor = neigh(neigh_local(i))
             if (neighbor>=0) then   ! If the neighbor exists...
                 findbis1 : do j = 0,1   ! Does the edge belong to this neighbor ?
                     num = corner(j)
                     ok = 0
                     findbis2 : do k = 0,7
                         if (elmnts(8*neighbor+k)==num) then
                             neighbor_corner(j) = k
                             ok = 1
                             exit findbis2
                         endif
                     enddo findbis2
                     if (ok==0) exit findbis1   ! NO, so let's see another neighbor
                     if (j==1) then   ! YES
                         call sort(neighbor_corner,2)
                         if (neighbor_corner(0)==0) then   ! So which edge of the neighbor is it ?
                          if (neighbor_corner(1)==1) neighbor_edge = 0
                          if (neighbor_corner(1)==3) neighbor_edge = 3
                          if (neighbor_corner(1)==4) neighbor_edge = 6
                         else if (neighbor_corner(0)==1) then
                          if (neighbor_corner(1)==2) neighbor_edge = 1
                          if (neighbor_corner(1)==5) neighbor_edge = 4
                         else if (neighbor_corner(0)==2) then
                          if (neighbor_corner(1)==3) neighbor_edge = 2
                          if (neighbor_corner(1)==6) neighbor_edge = 7
                         else if (neighbor_corner(0)==3) then
                          neighbor_edge = 10
                         else if (neighbor_corner(0)==4) then
                          if (neighbor_corner(1)==5) neighbor_edge = 5
                          if (neighbor_corner(1)==7) neighbor_edge = 11
                         else if (neighbor_corner(0)==5) then
                          neighbor_edge = 8
                         else if (neighbor_corner(0)==6) then
                          neighbor_edge = 9
                         else
                          print *,"Coherency Pb between edges and vertices of an element"
                         endif
                         num = part(neighbor)
                         if (num==proc) then   ! The neighbor is on the same processor than the element
                             if (neighbor<nel) then   ! The neighbor is an element we've ever seen
                                 do i_count = 0,n-1
                                     if (which_elem_in_proc(proc,i_count)==neighbor) then
                                         edges(n,ne) = edges(i_count,neighbor_edge)
                                         flag = .false.
                                         exit findbis0
                                     endif
                                 enddo
                             endif
                         else   ! The neighbor is not on the same processor than the element
                             do k = 0,n_other_proc-1
                                 if (num==other_proc(k,0))   exit findbis1 
                             enddo
                             other_proc(n_other_proc,0) = num
                             other_proc(n_other_proc,1) = neighbor
                             other_proc(n_other_proc,2) = neighbor_edge
                             n_other_proc = n_other_proc + 1 
                         endif
                     endif
                 enddo findbis1
             endif
         enddo findbis0
         if (flag) then   ! We've never seen this edge on this processor before
             edges(n,ne) = n_edges
             n_edges = n_edges + 1
             do i = 0,n_other_proc-1   ! If it is shared with other processors...
                 num = other_proc(i,0)
                 neighbor = other_proc(i,1)
                 neighbor_edge = other_proc(i,2)
                 !!! Ensuring the correspondence between the shared edges !!!
                 if (num<proc) then   ! It deals with a processor we've ever seen
                     edges_shared(num,memory(neighbor)%rank(proc)%Obj(neighbor_edge+6)) = edges(n,ne)
                 else   ! It deals with a processor we've never seen
                     edges_shared(num,ne_shared(num)) = edges(n,ne)
                     memory(nel)%rank(num)%Obj(ne+6) = ne_shared(num)
                 endif
                 ne_shared(num) = ne_shared(num) + 1
             enddo
         endif
     enddo
 enddo
 deallocate (corner, neighbor_corner)
 allocate (tmp2D(0:nparts-1,0:maxval(ne_shared)-1))
 tmp2D(0:nparts-1,0:maxval(ne_shared)-1) = edges_shared(0:nparts-1,0:maxval(ne_shared)-1)
 deallocate (edges_shared)
 allocate (edges_shared(0:nparts-1,0:maxval(ne_shared)-1))
 edges_shared(0:nparts-1,0:maxval(ne_shared)-1) = tmp2D(0:nparts-1,0:maxval(ne_shared)-1)
 deallocate (tmp2D)

 !!! Defining vertices shared with others processors !!!
 allocate (vertices_shared(0:nparts-1,0:7*nelem_in_proc(proc)-1))
 allocate (nv_shared(0:nparts-1));   nv_shared = 0
 allocate (L_Proc(0:nparts-1))
 do i = 0,n_vertices-1   ! i is the number of the considered vertex in the local numbering
     L_Proc = .true.;   L_Proc(proc) = .false.
     n = 0;   j = 0;   ok = 0
     do while (ok==0)   ! Looking for an element which contains the considered vertex
         if (elmnts_local(n,j)==i) then
             ok = 1
         else
             j = j+1
             if (j==8) then
                 n = n+1
                 j = 0
             endif
         endif
     enddo
     nel = which_elem_in_proc(proc,n)   ! nel is the number of the element in the global numbering
     !!! Sorting the 26 neighbors of the element (in the global numbering) !!!
     neigh(0) = nel-n_elem_xy-nx-1;   neigh(1) = neigh(0)+1;   neigh(2) = neigh(0)+2
     neigh(3) = nel-n_elem_xy-1;   neigh(4) = neigh(3)+1;   neigh(5) = neigh(3)+2
     neigh(6) = nel-n_elem_xy+nx-1;   neigh(7) = neigh(6)+1;   neigh(8) = neigh(6)+2
     neigh(9) = nel-nx-1;   neigh(10) = neigh(9)+1;   neigh(11) = neigh(9)+2
     neigh(12) = nel-1;   neigh(13) = nel+1
     neigh(14) = nel+nx-1;   neigh(15) = neigh(14)+1;   neigh(16) = neigh(14)+2
     neigh(17) = nel+n_elem_xy-nx-1;   neigh(18) = neigh(17)+1;   neigh(19) = neigh(17)+2
     neigh(20) = nel+n_elem_xy-1;   neigh(21) = neigh(20)+1;   neigh(22) = neigh(20)+2
     neigh(23) = nel+n_elem_xy+nx-1;   neigh(24) = neigh(23)+1;   neigh(25) = neigh(23)+2
     k = which_vertices(i)
     do j = 0,25   ! We look at a neighbor of the element
         neighbor = neigh(j)
         if (neighbor>=0 .and. neighbor<n_elem) then   ! If the neighbor exists...
             num = part(neighbor)
             if (L_Proc(num)) then
                 findter : do m = 0,7   ! Does the vertex belong to this neighbor ?
                     if (elmnts(8*neighbor+m)==k) then   ! YES
                         !!! The local numbering of the vertices follows the global numbering !!!
                         !!! Using memory to ensure the correspondence between the processors is therefore useless !!!
                         vertices_shared(num,nv_shared(num)) = i
                         nv_shared(num) = nv_shared(num) + 1
                         L_Proc(num) = .false.
                         exit findter
                     endif
                 enddo findter
             endif
         endif
     enddo
 enddo
 deallocate (L_Proc)
 allocate (tmp2D(0:nparts-1,0:maxval(nv_shared)-1))
 tmp2D(0:nparts-1,0:maxval(nv_shared)-1) = vertices_shared(0:nparts-1,0:maxval(nv_shared)-1)
 deallocate (vertices_shared)
 allocate (vertices_shared(0:nparts-1,0:maxval(nv_shared)-1))
 vertices_shared(0:nparts-1,0:maxval(nv_shared)-1) = tmp2D(0:nparts-1,0:maxval(nv_shared)-1)
 deallocate (tmp2D)

 !!! Defining the nodes !!!
 n_nodes = n_vertices
 if (nods_per_elem==20)   n_nodes = n_nodes + n_edges
 if (nods_per_elem==27)   n_nodes = n_nodes + n_edges + n_faces + nelem_in_proc(proc)
 allocate (coord_nodes(0:n_nodes-1,0:2))
 allocate (nodes(0:nelem_in_proc(proc)-1,0:nods_per_elem-1))
 allocate (L_Nodes(0:n_nodes-1));   L_Nodes = .true.
 do n = 0,nelem_in_proc(proc)-1
     do i = 0,7
         nv = elmnts_local(n,i)
         nodes(n,i) = nv
         j = which_vertices(nv)
         coord_vertices(i,0) = xco(j)
         coord_vertices(i,1) = yco(j)
         coord_vertices(i,2) = zco(j)
         if (L_Nodes(nv)) then
             coord_nodes(nv,:) = coord_vertices(i,:)
             L_Nodes(nv) = .false.
         endif
     enddo
     if ((nods_per_elem==20) .or. (nods_per_elem==27)) then
         do i = 0,11
             ne = edges(n,i)
             num = n_vertices + ne
             m = table_node_edge(i)
             nodes(n,m) = num
             if (L_Nodes(num)) then
                 select case (m)
                 case(8)
                     j = 0;   k = 1
                 case (9)
                     j = 1;   k = 2
                 case (10)
                     j = 2;   k = 3
                 case (11)
                     j = 0;   k = 3
                 case (12)
                     j = 0;   k = 4
                 case (13)
                     j = 1;   k = 5
                 case (14)
                     j = 2;   k = 6
                 case (15)
                     j = 3;   k = 7
                 case (16)
                     j = 4;   k = 5
                 case (17)
                     j = 5;   k = 6
                 case (18)
                     j = 6;   k = 7
                 case (19)
                     j = 4;   k = 7
                 end select
                 coord_nodes(num,:) = (coord_vertices(j,:)+coord_vertices(k,:))/2.d0
                 L_Nodes(num) = .false.
             endif
         enddo
         if (nods_per_elem==27) then
             do i = 0,5
                 nf = faces(n,i)
                 num = n_vertices + n_edges + nf
                 m = 20+i
                 nodes(n,m) = num
                 if (L_Nodes(num)) then
                     select case (m)
                     case (20)
                         j = 0;   k = 2
                     case (21)
                         j = 0;   k = 5
                     case (22)
                         j = 1;   k = 6
                     case (23)
                         j = 2;   k = 7
                     case (24)
                         j = 0;   k = 7
                     case (25)
                         j = 4;   k = 6
                     end select
                     coord_nodes(num,:) = (coord_vertices(j,:)+coord_vertices(k,:))/2.d0
                     L_Nodes(num) = .false.
                 endif
             enddo
             num = n_vertices + n_edges + n_faces + n
             nodes(n,26) = num
             coord_nodes(num,:) = (coord_vertices(0,:)+coord_vertices(6,:))/2.d0
         endif
     endif
 enddo
 deallocate (L_Nodes)

 !!! Calculating the Cartesian coordinates of each node !!!
 if (curve) then
     if (ellipticity)   allocate(Rsph(0:n_nodes-1))
     do n = 0,n_nodes-1
         R = coord_nodes(n,2)
         Y = dtan(deg2rad(coord_nodes(n,1)))
         X = dtan(deg2rad(coord_nodes(n,0)))
         D = dsqrt(1.d0 + Y**2 + X**2)
         if (ellipticity) then
             ! Passage en coordonnees cartesiennes
             xa = X/D;   ya = Y/D;   za = 1/D
             ! Rotation du chunk de ref vers le chunk reel
             xs = rot(0,0)*xa + rot(0,1)*ya + rot(0,2)*za
             ys = rot(1,0)*xa + rot(1,1)*ya + rot(1,2)*za
             zs = rot(2,0)*xa + rot(2,1)*ya + rot(2,2)*za
             ! Passage en spherique
             call cart2sph(xs,ys,zs,alpha,theta,phi)
             ! Calcul des rayons elliptiques
             P2 = (3.d0*dcos(theta)**2 - 1.d0) / 2.d0
             alpha = get_ellipticity(radius(n_layers-2))
             rad_loc(0) = radius(n_layers-2) * (1.d0 - 2.d0/3.d0 * alpha * P2)
             alpha = get_ellipticity(radius(n_layers-1))
             rad_loc(1) = radius(n_layers-1) * (1.d0 - 2.d0/3.d0 * alpha * P2)
             alpha = get_ellipticity(Rterre)
             rad_loc(2) = Rterre * (1.d0 - 2.d0/3.d0 * alpha * P2)
             alpha = get_ellipticity(R)
             rad_loc(3) = R * (1.d0 - 2.d0/3.d0 * alpha * P2)
         endif
         ! When a Moho and/or a surface topography is required, R is changed at some nodes
         if ((model==2) .or. (model==3)) then
             if (R>limit_inf) then
                 epsil_xi = dx/1000.d0
                 dist = dx/2.d0
                 i = 0
                 do while (egal(coord_nodes(n,0),x_len+i*dist,epsil_xi) .eqv. .false.)
                     i = i+1
                 enddo
                 epsil_eta = dy/1000.d0
                 dist = dy/2.d0
                 j = 0
                 do while (egal(coord_nodes(n,1),y_len+j*dist,epsil_eta) .eqv. .false.)
                     j = j+1
                 enddo
                 if (R<limit_med) then
                     epsil_r = dR(n_layers-2)/1000.d0
                     dist = dR(n_layers-2)/2.d0
                     k = 1
                     do while (egal(R,radius(n_layers-2)+k*dist,epsil_r) .eqv. .false.)
                         k = k+1
                     enddo
                     dist = (Rterre-radius(n_layers-2)-depth_moho(i,j)) / (2.d0*nz(n_layers-2))
                     R = radius(n_layers-2) + k*dist
                     if (ellipticity) then
                         Rsph(n) = R
                         dist = (rad_loc(2)-rad_loc(0)-depth_moho(i,j)) / (2.d0*nz(n_layers-2))
                         R = rad_loc(0) + k*dist
                     endif
                     if (2.d0*dist<dzmin(0,proc))   dzmin(0,proc) = 2.d0*dist
                     if (2.d0*dist>dzmax(0,proc))   dzmax(0,proc) = 2.d0*dist
                 else
                     epsil_r = dR(n_layers-1)/1000.d0
                     dist = dR(n_layers-1)/2.d0
                     k = 0
                     do while (egal(R,radius(n_layers-1)+k*dist,epsil_r) .eqv. .false.)
                         k = k+1
                     enddo
                     if (topo_log) then
                         dist = (depth_moho(i,j) + altitude_topo(i,j)) / (2.d0*nz(n_layers-1))
                     else
                         dist = depth_moho(i,j) / (2.d0*nz(n_layers-1))
                     endif
                     R = Rterre-depth_moho(i,j) + k*dist
                     if (ellipticity) then
                         Rsph(n) = R
                         R = rad_loc(2)-depth_moho(i,j) + k*dist
                     endif
                     if (2.d0*dist<dzmin(1,proc))   dzmin(1,proc) = 2.d0*dist
                     if (2.d0*dist>dzmax(1,proc))   dzmax(1,proc) = 2.d0*dist
                 endif
             else
                 if (ellipticity) then
                     Rsph(n) = R
                     R = rad_loc(3)
                 endif
             endif
         else
             if (topo_log) then
                 if (R>limit_med) then
                     epsil_xi = dx/1000.d0
                     dist = dx/2.d0
                     i = 0
                     do while (egal(coord_nodes(n,0),x_len+i*dist,epsil_xi) .eqv. .false.)
                         i = i+1
                     enddo
                     epsil_eta = dy/1000.d0
                     dist = dy/2.d0
                     j = 0
                     do while (egal(coord_nodes(n,1),y_len+j*dist,epsil_eta) .eqv. .false.)
                         j = j+1
                     enddo
                     epsil_r = dR(n_layers-1)/1000.d0
                     dist = dR(n_layers-1)/2.d0
                     k = 0
                     do while (egal(R,radius(n_layers-1)+k*dist,epsil_r) .eqv. .false.)
                         k = k+1
                     enddo
                     dist = (Rterre-radius(n_layers-1)+altitude_topo(i,j)) / (2.d0*nz(n_layers-1))
                     R = radius(n_layers-1) + k*dist
                     if (ellipticity) then
                         Rsph(n) = R
                         dist = (rad_loc(2)-rad_loc(1)+altitude_topo(i,j)) / (2.d0*nz(n_layers-1))
                         R = rad_loc(1) + k*dist
                     endif
                 else
                     if (ellipticity) then
                         Rsph(n) = R
                         R = rad_loc(3)
                     endif
                 endif
             else
                 if (ellipticity) then
                     Rsph(n) = R
                     R = rad_loc(3)
                 endif
             endif
         endif
         coord_nodes(n,0) = R*X/D
         coord_nodes(n,1) = R*Y/D
         coord_nodes(n,2) = R/D
     enddo
 endif

 if(write_vtk)then
     do n = 0,n_vertices-1
         coord_grid(which_vertices(n),:) = coord_nodes(n,:)/1000.d0   ! meters to kilometers
     enddo
 endif

 !!! Writing the meshfile !!!
 write (meshfilename(11:13), '(i3.3)') proc
 open (11, file = trim(meshfilename))
 write (11,"(1i3)") n_dim
 write (11,*)
 write (11,*) nx,ny,nz_tot
 write (11,*)
 write (11,*) i1,i2,j1,j2,k1,k2
 write (11,*)
 do n = 1,6
    write (11,*) cover_def(n,:)
 enddo
 write (11,*)
 write (11,"(l3,i3,l3)") curve, model, ellipticity
 if (curve) then
    do i = 0,2
       write (11,*) (rot(i,j),j=0,2)
    enddo
 endif
 write (11,*)
 write (11,"(1i7)") n_nodes
 do n = 0,n_nodes-1
     write (11,*) coord_nodes(n,:)
     if (ellipticity)   write (11,*) Rsph(n)
 enddo 
 write (11,*)
 write (11,"(1i7)") nelem_in_proc(proc)
 write (11,*)
 do n = 0,nelem_in_proc(proc)-1
     nel = which_elem_in_proc(proc,n)
     write (11,"(4i4,3i7)") Material(nel), Cover(nel), moho_position(nel), topo_position(nel), &
                            elmnts_indexes(nel,0),elmnts_indexes(nel,1),elmnts_indexes(nel,2)
 enddo
 write (11,*)
 write (11,"(1i3)") nods_per_elem
 if (nods_per_elem==8) then
     do n = 0,nelem_in_proc(proc)-1
         write (11,"(8i8)") (nodes(n,i),i=0,nods_per_elem-1)
     enddo
 else if (nods_per_elem==20) then
     do n = 0,nelem_in_proc(proc)-1
         write (11,"(20i8)") (nodes(n,i),i=0,nods_per_elem-1)
     enddo
 else if (nods_per_elem==27) then
     do n = 0,nelem_in_proc(proc)-1
         write (11,"(27i8)") (nodes(n,i),i=0,nods_per_elem-1)
     enddo
 else
     print *,"Bad number of nodes"
     stop
 endif
 write (11,*)
 write (11,"(1i7)") n_faces
 do n = 0,nelem_in_proc(proc)-1
     write (11,"(6i8)") (faces(n,i),i=0,5)
 enddo
 write (11,*)
 write (11,"(1i7)") n_edges
 do n = 0,nelem_in_proc(proc)-1
     write (11,"(12i8)") (edges(n,i),i=0,11)
 enddo
 write (11,*)
 write (11,"(1i7)") n_vertices
 do n = 0,nelem_in_proc(proc)-1
     write (11,"(8i8)") (elmnts_local(n,i),i=0,7)
 enddo
 write (11,*)
 write (11,"(1i6)") nparts
 do n = 0,nparts-1
     write (11,"(3i8)") nf_shared(n), ne_shared(n), nv_shared(n)
     do nf = 0,nf_shared(n)-1
         write (11,"(1i7)") faces_shared(n,nf)
     enddo
     do ne = 0,ne_shared(n)-1
         write (11,"(1i7)") edges_shared(n,ne)
     enddo
     do nv = 0,nv_shared(n)-1
         write (11,"(1i7)") vertices_shared(n,nv)
     enddo
 enddo
 close (11)

 !!! Writing chunk.out !!!
 if (output_chunk) then
     if (proc==0) then
         open(30,file="chunk_0.out")
         open(31,file="chunk_1.out")
         open(32,file="chunk_2.out")
         open(33,file="chunk_3.out")
         open(34,file="chunk_4.out")
         open(35,file="chunk_5.out")
     endif
     do k = 0,nz_tot-1
      do j = 0,ny-1
       do i = 0,nx-1
           if (k==0) then
               i_count = i + j*nx + k*n_elem_xy
               n_out = -1
               call find_numloc (nelem_in_proc(proc),which_elem_in_proc(proc,:),i_count,n_out)
               if (n_out>-1)   call write_mesh (nods_per_elem,coord_nodes(nodes(n_out,:),:),0)
           endif
           if (k==nz_tot-1) then
               i_count = i + j*nx + k*n_elem_xy
               n_out = -1
               call find_numloc (nelem_in_proc(proc),which_elem_in_proc(proc,:),i_count,n_out)
               if (n_out>-1)   call write_mesh (nods_per_elem,coord_nodes(nodes(n_out,:),:),5)
           endif
           if (j==0) then
               i_count = i + j*nx + k*n_elem_xy
               n_out = -1
               call find_numloc (nelem_in_proc(proc),which_elem_in_proc(proc,:),i_count,n_out)
               if (n_out>-1)   call write_mesh (nods_per_elem,coord_nodes(nodes(n_out,:),:),1)
           endif
           if (j==ny-1) then
               i_count = i + j*nx + k*n_elem_xy
               n_out = -1
               call find_numloc (nelem_in_proc(proc),which_elem_in_proc(proc,:),i_count,n_out)
               if (n_out>-1)   call write_mesh (nods_per_elem,coord_nodes(nodes(n_out,:),:),3)
           endif
           if (i==0) then
               i_count = i + j*nx + k*n_elem_xy
               n_out = -1
               call find_numloc (nelem_in_proc(proc),which_elem_in_proc(proc,:),i_count,n_out)
               if (n_out>-1)   call write_mesh (nods_per_elem,coord_nodes(nodes(n_out,:),:),4)
           endif
           if (i==nx-1) then
               i_count = i + j*nx + k*n_elem_xy
               n_out = -1
               call find_numloc (nelem_in_proc(proc),which_elem_in_proc(proc,:),i_count,n_out)
               if (n_out>-1)   call write_mesh (nods_per_elem,coord_nodes(nodes(n_out,:),:),2)
           endif
       enddo
      enddo
     enddo
     deallocate (nodes)
 endif

 !!! Writing the file for random forces !!!
 if (random_trac) then
     do n = 0,nelem_in_proc(proc)-1
         nel = which_elem_in_proc(proc,n)
         if (ellipticity) then
             call table_random_t (proc, n, Material(nel), coord_nodes(nodes(n,25),:), &
                                  moho_position(nel), topo_position(nel), &
                                  rot, Rterre, pi, ellipticity, Rsph(nodes(n,25)))
         else
             call table_random_t (proc, n, Material(nel), coord_nodes(nodes(n,25),:), &
                                  moho_position(nel), topo_position(nel), &
                                  rot, Rterre, pi, ellipticity)
         endif
     enddo 
 endif

 !!! Deallocation !!!
 deallocate (which_vertices, elmnts_local, faces, edges, coord_nodes, nodes)
 deallocate (nf_shared, ne_shared, nv_shared, faces_shared, edges_shared, vertices_shared)
 if (ellipticity)   deallocate(Rsph)
enddo


!!! write vtk grid
if(write_vtk)   call WRITE_VTS_GRID(nx+1,ny+1,nz_tot+1,coord_grid(:,0),coord_grid(:,1),coord_grid(:,2))


!!! Closing files !!!
if (output_chunk) then
    close (30)
    close (31)
    close (32)
    close (33)
    close (34)
    close (35)
endif
if (random_trac)   close (25)


!!! Returning the smallest and the largest vertical size among the elements above and below the Moho !!!
if ((model==2) .or. (model==3)) then
    do proc = 1,nparts-1
        if (dzmin(0,proc)<dzmin(0,0))   dzmin(0,0)=dzmin(0,proc)
        if (dzmin(1,proc)<dzmin(1,0))   dzmin(1,0)=dzmin(1,proc)
        if (dzmax(0,proc)>dzmax(0,0))   dzmax(0,0)=dzmax(0,proc)
        if (dzmax(1,proc)>dzmax(1,0))   dzmax(1,0)=dzmax(1,proc)
    enddo
    write(*,"(a,f5.1,a)") "The smallest vertical size in the crust is", dzmin(1,0)/1000, " km."
    write(*,"(a,f6.1,a)") "The smallest vertical size in the very upper mantle is", dzmin(0,0)/1000, " km."
    if (dzmax(1,0)>dRmax) then
        write(*,"(a,f5.1,a)") "WARNING: The largest element in the crust is larger than &
                               the horizontal step (dzmax=", dzmax(1,0)/1000, " km)."
    endif
    if (dzmax(0,0)>dRmax*vs(nb_couches_model-2)/vs(nb_couches_model-1)) then
        write(*,"(a,f6.1,a)") "WARNING: The largest element in the very upper mantle is larger than &
                               the horizontal step (dzmax=", dzmax(0,0)/1000, " km)."
    endif
    deallocate (dzmin, dzmax, depth_moho)
endif


!!! Deallocation !!!
do nel = 0,n_elem-1
    if (part(nel) /= nparts-1) then
        do proc = part(nel)+1,nparts-1
            deallocate (memory(nel)%rank(proc)%Obj)
        enddo
        deallocate (memory(nel)%rank)
    endif
enddo
if(write_vtk)   deallocate (coord_grid)
deallocate (elmnts_indexes)
deallocate (xco, yco, zco, moho_position, topo_position)
deallocate (elmnts, part, dxadj, dxadjncy)
deallocate (Material, Cover, nelem_in_proc, which_elem_in_proc, memory)
if (curve)   deallocate (radius, nz, dR)


contains


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
logical function egal(a,b,epsil)

doubleprecision :: a,b,epsil

egal = .false.
if (abs(a-b) < epsil) egal = .true.

end function


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
integer function table_node_edge(edge_num)

integer :: edge_num

select case (edge_num)
case (0)
    table_node_edge = 8
case (1)
    table_node_edge = 9
case (2)
    table_node_edge = 10
case (3)
    table_node_edge = 11
case (4)
    table_node_edge = 13
case (5)
    table_node_edge = 16
case (6)
    table_node_edge = 12
case (7)
    table_node_edge = 14
case (8)
    table_node_edge = 17
case (9)
    table_node_edge = 18
case (10)
    table_node_edge = 15
case (11)
    table_node_edge = 19
end select

end function


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
doubleprecision function deg2rad(val)

doubleprecision :: val

deg2rad = pi*val/180.d0

end function


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine bottom_remarks(rayon,rad,i,dx,Rterre,model)

implicit none

doubleprecision, intent(IN) :: dx, Rterre
doubleprecision, intent(INOUT) :: rayon
doubleprecision, dimension(0:11), intent(IN) :: rad
integer, intent(IN) :: i, model
doubleprecision :: dr
integer :: ok

dr = deg2rad(dx*Rterre)
ok = 0
if (i<2) then
    print *,"You are under the CMB. The elements of the core will be PML."
    rayon = rad(1)-dr
    ok = 1
else if (rad(i)-rayon<dr) then
    if (rad(i)-dr<rad(i-1)) then
        rayon = rad(i-1)
    else
        rayon = rad(i)-dr
    endif
    ok = 1
endif

if ((model==2 .or. model==3) .and. (rayon>6291000.d0)) then
    print *,"If you want the Moho, the bottom has to be deeper than 80 km."
    stop
endif

if (ok==1) write(*,"(a,f8.0)") "Optimization of the bottom: depth = ", Rterre-rayon

end subroutine bottom_remarks


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine sort(vector, length)

implicit none

integer, intent(IN) :: length
integer, dimension(0:length-1), intent(INOUT) :: vector
integer :: n, ok, i, tmp

n = length-1
ok = 1
do while (ok==1)
    ok = 0
    do i = 0,n-1
        if (vector(i)>vector(i+1)) then
            tmp = vector(i+1)
            vector(i+1) = vector(i)
            vector(i) = tmp
            ok = 1
        endif
    enddo
    n = n-1
enddo

end subroutine sort


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine define_model(num,nb_couches_model,rad,vs)

implicit none

integer, intent(IN) :: num
integer, intent(OUT) :: nb_couches_model
doubleprecision, dimension(0:11), intent(OUT) :: rad, vs

rad = -1.d0;   vs = -1.d0

if (num==1) then

    nb_couches_model = 12

    rad(0)=1221500; rad(1)=3480000; rad(2)=3630000; rad(3)=5600000
    rad(4)=5701000; rad(5)=5771000; rad(6)=5971000; rad(7)=6151000
    rad(8)=6291000; rad(9)=6346600; rad(10)=6356000

    vs(0)=3054.31; vs(1)=0.00; vs(2)=7265.97; vs(3)=6240.39
    vs(4)=5945.13; vs(5)=5516.02; vs(6)=4932.49; vs(7)=4643.90
    vs(8)=4470.52; vs(9)=4491.01; vs(10)=3900.00

else if (num==2) then

    nb_couches_model = 10

    rad(0)=1221500; rad(1)=3480000; rad(2)=3630000; rad(3)=5600000
    rad(4)=5701000; rad(5)=5771000; rad(6)=5971000; rad(7)=6151000
    rad(8)=6346600

    vs(0)=3054.31; vs(1)=0.00; vs(2)=7265.97; vs(3)=6240.39
    vs(4)=5945.13; vs(5)=5516.02; vs(6)=4932.49; vs(7)=4643.90
    vs(8)=4491.01

else if (num==3) then

    nb_couches_model = 9

    rad(0)=1221500; rad(1)=3480000; rad(2)=3630000; rad(3)=5600000
    rad(4)=5701000; rad(5)=5771000; rad(6)=5971000; rad(7)=6346600

    vs(0)=3054.31; vs(1)=0.00; vs(2)=7265.97; vs(3)=6240.39
    vs(4)=5945.13; vs(5)=5516.02; vs(6)=4932.49; vs(7)=4491.01

else if (num==5) then

    nb_couches_model = 7

    rad(0)=1221500; rad(1)=3480000; rad(2)=3630000; rad(3)=5701000
    rad(4)=5971000; rad(5)=6311000

    vs(0)=3054.31; vs(1)=0.00; vs(2)=7265.97; vs(3)=5920.61
    vs(4)=4877.09; vs(5)=4556.55

else if (num==6) then

    nb_couches_model = 11

    rad(0)=1221500; rad(1)=3480000; rad(2)=3630000; rad(3)=5600000
    rad(4)=5701000; rad(5)=5771000; rad(6)=5971000; rad(7)=6151000
    rad(8)=6291000; rad(9)=6346600

    vs(0)=3054.31; vs(1)=0.00; vs(2)=7265.97; vs(3)=6240.39
    vs(4)=5945.13; vs(5)=5516.02; vs(6)=4932.49; vs(7)=4643.90
    vs(8)=4470.52; vs(9)=4491.01

else

    print *,"ERROR: ",num, " is a wrong input"
    stop

endif

end subroutine define_model


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine read_moho(theta_deg,phi_deg,moho,m,n,latmin,latmax,dlat,lonmin,lonmax,dlong,depth)

implicit none

integer, intent(IN) :: m,n 
doubleprecision, intent(IN) :: latmin,latmax,dlat,lonmin,lonmax,dlong
doubleprecision, intent(OUT) :: depth
doubleprecision, intent(IN) :: theta_deg,phi_deg
doubleprecision, dimension(0:m-1,0:n-1), intent(IN) :: moho
integer :: i,j, k,l,kk,ll, u,v,uu,vv, ok, u1,u2,v1,v2
doubleprecision :: colatmax, phi_loc, sigma, norm, epsil, Wlong,Wlat
doubleprecision, dimension(0:20) :: theta_loc, pt_inter
doubleprecision, dimension(0:20,0:20) :: pt

epsil = dlong/1000.d0

j = 0
do while (lonmin+j*dlong < phi_deg-dlong/2.d0)
    j = j + 1
enddo
if (j<2 .or. j>n-3) then
    if (dabs((lonmax+dlong-360.d0)-lonmin) > epsil) then
        write (*,'(a,f6.1)') "PAS ASSEZ DE DONNEES POUR DETERMINER LA TOPO A LA LONGITUDE", phi_deg
        stop
    endif
endif
if (j==n) j = 0

colatmax = 90.d0 - latmin
i = 0
do while (colatmax-i*dlat > theta_deg+dlat/2.d0)
    i = i + 1
enddo
ok = 0
if (i<2 .or. i>m-3) then
    if (phi_deg >= lonmin+(n-1)*dlong) then
        v1 = n-1;   v2 = 0
        Wlong = (phi_deg - (lonmin+(n-1)*dlong)) / dlong
    else if (phi_deg<=lonmin) then
        v1 = n-1;   v2 = 0
        Wlong = (dlong - (lonmin-phi_deg)) / dlong
    else
        if (lonmin+j*dlong >= phi_deg) then
            v1 = j-1;   v2 = j
        else
            v1 = j;   v2 = j+1
        endif
        Wlong = (phi_deg - (lonmin+v1*dlong)) / dlong
    endif
    if (i==1) then
        if (colatmax-i*dlat >= theta_deg) then
            u1 = 1;   u2 = 2
        else
            u1 = 0;   u2 = 1
        endif
        Wlat = ((colatmax-u1*dlat) - theta_deg) / dlat
    else if (i==0) then
        if (colatmax-i*dlat >= theta_deg) then
            u1 = 0;   u2 = 1
            Wlat = ((colatmax-u1*dlat) - theta_deg) / dlat
        else
            if (colatmax+dlat > 180.d0) then
                u1 = 0;   u2 = 0
                Wlat = 0.d0
            else
                write (*,'(a,f5.1)') "PAS ASSEZ DE DONNEES POUR DETERMINER LA TOPO A LA COLATITUDE ", theta_deg
                stop
            endif
        endif
    else if (i==m-2) then
        if (colatmax-i*dlat <= theta_deg) then
            u1 = m-3;   u2 = m-2
        else
            u1 = m-2;   u2 = m-1
        endif
        Wlat = ((colatmax-u1*dlat) - theta_deg) / dlat
    else if (i==m-1) then
        if (colatmax-i*dlat <= theta_deg) then
            u1 = m-2;   u2 = m-1
            Wlat = ((colatmax-u1*dlat) - theta_deg) / dlat
        else
            if (colatmax-m*dlat < 0.d0) then
                u1 = m-1;   u2 = m-1
                Wlat = 0.d0
            else
                write (*,'(a,f5.1)') "PAS ASSEZ DE DONNEES POUR DETERMINER LA TOPO A LA COLATITUDE ", theta_deg
                stop
            endif
        endif
    endif
    depth = (1-Wlat) * ((1-Wlong)*moho(u1,v1) + Wlong*moho(u1,v2)) &
            + Wlat * ((1-Wlong)*moho(u2,v1) + Wlong*moho(u2,v2))
    ok = 1
endif


if (ok==0) then

 if (j-2<0) then ! j=0 ou j=1
    pt(0,0) = moho(i-2,n+(j-2))
    pt(5,0) = moho(i-1,n+(j-2))
    pt(10,0) = moho(i,n+(j-2))
    pt(15,0) = moho(i+1,n+(j-2))
    pt(20,0) = moho(i+2,n+(j-2))
 else ! j>1
    pt(0,0) = moho(i-2,j-2)
    pt(5,0) = moho(i-1,j-2)
    pt(10,0) = moho(i,j-2)
    pt(15,0) = moho(i+1,j-2)
    pt(20,0) = moho(i+2,j-2)
 endif

 if (j-1<0) then ! j=0
    pt(0,5) = moho(i-2,n-1)
    pt(5,5) = moho(i-1,n-1)
    pt(10,5) = moho(i,n-1)
    pt(15,5) = moho(i+1,n-1)
    pt(20,5) = moho(i+2,n-1)
 else ! j>0
    pt(0,5) = moho(i-2,j-1)
    pt(5,5) = moho(i-1,j-1)
    pt(10,5) = moho(i,j-1)
    pt(15,5) = moho(i+1,j-1)
    pt(20,5) = moho(i+2,j-1)
 endif

 pt(0,10) = moho(i-2,j)
 pt(5,10) = moho(i-1,j)
 pt(10,10) = moho(i,j)
 pt(15,10) = moho(i+1,j)
 pt(20,10) = moho(i+2,j)

 if (j+1>n-1) then ! j=n-1
    pt(0,15) = moho(i-2,0)
    pt(5,15) = moho(i-1,0)
    pt(10,15) = moho(i,0)
    pt(15,15) = moho(i+1,0)
    pt(20,15) = moho(i+2,0)
 else ! j<n-1
    pt(0,15) = moho(i-2,j+1)
    pt(5,15) = moho(i-1,j+1)
    pt(10,15) = moho(i,j+1)
    pt(15,15) = moho(i+1,j+1)
    pt(20,15) = moho(i+2,j+1)
 endif

 if (j+2>n-1) then ! j=n-1 ou j=n-2
    pt(0,20) = moho(i-2,(j+2)-n)
    pt(5,20) = moho(i-1,(j+2)-n)
    pt(10,20) = moho(i,(j+2)-n)
    pt(15,20) = moho(i+1,(j+2)-n)
    pt(20,20) = moho(i+2,(j+2)-n)
 else ! j<n-2
    pt(0,20) = moho(i-2,j+2)
    pt(5,20) = moho(i-1,j+2)
    pt(10,20) = moho(i,j+2)
    pt(15,20) = moho(i+1,j+2)
    pt(20,20) = moho(i+2,j+2)
 endif

 do l = 0,4
     u = 5*l
     do k = 0,4
         v = 5*k
         do ll = 0,4
             uu = u + (ll-2)
             if (uu>=0) then
                 do kk = 0,4
                     vv = v + (kk-2)
                     if (vv>=0) then
                         pt(uu,vv) = pt(u,v)
                     endif
                  enddo
             endif
         enddo
     enddo
 enddo

 sigma = dlong/2.d0
 if (j==0 .and. phi_deg>(n-1)*dlong)   j = n
 do l = 0,20
     theta_loc(l) = colatmax - (i*dlat + (l-10)*dlat/5.d0)
     do k = 0,20
         phi_loc = j*dlong + (k-10)*dlong/5.d0
         pt(l,k) = pt(l,k) * dexp(- (phi_loc-phi_deg)**2 / (2.d0*sigma**2))
     enddo
 enddo
 pt_inter(:) = 0
 do l = 0,20
     do k = 0,9
         pt_inter(l) = pt_inter(l) + (pt(l,2*k) + 4.d0*pt(l,2*k+1) + pt(l,2*k+2))
     enddo
     norm = dlong / (15.d0*sigma*dsqrt(2.d0*pi))
     pt_inter(l) = pt_inter(l) * norm
 enddo

 sigma = dlat/2.d0
 do l = 0,20
     pt_inter(l) = pt_inter(l) * dexp(- (theta_loc(l)-theta_deg)**2 / (2.d0*sigma**2))
 enddo
 depth = 0
 do k = 0,9
     depth = depth + (pt_inter(2*k) + 4.d0*pt_inter(2*k+1) + pt_inter(2*k+2))
 enddo
 norm = dlat / (15.d0*sigma*dsqrt(2.d0*pi))
 depth = depth * norm

endif

end subroutine read_moho


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE WRITE_VTS_GRID(nx,ny,nz,x,y,z)

integer :: E_IO,UNIT_VTK
integer :: nx,ny,nz,offset
real :: x(nx*ny*nz),y(nx*ny*nz),z(nx*ny*nz)
character(1), parameter:: end_rec = char(10)
character :: chnx*10,chny*10,chnz*10,choffset*10

write(chnx,'(I10)')nx-1
write(chny,'(I10)')ny-1
write(chnz,'(I10)')nz-1

UNIT_VTK = 40

open(unit       = Unit_VTK,        &
     file       = 'elements_grid.vts',      &
     form       = 'UNFORMATTED',   &
     access     = 'STREAM',        &
     action     = 'WRITE',         &
     status     = 'REPLACE',       &
     convert    = 'LITTLE_ENDIAN', &
     iostat     = E_IO)

offset = 0

write(unit=Unit_VTK,iostat=E_IO)'<?xml version = "1.0"?>'//end_rec
write(unit=Unit_VTK,iostat=E_IO)'<VTKFile type = "StructuredGrid" version="0.1" byte_order="LittleEndian">'//end_rec
write(unit=Unit_VTK,iostat=E_IO)'<StructuredGrid WholeExtent="'//'0 '//trim(adjustl(chnx))//' 0 '//trim(adjustl(chny))&
&//' 0 '//trim(adjustl(chnz))//'">'//end_rec
write(unit=Unit_VTK,iostat=E_IO)'<Piece Extent="'//'0 '//trim(adjustl(chnx))//' 0 '//trim(adjustl(chny))&
&//' 0 '//trim(adjustl(chnz))//'">'//end_rec
write(unit=Unit_VTK,iostat=E_IO)'<PointData>'//end_rec
write(unit=Unit_VTK,iostat=E_IO)'</PointData>'//end_rec
write(unit=Unit_VTK,iostat=E_IO)'<CellData>'//end_rec
write(unit=Unit_VTK,iostat=E_IO)'</CellData>'//end_rec
write(unit=Unit_VTK,iostat=E_IO)'<Points>'//end_rec
write(choffset,'(I10)')offset
write(unit=Unit_VTK,iostat=E_IO)'<DataArray type="Float32" NumberOfComponents="3"'&
&,' format="appended" offset="'//trim(adjustl(choffset))//'" />'//end_rec
write(unit=Unit_VTK,iostat=E_IO)'</Points>'//end_rec
write(unit=Unit_VTK,iostat=E_IO)'</Piece>'//end_rec
write(unit=Unit_VTK,iostat=E_IO)'</StructuredGrid>'//end_rec
write(unit=Unit_VTK,iostat=E_IO)'<AppendedData encoding="raw">'//end_rec
write(unit=Unit_VTK,iostat=E_IO)'_'
write(unit=Unit_VTK,iostat=E_IO)int(sizeof(x)+sizeof(y)+sizeof(z))
do n = 1,nx*ny*nz
 write(unit=Unit_VTK,iostat=E_IO)x(n),y(n),z(n)
enddo
write(unit=Unit_VTK,iostat=E_IO)end_rec
write(unit=Unit_VTK,iostat=E_IO)'</AppendedData>'//end_rec
write(unit=Unit_VTK,iostat=E_IO)'</VTKFile>'//end_rec

close(Unit=Unit_VTK)

RETURN

END SUBROUTINE WRITE_VTS_GRID


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
end program mesher
