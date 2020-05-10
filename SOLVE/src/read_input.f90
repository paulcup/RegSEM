subroutine read_input (Tdomain, rg)


use sdomains

implicit none

type (domain), intent (INOUT) :: Tdomain
integer, intent(IN) :: rg

logical, dimension(:), allocatable :: L_Face, L_Edge
integer :: length, i, j, mat, ok, nf, ne


! Read simulation features

open (11,file="input.spec",form="formatted",status="old")
read (11,*) Tdomain%sTimeParam%duration
read (11,*) Tdomain%mesh_file
length = len_trim(Tdomain%mesh_file) + 1
write (Tdomain%mesh_file(length:length+2),'(i3.3)') rg
read (11,*) Tdomain%material_file
read (11,*) Tdomain%aniso
read (11,*) Tdomain%ocean

read (11,*) Tdomain%save_trace
if (Tdomain%save_trace) then
   read (11,*) Tdomain%station_file
   read (11,*) Tdomain%comp_rot
   read (11,*) Tdomain%sTimeParam%samp_period
   Tdomain%sTimeParam%t_resamp = Tdomain%sTimeParam%samp_period
else 
   read (11,*)
   read (11,*)
   read (11,*)
endif

read (11,*) Tdomain%save_snapshot
if (Tdomain%save_snapshot) then
   read (11,*) Tdomain%sTimeParam%dnsnap
   read (11,*) Tdomain%nx_interp, Tdomain%ny_interp, Tdomain%nz_interp
else 
   read (11,*)
   read (11,*)
endif

read (11,*) Tdomain%n_sls
if (Tdomain%n_sls>0) then
   read (11,*) Tdomain%T1_att, Tdomain%T2_att
else
   read (11,*)
endif

read (11,*) Tdomain%adjoint
Tdomain%nb_simu = 1
if (Tdomain%adjoint)   Tdomain%nb_simu = 2

read (11,*) Tdomain%t_reversal_mirror

read(11,*) Tdomain%n_init

read (11,*) Tdomain%n_source
if (Tdomain%n_source>1 .and. Tdomain%save_trace) then
   if (Tdomain%comp_rot) then
      print *, "NO COMPONENT ROTATION WHEN USING SEVERAL SOURCES"
      Tdomain%comp_rot = .false.
   endif
endif
allocate (Tdomain%Ssource(0:Tdomain%n_source-1))
do i = 0, Tdomain%n_source - 1
   read (11,*) Tdomain%Ssource(i)%realcolat, Tdomain%Ssource(i)%reallong, Tdomain%Ssource(i)%depth
   read (11,*) Tdomain%Ssource(i)%i_type_source
   if (Tdomain%Ssource(i)%i_type_source == 1) then
      read (11,*) Tdomain%Ssource(i)%i_dir(0:2)
      read (11,*)
      read (11,*)
      read (11,*) Tdomain%Ssource(i)%i_time_function
      read (11,*) Tdomain%Ssource(i)%tau_b
      if (Tdomain%Ssource(i)%i_time_function == 3) then
         read (11,*) Tdomain%Ssource(i)%fh(0:3)
      else
         read (11,*)
      endif
   elseif (Tdomain%Ssource(i)%i_type_source == 2) then
      read (11,*) 
      read (11,*) Tdomain%Ssource(i)%Moment(0,0), Tdomain%Ssource(i)%Moment(1,1), Tdomain%Ssource(i)%Moment(2,2)
      read (11,*) Tdomain%Ssource(i)%Moment(0,1), Tdomain%Ssource(i)%Moment(0,2), Tdomain%Ssource(i)%Moment(1,2)
      Tdomain%Ssource(i)%Moment(1,0) = Tdomain%Ssource(i)%Moment(0,1)
      Tdomain%Ssource(i)%Moment(2,0) = Tdomain%Ssource(i)%Moment(0,2)
      Tdomain%Ssource(i)%Moment(2,1) = Tdomain%Ssource(i)%Moment(1,2)
      read (11,*) Tdomain%Ssource(i)%i_time_function
      read (11,*) Tdomain%Ssource(i)%tau_b
      if (Tdomain%Ssource(i)%i_time_function == 3) then
         read (11,*) Tdomain%Ssource(i)%fh(0:3)
      else
         read (11,*)
      endif
   elseif (Tdomain%Ssource(i)%i_type_source == 3) then
      do j = 1,6
         read (11,*)
      enddo
   elseif (Tdomain%Ssource(i)%i_type_source == 4) then
      read (11,*)
      read (11,*)
      read (11,*)
      read (11,*) Tdomain%Ssource(i)%i_time_function
      read (11,*) Tdomain%Ssource(i)%tau_b
      if (Tdomain%Ssource(i)%i_time_function == 3) then
         read (11,*) Tdomain%Ssource(i)%fh(0:3)
      else
         stop 'RANDOM TRACTIONS WANT SRC-TIME-FUNCTION = 3 !!!'
      endif
   else
      stop "WRONG SOURCE TYPE INPUT !!!"
   endif
!   read (11,*) Tdomain%MPML, Tdomain%MPML_coeff
   Tdomain%MPML = .false.
enddo
close (11)


! Read mesh properties

open (12, file=Tdomain%mesh_file, iostat=ok, status="old", form="formatted")
if (ok/=0) then
   write (*,*) "PROCESS ",rg, " CANNOT OPEN ITS MESH_FILE."
   stop
endif
read (12,*) Tdomain%n_dim
read (12,*) Tdomain%n_elem_x, Tdomain%n_elem_y, Tdomain%n_elem_z
read (12,*) Tdomain%i1_mirror, Tdomain%i2_mirror,Tdomain%j1_mirror, Tdomain%j2_mirror,Tdomain%k1_mirror, Tdomain%k2_mirror
do i = 1,6
   read (12,*) (Tdomain%cover(i,j),j=0,1)
enddo
read (12,*) Tdomain%curve, Tdomain%mesh, Tdomain%ellipticity
if (Tdomain%curve) then
   do i = 0,2
      read (12,*) (Tdomain%rot(i,j),j=0,2)
   enddo
endif
read (12,*) Tdomain%n_glob_nodes
allocate (Tdomain%Coord_nodes(0:Tdomain%n_dim-1,0:Tdomain%n_glob_nodes-1))
if (Tdomain%ellipticity)   allocate (Tdomain%Rsph_nodes(0:Tdomain%n_glob_nodes-1))
do i = 0,Tdomain%n_glob_nodes-1
   read (12,*) (Tdomain%Coord_nodes(j,i), j=0,Tdomain%n_dim-1)
   if (Tdomain%ellipticity)   read (12,*) Tdomain%Rsph_nodes(i)
enddo
read (12,*) Tdomain%n_elem
allocate (Tdomain%specel(0:Tdomain%n_elem-1))
do i = 0, Tdomain%n_elem - 1
   read(12,*) Tdomain%specel(i)%mat_index, Tdomain%specel(i)%cov_index, &
              Tdomain%specel(i)%moho_position, Tdomain%specel(i)%topo_position, &
              Tdomain%specel(i)%indx, Tdomain%specel(i)%indy, Tdomain%specel(i)%indz
enddo
if (Tdomain%n_dim == 3) then
   read (12,*) Tdomain%n_nodes
   do i = 0, Tdomain%n_elem - 1
      allocate (Tdomain%specel(i)%Control_Nodes(0:Tdomain%n_nodes-1))
      read(12,*) Tdomain%specel(i)%Control_Nodes(0:Tdomain%n_nodes-1)
   enddo
   read (12,*) Tdomain%n_face
   allocate (Tdomain%sFace(0:Tdomain%n_face-1))
   do i = 0, Tdomain%n_elem - 1
      read(12,*) Tdomain%specel(i)%Near_Faces(0:5)
   enddo
   read (12,*) Tdomain%n_edge
   allocate (Tdomain%sEdge(0:Tdomain%n_edge-1))
   do i = 0, Tdomain%n_elem - 1
      read(12,*) Tdomain%specel(i)%Near_Edges(0:11)
   enddo
   read (12,*) Tdomain%n_vertex
   allocate (Tdomain%sVertex(0:Tdomain%n_vertex-1))
   do i = 0, Tdomain%n_elem - 1
      read(12,*) Tdomain%specel(i)%Near_Vertices(0:7)
   enddo
   read (12,*) Tdomain%n_proc
   allocate (Tdomain%sComm(0:Tdomain%n_proc-1))
   do i = 0,Tdomain%n_proc-1
      read(12,*) Tdomain%sComm(i)%nb_faces, Tdomain%sComm(i)%nb_edges, Tdomain%sComm(i)%nb_vertices
      if (Tdomain%sComm(i)%nb_faces>0) then
          allocate (Tdomain%sComm(i)%faces(0:Tdomain%sComm(i)%nb_faces-1))
          do j = 0,Tdomain%sComm(i)%nb_faces-1
              read(12,*) Tdomain%sComm(i)%faces(j)
          enddo
      endif
      if (Tdomain%sComm(i)%nb_edges>0) then
          allocate (Tdomain%sComm(i)%edges(0:Tdomain%sComm(i)%nb_edges-1))
          do j = 0,Tdomain%sComm(i)%nb_edges-1
              read(12,*) Tdomain%sComm(i)%edges(j)
          enddo
      endif
      if (Tdomain%sComm(i)%nb_vertices>0) then
          allocate (Tdomain%sComm(i)%vertices(0:Tdomain%sComm(i)%nb_vertices-1))
          do j = 0,Tdomain%sComm(i)%nb_vertices-1
              read(12,*) Tdomain%sComm(i)%vertices(j)
          enddo
      endif
   enddo
else
   write (*,*) "A dimension different from 3 is not yet taken into account."
   stop
endif
close (12)


! Read material properties

open (13, file=Tdomain%material_file, status="old", form="formatted")
read (13,*) Tdomain%n_mat, Tdomain%special_z_crust
allocate(Tdomain%sSubdomain(0:Tdomain%n_mat-1))
do i = 0,Tdomain%n_mat-1
   read (13,*) Tdomain%sSubDomain(i)%material_type, Tdomain%sSubDomain(i)%NGLLx, &
               Tdomain%sSubDomain(i)%NGLLy, Tdomain%sSubDomain(i)%NGLLz
enddo
do i = 0,Tdomain%n_mat-1
   if (Tdomain%sSubdomain(i)%material_type == "P") then
      read (13,*) Tdomain%sSubdomain(i)%npow, Tdomain%sSubdomain(i)%Apow, &
                  Tdomain%sSubdomain(i)%Px, Tdomain%sSubdomain(i)%Left, &
                  Tdomain%sSubdomain(i)%Py, Tdomain%sSubdomain(i)%Forward, &
                  Tdomain%sSubdomain(i)%Pz, Tdomain%sSubdomain(i)%Down
   endif
enddo
close (13)
if (Tdomain%special_z_crust/=0) then
    Tdomain%sSubDomain(6)%material_type = Tdomain%sSubDomain(0)%material_type
    Tdomain%sSubDomain(6)%NGLLx = Tdomain%sSubDomain(0)%NGLLx
    Tdomain%sSubDomain(6)%NGLLy = Tdomain%sSubDomain(0)%NGLLy
    Tdomain%sSubDomain(6)%NGLLz = Tdomain%special_z_crust
    do i = 15, 18
        Tdomain%sSubDomain(i)%material_type = Tdomain%sSubDomain(i-13)%material_type
        Tdomain%sSubDomain(i)%NGLLx = Tdomain%sSubDomain(i-13)%NGLLx
        Tdomain%sSubDomain(i)%NGLLy = Tdomain%sSubDomain(i-13)%NGLLy
        Tdomain%sSubDomain(i)%NGLLz = Tdomain%special_z_crust
        Tdomain%sSubDomain(i)%npow = Tdomain%sSubDomain(i-13)%npow
        Tdomain%sSubDomain(i)%Apow = Tdomain%sSubDomain(i-13)%Apow
        Tdomain%sSubDomain(i)%Px = Tdomain%sSubDomain(i-13)%Px
        Tdomain%sSubDomain(i)%Py = Tdomain%sSubDomain(i-13)%Py
        Tdomain%sSubDomain(i)%Pz = Tdomain%sSubDomain(i-13)%Pz
        Tdomain%sSubdomain(i)%Left = Tdomain%sSubdomain(i-13)%Left
        Tdomain%sSubdomain(i)%Forward = Tdomain%sSubdomain(i-13)%Forward
        Tdomain%sSubdomain(i)%Down = Tdomain%sSubdomain(i-13)%Down
    enddo
    do i = 23, 26
        Tdomain%sSubDomain(i)%material_type = Tdomain%sSubDomain(i-12)%material_type
        Tdomain%sSubDomain(i)%NGLLx = Tdomain%sSubDomain(i-12)%NGLLx
        Tdomain%sSubDomain(i)%NGLLy = Tdomain%sSubDomain(i-12)%NGLLy
        Tdomain%sSubDomain(i)%NGLLz = Tdomain%special_z_crust
        Tdomain%sSubDomain(i)%npow = Tdomain%sSubDomain(i-12)%npow
        Tdomain%sSubDomain(i)%Apow = Tdomain%sSubDomain(i-12)%Apow
        Tdomain%sSubDomain(i)%Px = Tdomain%sSubDomain(i-12)%Px
        Tdomain%sSubDomain(i)%Py = Tdomain%sSubDomain(i-12)%Py
        Tdomain%sSubDomain(i)%Pz = Tdomain%sSubDomain(i-12)%Pz
        Tdomain%sSubdomain(i)%Left = Tdomain%sSubdomain(i-12)%Left
        Tdomain%sSubdomain(i)%Forward = Tdomain%sSubdomain(i-12)%Forward
        Tdomain%sSubdomain(i)%Down = Tdomain%sSubdomain(i-12)%Down
    enddo
    do i = 0, Tdomain%n_elem-1
        mat = Tdomain%specel(i)%mat_index
        if (mat==6 .or. &
            mat==15 .or. mat==16 .or. mat==17 .or. mat==18 .or. &
            mat==23 .or. mat==24 .or. mat==25 .or. mat==26) then
            print *,"YOU CANNOT USE A SPECIAL NGLLZ WHEN YOU HAVE PML AT THE SURFACE"
            stop
        endif
        if (Tdomain%specel(i)%moho_position==1) then
            select case (mat)
            case(0)
                Tdomain%specel(i)%mat_index = 6
            case(2)
                Tdomain%specel(i)%mat_index = 15
            case(3)
                Tdomain%specel(i)%mat_index = 16
            case(4)
                Tdomain%specel(i)%mat_index = 17
            case(5)
                Tdomain%specel(i)%mat_index = 18
            case(11)
                Tdomain%specel(i)%mat_index = 23
            case(12)
                Tdomain%specel(i)%mat_index = 24
            case(13)
                Tdomain%specel(i)%mat_index = 25
            case(14)
                Tdomain%specel(i)%mat_index = 26
            end select
        endif
    enddo
endif

allocate (L_Face(0:Tdomain%n_face-1))
L_Face = .true.
allocate (L_Edge(0:Tdomain%n_edge-1))
L_Edge = .true.
do i = 0,Tdomain%n_elem-1
    mat = Tdomain%specel(i)%mat_index
    Tdomain%specel(i)%ngllx = Tdomain%sSubDomain(mat)%NGLLx
    Tdomain%specel(i)%nglly = Tdomain%sSubDomain(mat)%NGLLy
    Tdomain%specel(i)%ngllz = Tdomain%sSubDomain(mat)%NGLLz
    do j = 0,5
        nf = Tdomain%specel(i)%Near_Faces(j)
        if (L_Face(nf)) then
            L_Face(nf) = .false.
            if (j==0 .or. j==5) then
                Tdomain%sFace(nf)%ngll1 = Tdomain%specel(i)%ngllx
                Tdomain%sFace(nf)%ngll2 = Tdomain%specel(i)%nglly
            else if (j==1 .or. j==3) then
                Tdomain%sFace(nf)%ngll1 = Tdomain%specel(i)%ngllx
                Tdomain%sFace(nf)%ngll2 = Tdomain%specel(i)%ngllz
            else
                Tdomain%sFace(nf)%ngll1 = Tdomain%specel(i)%nglly
                Tdomain%sFace(nf)%ngll2 = Tdomain%specel(i)%ngllz
            endif
        endif
    enddo
    do j = 0,11
        ne = Tdomain%specel(i)%Near_Edges(j)
        if (L_Edge(ne)) then
            L_Edge(ne) = .false.
            if (j==0 .or. j==2 .or. j==5 .or. j==9) then
                Tdomain%sEdge(ne)%ngll = Tdomain%specel(i)%ngllx
            else if (j==1 .or. j==3 .or. j==8 .or. j==11) then
                Tdomain%sEdge(ne)%ngll = Tdomain%specel(i)%nglly
            else
                Tdomain%sEdge(ne)%ngll = Tdomain%specel(i)%ngllz
            endif
        endif
    enddo
enddo
deallocate (L_Face,L_Edge)


! Read receivers features

if (Tdomain%save_trace) then
    open (14, file=Tdomain%station_file, status="old")
    read (14,*) Tdomain%n_receivers
    allocate (Tdomain%sReceiver(0:Tdomain%n_receivers-1))
    do i = 0, Tdomain%n_receivers-1
        read(14,*) Tdomain%sReceiver(i)%sta_name, & 
                   Tdomain%sReceiver(i)%realcolat, Tdomain%sReceiver(i)%reallong, Tdomain%sReceiver(i)%elevation
    enddo
    close (14)
endif


return
end subroutine read_Input

