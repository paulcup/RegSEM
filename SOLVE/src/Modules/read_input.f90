module read_input


implicit none

private
public :: param, voxet, mesh, homog, &
          load_param, load_voxet, load_mesh, load_homog

type :: param
 integer :: input_type   ! 1: a mesh, 2: a Cstar_rho.res file
end type param

type :: voxet
 integer :: voxettype   ! 1: adapts the upper bounds to delta, 2: adapts delta to the bounds
 integer, dimension(1:3) :: nb_pts
 doubleprecision, dimension(1:3) :: delta
 doubleprecision, dimension(1:3,1:2) :: bounds
end type voxet

type :: mesh
 integer :: nb_node, &   ! number of nodes in the whole domain
            nb_node_attrib, &   ! number of node attributes
            log_node_bmark, &   ! Is there a mark for the nodes (0:no, 1:yes) ?
            nb_ele, &   ! number of elements in the whole domain
            nb_ele_node, &   ! number of nodes per element (4 or 10)
            log_rmark, &   ! Is there a region mark -only for the elements- (0:no, 1:yes) ?
            node_valence_max   ! maximum valence among all the nodes
 integer, dimension(:), allocatable :: rmark, &   ! region mark of each element (default is -1)
                                       node_valence   ! valence of each node
 integer, dimension(:,:), allocatable :: ele_node, &   ! nodes of each element
                                         node_ele   ! elements of eahc node
 doubleprecision, dimension(:,:), allocatable :: node_coord(:,:)   ! coordinates of each node
 character(len=60) :: meshtype, &   ! tetgen
                      meshname   ! name of the .node, .ele and .neigh files
end type mesh

type :: homog
 doubleprecision, dimension(:), allocatable :: xseries,yseries,zseries
 doubleprecision, dimension(1:3,1:2) :: bounds
 doubleprecision, dimension(:,:,:), allocatable :: rho
 doubleprecision, dimension(:,:,:,:), allocatable :: Cstar
end type homog


contains

!------------------------------------------------------------------------------------
subroutine load_param(gparam,gmesh)


type(param), intent(INOUT) :: gparam
type(mesh), intent (INOUT) :: gmesh

integer :: ios, i
character(len=60) :: tmp


open (11,file="param.dat",status="old",iostat=ios)
if (ios>0)   stop 'NO param.dat FILE. I CANNOT WORK !!!'

read (11,*) gparam%input_type   ! 1: a mesh, 2: a Cstar_rho.res file

if (gparam%input_type==1) then
   read (11,'(a)') tmp   ! tetgen
   i = scan(tmp,"#")
   if (i==0)   i = len(tmp)+1
   write (gmesh%meshtype,'(a)') tmp(1:i-1)
   read (11,'(a)') tmp   ! name of the .node, .ele and .neigh files
   i = scan(tmp,"#")
   if (i==0)   i = len(tmp)+1
   write (gmesh%meshname,'(a)') tmp(1:i-1)
else if (gparam%input_type==2) then
   ! Line 2 and 3 of param.dat are useless in this case
else
   stop 'ERROR: input_type must be 1 or 2 !!!'
endif

close (11)


end subroutine load_param

!------------------------------------------------------------------------------------
subroutine load_voxet(gvoxet)


type(voxet), intent(INOUT) :: gvoxet

integer :: ios, i


open (11,file="voxet.dat",status="old",iostat=ios)
if (ios>0)   stop 'NO voxet.dat FILE. I CANNOT WORK !!!'

read (11,*) gvoxet%voxettype   ! 1: adapts the upper bounds to delta, 2: adapts delta to the bounds
if (gvoxet%voxettype==1) then
   do i = 1,3
      read (11,*) gvoxet%bounds(i,1), gvoxet%bounds(i,2), gvoxet%delta(i)
      gvoxet%nb_pts(i) = nint((gvoxet%bounds(i,2)-gvoxet%bounds(i,1))/gvoxet%delta(i)) + 1
      gvoxet%bounds(i,2) = gvoxet%bounds(i,1) + (gvoxet%nb_pts(i)-1)*gvoxet%delta(i)
   enddo
   print *, "Nb pts and new upper bounds :"
   do i = 1,3
      print *, gvoxet%nb_pts(i), gvoxet%bounds(i,2)
   enddo
else if (gvoxet%voxettype==2) then
   do i = 1,3
      read (11,*) gvoxet%bounds(i,1), gvoxet%bounds(i,2), gvoxet%delta(i)
      gvoxet%nb_pts(i) = nint((gvoxet%bounds(i,2)-gvoxet%bounds(i,1))/gvoxet%delta(i)) + 1
      gvoxet%delta(i) = (gvoxet%bounds(i,2)-gvoxet%bounds(i,1)) / (gvoxet%nb_pts(i)-1)
   enddo
   print *, "Nb pts and new delta :"
   do i = 1,3
      print *, gvoxet%nb_pts(i), gvoxet%delta(i)
   enddo
else
   stop 'ERROR: voxettype must be 1 or 2 !!!'
endif

close (11)


end subroutine load_voxet

!------------------------------------------------------------------------------------
subroutine load_mesh(gmesh)


type(mesh), intent (INOUT) :: gmesh

integer :: ios, d, i,j,k
integer, dimension(:,:), allocatable :: node_ele_table
doubleprecision, dimension(:,:), allocatable :: tmp_coord
character(len=60) :: node_file, ele_file, neigh_file
logical :: tetgenfiles_start_at_0 = .false.


if (trim(gmesh%meshtype)=='tetgen') then   ! reads the .node, .ele and .neigh files

   ! Reading the .node file
   ! It contains nb_node, nb_node_attrib, log_node_bmark, node_coord,
   ! and optionally node_attrib and node_bmark
   write (node_file,'(a,a)') trim(gmesh%meshname), ".node"
   open (21,file=trim(node_file),form="formatted",status="old",iostat=ios)
   if (ios>0)   stop 'CANNOT FIND THE .node FILE !!!'
   read (21,*) gmesh%nb_node, d, gmesh%nb_node_attrib, gmesh%log_node_bmark
   if (d/=3)   stop 'THE NODES ARE EXPECTED TO HAVE 3 COORDINATES !!!'
   if (gmesh%nb_node_attrib .neqv. 0)   print *, 'WARNING: CANNOT DEAL WITH NODE ATTRIBUTES'
   if (gmesh%log_node_bmark .eqv. 1)   print *, 'WARNING: CANNOT DEAL WITH NODE MARKS'
   allocate (gmesh%node_coord(1:3,1:gmesh%nb_node))
   do i = 1,gmesh%nb_node
      read (21,*) j, gmesh%node_coord(1,i), gmesh%node_coord(2,i), gmesh%node_coord(3,i)
   enddo
   close (21)

   ! Reading the .ele file
   ! It contains nb_ele, nb_ele_node, log_rmark, ele_node,
   ! and optionally rmark
   write (ele_file,'(a,a)') trim(gmesh%meshname), ".ele"
   open (22,file=trim(ele_file),form="formatted",status="old",iostat=ios)
   if (ios>0)   stop 'CANNOT FIND THE .ele FILE !!!'
   read (22,*) gmesh%nb_ele, gmesh%nb_ele_node, gmesh%log_rmark
   allocate (gmesh%ele_node(1:gmesh%nb_ele,1:gmesh%nb_ele_node))
   allocate (gmesh%rmark(1:gmesh%nb_ele));   gmesh%rmark(:) = -1
   if (gmesh%log_rmark==1) then
      read (22,*) j, gmesh%ele_node(1,:), gmesh%rmark(1)
      if (j==0)   tetgenfiles_start_at_0 = .true.
      do i = 2,gmesh%nb_ele
         read (22,*) j, gmesh%ele_node(i,:), gmesh%rmark(i)
      enddo
   else
      stop 'A REGION MARK FOR EACH ELEMENT IS EXPECTED !!!'
   endif
   close (22)
   if (tetgenfiles_start_at_0) then
      gmesh%ele_node = gmesh%ele_node + 1
   endif

   ! Determining the valence and the elements of each node
   allocate (gmesh%node_valence(1:gmesh%nb_node));   gmesh%node_valence(:) = 0
   allocate (node_ele_table(1:gmesh%nb_node,1:60))
   do i = 1,gmesh%nb_ele
    do j = 1,gmesh%nb_ele_node
       k = gmesh%ele_node(i,j)
       gmesh%node_valence(k) = gmesh%node_valence(k) + 1
       node_ele_table(k,gmesh%node_valence(k)) = i
    enddo
   enddo
   gmesh%node_valence_max = 0
   do i = 1,gmesh%nb_node
      if (gmesh%node_valence(i)>gmesh%node_valence_max) then
         gmesh%node_valence_max = gmesh%node_valence(i)
      endif
   enddo
   allocate (gmesh%node_ele(1:gmesh%nb_node,1:gmesh%node_valence_max))
   gmesh%node_ele(:,:) = -1
   do i = 1,gmesh%nb_node
    do j = 1,gmesh%node_valence(i)
       gmesh%node_ele(i,j) = node_ele_table(i,j)
    enddo
   enddo
   deallocate (node_ele_table)

else

   stop 'BAD MESHTYPE: ONLY TETGEN IS ENABLED'

endif


end subroutine load_mesh

!------------------------------------------------------------------------------------
subroutine load_homog(ghomog)


type(homog), intent (INOUT) :: ghomog

integer :: i,j,k, ios, i_count, length
integer, dimension(1:3) :: nb_pt_dir
doubleprecision, dimension(1:3) :: delta_dir
doubleprecision, dimension(:), allocatable :: tmp

integer, parameter :: nb_pt_val = 22   ! rho + C


allocate (tmp(1:nb_pt_val))
inquire(iolength=length) tmp
open (11,file="Cstar_rho.res",form="unformatted",access="direct",recl=length,status="old",iostat=ios)
if (ios>0)   stop 'BINARY FILE Cstar_rho.res IS NEEDED'
read (11,rec=1) tmp
ghomog%bounds(1,1) = tmp(1);   ghomog%bounds(1,2) = tmp(2);   nb_pt_dir(1) = int(tmp(3))
ghomog%bounds(2,1) = tmp(4);   ghomog%bounds(2,2) = tmp(5);   nb_pt_dir(2) = int(tmp(6))
ghomog%bounds(3,1) = tmp(7);   ghomog%bounds(3,2) = tmp(8);   nb_pt_dir(3) = int(tmp(9))
do i = 1,3
   delta_dir(i) = (ghomog%bounds(i,2)-ghomog%bounds(i,1)) / (nb_pt_dir(i)-1)
enddo
allocate (ghomog%rho(1:nb_pt_dir(1),1:nb_pt_dir(2),1:nb_pt_dir(3)))
allocate (ghomog%Cstar(1:nb_pt_dir(1),1:nb_pt_dir(2),1:nb_pt_dir(3),1:21))
do k = 1,nb_pt_dir(3)
 do j = 1,nb_pt_dir(2)
  do i = 1,nb_pt_dir(1)
     i_count = (k-1)*nb_pt_dir(1)*nb_pt_dir(2) + (j-1)*nb_pt_dir(1) + i + 1   ! +1 because of the header
     read (11,rec=i_count) tmp
     ghomog%rho(i,j,k) = tmp(1)
     ghomog%Cstar(i,j,k,1:21) = tmp(2:22)
  enddo
 enddo
enddo
close (11)

allocate (ghomog%xseries(1:nb_pt_dir(1)))
allocate (ghomog%yseries(1:nb_pt_dir(2)))
allocate (ghomog%zseries(1:nb_pt_dir(3)))
do i = 1,nb_pt_dir(1)
   ghomog%xseries(i) = ghomog%bounds(1,1) + (i-1)*delta_dir(1)
enddo
do j = 1,nb_pt_dir(2)
   ghomog%yseries(j) = ghomog%bounds(2,1) + (j-1)*delta_dir(2)
enddo
do k = 1,nb_pt_dir(3)
   ghomog%zseries(k) = ghomog%bounds(3,1) + (k-1)*delta_dir(3)
enddo


end subroutine load_homog


end module read_input
