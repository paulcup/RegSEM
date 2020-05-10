module sdomains


use selements
use sfaces
use sedges
use svertices
use scomms
use ssources
use stimeparams
use ssubdomains
use sreceivers
use sadjoints

type :: domain

type(timeparam) :: sTimeParam
type(source), dimension (:), pointer :: sSource
type(element), dimension(:), pointer :: specel
type(face), dimension (:), pointer :: sFace
type(comm), dimension (:), pointer :: sComm
type(edge), dimension (:), pointer :: sEdge
type(vertex), dimension (:), pointer :: sVertex
type(subdomain), dimension (:), pointer :: sSubdomain
type(receiver), dimension (:), pointer :: sReceiver
type(adjoint), dimension (:), pointer :: sAdjoint

logical :: curve, aniso, ocean, ellipticity, comp_rot, adjoint, save_trace, save_snapshot, save_model, MPML

integer :: n_source, n_dim, n_glob_nodes, n_mat, n_nodes, n_receivers, n_proc, &
           n_elem, n_face, n_edge, n_vertex, n_glob_points, n_sls, mesh, special_z_crust, &
           t_reversal_mirror, recl_mirror, nb_simu, n_source_adj, n_elem_x, n_elem_y, n_elem_z, &
           i1_mirror, i2_mirror, j1_mirror, j2_mirror, k1_mirror, k2_mirror, n_init, nx_interp, ny_interp, nz_interp

doubleprecision :: T1_att, T2_att, MPML_coeff
doubleprecision, dimension (0:2,0:2) :: rot
doubleprecision, dimension (1:6,0:1) :: cover
doubleprecision, dimension (:), pointer :: Rsph_nodes
doubleprecision, dimension (:,:), pointer :: Coord_nodes, GlobCoord

character (len=30) :: mesh_file, station_file, material_file

end type


contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
doubleprecision function Comp_shapefunc(which_nod,xi,eta,zeta)

integer :: which_nod
doubleprecision :: xi,eta,zeta

select case (which_nod)
 case (0)
    Comp_shapefunc = 0.125 * xi*(xi-1) * eta*(eta-1) * zeta*(zeta-1)
 case (1)
    Comp_shapefunc = 0.125 * xi*(xi+1) * eta*(eta-1) * zeta*(zeta-1)
 case (2)
    Comp_shapefunc = 0.125 * xi*(xi+1) * eta*(eta+1) * zeta*(zeta-1)
 case (3)
    Comp_shapefunc = 0.125 * xi*(xi-1) * eta*(eta+1) * zeta*(zeta-1)
 case (4)
    Comp_shapefunc = 0.125 * xi*(xi-1) * eta*(eta-1) * zeta*(zeta+1)
 case (5)
    Comp_shapefunc = 0.125 * xi*(xi+1) * eta*(eta-1) * zeta*(zeta+1)
 case (6)
    Comp_shapefunc = 0.125 * xi*(xi+1) * eta*(eta+1) * zeta*(zeta+1)
 case (7)
    Comp_shapefunc = 0.125 * xi*(xi-1) * eta*(eta+1) * zeta*(zeta+1)
 case (8)
    Comp_shapefunc = 0.25 * (1-xi**2) * eta*(eta-1) * zeta*(zeta-1)
 case (9)
    Comp_shapefunc = 0.25 * xi*(xi+1) * (1-eta**2) * zeta*(zeta-1)
 case (10)
    Comp_shapefunc = 0.25 * (1-xi**2) * eta*(eta+1) * zeta*(zeta-1)
 case (11)
    Comp_shapefunc = 0.25 * xi*(xi-1) * (1-eta**2) * zeta*(zeta-1)
 case (12)
    Comp_shapefunc = 0.25 * xi*(xi-1) * eta*(eta-1) * (1-zeta**2)
 case (13)
    Comp_shapefunc = 0.25 * xi*(xi+1) * eta*(eta-1) * (1-zeta**2)
 case (14)
    Comp_shapefunc = 0.25 * xi*(xi+1) * eta*(eta+1) * (1-zeta**2)
 case (15)
    Comp_shapefunc = 0.25 * xi*(xi-1) * eta*(eta+1) * (1-zeta**2)
 case (16)
    Comp_shapefunc = 0.25 * (1-xi**2) * eta*(eta-1) * zeta*(zeta+1)
 case (17)
    Comp_shapefunc = 0.25 * xi*(xi+1) * (1-eta**2) * zeta*(zeta+1)
 case (18)
    Comp_shapefunc = 0.25 * (1-xi**2) * eta*(eta+1) * zeta*(zeta+1)
 case (19)
    Comp_shapefunc = 0.25 * xi*(xi-1) * (1-eta**2) * zeta*(zeta+1)
 case(20)
    Comp_shapefunc = 0.5 * (1-xi**2) * (1-eta**2) * zeta*(zeta-1)
 case(21)
    Comp_shapefunc = 0.5 * (1-xi**2) * eta*(eta-1) * (1-zeta**2)
 case(22)
    Comp_shapefunc = 0.5 * xi*(xi+1) * (1-eta**2) * (1-zeta**2)
 case(23)
    Comp_shapefunc = 0.5 * (1-xi**2) * eta*(eta+1) * (1-zeta**2)
 case(24)
    Comp_shapefunc = 0.5 * xi*(xi-1) * (1-eta**2) * (1-zeta**2)
 case(25)
    Comp_shapefunc = 0.5 * (1-xi**2) * (1-eta**2) * zeta*(zeta+1)
 case(26)
    Comp_shapefunc = (1-xi**2) * (1-eta**2) * (1-zeta**2)
end select

return
end function

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
doubleprecision function Comp_derivshapefunc(which_nod,xi,eta,zeta,compo)

integer :: which_nod, compo
doubleprecision :: xi,eta,zeta

doubleprecision, dimension(0:2) :: df

select case (which_nod)
 case (0)
    df(0) = 0.125 * (2*xi-1) * eta*(eta-1) * zeta*(zeta-1)
    df(1) = 0.125 * xi*(xi-1) * (2*eta-1) * zeta*(zeta-1)
    df(2) = 0.125 * xi*(xi-1) * eta*(eta-1) * (2*zeta-1)
 case (1)
    df(0) = 0.125 * (2*xi+1) * eta*(eta-1) * zeta*(zeta-1)
    df(1) = 0.125 * xi*(xi+1) * (2*eta-1) * zeta*(zeta-1)
    df(2) = 0.125 * xi*(xi+1) * eta*(eta-1) * (2*zeta-1)
 case (2)
    df(0) = 0.125 * (2*xi+1) * eta*(eta+1) * zeta*(zeta-1)
    df(1) = 0.125 * xi*(xi+1) * (2*eta+1) * zeta*(zeta-1)
    df(2) = 0.125 * xi*(xi+1) * eta*(eta+1) * (2*zeta-1)
 case (3)
    df(0) = 0.125 * (2*xi-1) * eta*(eta+1) * zeta*(zeta-1)
    df(1) = 0.125 * xi*(xi-1) * (2*eta+1) * zeta*(zeta-1)
    df(2) = 0.125 * xi*(xi-1) * eta*(eta+1) * (2*zeta-1)
 case (4)
    df(0) = 0.125 * (2*xi-1) * eta*(eta-1) * zeta*(zeta+1)
    df(1) = 0.125 * xi*(xi-1) * (2*eta-1) * zeta*(zeta+1)
    df(2) = 0.125 * xi*(xi-1) * eta*(eta-1) * (2*zeta+1)
 case (5)
    df(0) = 0.125 * (2*xi+1) * eta*(eta-1) * zeta*(zeta+1)
    df(1) = 0.125 * xi*(xi+1) * (2*eta-1) * zeta*(zeta+1)
    df(2) = 0.125 * xi*(xi+1) * eta*(eta-1) * (2*zeta+1)
 case (6)
    df(0) = 0.125 * (2*xi+1) * eta*(eta+1) * zeta*(zeta+1)
    df(1) = 0.125 * xi*(xi+1) * (2*eta+1) * zeta*(zeta+1)
    df(2) = 0.125 * xi*(xi+1) * eta*(eta+1) * (2*zeta+1)
 case (7)
    df(0) = 0.125 * (2*xi-1) * eta*(eta+1) * zeta*(zeta+1)
    df(1) = 0.125 * xi*(xi-1) * (2*eta+1) * zeta*(zeta+1)
    df(2) = 0.125 * xi*(xi-1) * eta*(eta+1) * (2*zeta+1)
 case (8)
    df(0) = 0.25 * (-2*xi) * eta*(eta-1) * zeta*(zeta-1)
    df(1) = 0.25 * (1-xi**2) * (2*eta-1) * zeta*(zeta-1)
    df(2) = 0.25 * (1-xi**2) * eta*(eta-1) * (2*zeta-1)
 case (9)
    df(0) = 0.25 * (2*xi+1) * (1-eta**2) * zeta*(zeta-1)
    df(1) = 0.25 * xi*(xi+1) * (-2*eta) * zeta*(zeta-1)
    df(2) = 0.25 * xi*(xi+1) * (1-eta**2) * (2*zeta-1)
 case (10)
    df(0) = 0.25 * (-2*xi) * eta*(eta+1) * zeta*(zeta-1)
    df(1) = 0.25 * (1-xi**2) * (2*eta+1) * zeta*(zeta-1)
    df(2) = 0.25 * (1-xi**2) * eta*(eta+1) * (2*zeta-1)
 case (11)
    df(0) = 0.25 * (2*xi-1) * (1-eta**2) * zeta*(zeta-1)
    df(1) = 0.25 * xi*(xi-1) * (-2*eta) * zeta*(zeta-1)
    df(2) = 0.25 * xi*(xi-1) * (1-eta**2) * (2*zeta-1)
 case (12)
    df(0) = 0.25 * (2*xi-1) * eta*(eta-1) * (1-zeta**2)
    df(1) = 0.25 * xi*(xi-1) * (2*eta-1) * (1-zeta**2)
    df(2) = 0.25 * xi*(xi-1) * eta*(eta-1) * (-2*zeta)
 case (13)
    df(0) = 0.25 * (2*xi+1) * eta*(eta-1) * (1-zeta**2)
    df(1) = 0.25 * xi*(xi+1) * (2*eta-1) * (1-zeta**2)
    df(2) = 0.25 * xi*(xi+1) * eta*(eta-1) * (-2*zeta)
 case (14)
    df(0) = 0.25 * (2*xi+1) * eta*(eta+1) * (1-zeta**2)
    df(1) = 0.25 * xi*(xi+1) * (2*eta+1) * (1-zeta**2)
    df(2) = 0.25 * xi*(xi+1) * eta*(eta+1) * (-2*zeta)
 case (15)
    df(0) = 0.25 * (2*xi-1) * eta*(eta+1) * (1-zeta**2)
    df(1) = 0.25 * xi*(xi-1) * (2*eta+1) * (1-zeta**2)
    df(2) = 0.25 * xi*(xi-1) * eta*(eta+1) * (-2*zeta)
 case (16)
    df(0) = 0.25 * (-2*xi) * eta*(eta-1) * zeta*(zeta+1)
    df(1) = 0.25 * (1-xi**2) * (2*eta-1) * zeta*(zeta+1)
    df(2) = 0.25 * (1-xi**2) * eta*(eta-1) * (2*zeta+1)
 case (17)
    df(0) = 0.25 * (2*xi+1) * (1-eta**2) * zeta*(zeta+1)
    df(1) = 0.25 * xi*(xi+1) * (-2*eta) * zeta*(zeta+1)
    df(2) = 0.25 * xi*(xi+1) * (1-eta**2) * (2*zeta+1)
 case (18)
    df(0) = 0.25 * (-2*xi) * eta*(eta+1) * zeta*(zeta+1)
    df(1) = 0.25 * (1-xi**2) * (2*eta+1) * zeta*(zeta+1)
    df(2) = 0.25 * (1-xi**2) * eta*(eta+1) * (2*zeta+1)
 case (19)
    df(0) = 0.25 * (2*xi-1) * (1-eta**2) * zeta*(zeta+1)
    df(1) = 0.25 * xi*(xi-1) * (-2*eta) * zeta*(zeta+1)
    df(2) = 0.25 * xi*(xi-1) * (1-eta**2) * (2*zeta+1)
 case(20)
    df(0) = 0.5 * (-2*xi) * (1-eta**2) * zeta*(zeta-1)
    df(1) = 0.5 * (1-xi**2) * (-2*eta) * zeta*(zeta-1)
    df(2) = 0.5 * (1-xi**2) * (1-eta**2) * (2*zeta-1)
 case(21)
    df(0) = 0.5 * (-2*xi) * eta*(eta-1) * (1-zeta**2)
    df(1) = 0.5 * (1-xi**2) * (2*eta-1) * (1-zeta**2)
    df(2) = 0.5 * (1-xi**2) * eta*(eta-1) * (-2*zeta)
 case(22)
    df(0) = 0.5 * (2*xi+1) * (1-eta**2) * (1-zeta**2)
    df(1) = 0.5 * xi*(xi+1) * (-2*eta) * (1-zeta**2)
    df(2) = 0.5 * xi*(xi+1) * (1-eta**2) * (-2*zeta)
 case(23)
    df(0) = 0.5 * (-2*xi) * eta*(eta+1) * (1-zeta**2)
    df(1) = 0.5 * (1-xi**2) * (2*eta+1) * (1-zeta**2)
    df(2) = 0.5 * (1-xi**2) * eta*(eta+1) * (-2*zeta)
 case(24)
    df(0) = 0.5 * (2*xi-1) * (1-eta**2) * (1-zeta**2)
    df(1) = 0.5 * xi*(xi-1) * (-2*eta) * (1-zeta**2)
    df(2) = 0.5 * xi*(xi-1) * (1-eta**2) * (-2*zeta)
 case(25)
    df(0) = 0.5 * (-2*xi) * (1-eta**2) * zeta*(zeta+1)
    df(1) = 0.5 * (1-xi**2) * (-2*eta) * zeta*(zeta+1)
    df(2) = 0.5 * (1-xi**2) * (1-eta**2) * (2*zeta+1)
 case(26)
    df(0) = (-2*xi) * (1-eta**2) * (1-zeta**2)
    df(1) = (1-xi**2) * (-2*eta) * (1-zeta**2)
    df(2) = (1-xi**2) * (1-eta**2) * (-2*zeta)
end select

Comp_derivshapefunc = df(compo)

return
end function


end module sdomains
