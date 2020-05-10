module scomms


type :: comm

integer :: nb_faces, nb_edges, nb_vertices, ngll, ngllPML, ngllocean, ngllrandt
integer, dimension(:), pointer :: faces, edges, vertices
doubleprecision, dimension(:), pointer :: Give, Take
doubleprecision, dimension(:), pointer :: Giveocean, Takeocean
doubleprecision, dimension(:,:), pointer :: GivePML, TakePML
doubleprecision, dimension(:,:), pointer :: GiveForces, TakeForces
doubleprecision, dimension(:,:,:), pointer :: GiveForcesPML, TakeForcesPML
doubleprecision, dimension(:,:,:), pointer :: RandT

end type

end module scomms
