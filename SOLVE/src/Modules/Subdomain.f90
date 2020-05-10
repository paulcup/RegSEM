module ssubdomains


type Subdomain

logical :: Px, Py, Pz, Left, Forward, Down

integer :: NGLLx, NGLLy, NGLLz, npow

doubleprecision :: Apow
doubleprecision, dimension (:), pointer :: GLLcx, GLLpolx, GLLwx
doubleprecision, dimension (:,:), pointer :: hprimex, hTprimex
doubleprecision, dimension (:), pointer :: GLLcy, GLLpoly, GLLwy
doubleprecision, dimension (:,:), pointer :: hprimey, hTprimey
doubleprecision, dimension (:), pointer :: GLLcz, GLLpolz, GLLwz
doubleprecision, dimension (:,:), pointer :: hprimez, hTprimez

character(len=1) :: material_type

end type

end module ssubdomains
