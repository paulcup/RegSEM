subroutine get_Displ_Face2Elem (Tdomain, n, i_simu)


use sdomains

implicit none

type (Domain), intent (INOUT) :: Tdomain
integer, intent (IN) :: n, i_simu

integer :: nf,i,ngllx,nglly,ngllz
 

ngllx = Tdomain%specel(n)%ngllx; nglly = Tdomain%specel(n)%nglly;  ngllz = Tdomain%specel(n)%ngllz

nf = Tdomain%specel(n)%near_faces(0)
Tdomain%specel(n)%sSimu(i_simu)%Forces(1:ngllx-2,1:nglly-2,0,0:2) = Tdomain%sFace(nf)%sSimu(i_simu)%Forces(1:ngllx-2,1:nglly-2,0:2) 

nf = Tdomain%specel(n)%near_faces(1)
Tdomain%specel(n)%sSimu(i_simu)%Forces(1:ngllx-2,0,1:ngllz-2,0:2) = Tdomain%sFace(nf)%sSimu(i_simu)%Forces(1:ngllx-2,1:ngllz-2,0:2) 

nf = Tdomain%specel(n)%near_faces(2)
Tdomain%specel(n)%sSimu(i_simu)%Forces(ngllx-1,1:nglly-2,1:ngllz-2,0:2) = Tdomain%sFace(nf)%sSimu(i_simu)%Forces(1:nglly-2,1:ngllz-2,0:2) 

nf = Tdomain%specel(n)%near_faces(3)
Tdomain%specel(n)%sSimu(i_simu)%Forces(1:ngllx-2,nglly-1,1:ngllz-2,0:2) = Tdomain%sFace(nf)%sSimu(i_simu)%Forces(1:ngllx-2,1:ngllz-2,0:2) 

nf = Tdomain%specel(n)%near_faces(4)
Tdomain%specel(n)%sSimu(i_simu)%Forces(0,1:nglly-2,1:ngllz-2,0:2) = Tdomain%sFace(nf)%sSimu(i_simu)%Forces(1:nglly-2,1:ngllz-2,0:2) 

nf = Tdomain%specel(n)%near_faces(5)
Tdomain%specel(n)%sSimu(i_simu)%Forces(1:ngllx-2,1:nglly-2,ngllz-1,0:2) = Tdomain%sFace(nf)%sSimu(i_simu)%Forces(1:ngllx-2,1:nglly-2,0:2) 

return
end subroutine get_Displ_Face2Elem
