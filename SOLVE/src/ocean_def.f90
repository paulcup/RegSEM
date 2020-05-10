subroutine ocean_definition (Tdomain)


use sdomains

implicit none

type (Domain), intent (INOUT) :: Tdomain

integer :: n, nf, ne, nv


do n = 0,Tdomain%n_face-1
   Tdomain%sFace(n)%ocean = .false.
enddo
do n = 0,Tdomain%n_edge-1
   Tdomain%sEdge(n)%ocean = .false.
enddo
do n = 0,Tdomain%n_vertex-1
   Tdomain%sVertex(n)%ocean = .false.
enddo

do n = 0,Tdomain%n_elem-1
    Tdomain%specel(n)%ocean = .false.
    if ((Tdomain%ocean .eqv. .true.) .and. (Tdomain%specel(n)%topo_position==1)) then
        Tdomain%specel(n)%ocean = .true.
        nf = Tdomain%specel(n)%Near_Faces(5);   Tdomain%sFace(nf)%ocean = .true.
        ne = Tdomain%specel(n)%Near_Edges(5);   Tdomain%sEdge(ne)%ocean = .true.
        ne = Tdomain%specel(n)%Near_Edges(8);   Tdomain%sEdge(ne)%ocean = .true.
        ne = Tdomain%specel(n)%Near_Edges(9);   Tdomain%sEdge(ne)%ocean = .true.
        ne = Tdomain%specel(n)%Near_Edges(11);   Tdomain%sEdge(ne)%ocean = .true.
        nv = Tdomain%specel(n)%Near_Vertices(4);   Tdomain%sVertex(nv)%ocean = .true.
        nv = Tdomain%specel(n)%Near_Vertices(5);   Tdomain%sVertex(nv)%ocean = .true.
        nv = Tdomain%specel(n)%Near_Vertices(6);   Tdomain%sVertex(nv)%ocean = .true.
        nv = Tdomain%specel(n)%Near_Vertices(7);   Tdomain%sVertex(nv)%ocean = .true.
    endif
enddo

return
end subroutine ocean_definition

