subroutine global_numbering (Tdomain)


use sdomains

implicit none

type(domain), target, intent (INOUT) :: Tdomain

integer :: n, icount, i,j,k, ngllx,nglly,ngllz, nf,ne,nv


icount = 0

do n = 0,Tdomain%n_elem-1
   ngllx = Tdomain%specel(n)%ngllx
   nglly = Tdomain%specel(n)%nglly
   ngllz = Tdomain%specel(n)%ngllz
   allocate (Tdomain%specel(n)%Iglobnum(0:ngllx-1, 0:nglly-1, 0:ngllz-1))
   Tdomain%specel(n)%Iglobnum = -1
   do k = 1,ngllz-2
       do j = 1,nglly-2
           do i = 1,ngllx-2
               Tdomain%specel(n)%Iglobnum(i,j,k) = icount
               icount = icount + 1
           enddo
       enddo
   enddo
enddo

do n = 0,Tdomain%n_face-1
   ngllx = Tdomain%sFace(n)%ngll1
   nglly = Tdomain%sFace(n)%ngll2
   allocate (Tdomain%sFace(n)%Iglobnum_Face(1:ngllx-2,1:nglly-2))
   Tdomain%sFace(n)%Iglobnum_Face = -1
   do j = 1,nglly-2
       do i = 1,ngllx-2
           Tdomain%sFace(n)%Iglobnum_Face(i,j) = icount
           icount = icount + 1
       enddo
   enddo
enddo

do n = 0,Tdomain%n_edge-1
    ngllx = Tdomain%sEdge(n)%ngll
    allocate (Tdomain%sEdge(n)%Iglobnum_Edge(1:ngllx-2))
    Tdomain%sEdge(n)%Iglobnum_Edge = -1
    do i = 1,ngllx-2
        Tdomain%sEdge(n)%Iglobnum_Edge(i) = icount
        icount = icount + 1
    enddo
enddo

do n = 0,Tdomain%n_vertex-1
    Tdomain%sVertex(n)%Iglobnum_Vertex = icount
    icount = icount + 1
enddo

Tdomain%n_glob_points = icount 

do n = 0,Tdomain%n_elem - 1
    ngllx = Tdomain%specel(n)%ngllx
    nglly = Tdomain%specel(n)%nglly
    ngllz = Tdomain%specel(n)%ngllz
    ! Taking information from faces
    nf = Tdomain%specel(n)%Near_Faces(0)
    Tdomain%specel(n)%Iglobnum (1:ngllx-2, 1:nglly-2, 0) = Tdomain%sFace(nf)%Iglobnum_Face(1:ngllx-2, 1:nglly-2)
    nf = Tdomain%specel(n)%Near_Faces(5)
    Tdomain%specel(n)%Iglobnum (1:ngllx-2, 1:nglly-2, ngllz-1) = Tdomain%sFace(nf)%Iglobnum_Face(1:ngllx-2, 1:nglly-2)
    nf = Tdomain%specel(n)%Near_Faces(1)
    Tdomain%specel(n)%Iglobnum (1:ngllx-2, 0, 1:ngllz-2) = Tdomain%sFace(nf)%Iglobnum_Face(1:ngllx-2, 1:ngllz-2)
    nf = Tdomain%specel(n)%Near_Faces(3)
    Tdomain%specel(n)%Iglobnum (1:ngllx-2, nglly-1, 1:ngllz-2) = Tdomain%sFace(nf)%Iglobnum_Face(1:ngllx-2, 1:ngllz-2)
    nf = Tdomain%specel(n)%Near_Faces(2)
    Tdomain%specel(n)%Iglobnum (ngllx-1, 1:nglly-2, 1:ngllz-2) = Tdomain%sFace(nf)%Iglobnum_Face(1:nglly-2, 1:ngllz-2)
    nf = Tdomain%specel(n)%Near_Faces(4)
    Tdomain%specel(n)%Iglobnum (0, 1:nglly-2, 1:ngllz-2) = Tdomain%sFace(nf)%Iglobnum_Face(1:nglly-2, 1:ngllz-2)
    ! Taking information from edges
    ne = Tdomain%specel(n)%Near_Edges(0)
    Tdomain%specel(n)%Iglobnum(1:ngllx-2,0,0) = Tdomain%sEdge(ne)%Iglobnum_Edge(1:ngllx-2)
    ne = Tdomain%specel(n)%Near_Edges(1)
    Tdomain%specel(n)%Iglobnum(ngllx-1,1:nglly-2,0) = Tdomain%sEdge(ne)%Iglobnum_Edge(1:nglly-2)
    ne = Tdomain%specel(n)%Near_Edges(2)
    Tdomain%specel(n)%Iglobnum(1:ngllx-2,nglly-1,0) = Tdomain%sEdge(ne)%Iglobnum_Edge(1:ngllx-2)
    ne = Tdomain%specel(n)%Near_Edges(3)
    Tdomain%specel(n)%Iglobnum(0,1:nglly-2,0) = Tdomain%sEdge(ne)%Iglobnum_Edge(1:nglly-2)
    ne = Tdomain%specel(n)%Near_Edges(4)
    Tdomain%specel(n)%Iglobnum(ngllx-1,0,1:ngllz-2) = Tdomain%sEdge(ne)%Iglobnum_Edge(1:ngllz-2)
    ne = Tdomain%specel(n)%Near_Edges(5)
    Tdomain%specel(n)%Iglobnum(1:ngllx-2,0,ngllz-1) = Tdomain%sEdge(ne)%Iglobnum_Edge(1:ngllx-2)
    ne = Tdomain%specel(n)%Near_Edges(6)
    Tdomain%specel(n)%Iglobnum(0,0,1:ngllz-2) = Tdomain%sEdge(ne)%Iglobnum_Edge(1:ngllz-2)
    ne = Tdomain%specel(n)%Near_Edges(7)
    Tdomain%specel(n)%Iglobnum(ngllx-1,nglly-1,1:ngllz-2) = Tdomain%sEdge(ne)%Iglobnum_Edge(1:ngllz-2)
    ne = Tdomain%specel(n)%Near_Edges(8)
    Tdomain%specel(n)%Iglobnum(ngllx-1,1:nglly-2,ngllz-1) = Tdomain%sEdge(ne)%Iglobnum_Edge(1:nglly-2)
    ne = Tdomain%specel(n)%Near_Edges(9)
    Tdomain%specel(n)%Iglobnum(1:ngllx-2,nglly-1,ngllz-1) = Tdomain%sEdge(ne)%Iglobnum_Edge(1:ngllx-2)
    ne = Tdomain%specel(n)%Near_Edges(10)
    Tdomain%specel(n)%Iglobnum(0,nglly-1,1:ngllz-2) = Tdomain%sEdge(ne)%Iglobnum_Edge(1:ngllz-2)
    ne = Tdomain%specel(n)%Near_Edges(11)
    Tdomain%specel(n)%Iglobnum(0,1:nglly-2,ngllz-1) = Tdomain%sEdge(ne)%Iglobnum_Edge(1:nglly-2)
    ! Taking information from vertices
    nv = Tdomain%specel(n)%Near_Vertices(0)
    Tdomain%specel(n)%Iglobnum(0,0,0) = Tdomain%sVertex(nv)%Iglobnum_Vertex
    nv = Tdomain%specel(n)%Near_Vertices(1)
    Tdomain%specel(n)%Iglobnum(ngllx-1,0,0) = Tdomain%sVertex(nv)%Iglobnum_Vertex
    nv = Tdomain%specel(n)%Near_Vertices(2)
    Tdomain%specel(n)%Iglobnum(ngllx-1,nglly-1,0) = Tdomain%sVertex(nv)%Iglobnum_Vertex
    nv = Tdomain%specel(n)%Near_Vertices(3)
    Tdomain%specel(n)%Iglobnum(0,nglly-1,0) = Tdomain%sVertex(nv)%Iglobnum_Vertex
    nv = Tdomain%specel(n)%Near_Vertices(4)
    Tdomain%specel(n)%Iglobnum(0,0,ngllz-1) = Tdomain%sVertex(nv)%Iglobnum_Vertex
    nv = Tdomain%specel(n)%Near_Vertices(5)
    Tdomain%specel(n)%Iglobnum(ngllx-1,0,ngllz-1) = Tdomain%sVertex(nv)%Iglobnum_Vertex
    nv = Tdomain%specel(n)%Near_Vertices(6)
    Tdomain%specel(n)%Iglobnum(ngllx-1,nglly-1,ngllz-1) = Tdomain%sVertex(nv)%Iglobnum_Vertex
    nv = Tdomain%specel(n)%Near_Vertices(7)
    Tdomain%specel(n)%Iglobnum(0,nglly-1,ngllz-1) = Tdomain%sVertex(nv)%Iglobnum_Vertex
enddo

do n = 0,Tdomain%n_face-1
    deallocate (Tdomain%sFace(n)%Iglobnum_Face)
enddo
do n = 0,Tdomain%n_edge-1
    deallocate (Tdomain%sEdge(n)%Iglobnum_Edge)
enddo

return
end
