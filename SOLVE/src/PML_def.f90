subroutine PML_definition (Tdomain)


use sdomains

implicit none

type (Domain), intent (INOUT) :: Tdomain

integer :: n, i, nf, ne, nv, mat


do n = 0,Tdomain%n_elem-1 
   mat = Tdomain%specel(n)%mat_index
   Tdomain%specel(n)%PML = .false.
   if (Tdomain%sSubDomain(mat)%material_type == "P") Tdomain%specel(n)%PML = .true.
enddo

do n = 0,Tdomain%n_face-1
   Tdomain%sFace(n)%PML = .true.
   Tdomain%sFace(n)%Abs = .false.
enddo
do n = 0,Tdomain%n_edge-1
   Tdomain%sEdge(n)%PML = .true.
   Tdomain%sEdge(n)%Abs = .false.
enddo
do n = 0,Tdomain%n_vertex-1
   Tdomain%sVertex(n)%PML = .true.
   Tdomain%sVertex(n)%Abs = .false.
enddo

do n = 0,Tdomain%n_elem-1
 if (Tdomain%specel(n)%PML) then
    mat = Tdomain%specel(n)%mat_index
    if (Tdomain%sSubdomain(mat)%Px) then
       if (Tdomain%sSubdomain(mat)%Py) then
          if (Tdomain%sSubdomain(mat)%Pz) then
             if (Tdomain%sSubdomain(mat)%Left) then
                if (Tdomain%sSubdomain(mat)%Forward) then
                   if (Tdomain%sSubdomain(mat)%Down) then
                      nv = Tdomain%specel(n)%Near_Vertices(6);   Tdomain%sVertex(nv)%PML = .false.
                   else
                      nv = Tdomain%specel(n)%Near_Vertices(2);   Tdomain%sVertex(nv)%PML = .false.
                   endif
                else
                   if (Tdomain%sSubdomain(mat)%Down) then
                      nv = Tdomain%specel(n)%Near_Vertices(5);   Tdomain%sVertex(nv)%PML = .false.
                   else
                      nv = Tdomain%specel(n)%Near_Vertices(1);   Tdomain%sVertex(nv)%PML = .false.
                   endif
                endif
             else
                if (Tdomain%sSubdomain(mat)%Forward) then
                   if (Tdomain%sSubdomain(mat)%Down) then
                      nv = Tdomain%specel(n)%Near_Vertices(7);   Tdomain%sVertex(nv)%PML = .false.
                   else
                      nv = Tdomain%specel(n)%Near_Vertices(3);   Tdomain%sVertex(nv)%PML = .false.
                   endif
                else
                   if (Tdomain%sSubdomain(mat)%Down) then
                      nv = Tdomain%specel(n)%Near_Vertices(4);   Tdomain%sVertex(nv)%PML = .false.
                   else
                      nv = Tdomain%specel(n)%Near_Vertices(0);   Tdomain%sVertex(nv)%PML = .false.
                   endif
                endif
             endif
          else
             if (Tdomain%sSubdomain(mat)%Left) then
                if (Tdomain%sSubdomain(mat)%Forward) then
                   ne = Tdomain%specel(n)%Near_Edges(7);   Tdomain%sEdge(ne)%PML = .false.
                   nv = Tdomain%specel(n)%Near_Vertices(2);   Tdomain%sVertex(nv)%PML = .false.
                   nv = Tdomain%specel(n)%Near_Vertices(6);   Tdomain%sVertex(nv)%PML = .false.
                else
                   ne = Tdomain%specel(n)%Near_Edges(4);   Tdomain%sEdge(ne)%PML = .false.
                   nv = Tdomain%specel(n)%Near_Vertices(1);   Tdomain%sVertex(nv)%PML = .false.
                   nv = Tdomain%specel(n)%Near_Vertices(5);   Tdomain%sVertex(nv)%PML = .false.
                endif
             else
                if (Tdomain%sSubdomain(mat)%Forward) then
                   ne = Tdomain%specel(n)%Near_Edges(10);   Tdomain%sEdge(ne)%PML = .false.
                   nv = Tdomain%specel(n)%Near_Vertices(3);   Tdomain%sVertex(nv)%PML = .false.
                   nv = Tdomain%specel(n)%Near_Vertices(7);   Tdomain%sVertex(nv)%PML = .false.
                else
                   ne = Tdomain%specel(n)%Near_Edges(6);   Tdomain%sEdge(ne)%PML = .false.
                   nv = Tdomain%specel(n)%Near_Vertices(0);   Tdomain%sVertex(nv)%PML = .false.
                   nv = Tdomain%specel(n)%Near_Vertices(4);   Tdomain%sVertex(nv)%PML = .false.
                endif
             endif
          endif
       else
          if (Tdomain%sSubdomain(mat)%Pz) then
             if (Tdomain%sSubdomain(mat)%Left) then
                if (Tdomain%sSubdomain(mat)%Down) then
                   ne = Tdomain%specel(n)%Near_Edges(8);   Tdomain%sEdge(ne)%PML = .false.
                   nv = Tdomain%specel(n)%Near_Vertices(5);   Tdomain%sVertex(nv)%PML = .false.
                   nv = Tdomain%specel(n)%Near_Vertices(6);   Tdomain%sVertex(nv)%PML = .false.
                else
                   ne = Tdomain%specel(n)%Near_Edges(1);   Tdomain%sEdge(ne)%PML = .false.
                   nv = Tdomain%specel(n)%Near_Vertices(1);   Tdomain%sVertex(nv)%PML = .false.
                   nv = Tdomain%specel(n)%Near_Vertices(2);   Tdomain%sVertex(nv)%PML = .false.
                endif
             else
                if (Tdomain%sSubdomain(mat)%Down) then
                   ne = Tdomain%specel(n)%Near_Edges(11);   Tdomain%sEdge(ne)%PML = .false.
                   nv = Tdomain%specel(n)%Near_Vertices(4);   Tdomain%sVertex(nv)%PML = .false.
                   nv = Tdomain%specel(n)%Near_Vertices(7);   Tdomain%sVertex(nv)%PML = .false.
                else
                   ne = Tdomain%specel(n)%Near_Edges(3);   Tdomain%sEdge(ne)%PML = .false.
                   nv = Tdomain%specel(n)%Near_Vertices(0);   Tdomain%sVertex(nv)%PML = .false.
                   nv = Tdomain%specel(n)%Near_Vertices(3);   Tdomain%sVertex(nv)%PML = .false.
                endif
             endif
          else
            if (Tdomain%sSubdomain(mat)%Left) then
               nf = Tdomain%specel(n)%Near_Faces(2);   Tdomain%sFace(nf)%PML = .false.
               ne = Tdomain%specel(n)%Near_Edges(1);   Tdomain%sEdge(ne)%PML = .false.
               ne = Tdomain%specel(n)%Near_Edges(4);   Tdomain%sEdge(ne)%PML = .false.
               ne = Tdomain%specel(n)%Near_Edges(7);   Tdomain%sEdge(ne)%PML = .false.
               ne = Tdomain%specel(n)%Near_Edges(8);   Tdomain%sEdge(ne)%PML = .false.
               nv = Tdomain%specel(n)%Near_Vertices(1);   Tdomain%sVertex(nv)%PML = .false.
               nv = Tdomain%specel(n)%Near_Vertices(2);   Tdomain%sVertex(nv)%PML = .false.
               nv = Tdomain%specel(n)%Near_Vertices(5);   Tdomain%sVertex(nv)%PML = .false.
               nv = Tdomain%specel(n)%Near_Vertices(6);   Tdomain%sVertex(nv)%PML = .false.
            else
               nf = Tdomain%specel(n)%Near_Faces(4);   Tdomain%sFace(nf)%PML = .false.
               ne = Tdomain%specel(n)%Near_Edges(3);   Tdomain%sEdge(ne)%PML = .false.
               ne = Tdomain%specel(n)%Near_Edges(6);   Tdomain%sEdge(ne)%PML = .false.
               ne = Tdomain%specel(n)%Near_Edges(10);   Tdomain%sEdge(ne)%PML = .false.
               ne = Tdomain%specel(n)%Near_Edges(11);   Tdomain%sEdge(ne)%PML = .false.
               nv = Tdomain%specel(n)%Near_Vertices(0);   Tdomain%sVertex(nv)%PML = .false.
               nv = Tdomain%specel(n)%Near_Vertices(3);   Tdomain%sVertex(nv)%PML = .false.
               nv = Tdomain%specel(n)%Near_Vertices(4);   Tdomain%sVertex(nv)%PML = .false.
               nv = Tdomain%specel(n)%Near_Vertices(7);   Tdomain%sVertex(nv)%PML = .false.
            endif
          endif
       endif
    else
       if (Tdomain%sSubdomain(mat)%Py) then
          if (Tdomain%sSubdomain(mat)%Pz) then
             if (Tdomain%sSubdomain(mat)%Forward) then
                if (Tdomain%sSubdomain(mat)%Down) then
                   ne = Tdomain%specel(n)%Near_Edges(9);   Tdomain%sEdge(ne)%PML = .false.
                   nv = Tdomain%specel(n)%Near_Vertices(6);   Tdomain%sVertex(nv)%PML = .false.
                   nv = Tdomain%specel(n)%Near_Vertices(7);   Tdomain%sVertex(nv)%PML = .false.
                else
                   ne = Tdomain%specel(n)%Near_Edges(2);   Tdomain%sEdge(ne)%PML = .false.
                   nv = Tdomain%specel(n)%Near_Vertices(2);   Tdomain%sVertex(nv)%PML = .false.
                   nv = Tdomain%specel(n)%Near_Vertices(3);   Tdomain%sVertex(nv)%PML = .false.
                endif
             else
                if (Tdomain%sSubdomain(mat)%Down) then
                   ne = Tdomain%specel(n)%Near_Edges(5);   Tdomain%sEdge(ne)%PML = .false.
                   nv = Tdomain%specel(n)%Near_Vertices(4);   Tdomain%sVertex(nv)%PML = .false.
                   nv = Tdomain%specel(n)%Near_Vertices(5);   Tdomain%sVertex(nv)%PML = .false.
                else
                   ne = Tdomain%specel(n)%Near_Edges(0);   Tdomain%sEdge(ne)%PML = .false.
                   nv = Tdomain%specel(n)%Near_Vertices(0);   Tdomain%sVertex(nv)%PML = .false.
                   nv = Tdomain%specel(n)%Near_Vertices(1);   Tdomain%sVertex(nv)%PML = .false.
                endif
             endif
          else
             if (Tdomain%sSubdomain(mat)%Forward) then
                nf = Tdomain%specel(n)%Near_Faces(3);   Tdomain%sFace(nf)%PML = .false.
                ne = Tdomain%specel(n)%Near_Edges(2);   Tdomain%sEdge(ne)%PML = .false.
                ne = Tdomain%specel(n)%Near_Edges(7);   Tdomain%sEdge(ne)%PML = .false.
                ne = Tdomain%specel(n)%Near_Edges(9);   Tdomain%sEdge(ne)%PML = .false.
                ne = Tdomain%specel(n)%Near_Edges(10);   Tdomain%sEdge(ne)%PML = .false.
                nv = Tdomain%specel(n)%Near_Vertices(2);   Tdomain%sVertex(nv)%PML = .false.
                nv = Tdomain%specel(n)%Near_Vertices(3);   Tdomain%sVertex(nv)%PML = .false.
                nv = Tdomain%specel(n)%Near_Vertices(6);   Tdomain%sVertex(nv)%PML = .false.
                nv = Tdomain%specel(n)%Near_Vertices(7);   Tdomain%sVertex(nv)%PML = .false.
             else
                nf = Tdomain%specel(n)%Near_Faces(1);   Tdomain%sFace(nf)%PML = .false.
                ne = Tdomain%specel(n)%Near_Edges(0);   Tdomain%sEdge(ne)%PML = .false.
                ne = Tdomain%specel(n)%Near_Edges(4);   Tdomain%sEdge(ne)%PML = .false.
                ne = Tdomain%specel(n)%Near_Edges(5);   Tdomain%sEdge(ne)%PML = .false.
                ne = Tdomain%specel(n)%Near_Edges(6);   Tdomain%sEdge(ne)%PML = .false.
                nv = Tdomain%specel(n)%Near_Vertices(0);   Tdomain%sVertex(nv)%PML = .false.
                nv = Tdomain%specel(n)%Near_Vertices(1);   Tdomain%sVertex(nv)%PML = .false.
                nv = Tdomain%specel(n)%Near_Vertices(4);   Tdomain%sVertex(nv)%PML = .false.
                nv = Tdomain%specel(n)%Near_Vertices(5);   Tdomain%sVertex(nv)%PML = .false.
             endif
          endif
       else
          if (Tdomain%sSubdomain(mat)%Pz) then
             if (Tdomain%sSubdomain(mat)%Down) then
                nf = Tdomain%specel(n)%Near_Faces(5);   Tdomain%sFace(nf)%PML = .false.
                ne = Tdomain%specel(n)%Near_Edges(5);   Tdomain%sEdge(ne)%PML = .false.
                ne = Tdomain%specel(n)%Near_Edges(8);   Tdomain%sEdge(ne)%PML = .false.
                ne = Tdomain%specel(n)%Near_Edges(9);   Tdomain%sEdge(ne)%PML = .false.
                ne = Tdomain%specel(n)%Near_Edges(11);   Tdomain%sEdge(ne)%PML = .false.
                nv = Tdomain%specel(n)%Near_Vertices(4);   Tdomain%sVertex(nv)%PML = .false.
                nv = Tdomain%specel(n)%Near_Vertices(5);   Tdomain%sVertex(nv)%PML = .false.
                nv = Tdomain%specel(n)%Near_Vertices(6);   Tdomain%sVertex(nv)%PML = .false.
                nv = Tdomain%specel(n)%Near_Vertices(7);   Tdomain%sVertex(nv)%PML = .false.
             else
                nf = Tdomain%specel(n)%Near_Faces(0);   Tdomain%sFace(nf)%PML = .false.
                ne = Tdomain%specel(n)%Near_Edges(0);   Tdomain%sEdge(ne)%PML = .false.
                ne = Tdomain%specel(n)%Near_Edges(1);   Tdomain%sEdge(ne)%PML = .false.
                ne = Tdomain%specel(n)%Near_Edges(2);   Tdomain%sEdge(ne)%PML = .false.
                ne = Tdomain%specel(n)%Near_Edges(3);   Tdomain%sEdge(ne)%PML = .false.
                nv = Tdomain%specel(n)%Near_Vertices(0);   Tdomain%sVertex(nv)%PML = .false.
                nv = Tdomain%specel(n)%Near_Vertices(1);   Tdomain%sVertex(nv)%PML = .false.
                nv = Tdomain%specel(n)%Near_Vertices(2);   Tdomain%sVertex(nv)%PML = .false.
                nv = Tdomain%specel(n)%Near_Vertices(3);   Tdomain%sVertex(nv)%PML = .false.
             endif
          else
             print *,"Pb: There is a PML element which is neither Px nor Py nor Pz !?!"
          endif
       endif
    endif
 else
    do i = 0,5
        nf = Tdomain%specel(n)%Near_Faces(i)
        Tdomain%sFace(nf)%PML = .false.
    enddo
    do i = 0,11
        ne = Tdomain%specel(n)%Near_Edges(i)
        Tdomain%sEdge(ne)%PML = .false.
    enddo
    do i = 0,7
        nv = Tdomain%specel(n)%Near_Vertices(i)
        Tdomain%sVertex(nv)%PML = .false.
    enddo
 endif
enddo

do n = 0,Tdomain%n_elem-1
    mat = Tdomain%specel(n)%mat_index
    if (Tdomain%sSubdomain(mat)%Px) then
        if (Tdomain%sSubdomain(mat)%Left) then
            nf = Tdomain%specel(n)%Near_Faces(4);   Tdomain%sFace(nf)%Abs = .true.
            ne = Tdomain%specel(n)%Near_Edges(3);   Tdomain%sEdge(ne)%Abs = .true.
            ne = Tdomain%specel(n)%Near_Edges(6);   Tdomain%sEdge(ne)%Abs = .true.
            ne = Tdomain%specel(n)%Near_Edges(10);   Tdomain%sEdge(ne)%Abs = .true.
            ne = Tdomain%specel(n)%Near_Edges(11);   Tdomain%sEdge(ne)%Abs = .true.
            nv = Tdomain%specel(n)%Near_Vertices(0);   Tdomain%sVertex(nv)%Abs = .true.
            nv = Tdomain%specel(n)%Near_Vertices(3);   Tdomain%sVertex(nv)%Abs = .true.
            nv = Tdomain%specel(n)%Near_Vertices(4);   Tdomain%sVertex(nv)%Abs = .true.
            nv = Tdomain%specel(n)%Near_Vertices(7);   Tdomain%sVertex(nv)%Abs = .true.
        else
            nf = Tdomain%specel(n)%Near_Faces(2);   Tdomain%sFace(nf)%Abs = .true.
            ne = Tdomain%specel(n)%Near_Edges(1);   Tdomain%sEdge(ne)%Abs = .true.
            ne = Tdomain%specel(n)%Near_Edges(4);   Tdomain%sEdge(ne)%Abs = .true.
            ne = Tdomain%specel(n)%Near_Edges(7);   Tdomain%sEdge(ne)%Abs = .true.
            ne = Tdomain%specel(n)%Near_Edges(8);   Tdomain%sEdge(ne)%Abs = .true.
            nv = Tdomain%specel(n)%Near_Vertices(1);   Tdomain%sVertex(nv)%Abs = .true.
            nv = Tdomain%specel(n)%Near_Vertices(2);   Tdomain%sVertex(nv)%Abs = .true.
            nv = Tdomain%specel(n)%Near_Vertices(5);   Tdomain%sVertex(nv)%Abs = .true.
            nv = Tdomain%specel(n)%Near_Vertices(6);   Tdomain%sVertex(nv)%Abs = .true.
        endif
    endif
    if (Tdomain%sSubdomain(mat)%Py) then
        if (Tdomain%sSubdomain(mat)%Forward) then
            nf = Tdomain%specel(n)%Near_Faces(1);   Tdomain%sFace(nf)%Abs = .true.
            ne = Tdomain%specel(n)%Near_Edges(0);   Tdomain%sEdge(ne)%Abs = .true.
            ne = Tdomain%specel(n)%Near_Edges(4);   Tdomain%sEdge(ne)%Abs = .true.
            ne = Tdomain%specel(n)%Near_Edges(5);   Tdomain%sEdge(ne)%Abs = .true.
            ne = Tdomain%specel(n)%Near_Edges(6);   Tdomain%sEdge(ne)%Abs = .true.
            nv = Tdomain%specel(n)%Near_Vertices(0);   Tdomain%sVertex(nv)%Abs = .true.
            nv = Tdomain%specel(n)%Near_Vertices(1);   Tdomain%sVertex(nv)%Abs = .true.
            nv = Tdomain%specel(n)%Near_Vertices(4);   Tdomain%sVertex(nv)%Abs = .true.
            nv = Tdomain%specel(n)%Near_Vertices(5);   Tdomain%sVertex(nv)%Abs = .true.
        else
            nf = Tdomain%specel(n)%Near_Faces(3);   Tdomain%sFace(nf)%Abs = .true.
            ne = Tdomain%specel(n)%Near_Edges(2);   Tdomain%sEdge(ne)%Abs = .true.
            ne = Tdomain%specel(n)%Near_Edges(7);   Tdomain%sEdge(ne)%Abs = .true.
            ne = Tdomain%specel(n)%Near_Edges(9);   Tdomain%sEdge(ne)%Abs = .true.
            ne = Tdomain%specel(n)%Near_Edges(10);   Tdomain%sEdge(ne)%Abs = .true.
            nv = Tdomain%specel(n)%Near_Vertices(2);   Tdomain%sVertex(nv)%Abs = .true.
            nv = Tdomain%specel(n)%Near_Vertices(3);   Tdomain%sVertex(nv)%Abs = .true.
            nv = Tdomain%specel(n)%Near_Vertices(6);   Tdomain%sVertex(nv)%Abs = .true.
            nv = Tdomain%specel(n)%Near_Vertices(7);   Tdomain%sVertex(nv)%Abs = .true.
        endif
    endif
    if (Tdomain%sSubdomain(mat)%Pz) then
        if (Tdomain%sSubdomain(mat)%Down) then
            nf = Tdomain%specel(n)%Near_Faces(0);   Tdomain%sFace(nf)%Abs = .true.
            ne = Tdomain%specel(n)%Near_Edges(0);   Tdomain%sEdge(ne)%Abs = .true.
            ne = Tdomain%specel(n)%Near_Edges(1);   Tdomain%sEdge(ne)%Abs = .true.
            ne = Tdomain%specel(n)%Near_Edges(2);   Tdomain%sEdge(ne)%Abs = .true.
            ne = Tdomain%specel(n)%Near_Edges(3);   Tdomain%sEdge(ne)%Abs = .true.
            nv = Tdomain%specel(n)%Near_Vertices(0);   Tdomain%sVertex(nv)%Abs = .true.
            nv = Tdomain%specel(n)%Near_Vertices(1);   Tdomain%sVertex(nv)%Abs = .true.
            nv = Tdomain%specel(n)%Near_Vertices(2);   Tdomain%sVertex(nv)%Abs = .true.
            nv = Tdomain%specel(n)%Near_Vertices(3);   Tdomain%sVertex(nv)%Abs = .true.
        else
            nf = Tdomain%specel(n)%Near_Faces(5);   Tdomain%sFace(nf)%Abs = .true.
            ne = Tdomain%specel(n)%Near_Edges(5);   Tdomain%sEdge(ne)%Abs = .true.
            ne = Tdomain%specel(n)%Near_Edges(8);   Tdomain%sEdge(ne)%Abs = .true.
            ne = Tdomain%specel(n)%Near_Edges(9);   Tdomain%sEdge(ne)%Abs = .true.
            ne = Tdomain%specel(n)%Near_Edges(11);   Tdomain%sEdge(ne)%Abs = .true.
            nv = Tdomain%specel(n)%Near_Vertices(4);   Tdomain%sVertex(nv)%Abs = .true.
            nv = Tdomain%specel(n)%Near_Vertices(5);   Tdomain%sVertex(nv)%Abs = .true.
            nv = Tdomain%specel(n)%Near_Vertices(6);   Tdomain%sVertex(nv)%Abs = .true.
            nv = Tdomain%specel(n)%Near_Vertices(7);   Tdomain%sVertex(nv)%Abs = .true.
        endif
    endif
enddo


return
end subroutine PML_definition

