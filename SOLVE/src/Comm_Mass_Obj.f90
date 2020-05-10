subroutine Comm_Mass_Obj (Tdomain, n)


use sdomains

implicit none

type(domain), intent(INOUT) :: Tdomain
integer, intent(IN) :: n

integer :: i,j,k, nf,ne,nv, ngll,ngllPML,ngllocean


    ngll = 0
    ngllPML = 0
    ngllocean = 0

    ! Faces
    do i = 0,Tdomain%sComm(n)%nb_faces-1
        nf = Tdomain%sComm(n)%faces(i)
        do j = 1,Tdomain%sFace(nf)%ngll2-2
            do k = 1,Tdomain%sFace(nf)%ngll1-2
                Tdomain%sFace(nf)%MassMat(k,j) = Tdomain%sFace(nf)%MassMat(k,j) + Tdomain%sComm(n)%Take(ngll)
                ngll = ngll + 1
            enddo
        enddo
        if (Tdomain%sFace(nf)%PML) then
            do j = 1,Tdomain%sFace(nf)%ngll2-2
                do k = 1,Tdomain%sFace(nf)%ngll1-2
                    Tdomain%sFace(nf)%DumpMass(k,j,0:2) = Tdomain%sFace(nf)%DumpMass(k,j,0:2) + Tdomain%sComm(n)%TakePML(ngllPML,0:2)
                    ngllPML = ngllPML + 1
                enddo
            enddo
        endif
    enddo

    ! Edges
    do i = 0,Tdomain%sComm(n)%nb_edges-1
        ne = Tdomain%sComm(n)%edges(i)
        do j = 1,Tdomain%sEdge(ne)%ngll-2
            Tdomain%sEdge(ne)%MassMat(j) = Tdomain%sEdge(ne)%MassMat(j) + Tdomain%sComm(n)%Take(ngll)
            ngll = ngll + 1
        enddo
        if (Tdomain%sEdge(ne)%ocean) then
            do j = 1,Tdomain%sEdge(ne)%ngll-2
                Tdomain%sEdge(ne)%Mocean(j) = Tdomain%sEdge(ne)%Mocean(j) + Tdomain%sComm(n)%Takeocean(ngllocean)
                ngllocean = ngllocean + 1
            enddo
        endif
        if (Tdomain%sEdge(ne)%PML) then
            do j = 1,Tdomain%sEdge(ne)%ngll-2
                Tdomain%sEdge(ne)%DumpMass(j,0:2) = Tdomain%sEdge(ne)%DumpMass(j,0:2) + Tdomain%sComm(n)%TakePML(ngllPML,0:2)
                ngllPML = ngllPML + 1
            enddo
        endif
    enddo

    ! Vertices
    do i = 0,Tdomain%sComm(n)%nb_vertices-1
        nv = Tdomain%sComm(n)%vertices(i)
        Tdomain%sVertex(nv)%MassMat = Tdomain%sVertex(nv)%MassMat + Tdomain%sComm(n)%Take(ngll)
        ngll = ngll + 1
        if (Tdomain%sVertex(nv)%ocean) then
            Tdomain%sVertex(nv)%Mocean = Tdomain%sVertex(nv)%Mocean + Tdomain%sComm(n)%Takeocean(ngllocean)
            ngllocean = ngllocean + 1
        endif
        if (Tdomain%sVertex(nv)%PML) then
            Tdomain%sVertex(nv)%DumpMass(0:2) = Tdomain%sVertex(nv)%DumpMass(0:2) + Tdomain%sComm(n)%TakePML(ngllPML,0:2)
            ngllPML = ngllPML + 1
        endif
    enddo


end subroutine Comm_Mass_Obj
