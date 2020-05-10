subroutine get_Mass_Elem2Obj (Tdomain, n)


use sdomains

implicit none

type(domain), intent(INOUT) :: Tdomain
integer, intent(IN) :: n

integer :: ngllx,nglly,ngllz, ngll1,ngll2, ngll, i, nf,ne,nv


    ngllx = Tdomain%specel(n)%ngllx
    nglly = Tdomain%specel(n)%nglly
    ngllz = Tdomain%specel(n)%ngllz

    ! Faces
    do i = 0,5
        nf = Tdomain%specel(n)%Near_Faces(i)
        ngll1 = Tdomain%sFace(nf)%ngll1
        ngll2 = Tdomain%sFace(nf)%ngll2
        select case (i)
         case (0)
            Tdomain%sFace(nf)%MassMat(:,:) = Tdomain%sFace(nf)%MassMat(:,:) + &
                                             Tdomain%specel(n)%MassMat(1:ngll1-2,1:ngll2-2,0)
         case (1)
            Tdomain%sFace(nf)%MassMat(:,:) = Tdomain%sFace(nf)%MassMat(:,:) + &
                                             Tdomain%specel(n)%MassMat(1:ngll1-2,0,1:ngll2-2)
         case (2)
            Tdomain%sFace(nf)%MassMat(:,:) = Tdomain%sFace(nf)%MassMat(:,:) + &
                                             Tdomain%specel(n)%MassMat(ngllx-1,1:ngll1-2,1:ngll2-2)
         case (3)
            Tdomain%sFace(nf)%MassMat(:,:) = Tdomain%sFace(nf)%MassMat(:,:) + &
                                             Tdomain%specel(n)%MassMat(1:ngll1-2,nglly-1,1:ngll2-2)
         case (4)
            Tdomain%sFace(nf)%MassMat(:,:) = Tdomain%sFace(nf)%MassMat(:,:) + &
                                             Tdomain%specel(n)%MassMat(0,1:ngll1-2,1:ngll2-2)
         case (5)
            Tdomain%sFace(nf)%MassMat(:,:) = Tdomain%sFace(nf)%MassMat(:,:) + &
                                             Tdomain%specel(n)%MassMat(1:ngll1-2,1:ngll2-2,ngllz-1)
            if (Tdomain%sFace(nf)%ocean) then
                Tdomain%sFace(nf)%Mocean(:,:) = Tdomain%specel(n)%Mocean(1:ngll1-2,1:ngll2-2)
            endif
        end select
        if (Tdomain%sFace(nf)%PML) then
            select case (i)
             case (0)
                Tdomain%sFace(nf)%DumpMass(:,:,:) = Tdomain%sFace(nf)%DumpMass(:,:,:) + &
                                                  Tdomain%specel(n)%DumpMass(1:ngll1-2,1:ngll2-2,0,:)
             case (1)
                Tdomain%sFace(nf)%DumpMass(:,:,:) = Tdomain%sFace(nf)%DumpMass(:,:,:) + &
                                                  Tdomain%specel(n)%DumpMass(1:ngll1-2,0,1:ngll2-2,:)
             case (2)
                Tdomain%sFace(nf)%DumpMass(:,:,:) = Tdomain%sFace(nf)%DumpMass(:,:,:) + &
                                                  Tdomain%specel(n)%DumpMass(ngllx-1,1:ngll1-2,1:ngll2-2,:)
             case (3)
                Tdomain%sFace(nf)%DumpMass(:,:,:) = Tdomain%sFace(nf)%DumpMass(:,:,:) + &
                                                  Tdomain%specel(n)%DumpMass(1:ngll1-2,nglly-1,1:ngll2-2,:)
             case (4)
                Tdomain%sFace(nf)%DumpMass(:,:,:) = Tdomain%sFace(nf)%DumpMass(:,:,:) + &
                                                  Tdomain%specel(n)%DumpMass(0,1:ngll1-2,1:ngll2-2,:)
             case (5)
                Tdomain%sFace(nf)%DumpMass(:,:,:) = Tdomain%sFace(nf)%DumpMass(:,:,:) + &
                                                  Tdomain%specel(n)%DumpMass(1:ngll1-2,1:ngll2-2,ngllz-1,:)
            end select
        endif
    enddo

    ! Edges
    do i = 0,11
        ne = Tdomain%specel(n)%Near_Edges(i)
        ngll = Tdomain%sEdge(ne)%ngll
        select case (i)
         case (0)
            Tdomain%sEdge(ne)%MassMat(:) = Tdomain%sEdge(ne)%MassMat(:) + &
                                           Tdomain%specel(n)%MassMat(1:ngll-2,0,0)
         case (1)
            Tdomain%sEdge(ne)%MassMat(:) = Tdomain%sEdge(ne)%MassMat(:) + &
                                           Tdomain%specel(n)%MassMat(ngllx-1,1:ngll-2,0)
         case (2)
            Tdomain%sEdge(ne)%MassMat(:) = Tdomain%sEdge(ne)%MassMat(:) + &
                                           Tdomain%specel(n)%MassMat(1:ngll-2,nglly-1,0)
         case (3)
            Tdomain%sEdge(ne)%MassMat(:) = Tdomain%sEdge(ne)%MassMat(:) + &
                                           Tdomain%specel(n)%MassMat(0,1:ngll-2,0)
         case (4)
            Tdomain%sEdge(ne)%MassMat(:) = Tdomain%sEdge(ne)%MassMat(:) + &
                                           Tdomain%specel(n)%MassMat(ngllx-1,0,1:ngll-2)
         case (5)
            Tdomain%sEdge(ne)%MassMat(:) = Tdomain%sEdge(ne)%MassMat(:) + &
                                           Tdomain%specel(n)%MassMat(1:ngll-2,0,ngllz-1)
            if (Tdomain%sEdge(ne)%ocean) Tdomain%sEdge(ne)%Mocean(:) = Tdomain%sEdge(ne)%Mocean(:) + &
                                         Tdomain%specel(n)%Mocean(1:ngll-2,0)
         case (6)
            Tdomain%sEdge(ne)%MassMat(:) = Tdomain%sEdge(ne)%MassMat(:) + &
                                           Tdomain%specel(n)%MassMat(0,0,1:ngll-2)
         case (7)
            Tdomain%sEdge(ne)%MassMat(:) = Tdomain%sEdge(ne)%MassMat(:) + &
                                           Tdomain%specel(n)%MassMat(ngllx-1,nglly-1,1:ngll-2)
         case (8)
            Tdomain%sEdge(ne)%MassMat(:) = Tdomain%sEdge(ne)%MassMat(:) + &
                                           Tdomain%specel(n)%MassMat(ngllx-1,1:ngll-2,ngllz-1)
            if (Tdomain%sEdge(ne)%ocean) Tdomain%sEdge(ne)%Mocean(:) = Tdomain%sEdge(ne)%Mocean(:) + &
                                         Tdomain%specel(n)%Mocean(ngllx-1,1:ngll-2)
         case (9)
            Tdomain%sEdge(ne)%MassMat(:) = Tdomain%sEdge(ne)%MassMat(:) + &
                                           Tdomain%specel(n)%MassMat(1:ngll-2,nglly-1,ngllz-1)
            if (Tdomain%sEdge(ne)%ocean) Tdomain%sEdge(ne)%Mocean(:) = Tdomain%sEdge(ne)%Mocean(:) + &
                                         Tdomain%specel(n)%Mocean(1:ngll-2,nglly-1)
         case (10)
            Tdomain%sEdge(ne)%MassMat(:) = Tdomain%sEdge(ne)%MassMat(:) + &
                                           Tdomain%specel(n)%MassMat(0,nglly-1,1:ngll-2)
         case (11)
            Tdomain%sEdge(ne)%MassMat(:) = Tdomain%sEdge(ne)%MassMat(:) + &
                                           Tdomain%specel(n)%MassMat(0,1:ngll-2,ngllz-1)
            if (Tdomain%sEdge(ne)%ocean) Tdomain%sEdge(ne)%Mocean(:) = Tdomain%sEdge(ne)%Mocean(:) + &
                                         Tdomain%specel(n)%Mocean(0,1:ngll-2)
        end select
        if (Tdomain%sEdge(ne)%PML) then
            select case (i)
             case (0)
                Tdomain%sEdge(ne)%DumpMass(:,:) = Tdomain%sEdge(ne)%DumpMass(:,:) + &
                                                  Tdomain%specel(n)%DumpMass(1:ngll-2,0,0,:)
             case (1)
                Tdomain%sEdge(ne)%DumpMass(:,:) = Tdomain%sEdge(ne)%DumpMass(:,:) + &
                                                  Tdomain%specel(n)%DumpMass(ngllx-1,1:ngll-2,0,:)
             case (2)
                Tdomain%sEdge(ne)%DumpMass(:,:) = Tdomain%sEdge(ne)%DumpMass(:,:) + &
                                                  Tdomain%specel(n)%DumpMass(1:ngll-2,nglly-1,0,:)
             case (3)
                Tdomain%sEdge(ne)%DumpMass(:,:) = Tdomain%sEdge(ne)%DumpMass(:,:) + &
                                                  Tdomain%specel(n)%DumpMass(0,1:ngll-2,0,:)
             case (4)
                Tdomain%sEdge(ne)%DumpMass(:,:) = Tdomain%sEdge(ne)%DumpMass(:,:) + &
                                                  Tdomain%specel(n)%DumpMass(ngllx-1,0,1:ngll-2,:)
             case (5)
                Tdomain%sEdge(ne)%DumpMass(:,:) = Tdomain%sEdge(ne)%DumpMass(:,:) + &
                                                  Tdomain%specel(n)%DumpMass(1:ngll-2,0,ngllz-1,:)
             case (6)
                Tdomain%sEdge(ne)%DumpMass(:,:) = Tdomain%sEdge(ne)%DumpMass(:,:) + &
                                                  Tdomain%specel(n)%DumpMass(0,0,1:ngll-2,:)
             case (7)
                Tdomain%sEdge(ne)%DumpMass(:,:) = Tdomain%sEdge(ne)%DumpMass(:,:) + &
                                                  Tdomain%specel(n)%DumpMass(ngllx-1,nglly-1,1:ngll-2,:)
             case (8)
                Tdomain%sEdge(ne)%DumpMass(:,:) = Tdomain%sEdge(ne)%DumpMass(:,:) + &
                                                  Tdomain%specel(n)%DumpMass(ngllx-1,1:ngll-2,ngllz-1,:)
             case (9)
                Tdomain%sEdge(ne)%DumpMass(:,:) = Tdomain%sEdge(ne)%DumpMass(:,:) + &
                                                  Tdomain%specel(n)%DumpMass(1:ngll-2,nglly-1,ngllz-1,:)
             case (10)
                Tdomain%sEdge(ne)%DumpMass(:,:) = Tdomain%sEdge(ne)%DumpMass(:,:) + &
                                                  Tdomain%specel(n)%DumpMass(0,nglly-1,1:ngll-2,:)
             case (11)
                Tdomain%sEdge(ne)%DumpMass(:,:) = Tdomain%sEdge(ne)%DumpMass(:,:) + &
                                                  Tdomain%specel(n)%DumpMass(0,1:ngll-2,ngllz-1,:)
            end select
        endif
    enddo

    ! Vertices
    do i = 0,7
        nv = Tdomain%specel(n)%Near_Vertices(i)
        select case (i)
         case (0)
            Tdomain%sVertex(nv)%MassMat = Tdomain%sVertex(nv)%MassMat + &
                                          Tdomain%specel(n)%MassMat(0,0,0)
         case (1)
            Tdomain%sVertex(nv)%MassMat = Tdomain%sVertex(nv)%MassMat + &
                                          Tdomain%specel(n)%MassMat(ngllx-1,0,0)
         case (2)
            Tdomain%sVertex(nv)%MassMat = Tdomain%sVertex(nv)%MassMat + &
                                          Tdomain%specel(n)%MassMat(ngllx-1,nglly-1,0)
         case (3)
            Tdomain%sVertex(nv)%MassMat = Tdomain%sVertex(nv)%MassMat + &
                                          Tdomain%specel(n)%MassMat(0,nglly-1,0)
         case (4)
            Tdomain%sVertex(nv)%MassMat = Tdomain%sVertex(nv)%MassMat + &
                                          Tdomain%specel(n)%MassMat(0,0,ngllz-1)
            if (Tdomain%sVertex(nv)%ocean) Tdomain%sVertex(nv)%Mocean = Tdomain%sVertex(nv)%Mocean + &
                                           Tdomain%specel(n)%Mocean(0,0)
         case (5)
            Tdomain%sVertex(nv)%MassMat = Tdomain%sVertex(nv)%MassMat + &
                                          Tdomain%specel(n)%MassMat(ngllx-1,0,ngllz-1)
            if (Tdomain%sVertex(nv)%ocean) Tdomain%sVertex(nv)%Mocean = Tdomain%sVertex(nv)%Mocean + &
                                            Tdomain%specel(n)%Mocean(ngllx-1,0)
         case (6)
            Tdomain%sVertex(nv)%MassMat = Tdomain%sVertex(nv)%MassMat + &
                                          Tdomain%specel(n)%MassMat(ngllx-1,nglly-1,ngllz-1)
            if (Tdomain%sVertex(nv)%ocean) Tdomain%sVertex(nv)%Mocean = Tdomain%sVertex(nv)%Mocean + &
                                            Tdomain%specel(n)%Mocean(ngllx-1,nglly-1)
         case (7)
            Tdomain%sVertex(nv)%MassMat = Tdomain%sVertex(nv)%MassMat + &
                                          Tdomain%specel(n)%MassMat(0,nglly-1,ngllz-1)
            if (Tdomain%sVertex(nv)%ocean) Tdomain%sVertex(nv)%Mocean = Tdomain%sVertex(nv)%Mocean + &
                                           Tdomain%specel(n)%Mocean(0,nglly-1)
        end select
        if (Tdomain%sVertex(nv)%PML) then
            select case (i)
             case (0)
                Tdomain%sVertex(nv)%DumpMass(:) = Tdomain%sVertex(nv)%DumpMass(:) + &
                                              Tdomain%specel(n)%DumpMass(0,0,0,:)
             case (1)
                Tdomain%sVertex(nv)%DumpMass(:) = Tdomain%sVertex(nv)%DumpMass(:) + &
                                              Tdomain%specel(n)%DumpMass(ngllx-1,0,0,:)
             case (2)
                Tdomain%sVertex(nv)%DumpMass(:) = Tdomain%sVertex(nv)%DumpMass(:) + &
                                              Tdomain%specel(n)%DumpMass(ngllx-1,nglly-1,0,:)
             case (3)
                Tdomain%sVertex(nv)%DumpMass(:) = Tdomain%sVertex(nv)%DumpMass(:) + &
                                              Tdomain%specel(n)%DumpMass(0,nglly-1,0,:)
             case (4)
                Tdomain%sVertex(nv)%DumpMass(:) = Tdomain%sVertex(nv)%DumpMass(:) + &
                                              Tdomain%specel(n)%DumpMass(0,0,ngllz-1,:)
             case (5)
                Tdomain%sVertex(nv)%DumpMass(:) = Tdomain%sVertex(nv)%DumpMass(:) + &
                                              Tdomain%specel(n)%DumpMass(ngllx-1,0,ngllz-1,:)
             case (6)
                Tdomain%sVertex(nv)%DumpMass(:) = Tdomain%sVertex(nv)%DumpMass(:) + &
                                              Tdomain%specel(n)%DumpMass(ngllx-1,nglly-1,ngllz-1,:)
             case (7)
                Tdomain%sVertex(nv)%DumpMass(:) = Tdomain%sVertex(nv)%DumpMass(:) + &
                                              Tdomain%specel(n)%DumpMass(0,nglly-1,ngllz-1,:)
            end select
        endif
    enddo

    if (Tdomain%specel(n)%ocean .eqv. .true.)   deallocate (Tdomain%specel(n)%Mocean)


end subroutine get_Mass_Elem2Obj
