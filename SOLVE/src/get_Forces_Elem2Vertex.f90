subroutine get_Forces_Elem2Vertex (Elem, Vert, i, i_simu)


use selements
use svertices

implicit none

type (Element), intent (INOUT) :: Elem
type (Vertex), intent (INOUT) :: Vert
integer, intent (IN) :: i, i_simu

integer :: ngllx,nglly,ngllz


       ngllx = Elem%ngllx
       nglly = Elem%nglly
       ngllz = Elem%ngllz

       select case (i)
        case (0)
           Vert%sSimu(i_simu)%Forces(0:2) = Vert%sSimu(i_simu)%Forces(0:2) + &
                                            Elem%sSimu(i_simu)%Forces(0,0,0,0:2)
           if (Vert%PML .and. Elem%PML) then
            Vert%sSimu(i_simu)%Forces1(0:2) = Vert%sSimu(i_simu)%Forces1(0:2) + &
                                              Elem%sSimu(i_simu)%Forces1(0,0,0,0:2)
            Vert%sSimu(i_simu)%Forces2(0:2) = Vert%sSimu(i_simu)%Forces2(0:2) + &
                                              Elem%sSimu(i_simu)%Forces2(0,0,0,0:2)
            Vert%sSimu(i_simu)%Forces3(0:2) = Vert%sSimu(i_simu)%Forces3(0:2) + &
                                              Elem%sSimu(i_simu)%Forces3(0,0,0,0:2)
           endif
        case (1)
           Vert%sSimu(i_simu)%Forces(0:2) = Vert%sSimu(i_simu)%Forces(0:2) + &
                                            Elem%sSimu(i_simu)%Forces(ngllx-1,0,0,0:2)
           if (Vert%PML .and. Elem%PML) then
            Vert%sSimu(i_simu)%Forces1(0:2) = Vert%sSimu(i_simu)%Forces1(0:2) + &
                                              Elem%sSimu(i_simu)%Forces1(ngllx-1,0,0,0:2)
            Vert%sSimu(i_simu)%Forces2(0:2) = Vert%sSimu(i_simu)%Forces2(0:2) + &
                                              Elem%sSimu(i_simu)%Forces2(ngllx-1,0,0,0:2)
            Vert%sSimu(i_simu)%Forces3(0:2) = Vert%sSimu(i_simu)%Forces3(0:2) + &
                                              Elem%sSimu(i_simu)%Forces3(ngllx-1,0,0,0:2)
           endif
        case (2)
           Vert%sSimu(i_simu)%Forces(0:2) = Vert%sSimu(i_simu)%Forces(0:2) + &
                                            Elem%sSimu(i_simu)%Forces(ngllx-1,nglly-1,0,0:2)
           if (Vert%PML .and. Elem%PML) then
            Vert%sSimu(i_simu)%Forces1(0:2) = Vert%sSimu(i_simu)%Forces1(0:2) + &
                                              Elem%sSimu(i_simu)%Forces1(ngllx-1,nglly-1,0,0:2)
            Vert%sSimu(i_simu)%Forces2(0:2) = Vert%sSimu(i_simu)%Forces2(0:2) + &
                                              Elem%sSimu(i_simu)%Forces2(ngllx-1,nglly-1,0,0:2)
            Vert%sSimu(i_simu)%Forces3(0:2) = Vert%sSimu(i_simu)%Forces3(0:2) + &
                                              Elem%sSimu(i_simu)%Forces3(ngllx-1,nglly-1,0,0:2)
           endif
        case (3)
           Vert%sSimu(i_simu)%Forces(0:2) = Vert%sSimu(i_simu)%Forces(0:2) + &
                                            Elem%sSimu(i_simu)%Forces(0,nglly-1,0,0:2)
           if (Vert%PML .and. Elem%PML) then
            Vert%sSimu(i_simu)%Forces1(0:2) = Vert%sSimu(i_simu)%Forces1(0:2) + &
                                              Elem%sSimu(i_simu)%Forces1(0,nglly-1,0,0:2)
            Vert%sSimu(i_simu)%Forces2(0:2) = Vert%sSimu(i_simu)%Forces2(0:2) + &
                                              Elem%sSimu(i_simu)%Forces2(0,nglly-1,0,0:2)
            Vert%sSimu(i_simu)%Forces3(0:2) = Vert%sSimu(i_simu)%Forces3(0:2) + &
                                              Elem%sSimu(i_simu)%Forces3(0,nglly-1,0,0:2)
           endif
        case (4)
           Vert%sSimu(i_simu)%Forces(0:2) = Vert%sSimu(i_simu)%Forces(0:2) + &
                                            Elem%sSimu(i_simu)%Forces(0,0,ngllz-1,0:2)
           if (Vert%PML .and. Elem%PML) then
            Vert%sSimu(i_simu)%Forces1(0:2) = Vert%sSimu(i_simu)%Forces1(0:2) + &
                                              Elem%sSimu(i_simu)%Forces1(0,0,ngllz-1,0:2)
            Vert%sSimu(i_simu)%Forces2(0:2) = Vert%sSimu(i_simu)%Forces2(0:2) + &
                                              Elem%sSimu(i_simu)%Forces2(0,0,ngllz-1,0:2)
            Vert%sSimu(i_simu)%Forces3(0:2) = Vert%sSimu(i_simu)%Forces3(0:2) + &
                                              Elem%sSimu(i_simu)%Forces3(0,0,ngllz-1,0:2)
           endif
        case (5)
           Vert%sSimu(i_simu)%Forces(0:2) = Vert%sSimu(i_simu)%Forces(0:2) + &
                                            Elem%sSimu(i_simu)%Forces(ngllx-1,0,ngllz-1,0:2)
           if (Vert%PML .and. Elem%PML) then
            Vert%sSimu(i_simu)%Forces1(0:2) = Vert%sSimu(i_simu)%Forces1(0:2) + &
                                              Elem%sSimu(i_simu)%Forces1(ngllx-1,0,ngllz-1,0:2)
            Vert%sSimu(i_simu)%Forces2(0:2) = Vert%sSimu(i_simu)%Forces2(0:2) + &
                                              Elem%sSimu(i_simu)%Forces2(ngllx-1,0,ngllz-1,0:2)
            Vert%sSimu(i_simu)%Forces3(0:2) = Vert%sSimu(i_simu)%Forces3(0:2) + &
                                              Elem%sSimu(i_simu)%Forces3(ngllx-1,0,ngllz-1,0:2)
           endif
        case (6)
           Vert%sSimu(i_simu)%Forces(0:2) = Vert%sSimu(i_simu)%Forces(0:2) + &
                                            Elem%sSimu(i_simu)%Forces(ngllx-1,nglly-1,ngllz-1,0:2)
           if (Vert%PML .and. Elem%PML) then
            Vert%sSimu(i_simu)%Forces1(0:2) = Vert%sSimu(i_simu)%Forces1(0:2) + &
                                              Elem%sSimu(i_simu)%Forces1(ngllx-1,nglly-1,ngllz-1,0:2)
            Vert%sSimu(i_simu)%Forces2(0:2) = Vert%sSimu(i_simu)%Forces2(0:2) + &
                                              Elem%sSimu(i_simu)%Forces2(ngllx-1,nglly-1,ngllz-1,0:2)
            Vert%sSimu(i_simu)%Forces3(0:2) = Vert%sSimu(i_simu)%Forces3(0:2) + &
                                              Elem%sSimu(i_simu)%Forces3(ngllx-1,nglly-1,ngllz-1,0:2)
           endif
        case (7)
           Vert%sSimu(i_simu)%Forces(0:2) = Vert%sSimu(i_simu)%Forces(0:2) + &
                                            Elem%sSimu(i_simu)%Forces(0,nglly-1,ngllz-1,0:2)
           if (Vert%PML .and. Elem%PML) then
            Vert%sSimu(i_simu)%Forces1(0:2) = Vert%sSimu(i_simu)%Forces1(0:2) + &
                                              Elem%sSimu(i_simu)%Forces1(0,nglly-1,ngllz-1,0:2)
            Vert%sSimu(i_simu)%Forces2(0:2) = Vert%sSimu(i_simu)%Forces2(0:2) + &
                                              Elem%sSimu(i_simu)%Forces2(0,nglly-1,ngllz-1,0:2)
            Vert%sSimu(i_simu)%Forces3(0:2) = Vert%sSimu(i_simu)%Forces3(0:2) + &
                                              Elem%sSimu(i_simu)%Forces3(0,nglly-1,ngllz-1,0:2)
           endif
       end select

       return


end subroutine get_Forces_Elem2Vertex
