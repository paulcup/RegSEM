subroutine get_Forces_Elem2Edge (Elem, Edg, i, i_simu)


use selements
use sedges

implicit none

type (Element), intent (INOUT) :: Elem
type (Edge), intent (INOUT) :: Edg
integer, intent (IN) :: i, i_simu

integer :: ngllx,nglly,ngllz, ngll


       ngllx = Elem%ngllx
       nglly = Elem%nglly
       ngllz = Elem%ngllz

       ngll = Edg%ngll

       select case (i)
        case (0)
           Edg%sSimu(i_simu)%Forces(1:ngll-2,0:2) = Edg%sSimu(i_simu)%Forces(1:ngll-2,0:2) + &
                                                    Elem%sSimu(i_simu)%Forces(1:ngll-2,0,0,0:2)
           if (Edg%PML .and. Elem%PML) then
            Edg%sSimu(i_simu)%Forces1(1:ngll-2,0:2) = Edg%sSimu(i_simu)%Forces1(1:ngll-2,0:2) + &
                                                      Elem%sSimu(i_simu)%Forces1(1:ngll-2,0,0,0:2)
            Edg%sSimu(i_simu)%Forces2(1:ngll-2,0:2) = Edg%sSimu(i_simu)%Forces2(1:ngll-2,0:2) + &
                                                      Elem%sSimu(i_simu)%Forces2(1:ngll-2,0,0,0:2)
            Edg%sSimu(i_simu)%Forces3(1:ngll-2,0:2) = Edg%sSimu(i_simu)%Forces3(1:ngll-2,0:2) + &
                                                      Elem%sSimu(i_simu)%Forces3(1:ngll-2,0,0,0:2)
           endif
        case (1)
           Edg%sSimu(i_simu)%Forces(1:ngll-2,0:2) = Edg%sSimu(i_simu)%Forces(1:ngll-2,0:2) + &
                                                    Elem%sSimu(i_simu)%Forces(ngllx-1,1:ngll-2,0,0:2)
           if (Edg%PML .and. Elem%PML) then
            Edg%sSimu(i_simu)%Forces1(1:ngll-2,0:2) = Edg%sSimu(i_simu)%Forces1(1:ngll-2,0:2) + &
                                                      Elem%sSimu(i_simu)%Forces1(ngllx-1,1:ngll-2,0,0:2)
            Edg%sSimu(i_simu)%Forces2(1:ngll-2,0:2) = Edg%sSimu(i_simu)%Forces2(1:ngll-2,0:2) + &
                                                      Elem%sSimu(i_simu)%Forces2(ngllx-1,1:ngll-2,0,0:2)
            Edg%sSimu(i_simu)%Forces3(1:ngll-2,0:2) = Edg%sSimu(i_simu)%Forces3(1:ngll-2,0:2) + &
                                                      Elem%sSimu(i_simu)%Forces3(ngllx-1,1:ngll-2,0,0:2)
           endif
        case (2)
           Edg%sSimu(i_simu)%Forces(1:ngll-2,0:2) = Edg%sSimu(i_simu)%Forces(1:ngll-2,0:2) + &
                                                    Elem%sSimu(i_simu)%Forces(1:ngll-2,nglly-1,0,0:2)
           if (Edg%PML .and. Elem%PML) then
            Edg%sSimu(i_simu)%Forces1(1:ngll-2,0:2) = Edg%sSimu(i_simu)%Forces1(1:ngll-2,0:2) + &
                                                      Elem%sSimu(i_simu)%Forces1(1:ngll-2,nglly-1,0,0:2)
            Edg%sSimu(i_simu)%Forces2(1:ngll-2,0:2) = Edg%sSimu(i_simu)%Forces2(1:ngll-2,0:2) + &
                                                      Elem%sSimu(i_simu)%Forces2(1:ngll-2,nglly-1,0,0:2)
            Edg%sSimu(i_simu)%Forces3(1:ngll-2,0:2) = Edg%sSimu(i_simu)%Forces3(1:ngll-2,0:2) + &
                                                      Elem%sSimu(i_simu)%Forces3(1:ngll-2,nglly-1,0,0:2)
           endif
        case (3)
           Edg%sSimu(i_simu)%Forces(1:ngll-2,0:2) = Edg%sSimu(i_simu)%Forces(1:ngll-2,0:2) + &
                                                    Elem%sSimu(i_simu)%Forces(0,1:ngll-2,0,0:2)
           if (Edg%PML .and. Elem%PML) then
            Edg%sSimu(i_simu)%Forces1(1:ngll-2,0:2) = Edg%sSimu(i_simu)%Forces1(1:ngll-2,0:2) + &
                                                      Elem%sSimu(i_simu)%Forces1(0,1:ngll-2,0,0:2)
            Edg%sSimu(i_simu)%Forces2(1:ngll-2,0:2) = Edg%sSimu(i_simu)%Forces2(1:ngll-2,0:2) + &
                                                      Elem%sSimu(i_simu)%Forces2(0,1:ngll-2,0,0:2)
            Edg%sSimu(i_simu)%Forces3(1:ngll-2,0:2) = Edg%sSimu(i_simu)%Forces3(1:ngll-2,0:2) + &
                                                      Elem%sSimu(i_simu)%Forces3(0,1:ngll-2,0,0:2)
           endif
        case (4)
           Edg%sSimu(i_simu)%Forces(1:ngll-2,0:2) = Edg%sSimu(i_simu)%Forces(1:ngll-2,0:2) + &
                                                    Elem%sSimu(i_simu)%Forces(ngllx-1,0,1:ngll-2,0:2)
           if (Edg%PML .and. Elem%PML) then
            Edg%sSimu(i_simu)%Forces1(1:ngll-2,0:2) = Edg%sSimu(i_simu)%Forces1(1:ngll-2,0:2) + &
                                                      Elem%sSimu(i_simu)%Forces1(ngllx-1,0,1:ngll-2,0:2)
            Edg%sSimu(i_simu)%Forces2(1:ngll-2,0:2) = Edg%sSimu(i_simu)%Forces2(1:ngll-2,0:2) + &
                                                      Elem%sSimu(i_simu)%Forces2(ngllx-1,0,1:ngll-2,0:2)
            Edg%sSimu(i_simu)%Forces3(1:ngll-2,0:2) = Edg%sSimu(i_simu)%Forces3(1:ngll-2,0:2) + &
                                                      Elem%sSimu(i_simu)%Forces3(ngllx-1,0,1:ngll-2,0:2)
           endif
        case (5)
           Edg%sSimu(i_simu)%Forces(1:ngll-2,0:2) = Edg%sSimu(i_simu)%Forces(1:ngll-2,0:2) + &
                                                    Elem%sSimu(i_simu)%Forces(1:ngll-2,0,ngllz-1,0:2)
           if (Edg%PML .and. Elem%PML) then
            Edg%sSimu(i_simu)%Forces1(1:ngll-2,0:2) = Edg%sSimu(i_simu)%Forces1(1:ngll-2,0:2) + &
                                                      Elem%sSimu(i_simu)%Forces1(1:ngll-2,0,ngllz-1,0:2)
            Edg%sSimu(i_simu)%Forces2(1:ngll-2,0:2) = Edg%sSimu(i_simu)%Forces2(1:ngll-2,0:2) + &
                                                      Elem%sSimu(i_simu)%Forces2(1:ngll-2,0,ngllz-1,0:2)
            Edg%sSimu(i_simu)%Forces3(1:ngll-2,0:2) = Edg%sSimu(i_simu)%Forces3(1:ngll-2,0:2) + &
                                                      Elem%sSimu(i_simu)%Forces3(1:ngll-2,0,ngllz-1,0:2)
           endif
        case (6)
           Edg%sSimu(i_simu)%Forces(1:ngll-2,0:2) = Edg%sSimu(i_simu)%Forces(1:ngll-2,0:2) + &
                                                    Elem%sSimu(i_simu)%Forces(0,0,1:ngll-2,0:2)
           if (Edg%PML .and. Elem%PML) then
            Edg%sSimu(i_simu)%Forces1(1:ngll-2,0:2) = Edg%sSimu(i_simu)%Forces1(1:ngll-2,0:2) + &
                                                      Elem%sSimu(i_simu)%Forces1(0,0,1:ngll-2,0:2)
            Edg%sSimu(i_simu)%Forces2(1:ngll-2,0:2) = Edg%sSimu(i_simu)%Forces2(1:ngll-2,0:2) + &
                                                      Elem%sSimu(i_simu)%Forces2(0,0,1:ngll-2,0:2)
            Edg%sSimu(i_simu)%Forces3(1:ngll-2,0:2) = Edg%sSimu(i_simu)%Forces3(1:ngll-2,0:2) + &
                                                      Elem%sSimu(i_simu)%Forces3(0,0,1:ngll-2,0:2)
           endif
        case (7)
           Edg%sSimu(i_simu)%Forces(1:ngll-2,0:2) = Edg%sSimu(i_simu)%Forces(1:ngll-2,0:2) + &
                                                    Elem%sSimu(i_simu)%Forces(ngllx-1,nglly-1,1:ngll-2,0:2)
           if (Edg%PML .and. Elem%PML) then
            Edg%sSimu(i_simu)%Forces1(1:ngll-2,0:2) = Edg%sSimu(i_simu)%Forces1(1:ngll-2,0:2) + &
                                                      Elem%sSimu(i_simu)%Forces1(ngllx-1,nglly-1,1:ngll-2,0:2)
            Edg%sSimu(i_simu)%Forces2(1:ngll-2,0:2) = Edg%sSimu(i_simu)%Forces2(1:ngll-2,0:2) + &
                                                      Elem%sSimu(i_simu)%Forces2(ngllx-1,nglly-1,1:ngll-2,0:2)
            Edg%sSimu(i_simu)%Forces3(1:ngll-2,0:2) = Edg%sSimu(i_simu)%Forces3(1:ngll-2,0:2) + &
                                                      Elem%sSimu(i_simu)%Forces3(ngllx-1,nglly-1,1:ngll-2,0:2)
           endif
        case (8)
           Edg%sSimu(i_simu)%Forces(1:ngll-2,0:2) = Edg%sSimu(i_simu)%Forces(1:ngll-2,0:2) + &
                                                    Elem%sSimu(i_simu)%Forces(ngllx-1,1:ngll-2,ngllz-1,0:2)
           if (Edg%PML .and. Elem%PML) then
            Edg%sSimu(i_simu)%Forces1(1:ngll-2,0:2) = Edg%sSimu(i_simu)%Forces1(1:ngll-2,0:2) + &
                                                      Elem%sSimu(i_simu)%Forces1(ngllx-1,1:ngll-2,ngllz-1,0:2)
            Edg%sSimu(i_simu)%Forces2(1:ngll-2,0:2) = Edg%sSimu(i_simu)%Forces2(1:ngll-2,0:2) + &
                                                      Elem%sSimu(i_simu)%Forces2(ngllx-1,1:ngll-2,ngllz-1,0:2)
            Edg%sSimu(i_simu)%Forces3(1:ngll-2,0:2) = Edg%sSimu(i_simu)%Forces3(1:ngll-2,0:2) + &
                                                      Elem%sSimu(i_simu)%Forces3(ngllx-1,1:ngll-2,ngllz-1,0:2)
           endif
        case (9)
           Edg%sSimu(i_simu)%Forces(1:ngll-2,0:2) = Edg%sSimu(i_simu)%Forces(1:ngll-2,0:2) + &
                                                    Elem%sSimu(i_simu)%Forces(1:ngll-2,nglly-1,ngllz-1,0:2)
           if (Edg%PML .and. Elem%PML) then
            Edg%sSimu(i_simu)%Forces1(1:ngll-2,0:2) = Edg%sSimu(i_simu)%Forces1(1:ngll-2,0:2) + &
                                                      Elem%sSimu(i_simu)%Forces1(1:ngll-2,nglly-1,ngllz-1,0:2)
            Edg%sSimu(i_simu)%Forces2(1:ngll-2,0:2) = Edg%sSimu(i_simu)%Forces2(1:ngll-2,0:2) + &
                                                      Elem%sSimu(i_simu)%Forces2(1:ngll-2,nglly-1,ngllz-1,0:2)
            Edg%sSimu(i_simu)%Forces3(1:ngll-2,0:2) = Edg%sSimu(i_simu)%Forces3(1:ngll-2,0:2) + &
                                                      Elem%sSimu(i_simu)%Forces3(1:ngll-2,nglly-1,ngllz-1,0:2)
           endif
        case (10)
           Edg%sSimu(i_simu)%Forces(1:ngll-2,0:2) = Edg%sSimu(i_simu)%Forces(1:ngll-2,0:2) + &
                                                    Elem%sSimu(i_simu)%Forces(0,nglly-1,1:ngll-2,0:2)
           if (Edg%PML .and. Elem%PML) then
            Edg%sSimu(i_simu)%Forces1(1:ngll-2,0:2) = Edg%sSimu(i_simu)%Forces1(1:ngll-2,0:2) + &
                                                      Elem%sSimu(i_simu)%Forces1(0,nglly-1,1:ngll-2,0:2)
            Edg%sSimu(i_simu)%Forces2(1:ngll-2,0:2) = Edg%sSimu(i_simu)%Forces2(1:ngll-2,0:2) + &
                                                      Elem%sSimu(i_simu)%Forces2(0,nglly-1,1:ngll-2,0:2)
            Edg%sSimu(i_simu)%Forces3(1:ngll-2,0:2) = Edg%sSimu(i_simu)%Forces3(1:ngll-2,0:2) + &
                                                      Elem%sSimu(i_simu)%Forces3(0,nglly-1,1:ngll-2,0:2)
           endif
        case (11)
           Edg%sSimu(i_simu)%Forces(1:ngll-2,0:2) = Edg%sSimu(i_simu)%Forces(1:ngll-2,0:2) + &
                                                    Elem%sSimu(i_simu)%Forces(0,1:ngll-2,ngllz-1,0:2)
           if (Edg%PML .and. Elem%PML) then
            Edg%sSimu(i_simu)%Forces1(1:ngll-2,0:2) = Edg%sSimu(i_simu)%Forces1(1:ngll-2,0:2) + &
                                                      Elem%sSimu(i_simu)%Forces1(0,1:ngll-2,ngllz-1,0:2)
            Edg%sSimu(i_simu)%Forces2(1:ngll-2,0:2) = Edg%sSimu(i_simu)%Forces2(1:ngll-2,0:2) + &
                                                      Elem%sSimu(i_simu)%Forces2(0,1:ngll-2,ngllz-1,0:2)
            Edg%sSimu(i_simu)%Forces3(1:ngll-2,0:2) = Edg%sSimu(i_simu)%Forces3(1:ngll-2,0:2) + &
                                                      Elem%sSimu(i_simu)%Forces3(0,1:ngll-2,ngllz-1,0:2)
           endif
       end select

       return


end subroutine get_Forces_Elem2Edge
