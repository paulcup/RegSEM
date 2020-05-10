subroutine get_Forces_Elem2Face (Elem, Fac, i, i_simu)


use selements
use sfaces

implicit none

type (Element), intent (INOUT) :: Elem
type (Face), intent (INOUT) :: Fac
integer, intent (IN) :: i, i_simu

integer :: ngllx,nglly,ngllz, ngll1,ngll2


       ngllx = Elem%ngllx
       nglly = Elem%nglly
       ngllz = Elem%ngllz

       ngll1 = Fac%ngll1
       ngll2 = Fac%ngll2

       select case (i)
        case (0)
           Fac%sSimu(i_simu)%Forces(1:ngll1-2,1:ngll2-2,0:2) = Fac%sSimu(i_simu)%Forces(1:ngll1-2,1:ngll2-2,0:2) + &
                                                               Elem%sSimu(i_simu)%Forces(1:ngll1-2,1:ngll2-2,0,0:2)
           if (Fac%PML) then
            Fac%sSimu(i_simu)%Forces1(1:ngll1-2,1:ngll2-2,0:2) = Fac%sSimu(i_simu)%Forces1(1:ngll1-2,1:ngll2-2,0:2) + &
                                                                 Elem%sSimu(i_simu)%Forces1(1:ngll1-2,1:ngll2-2,0,0:2)
            Fac%sSimu(i_simu)%Forces2(1:ngll1-2,1:ngll2-2,0:2) = Fac%sSimu(i_simu)%Forces2(1:ngll1-2,1:ngll2-2,0:2) + &
                                                                 Elem%sSimu(i_simu)%Forces2(1:ngll1-2,1:ngll2-2,0,0:2)
            Fac%sSimu(i_simu)%Forces3(1:ngll1-2,1:ngll2-2,0:2) = Fac%sSimu(i_simu)%Forces3(1:ngll1-2,1:ngll2-2,0:2) + &
                                                                 Elem%sSimu(i_simu)%Forces3(1:ngll1-2,1:ngll2-2,0,0:2)
           endif
        case (1)
           Fac%sSimu(i_simu)%Forces(1:ngll1-2,1:ngll2-2,0:2) = Fac%sSimu(i_simu)%Forces(1:ngll1-2,1:ngll2-2,0:2) + &
                                                               Elem%sSimu(i_simu)%Forces(1:ngll1-2,0,1:ngll2-2,0:2)
           if (Fac%PML) then
            Fac%sSimu(i_simu)%Forces1(1:ngll1-2,1:ngll2-2,0:2) = Fac%sSimu(i_simu)%Forces1(1:ngll1-2,1:ngll2-2,0:2) + &
                                                                 Elem%sSimu(i_simu)%Forces1(1:ngll1-2,0,1:ngll2-2,0:2)
            Fac%sSimu(i_simu)%Forces2(1:ngll1-2,1:ngll2-2,0:2) = Fac%sSimu(i_simu)%Forces2(1:ngll1-2,1:ngll2-2,0:2) + &
                                                                 Elem%sSimu(i_simu)%Forces2(1:ngll1-2,0,1:ngll2-2,0:2)
            Fac%sSimu(i_simu)%Forces3(1:ngll1-2,1:ngll2-2,0:2) = Fac%sSimu(i_simu)%Forces3(1:ngll1-2,1:ngll2-2,0:2) + &
                                                                 Elem%sSimu(i_simu)%Forces3(1:ngll1-2,0,1:ngll2-2,0:2)
           endif
        case (2)
           Fac%sSimu(i_simu)%Forces(1:ngll1-2,1:ngll2-2,0:2) = Fac%sSimu(i_simu)%Forces(1:ngll1-2,1:ngll2-2,0:2) + &
                                                               Elem%sSimu(i_simu)%Forces(ngllx-1,1:ngll1-2,1:ngll2-2,0:2)
           if (Fac%PML) then
            Fac%sSimu(i_simu)%Forces1(1:ngll1-2,1:ngll2-2,0:2) = Fac%sSimu(i_simu)%Forces1(1:ngll1-2,1:ngll2-2,0:2) + &
                                                                 Elem%sSimu(i_simu)%Forces1(ngllx-1,1:ngll1-2,1:ngll2-2,0:2)
            Fac%sSimu(i_simu)%Forces2(1:ngll1-2,1:ngll2-2,0:2) = Fac%sSimu(i_simu)%Forces2(1:ngll1-2,1:ngll2-2,0:2) + &
                                                                 Elem%sSimu(i_simu)%Forces2(ngllx-1,1:ngll1-2,1:ngll2-2,0:2)
            Fac%sSimu(i_simu)%Forces3(1:ngll1-2,1:ngll2-2,0:2) = Fac%sSimu(i_simu)%Forces3(1:ngll1-2,1:ngll2-2,0:2) + &
                                                                 Elem%sSimu(i_simu)%Forces3(ngllx-1,1:ngll1-2,1:ngll2-2,0:2)
           endif
        case (3)
          Fac%sSimu(i_simu)%Forces(1:ngll1-2,1:ngll2-2,0:2) = Fac%sSimu(i_simu)%Forces(1:ngll1-2,1:ngll2-2,0:2) + &
                                                              Elem%sSimu(i_simu)%Forces(1:ngll1-2,nglly-1,1:ngll2-2,0:2)
           if (Fac%PML) then
            Fac%sSimu(i_simu)%Forces1(1:ngll1-2,1:ngll2-2,0:2) = Fac%sSimu(i_simu)%Forces1(1:ngll1-2,1:ngll2-2,0:2) + &
                                                                 Elem%sSimu(i_simu)%Forces1(1:ngll1-2,nglly-1,1:ngll2-2,0:2)
            Fac%sSimu(i_simu)%Forces2(1:ngll1-2,1:ngll2-2,0:2) = Fac%sSimu(i_simu)%Forces2(1:ngll1-2,1:ngll2-2,0:2) + &
                                                                 Elem%sSimu(i_simu)%Forces2(1:ngll1-2,nglly-1,1:ngll2-2,0:2)
            Fac%sSimu(i_simu)%Forces3(1:ngll1-2,1:ngll2-2,0:2) = Fac%sSimu(i_simu)%Forces3(1:ngll1-2,1:ngll2-2,0:2) + &
                                                                 Elem%sSimu(i_simu)%Forces3(1:ngll1-2,nglly-1,1:ngll2-2,0:2)
           endif
        case (4)
           Fac%sSimu(i_simu)%Forces(1:ngll1-2,1:ngll2-2,0:2) = Fac%sSimu(i_simu)%Forces(1:ngll1-2,1:ngll2-2,0:2) + &
                                                               Elem%sSimu(i_simu)%Forces(0,1:ngll1-2,1:ngll2-2,0:2)
           if (Fac%PML) then
            Fac%sSimu(i_simu)%Forces1(1:ngll1-2,1:ngll2-2,0:2) = Fac%sSimu(i_simu)%Forces1(1:ngll1-2,1:ngll2-2,0:2) + &
                                                                 Elem%sSimu(i_simu)%Forces1(0,1:ngll1-2,1:ngll2-2,0:2)
            Fac%sSimu(i_simu)%Forces2(1:ngll1-2,1:ngll2-2,0:2) = Fac%sSimu(i_simu)%Forces2(1:ngll1-2,1:ngll2-2,0:2) + &
                                                                 Elem%sSimu(i_simu)%Forces2(0,1:ngll1-2,1:ngll2-2,0:2)
            Fac%sSimu(i_simu)%Forces3(1:ngll1-2,1:ngll2-2,0:2) = Fac%sSimu(i_simu)%Forces3(1:ngll1-2,1:ngll2-2,0:2) + &
                                                                 Elem%sSimu(i_simu)%Forces3(0,1:ngll1-2,1:ngll2-2,0:2)
           endif
        case (5)
           Fac%sSimu(i_simu)%Forces(1:ngll1-2,1:ngll2-2,0:2) = Fac%sSimu(i_simu)%Forces(1:ngll1-2,1:ngll2-2,0:2) + &
                                                               Elem%sSimu(i_simu)%Forces(1:ngll1-2,1:ngll2-2,ngllz-1,0:2)
           if (Fac%PML) then
            Fac%sSimu(i_simu)%Forces1(1:ngll1-2,1:ngll2-2,0:2) = Fac%sSimu(i_simu)%Forces1(1:ngll1-2,1:ngll2-2,0:2) + &
                                                                 Elem%sSimu(i_simu)%Forces1(1:ngll1-2,1:ngll2-2,ngllz-1,0:2)
            Fac%sSimu(i_simu)%Forces2(1:ngll1-2,1:ngll2-2,0:2) = Fac%sSimu(i_simu)%Forces2(1:ngll1-2,1:ngll2-2,0:2) + &
                                                                 Elem%sSimu(i_simu)%Forces2(1:ngll1-2,1:ngll2-2,ngllz-1,0:2)
            Fac%sSimu(i_simu)%Forces3(1:ngll1-2,1:ngll2-2,0:2) = Fac%sSimu(i_simu)%Forces3(1:ngll1-2,1:ngll2-2,0:2) + &
                                                                 Elem%sSimu(i_simu)%Forces3(1:ngll1-2,1:ngll2-2,ngllz-1,0:2)
           endif
       end select

       return


end subroutine get_Forces_Elem2Face
