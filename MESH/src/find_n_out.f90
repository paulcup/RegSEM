subroutine find_numloc (nb_elem,elem,ii,n_out)


implicit none

integer, intent(IN) :: nb_elem, ii
integer, intent(INOUT) :: n_out
integer, dimension(0:nb_elem-1), intent(IN) :: elem
integer :: n


boucle : do n = 0,nb_elem-1
    if (elem(n)==ii) then
        n_out = n
        exit boucle
    endif
enddo boucle


end subroutine find_numloc
