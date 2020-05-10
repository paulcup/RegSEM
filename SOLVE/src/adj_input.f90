subroutine adj_input (Tdomain)


use sdomains

implicit none

type (domain), intent (INOUT) :: Tdomain

integer :: n_src, n, i
doubleprecision :: t


open (60,file="src_adj",form="formatted",status="old")
read (60,*) Tdomain%n_source_adj
allocate (Tdomain%sAdjoint(0:Tdomain%n_source_adj-1))
do n_src = 0, Tdomain%n_source_adj-1
    read (60,*) Tdomain%sAdjoint(n_src)%realcolat, Tdomain%sAdjoint(n_src)%reallong, Tdomain%sAdjoint(n_src)%depth
    read (60,*) n
    if (n/=Tdomain%sTimeParam%ntime) then
        print *, "THE NUMBER OF SAMPLES IN THE ADJOINT SOURCE IS DIFFERENT FROM THE NUMBER OF ITERATION !!!"
        stop
    endif
    allocate (Tdomain%sAdjoint(n_src)%timefunc(0:2,0:n-1))
    do i = 0,n-1
        read (60,*) Tdomain%sAdjoint(n_src)%timefunc(0,i), &
                    Tdomain%sAdjoint(n_src)%timefunc(1,i), &
                    Tdomain%sAdjoint(n_src)%timefunc(2,i)
    enddo
enddo
close (60)


return
end subroutine adj_input

