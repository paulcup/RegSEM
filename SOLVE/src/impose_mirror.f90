subroutine impose_mirror (opt, Tdomain, ntime, rg)


use sdomains

implicit none

type (Domain), intent (INOUT) :: Tdomain
integer, intent (IN) :: ntime, rg

integer :: len_DP, n, i, x,y,z, ngllx,nglly,ngllz, comp
doubleprecision :: current_t, t_max, wtaper
doubleprecision, dimension(0:Tdomain%recl_mirror-1) :: tmp
character*60 :: mirror_file
character :: opt*5
doubleprecision, parameter :: taperwidth = 200.d0   ! Mettre 0.d0 pour ne pas avoir de taper


if (ntime==0) then
    inquire(iolength=len_DP) tmp
    if(opt=='displ')then
    if (len_DP/=0) then
     write (mirror_file,"(a,I3.3)") "mirror.displ.",rg
     open (44,file=mirror_file,access='direct',form="unformatted",recl=len_DP,status='old')
    endif
    elseif(opt=='force')then
    if (len_DP/=0) then
     write (mirror_file,"(a,I3.3)") "mirror.force.",rg
     open (43,file=mirror_file,access='direct',form="unformatted",recl=len_DP,status='old')
    endif
    else
    write(*,*)'WRONG ARGUMENT IN rec_mirror.f90 op must be "displ" or "force"'
    endif
endif

if (Tdomain%recl_mirror/=0.and.opt=='displ')   read(44,rec=Tdomain%sTimeParam%ntime-ntime) tmp
if (Tdomain%recl_mirror/=0.and.opt=='force')   read(43,rec=Tdomain%sTimeParam%ntime-ntime) tmp

! Taper applique au debut du mirroir pour eviter une arrivee trop brusque de l'energie.
! Ca peut arriver si le mirroir n'est pas bien a zero.
!current_t = Tdomain%sTimeParam%rtime
!t_max = Tdomain%sTimeParam%duration + Tdomain%sTimeParam%dt   ! On ajoute un dt pour etre bien sur que t est trjs < t_max
!call wtcoef(current_t,0.d0,taperwidth,t_max,t_max,wtaper)
!tmp = wtaper*tmp

!i = 0
!do n = 0, Tdomain%n_face-1
!    if (Tdomain%sFace(n)%mirror) then
!        do x = 1, Tdomain%sFace(n)%ngll1-2
!         do y = 1, Tdomain%sFace(n)%ngll2-2
!          do comp = 0,2
!              Tdomain%sFace(n)%sSimu(0)%Forces(x,y,comp) = tmp(i)
!              i = i + 1
!          enddo
!         enddo
!        enddo
!    endif
!enddo
!do n = 0,Tdomain%n_edge-1
!    if (Tdomain%sEdge(n)%mirror) then
!        do x = 1, Tdomain%sEdge(n)%ngll-2
!         do comp = 0,2
!             Tdomain%sEdge(n)%sSimu(0)%Forces(x,comp) = tmp(i)
!             i = i + 1
!         enddo
!        enddo
!    endif
!enddo
!do n = 0,Tdomain%n_vertex-1
!    if (Tdomain%sVertex(n)%mirror) then
!        do comp = 0,2
!            Tdomain%sVertex(n)%sSimu(0)%Forces(comp) = tmp(i)
!            i = i + 1
!        enddo
!    endif
!enddo

i = 0
do n = 0,Tdomain%n_elem-1
   if(Tdomain%specel(n)%mirror_position == 1)then
      ngllx = Tdomain%specel(n)%ngllx
      nglly = Tdomain%specel(n)%nglly
      ngllz = Tdomain%specel(n)%ngllz
      do z = 0,ngllz-1
         do y = 0,nglly-1
            do x = 0,ngllx-1
               if(Tdomain%specel(n)%win_mirror(x,y,z)==1)then
		 do comp = 0,2
		    if (opt=='displ')  Tdomain%specel(n)%sSimu(0)%Forces(x,y,z,comp) = Tdomain%specel(n)%sSimu(0)%Forces(x,y,z,comp)-tmp(i)
		    if (opt=='force')  Tdomain%specel(n)%sSimu(0)%Forces(x,y,z,comp) = Tdomain%specel(n)%sSimu(0)%Forces(x,y,z,comp)+tmp(i)
		    i = i+1
	         enddo
	       endif
            enddo
         enddo
      enddo
   endif
enddo

if (i/=Tdomain%recl_mirror)   stop 'A TMP BUFFER HAS A LENGTH DIFFERENT FROM RECL_MIRROR !!!'

if (ntime==Tdomain%sTimeParam%ntime-1 .and. Tdomain%recl_mirror/=0 .and. opt=='displ')   close (44)
if (ntime==Tdomain%sTimeParam%ntime-1 .and. Tdomain%recl_mirror/=0 .and. opt=='force')   close (43)


return
end subroutine impose_mirror
