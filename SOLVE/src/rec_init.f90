!==============================================
subroutine rec_init (Tdomain, ntime, ninit, rg)
!==============================================
!
! Written by Yder MASSON (masson@ipgp.fr) june 20 2012  
! 
! record fields variables for point matching when running the time reversal mirror simulation with attenation 
!
use sdomains
!
implicit none
!
type (Domain), intent (INOUT) :: Tdomain
integer, intent (IN) :: ntime, ninit, rg
!
integer :: len_DP, n
character*60 :: init_file
!
if (ntime==0) then
!
   inquire(iolength=len_DP)&
        &(Tdomain%specel(n)%sSimu(0)%Displ,n=0,Tdomain%n_elem-1),         &
        &(Tdomain%sFace(n)%sSimu(0)%Displ,n=0,Tdomain%n_face-1),          &
        &(Tdomain%sEdge(n)%sSimu(0)%Displ,n=0,Tdomain%n_edge-1),          &
        &(Tdomain%sVertex(n)%sSimu(0)%Displ,n=0,Tdomain%n_vertex-1),      &
!        &                                                                 &
        &(Tdomain%specel(n)%sSimu(0)%Veloc,n=0,Tdomain%n_elem-1),         &
        &(Tdomain%sFace(n)%sSimu(0)%Veloc,n=0,Tdomain%n_face-1),          &
        &(Tdomain%sEdge(n)%sSimu(0)%Veloc,n=0,Tdomain%n_edge-1),          &
        &(Tdomain%sVertex(n)%sSimu(0)%Veloc,n=0,Tdomain%n_vertex-1),      &
!        &                                                                 &
        &(Tdomain%specel(n)%sSimu(0)%R_xx_,n=0,Tdomain%n_elem-1),         &
        &(Tdomain%specel(n)%sSimu(0)%R_yy_,n=0,Tdomain%n_elem-1),         &
        &(Tdomain%specel(n)%sSimu(0)%R_xy_,n=0,Tdomain%n_elem-1),         &
        &(Tdomain%specel(n)%sSimu(0)%R_xz_,n=0,Tdomain%n_elem-1),         &
        &(Tdomain%specel(n)%sSimu(0)%R_yz_,n=0,Tdomain%n_elem-1),         &
!        &                                                                 &
        &(Tdomain%specel(n)%sSimu(0)%epsilondev_xx_,n=0,Tdomain%n_elem-1),&
        &(Tdomain%specel(n)%sSimu(0)%epsilondev_yy_,n=0,Tdomain%n_elem-1),&
        &(Tdomain%specel(n)%sSimu(0)%epsilondev_xy_,n=0,Tdomain%n_elem-1),&
        &(Tdomain%specel(n)%sSimu(0)%epsilondev_xz_,n=0,Tdomain%n_elem-1),&
        &(Tdomain%specel(n)%sSimu(0)%epsilondev_yz_,n=0,Tdomain%n_elem-1)
!   
   if (len_DP/=0) then
      write (init_file,"(a,I3.3)") "init",rg
      open (45,file=init_file,access='direct',form="unformatted",recl=len_DP,status='replace')
   endif
!
endif
!
if(mod(Tdomain%sTimeParam%ntime-ntime-1,ninit)==0)then
!
   write (45,rec=(Tdomain%sTimeParam%ntime-ntime-1)/ninit+1)              &
        &(Tdomain%specel(n)%sSimu(0)%Displ,n=0,Tdomain%n_elem-1),         &
        &(Tdomain%sFace(n)%sSimu(0)%Displ,n=0,Tdomain%n_face-1),          &
        &(Tdomain%sEdge(n)%sSimu(0)%Displ,n=0,Tdomain%n_edge-1),          &
        &(Tdomain%sVertex(n)%sSimu(0)%Displ,n=0,Tdomain%n_vertex-1),      &
!        &                                                                 &
        &(-Tdomain%specel(n)%sSimu(0)%Veloc,n=0,Tdomain%n_elem-1),        &
        &(-Tdomain%sFace(n)%sSimu(0)%Veloc,n=0,Tdomain%n_face-1),         &
        &(-Tdomain%sEdge(n)%sSimu(0)%Veloc,n=0,Tdomain%n_edge-1),         &
        &(-Tdomain%sVertex(n)%sSimu(0)%Veloc,n=0,Tdomain%n_vertex-1),     &
!        &                                                                 &
        &(Tdomain%specel(n)%sSimu(0)%R_xx_,n=0,Tdomain%n_elem-1),         &
        &(Tdomain%specel(n)%sSimu(0)%R_yy_,n=0,Tdomain%n_elem-1),         &
        &(Tdomain%specel(n)%sSimu(0)%R_xy_,n=0,Tdomain%n_elem-1),         &
        &(Tdomain%specel(n)%sSimu(0)%R_xz_,n=0,Tdomain%n_elem-1),         &
        &(Tdomain%specel(n)%sSimu(0)%R_yz_,n=0,Tdomain%n_elem-1),         & 
!        &                                                                 &
        &(Tdomain%specel(n)%sSimu(0)%epsilondev_xx_,n=0,Tdomain%n_elem-1),&
        &(Tdomain%specel(n)%sSimu(0)%epsilondev_yy_,n=0,Tdomain%n_elem-1),&
        &(Tdomain%specel(n)%sSimu(0)%epsilondev_xy_,n=0,Tdomain%n_elem-1),&
        &(Tdomain%specel(n)%sSimu(0)%epsilondev_xz_,n=0,Tdomain%n_elem-1),&
        &(Tdomain%specel(n)%sSimu(0)%epsilondev_yz_,n=0,Tdomain%n_elem-1)
!
endif
!
if (ntime==Tdomain%sTimeParam%ntime-1)   close (45)
!
!======================
end subroutine rec_init
!======================
