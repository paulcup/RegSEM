!=================================================
subroutine impose_init (Tdomain, ntime, ninit, rg)
!=================================================
!
! Written by Yder MASSON (masson@ipgp.fr) june 20 2012  
! 
! reset fields variables (point matching) when running the time reversal mirror simulation with attenation 
!
use sdomains
!
implicit none
!
type (Domain), intent (INOUT) :: Tdomain
integer, intent (IN) :: ntime, ninit, rg
!
integer :: len_DP, n, ngllx, nglly, ngllz, comp, nf, ne, nv, n_solid, i, ngll, ngll1, ngll2
logical :: pml
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
      open (45,file=init_file,access='direct',form="unformatted",recl=len_DP,status='old')
   endif
!
endif
!
if(mod(ntime,ninit)==0)then
!
! Zero fields
!
n_solid = Tdomain%n_sls
!
do n = 0,Tdomain%n_elem-1
    ngllx = Tdomain%specel(n)%ngllx
    nglly = Tdomain%specel(n)%nglly
    ngllz = Tdomain%specel(n)%ngllz
    Tdomain%specel(n)%ssimu(0)%Veloc = 0.d0
    Tdomain%specel(n)%ssimu(0)%Accel = 0.d0
    if (Tdomain%specel(n)%PML) then
            Tdomain%specel(n)%ssimu(0)%Diagonal_Stress = 0.d0
            Tdomain%specel(n)%ssimu(0)%Diagonal_Stress1 = 0.d0
            Tdomain%specel(n)%ssimu(0)%Diagonal_Stress2 = 0.d0
            Tdomain%specel(n)%ssimu(0)%Diagonal_Stress3 = 0.d0
            Tdomain%specel(n)%ssimu(0)%Residual_Stress = 0.d0
            Tdomain%specel(n)%ssimu(0)%Residual_Stress1 = 0.d0
            Tdomain%specel(n)%ssimu(0)%Residual_Stress2 = 0.d0
            Tdomain%specel(n)%ssimu(0)%Veloc1 = 0.d0
            Tdomain%specel(n)%ssimu(0)%Veloc2 = 0.d0
            Tdomain%specel(n)%ssimu(0)%Veloc3 = 0.d0
        if (Tdomain%curve) then
                Tdomain%specel(n)%ssimu(0)%Residual_Stress3 = 0.d0
        endif
    else
            Tdomain%specel(n)%ssimu(0)%Displ = 0.d0
            Tdomain%specel(n)%ssimu(0)%V0 = 0.d0
        if (n_solid>0) then
                Tdomain%specel(n)%ssimu(0)%epsilondev_xx_ = 0.d0
                Tdomain%specel(n)%ssimu(0)%epsilondev_yy_ = 0.d0
                Tdomain%specel(n)%ssimu(0)%epsilondev_xy_ = 0.d0
                Tdomain%specel(n)%ssimu(0)%epsilondev_xz_ = 0.d0
                Tdomain%specel(n)%ssimu(0)%epsilondev_yz_ = 0.d0
                Tdomain%specel(n)%ssimu(0)%R_xx_ = 0.d0
                Tdomain%specel(n)%ssimu(0)%R_yy_ = 0.d0
                Tdomain%specel(n)%ssimu(0)%R_xy_ = 0.d0
                Tdomain%specel(n)%ssimu(0)%R_xz_ = 0.d0
                Tdomain%specel(n)%ssimu(0)%R_yz_ = 0.d0
        endif
    endif
enddo

do n = 0,Tdomain%n_face-1
    ngll1 = Tdomain%sFace(n)%ngll1
    ngll2 = Tdomain%sFace(n)%ngll2
        Tdomain%sFace(n)%ssimu(0)%Veloc = 0.d0
        Tdomain%sFace(n)%ssimu(0)%Accel = 0.d0
        Tdomain%sFace(n)%ssimu(0)%Forces = 0.d0
    if (Tdomain%sFace(n)%PML) then 
            Tdomain%sFace(n)%ssimu(0)%Veloc1 = 0.d0
            Tdomain%sFace(n)%ssimu(0)%Veloc2 = 0.d0
            Tdomain%sFace(n)%ssimu(0)%Veloc3 = 0.d0
    else
            Tdomain%sFace(n)%ssimu(0)%Displ = 0.d0
            Tdomain%sFace(n)%ssimu(0)%V0 = 0.d0
    endif
enddo

do n = 0,Tdomain%n_edge-1
    ngll = Tdomain%sEdge(n)%ngll
        Tdomain%sEdge(n)%ssimu(0)%Veloc = 0.d0
        Tdomain%sEdge(n)%ssimu(0)%Accel = 0.d0
        Tdomain%sEdge(n)%ssimu(0)%Forces = 0.d0
    if (Tdomain%sEdge(n)%PML) then
            Tdomain%sEdge(n)%ssimu(0)%Veloc1 = 0.d0
            Tdomain%sEdge(n)%ssimu(0)%Veloc2 = 0.d0
            Tdomain%sEdge(n)%ssimu(0)%Veloc3 = 0.d0
    else
            Tdomain%sEdge(n)%ssimu(0)%Displ = 0.d0
            Tdomain%sEdge(n)%ssimu(0)%V0 = 0.d0
    endif
enddo

do n = 0,Tdomain%n_vertex-1
        Tdomain%sVertex(n)%ssimu(0)%Veloc = 0.d0
        Tdomain%sVertex(n)%ssimu(0)%Accel = 0.d0
        Tdomain%sVertex(n)%ssimu(0)%Forces = 0.d0
    if (Tdomain%sVertex(n)%PML) then
            Tdomain%sVertex(n)%ssimu(0)%Veloc1 = 0.d0
            Tdomain%sVertex(n)%ssimu(0)%Veloc2 = 0.d0
            Tdomain%sVertex(n)%ssimu(0)%Veloc3 = 0.d0
    else
            Tdomain%sVertex(n)%ssimu(0)%Displ = 0.d0
            Tdomain%sVertex(n)%ssimu(0)%V0 = 0.d0
    endif
enddo
!
   read (45,rec=ntime/ninit+1) &
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
! Apply mirror to initial conditions
!
do n = 0,Tdomain%n_elem-1

    ngllx = Tdomain%specel(n)%ngllx
    nglly = Tdomain%specel(n)%nglly
    ngllz = Tdomain%specel(n)%ngllz

    Tdomain%specel(n)%sSimu(0)%Displ(:,:,:,0) = Tdomain%specel(n)%sSimu(0)%Displ(:,:,:,0) * Tdomain%specel(n)%win_mirror(1:ngllx-2, 1:nglly-2, 1:ngllz-2)
    Tdomain%specel(n)%sSimu(0)%Displ(:,:,:,1) = Tdomain%specel(n)%sSimu(0)%Displ(:,:,:,1) * Tdomain%specel(n)%win_mirror(1:ngllx-2, 1:nglly-2, 1:ngllz-2)
    Tdomain%specel(n)%sSimu(0)%Displ(:,:,:,2) = Tdomain%specel(n)%sSimu(0)%Displ(:,:,:,2) * Tdomain%specel(n)%win_mirror(1:ngllx-2, 1:nglly-2, 1:ngllz-2)

    Tdomain%specel(n)%sSimu(0)%Veloc(:,:,:,0) = Tdomain%specel(n)%sSimu(0)%Veloc(:,:,:,0) * Tdomain%specel(n)%win_mirror(1:ngllx-2, 1:nglly-2, 1:ngllz-2)
    Tdomain%specel(n)%sSimu(0)%Veloc(:,:,:,1) = Tdomain%specel(n)%sSimu(0)%Veloc(:,:,:,1) * Tdomain%specel(n)%win_mirror(1:ngllx-2, 1:nglly-2, 1:ngllz-2)
    Tdomain%specel(n)%sSimu(0)%Veloc(:,:,:,2) = Tdomain%specel(n)%sSimu(0)%Veloc(:,:,:,2) * Tdomain%specel(n)%win_mirror(1:ngllx-2, 1:nglly-2, 1:ngllz-2)


    if(Tdomain%specel(n)%mirror_position<=1)then

       Tdomain%specel(n)%sSimu(0)%R_xx_ = 0
       Tdomain%specel(n)%sSimu(0)%R_yy_ = 0
       Tdomain%specel(n)%sSimu(0)%R_xy_ = 0
       Tdomain%specel(n)%sSimu(0)%R_xz_ = 0
       Tdomain%specel(n)%sSimu(0)%R_yz_ = 0

       Tdomain%specel(n)%sSimu(0)%epsilondev_xx_ = 0
       Tdomain%specel(n)%sSimu(0)%epsilondev_yy_ = 0
       Tdomain%specel(n)%sSimu(0)%epsilondev_xy_ = 0
       Tdomain%specel(n)%sSimu(0)%epsilondev_xz_ = 0
       Tdomain%specel(n)%sSimu(0)%epsilondev_yz_ = 0

    endif
!
do comp = 0,2
!
nf = Tdomain%specel(n)%near_faces(0)
pml = .not.Tdomain%sFace(nf)%PML
if(pml)Tdomain%sFace(nf)%sSimu(0)%Displ(1:ngllx-2,1:nglly-2,comp) = Tdomain%sFace(nf)%sSimu(0)%Displ(1:ngllx-2,1:nglly-2,comp)&
                                                                  *Tdomain%specel(n)%win_mirror(1:ngllx-2,1:nglly-2,0)

nf = Tdomain%specel(n)%near_faces(1)
pml = .not.Tdomain%sFace(nf)%PML
if(pml)Tdomain%sFace(nf)%sSimu(0)%Displ(1:ngllx-2,1:ngllz-2,comp) = Tdomain%sFace(nf)%sSimu(0)%Displ(1:ngllx-2,1:ngllz-2,comp)&
                                                                  *Tdomain%specel(n)%win_mirror(1:ngllx-2,0,1:ngllz-2)

nf = Tdomain%specel(n)%near_faces(2)
pml = .not.Tdomain%sFace(nf)%PML
if(pml)Tdomain%sFace(nf)%sSimu(0)%Displ(1:nglly-2,1:ngllz-2,comp) = Tdomain%sFace(nf)%sSimu(0)%Displ(1:nglly-2,1:ngllz-2,comp)&
                                                                  *Tdomain%specel(n)%win_mirror(ngllx-1,1:nglly-2,1:ngllz-2)

nf = Tdomain%specel(n)%near_faces(3)
pml = .not.Tdomain%sFace(nf)%PML
if(pml)Tdomain%sFace(nf)%sSimu(0)%Displ(1:ngllx-2,1:ngllz-2,comp) = Tdomain%sFace(nf)%sSimu(0)%Displ(1:ngllx-2,1:ngllz-2,comp)&
                                                                  *Tdomain%specel(n)%win_mirror(1:ngllx-2,nglly-1,1:ngllz-2)

nf = Tdomain%specel(n)%near_faces(4)
pml = .not.Tdomain%sFace(nf)%PML
if(pml)Tdomain%sFace(nf)%sSimu(0)%Displ(1:nglly-2,1:ngllz-2,comp) = Tdomain%sFace(nf)%sSimu(0)%Displ(1:nglly-2,1:ngllz-2,comp)&
                                                                  *Tdomain%specel(n)%win_mirror(0,1:nglly-2,1:ngllz-2)

nf = Tdomain%specel(n)%near_faces(5)
pml = .not.Tdomain%sFace(nf)%PML
if(pml)Tdomain%sFace(nf)%sSimu(0)%Displ(1:ngllx-2,1:nglly-2,comp) = Tdomain%sFace(nf)%sSimu(0)%Displ(1:ngllx-2,1:nglly-2,comp)&
                                                                  *Tdomain%specel(n)%win_mirror(1:ngllx-2,1:nglly-2,ngllz-1)
!

ne = Tdomain%specel(n)%near_edges(0)
pml = .not.Tdomain%sEdge(ne)%PML
if(pml)Tdomain%sEdge(ne)%sSimu(0)%Displ(1:ngllx-2,comp) = Tdomain%sEdge(ne)%sSimu(0)%Displ(1:ngllx-2,comp)&
                                                        *Tdomain%specel(n)%win_mirror(1:ngllx-2,0,0) 

ne = Tdomain%specel(n)%near_edges(1)
pml = .not.Tdomain%sEdge(ne)%PML
if(pml)Tdomain%sEdge(ne)%sSimu(0)%Displ(1:nglly-2,comp) = Tdomain%sEdge(ne)%sSimu(0)%Displ(1:nglly-2,comp)&
                                                         *Tdomain%specel(n)%win_mirror(ngllx-1,1:nglly-2,0)

ne = Tdomain%specel(n)%near_edges(2)
pml = .not.Tdomain%sEdge(ne)%PML
if(pml)Tdomain%sEdge(ne)%sSimu(0)%Displ(1:ngllx-2,comp) = Tdomain%sEdge(ne)%sSimu(0)%Displ(1:ngllx-2,comp)&
                                                        *Tdomain%specel(n)%win_mirror(1:ngllx-2,nglly-1,0) 

ne = Tdomain%specel(n)%near_edges(3)
pml = .not.Tdomain%sEdge(ne)%PML
if(pml)Tdomain%sEdge(ne)%sSimu(0)%Displ(1:nglly-2,comp) = Tdomain%sEdge(ne)%sSimu(0)%Displ(1:nglly-2,comp)&
                                                         *Tdomain%specel(n)%win_mirror(0,1:nglly-2,0) 

ne = Tdomain%specel(n)%near_edges(4)
pml = .not.Tdomain%sEdge(ne)%PML
if(pml)Tdomain%sEdge(ne)%sSimu(0)%Displ(1:ngllz-2,comp) = Tdomain%sEdge(ne)%sSimu(0)%Displ(1:ngllz-2,comp)&
                                                        *Tdomain%specel(n)%win_mirror(ngllx-1,0,1:ngllz-2) 

ne = Tdomain%specel(n)%near_edges(5)
pml = .not.Tdomain%sEdge(ne)%PML
if(pml)Tdomain%sEdge(ne)%sSimu(0)%Displ(1:ngllx-2,comp) = Tdomain%sEdge(ne)%sSimu(0)%Displ(1:ngllx-2,comp)&
                                                        *Tdomain%specel(n)%win_mirror(1:ngllx-2,0,ngllz-1) 

ne = Tdomain%specel(n)%near_edges(6)
pml = .not.Tdomain%sEdge(ne)%PML
if(pml)Tdomain%sEdge(ne)%sSimu(0)%Displ(1:ngllz-2,comp) = Tdomain%sEdge(ne)%sSimu(0)%Displ(1:ngllz-2,comp)&
                                                        *Tdomain%specel(n)%win_mirror(0,0,1:ngllz-2) 

ne = Tdomain%specel(n)%near_edges(7)
pml = .not.Tdomain%sEdge(ne)%PML
if(pml)Tdomain%sEdge(ne)%sSimu(0)%Displ(1:ngllz-2,comp) = Tdomain%sEdge(ne)%sSimu(0)%Displ(1:ngllz-2,comp)&
                                                         * Tdomain%specel(n)%win_mirror(ngllx-1,nglly-1,1:ngllz-2)

ne = Tdomain%specel(n)%near_edges(8)
pml = .not.Tdomain%sEdge(ne)%PML
if(pml)Tdomain%sEdge(ne)%sSimu(0)%Displ(1:nglly-2,comp) = Tdomain%sEdge(ne)%sSimu(0)%Displ(1:nglly-2,comp)&
                                                         *Tdomain%specel(n)%win_mirror(ngllx-1,1:nglly-2,ngllz-1)

ne = Tdomain%specel(n)%near_edges(9)
pml = .not.Tdomain%sEdge(ne)%PML
if(pml)Tdomain%sEdge(ne)%sSimu(0)%Displ(1:ngllx-2,comp) = Tdomain%sEdge(ne)%sSimu(0)%Displ(1:ngllx-2,comp)&
                                                        *Tdomain%specel(n)%win_mirror(1:ngllx-2,nglly-1,ngllz-1) 

ne = Tdomain%specel(n)%near_edges(10)
pml = .not.Tdomain%sEdge(ne)%PML
if(pml)Tdomain%sEdge(ne)%sSimu(0)%Displ(1:ngllz-2,comp) = Tdomain%sEdge(ne)%sSimu(0)%Displ(1:ngllz-2,comp)&
                                                        *Tdomain%specel(n)%win_mirror(0,nglly-1,1:ngllz-2) 

ne = Tdomain%specel(n)%near_edges(11)
pml = .not.Tdomain%sEdge(ne)%PML
if(pml)Tdomain%sEdge(ne)%sSimu(0)%Displ(1:nglly-2,comp) = Tdomain%sEdge(ne)%sSimu(0)%Displ(1:nglly-2,comp)&
                                                        *Tdomain%specel(n)%win_mirror(0,1:nglly-2,ngllz-1) 
!

nv = Tdomain%specel(n)%near_vertices(0)
pml = .not.Tdomain%sVertex(nv)%PML
if(pml)Tdomain%sVertex(nv)%sSimu(0)%Displ(comp)= Tdomain%sVertex(nv)%sSimu(0)%Displ(comp)&
                                                *Tdomain%specel(n)%win_mirror(0,0,0) 

nv = Tdomain%specel(n)%near_vertices(1)
pml = .not.Tdomain%sVertex(nv)%PML
if(pml)Tdomain%sVertex(nv)%sSimu(0)%Displ(comp) = Tdomain%sVertex(nv)%sSimu(0)%Displ(comp)&
                                                *Tdomain%specel(n)%win_mirror(ngllx-1,0,0)

nv = Tdomain%specel(n)%near_vertices(2)
pml = .not.Tdomain%sVertex(nv)%PML
if(pml)Tdomain%sVertex(nv)%sSimu(0)%Displ(comp) = Tdomain%sVertex(nv)%sSimu(0)%Displ(comp)&
                                                *Tdomain%specel(n)%win_mirror(ngllx-1,nglly-1,0) 

nv = Tdomain%specel(n)%near_vertices(3)
pml = .not.Tdomain%sVertex(nv)%PML
if(pml)Tdomain%sVertex(nv)%sSimu(0)%Displ(comp) = Tdomain%sVertex(nv)%sSimu(0)%Displ(comp)&
                                                *Tdomain%specel(n)%win_mirror(0,nglly-1,0)

nv = Tdomain%specel(n)%near_vertices(4)
pml = .not.Tdomain%sVertex(nv)%PML
if(pml)Tdomain%sVertex(nv)%sSimu(0)%Displ(comp) = Tdomain%sVertex(nv)%sSimu(0)%Displ(comp)&
                                                 *Tdomain%specel(n)%win_mirror(0,0,ngllz-1)

nv = Tdomain%specel(n)%near_vertices(5)
pml = .not.Tdomain%sVertex(nv)%PML
if(pml)Tdomain%sVertex(nv)%sSimu(0)%Displ(comp) = Tdomain%sVertex(nv)%sSimu(0)%Displ(comp)&
                                                 *Tdomain%specel(n)%win_mirror(ngllx-1,0,ngllz-1)

nv = Tdomain%specel(n)%near_vertices(6)
pml = .not.Tdomain%sVertex(nv)%PML
if(pml)Tdomain%sVertex(nv)%sSimu(0)%Displ(comp) = Tdomain%sVertex(nv)%sSimu(0)%Displ(comp)&
                                                 *Tdomain%specel(n)%win_mirror(ngllx-1,nglly-1,ngllz-1)

nv = Tdomain%specel(n)%near_vertices(7)
pml = .not.Tdomain%sVertex(nv)%PML
if(pml)Tdomain%sVertex(nv)%sSimu(0)%Displ(comp) = Tdomain%sVertex(nv)%sSimu(0)%Displ(comp)&
                                                *Tdomain%specel(n)%win_mirror(0,nglly-1,ngllz-1)
!
! veloc
!
nf = Tdomain%specel(n)%near_faces(0)
Tdomain%sFace(nf)%sSimu(0)%Veloc(1:ngllx-2,1:nglly-2,comp) = Tdomain%sFace(nf)%sSimu(0)%Veloc(1:ngllx-2,1:nglly-2,comp)&
                                                                  *Tdomain%specel(n)%win_mirror(1:ngllx-2,1:nglly-2,0)

nf = Tdomain%specel(n)%near_faces(1)
Tdomain%sFace(nf)%sSimu(0)%Veloc(1:ngllx-2,1:ngllz-2,comp) = Tdomain%sFace(nf)%sSimu(0)%Veloc(1:ngllx-2,1:ngllz-2,comp)&
                                                                  *Tdomain%specel(n)%win_mirror(1:ngllx-2,0,1:ngllz-2)

nf = Tdomain%specel(n)%near_faces(2)
Tdomain%sFace(nf)%sSimu(0)%Veloc(1:nglly-2,1:ngllz-2,comp) = Tdomain%sFace(nf)%sSimu(0)%Veloc(1:nglly-2,1:ngllz-2,comp)&
                                                                  *Tdomain%specel(n)%win_mirror(ngllx-1,1:nglly-2,1:ngllz-2)

nf = Tdomain%specel(n)%near_faces(3)
Tdomain%sFace(nf)%sSimu(0)%Veloc(1:ngllx-2,1:ngllz-2,comp) = Tdomain%sFace(nf)%sSimu(0)%Veloc(1:ngllx-2,1:ngllz-2,comp)&
                                                                  *Tdomain%specel(n)%win_mirror(1:ngllx-2,nglly-1,1:ngllz-2)

nf = Tdomain%specel(n)%near_faces(4)
Tdomain%sFace(nf)%sSimu(0)%Veloc(1:nglly-2,1:ngllz-2,comp) = Tdomain%sFace(nf)%sSimu(0)%Veloc(1:nglly-2,1:ngllz-2,comp)&
                                                                  *Tdomain%specel(n)%win_mirror(0,1:nglly-2,1:ngllz-2)

nf = Tdomain%specel(n)%near_faces(5)
Tdomain%sFace(nf)%sSimu(0)%Veloc(1:ngllx-2,1:nglly-2,comp) = Tdomain%sFace(nf)%sSimu(0)%Veloc(1:ngllx-2,1:nglly-2,comp)&
                                                                  *Tdomain%specel(n)%win_mirror(1:ngllx-2,1:nglly-2,ngllz-1)
!

ne = Tdomain%specel(n)%near_edges(0)
Tdomain%sEdge(ne)%sSimu(0)%Veloc(1:ngllx-2,comp) = Tdomain%sEdge(ne)%sSimu(0)%Veloc(1:ngllx-2,comp)&
                                                        *Tdomain%specel(n)%win_mirror(1:ngllx-2,0,0) 

ne = Tdomain%specel(n)%near_edges(1)
 Tdomain%sEdge(ne)%sSimu(0)%Veloc(1:nglly-2,comp) = Tdomain%sEdge(ne)%sSimu(0)%Veloc(1:nglly-2,comp)&
                                                         *Tdomain%specel(n)%win_mirror(ngllx-1,1:nglly-2,0)

ne = Tdomain%specel(n)%near_edges(2)
Tdomain%sEdge(ne)%sSimu(0)%Veloc(1:ngllx-2,comp) = Tdomain%sEdge(ne)%sSimu(0)%Veloc(1:ngllx-2,comp)&
                                                        *Tdomain%specel(n)%win_mirror(1:ngllx-2,nglly-1,0) 

ne = Tdomain%specel(n)%near_edges(3)
Tdomain%sEdge(ne)%sSimu(0)%Veloc(1:nglly-2,comp) = Tdomain%sEdge(ne)%sSimu(0)%Veloc(1:nglly-2,comp)&
                                                         *Tdomain%specel(n)%win_mirror(0,1:nglly-2,0) 

ne = Tdomain%specel(n)%near_edges(4)
Tdomain%sEdge(ne)%sSimu(0)%Veloc(1:ngllz-2,comp) = Tdomain%sEdge(ne)%sSimu(0)%Veloc(1:ngllz-2,comp)&
                                                        *Tdomain%specel(n)%win_mirror(ngllx-1,0,1:ngllz-2) 

ne = Tdomain%specel(n)%near_edges(5)
Tdomain%sEdge(ne)%sSimu(0)%Veloc(1:ngllx-2,comp) = Tdomain%sEdge(ne)%sSimu(0)%Veloc(1:ngllx-2,comp)&
                                                        *Tdomain%specel(n)%win_mirror(1:ngllx-2,0,ngllz-1) 

ne = Tdomain%specel(n)%near_edges(6)
Tdomain%sEdge(ne)%sSimu(0)%Veloc(1:ngllz-2,comp) = Tdomain%sEdge(ne)%sSimu(0)%Veloc(1:ngllz-2,comp)&
                                                        *Tdomain%specel(n)%win_mirror(0,0,1:ngllz-2) 

ne = Tdomain%specel(n)%near_edges(7)
 Tdomain%sEdge(ne)%sSimu(0)%Veloc(1:ngllz-2,comp) = Tdomain%sEdge(ne)%sSimu(0)%Veloc(1:ngllz-2,comp)&
                                                         * Tdomain%specel(n)%win_mirror(ngllx-1,nglly-1,1:ngllz-2)

ne = Tdomain%specel(n)%near_edges(8)
 Tdomain%sEdge(ne)%sSimu(0)%Veloc(1:nglly-2,comp) = Tdomain%sEdge(ne)%sSimu(0)%Veloc(1:nglly-2,comp)&
                                                         *Tdomain%specel(n)%win_mirror(ngllx-1,1:nglly-2,ngllz-1)

ne = Tdomain%specel(n)%near_edges(9)
Tdomain%sEdge(ne)%sSimu(0)%Veloc(1:ngllx-2,comp) = Tdomain%sEdge(ne)%sSimu(0)%Veloc(1:ngllx-2,comp)&
                                                        *Tdomain%specel(n)%win_mirror(1:ngllx-2,nglly-1,ngllz-1) 

ne = Tdomain%specel(n)%near_edges(10)
Tdomain%sEdge(ne)%sSimu(0)%Veloc(1:ngllz-2,comp) = Tdomain%sEdge(ne)%sSimu(0)%Veloc(1:ngllz-2,comp)&
                                                        *Tdomain%specel(n)%win_mirror(0,nglly-1,1:ngllz-2) 

ne = Tdomain%specel(n)%near_edges(11)
Tdomain%sEdge(ne)%sSimu(0)%Veloc(1:nglly-2,comp) = Tdomain%sEdge(ne)%sSimu(0)%Veloc(1:nglly-2,comp)&
                                                        *Tdomain%specel(n)%win_mirror(0,1:nglly-2,ngllz-1) 
!

nv = Tdomain%specel(n)%near_vertices(0)
 Tdomain%sVertex(nv)%sSimu(0)%Veloc(comp)= Tdomain%sVertex(nv)%sSimu(0)%Veloc(comp)&
                                                *Tdomain%specel(n)%win_mirror(0,0,0) 

nv = Tdomain%specel(n)%near_vertices(1)
Tdomain%sVertex(nv)%sSimu(0)%Veloc(comp) = Tdomain%sVertex(nv)%sSimu(0)%Veloc(comp)&
                                                *Tdomain%specel(n)%win_mirror(ngllx-1,0,0)

nv = Tdomain%specel(n)%near_vertices(2)
Tdomain%sVertex(nv)%sSimu(0)%Veloc(comp) = Tdomain%sVertex(nv)%sSimu(0)%Veloc(comp)&
                                                *Tdomain%specel(n)%win_mirror(ngllx-1,nglly-1,0) 

nv = Tdomain%specel(n)%near_vertices(3)
Tdomain%sVertex(nv)%sSimu(0)%Veloc(comp) = Tdomain%sVertex(nv)%sSimu(0)%Veloc(comp)&
                                                *Tdomain%specel(n)%win_mirror(0,nglly-1,0)

nv = Tdomain%specel(n)%near_vertices(4)
 Tdomain%sVertex(nv)%sSimu(0)%Veloc(comp) = Tdomain%sVertex(nv)%sSimu(0)%Veloc(comp)&
                                                 *Tdomain%specel(n)%win_mirror(0,0,ngllz-1)

nv = Tdomain%specel(n)%near_vertices(5)
 Tdomain%sVertex(nv)%sSimu(0)%Veloc(comp) = Tdomain%sVertex(nv)%sSimu(0)%Veloc(comp)&
                                                 *Tdomain%specel(n)%win_mirror(ngllx-1,0,ngllz-1)

nv = Tdomain%specel(n)%near_vertices(6)
 Tdomain%sVertex(nv)%sSimu(0)%Veloc(comp) = Tdomain%sVertex(nv)%sSimu(0)%Veloc(comp)&
                                                 *Tdomain%specel(n)%win_mirror(ngllx-1,nglly-1,ngllz-1)

nv = Tdomain%specel(n)%near_vertices(7)
Tdomain%sVertex(nv)%sSimu(0)%Veloc(comp) = Tdomain%sVertex(nv)%sSimu(0)%Veloc(comp)&
                                                *Tdomain%specel(n)%win_mirror(0,nglly-1,ngllz-1)

!
enddo
!
enddo
!
endif
!
if (ntime==Tdomain%sTimeParam%ntime-1)   close (45)
!
!=========================
end subroutine impose_init
!=========================
