!--------------------------------------------------------------
subroutine put_init_cond (Tdomain)
!--------------------------------------------------------------

  use sdomains
  use angles
  use init_cond

  implicit none

  type (domain), intent (INOUT) :: Tdomain

  logical, dimension(:), allocatable :: L_Face, L_Edge, L_Vertex
  integer :: i,j,k, n, idef, ngllx,nglly,ngllz, mat, a,b, nf,ne,nv, ngll1,ngll2, ngll
  doubleprecision :: x,y,z, xa,ya,za, r,theta,phi, phi_bis, ct,st,cp,sp, colat_aigle
  doubleprecision, dimension(3) :: U,V
  doubleprecision, dimension(0:2) :: tmpU,tmpV
  doubleprecision, dimension(0:2,0:2) :: Pcs, Rot,tRot, rot00,tRot00
  doubleprecision, dimension(:,:,:), allocatable :: Rsph
  doubleprecision, dimension(:,:,:,:), allocatable :: displ_loc,veloc_loc


Rot = Tdomain%rot
tRot = transpose(Rot)


if (Tdomain%curve) then
    call init_init_cond()
else
    write (*,*) "Initial conditions cannot be used in a cuboid so far."
    stop
endif


allocate (L_Face(0:Tdomain%n_face-1))
L_Face = .true.
allocate (L_Edge(0:Tdomain%n_edge-1))
L_Edge = .true.
allocate (L_Vertex(0:Tdomain%n_vertex-1))
L_Vertex = .true.
do n = 0,Tdomain%n_elem-1

   mat = Tdomain%specel(n)%mat_index
   ngllx = Tdomain%specel(n)%ngllx
   nglly = Tdomain%specel(n)%nglly
   ngllz = Tdomain%specel(n)%ngllz

   if (Tdomain%ellipticity) then
       allocate (Rsph(0:ngllx-1,0:nglly-1,0:ngllz-1))
       call r_spheric(ngllx,nglly,ngllz,Tdomain,mat,n,Rsph)
   endif

   allocate (displ_loc(0:ngllx-1,0:nglly-1,0:ngllz-1,0:2))
   allocate (veloc_loc(0:ngllx-1,0:nglly-1,0:ngllz-1,0:2))
   do k = 0,ngllz-1
    do j = 0,nglly-1
     do i = 0,ngllx-1
        ! Taking cartesian coordinates of the GLL
        idef = Tdomain%specel(n)%Iglobnum(i,j,k)
        x = Tdomain%Globcoord(0,idef)
        y = Tdomain%Globcoord(1,idef)
        z = Tdomain%Globcoord(2,idef)
        ! Coordinates of the GLL in the real chunk
        xa = x;   ya = y;   za = z
        x = Rot(0,0)*xa + Rot(0,1)*ya + Rot(0,2)*za
        y = Rot(1,0)*xa + Rot(1,1)*ya + Rot(1,2)*za
        z = Rot(2,0)*xa + Rot(2,1)*ya + Rot(2,2)*za
        ! Convert cartesian coordinates of the GLL into spherical coordinates
        call cart2sph(x,y,z,r,theta,phi)   !!! phi dans [0;360]
        phi_bis = phi;   if (phi_bis>pi) phi_bis = phi_bis - 2*pi   !!! phi_bis dans [-180;180]
        ! Getting the initial conditions at the GLL
        if (Tdomain%ellipticity) r = Rsph(i,j,k)
        call get_init_cond(r,theta,phi_bis,U,V)
        ! Convert spherical coordinates of the vectors U and V into cartesian coordinates
        ct=dcos(theta)
        st=dsin(theta)
        cp=dcos(phi)
        sp=dsin(phi)
        Pcs(0,0) = st*cp; Pcs(0,1) = ct*cp; Pcs(0,2) = -sp
        Pcs(1,0) = st*sp; Pcs(1,1) = ct*sp; Pcs(1,2) = cp
        Pcs(2,0) = ct   ; Pcs(2,1) = -st  ; Pcs(2,2) = 0.0d0
        do a = 0,2
            tmpU(a) = 0.0d0
            tmpV(a) = 0.0d0
            do b = 0,2
                tmpU(a) = tmpU(a) + Pcs(a,b)*U(b+1)
                tmpV(a) = tmpV(a) + Pcs(a,b)*V(b+1)
            enddo
        enddo
        ! Rotation (du chunk reel au chunk de reference) des vecteurs U et V
        U(1) = tRot(0,0)*tmpU(0) + tRot(0,1)*tmpU(1) + tRot(0,2)*tmpU(2)
        U(2) = tRot(1,0)*tmpU(0) + tRot(1,1)*tmpU(1) + tRot(1,2)*tmpU(2)
        U(3) = tRot(2,0)*tmpU(0) + tRot(2,1)*tmpU(1) + tRot(2,2)*tmpU(2)
        V(1) = tRot(0,0)*tmpV(0) + tRot(0,1)*tmpV(1) + tRot(0,2)*tmpV(2)
        V(2) = tRot(1,0)*tmpV(0) + tRot(1,1)*tmpV(1) + tRot(1,2)*tmpV(2)
        V(3) = tRot(2,0)*tmpV(0) + tRot(2,1)*tmpV(1) + tRot(2,2)*tmpV(2)
        ! Stockage du resultat
        displ_loc(i,j,k,0:2) = U(1:3)
        veloc_loc(i,j,k,0:2) = V(1:3)
     enddo
    enddo
   enddo

   if (.not. Tdomain%specel(n)%PML) then
      Tdomain%specel(n)%sSimu(0)%Displ(1:ngllx-2,1:nglly-2,1:ngllz-2,0:2) = &
       Tdomain%specel(n)%sSimu(0)%Displ(1:ngllx-2,1:nglly-2,1:ngllz-2,0:2) + displ_loc(1:ngllx-2,1:nglly-2,1:ngllz-2,0:2)
      Tdomain%specel(n)%sSimu(0)%Veloc(1:ngllx-2,1:nglly-2,1:ngllz-2,0:2) = &
       Tdomain%specel(n)%sSimu(0)%Veloc(1:ngllx-2,1:nglly-2,1:ngllz-2,0:2) + veloc_loc(1:ngllx-2,1:nglly-2,1:ngllz-2,0:2)
   endif

   do i = 0,5
       nf = Tdomain%specel(n)%Near_Faces(i)
       if (L_Face(nf)) then
           L_Face(nf) = .false.
           if (.not.Tdomain%sface(nf)%PML) then
               ngll1 = Tdomain%sface(nf)%ngll1
               ngll2 = Tdomain%sface(nf)%ngll2
               select case (i)
                case (0)
                   Tdomain%sface(nf)%sSimu(0)%Displ(1:ngll1-2,1:ngll2-2,0:2) = &
                    Tdomain%sFace(nf)%sSimu(0)%Displ(1:ngll1-2,1:ngll2-2,0:2) + displ_loc(1:ngll1-2,1:ngll2-2,0,0:2)
                   Tdomain%sface(nf)%sSimu(0)%Veloc(1:ngll1-2,1:ngll2-2,0:2) = &
                    Tdomain%sFace(nf)%sSimu(0)%Veloc(1:ngll1-2,1:ngll2-2,0:2) + veloc_loc(1:ngll1-2,1:ngll2-2,0,0:2)
                case (1)
                   Tdomain%sFace(nf)%sSimu(0)%Displ(1:ngll1-2,1:ngll2-2,0:2) = &
                    Tdomain%sFace(nf)%sSimu(0)%Displ(1:ngll1-2,1:ngll2-2,0:2) + displ_loc(1:ngll1-2,0,1:ngll2-2,0:2)
                   Tdomain%sFace(nf)%sSimu(0)%Veloc(1:ngll1-2,1:ngll2-2,0:2) = &
                    Tdomain%sFace(nf)%sSimu(0)%Veloc(1:ngll1-2,1:ngll2-2,0:2) + veloc_loc(1:ngll1-2,0,1:ngll2-2,0:2)
                case (2)
                   Tdomain%sFace(nf)%sSimu(0)%Displ(1:ngll1-2,1:ngll2-2,0:2) = &
                    Tdomain%sFace(nf)%sSimu(0)%Displ(1:ngll1-2,1:ngll2-2,0:2) + displ_loc(ngllx-1,1:ngll1-2,1:ngll2-2,0:2)
                   Tdomain%sFace(nf)%sSimu(0)%Veloc(1:ngll1-2,1:ngll2-2,0:2) = &
                    Tdomain%sFace(nf)%sSimu(0)%Veloc(1:ngll1-2,1:ngll2-2,0:2) + veloc_loc(ngllx-1,1:ngll1-2,1:ngll2-2,0:2)
                case (3)
                   Tdomain%sFace(nf)%sSimu(0)%Displ(1:ngll1-2,1:ngll2-2,0:2) = &
                    Tdomain%sFace(nf)%sSimu(0)%Displ(1:ngll1-2,1:ngll2-2,0:2) + displ_loc(1:ngll1-2,nglly-1,1:ngll2-2,0:2)
                   Tdomain%sFace(nf)%sSimu(0)%Veloc(1:ngll1-2,1:ngll2-2,0:2) = &
                    Tdomain%sFace(nf)%sSimu(0)%Veloc(1:ngll1-2,1:ngll2-2,0:2) + veloc_loc(1:ngll1-2,nglly-1,1:ngll2-2,0:2)
                case (4)
                   Tdomain%sFace(nf)%sSimu(0)%Displ(1:ngll1-2,1:ngll2-2,0:2) = &
                    Tdomain%sFace(nf)%sSimu(0)%Displ(1:ngll1-2,1:ngll2-2,0:2) + displ_loc(0,1:ngll1-2,1:ngll2-2,0:2)
                   Tdomain%sFace(nf)%sSimu(0)%Veloc(1:ngll1-2,1:ngll2-2,0:2) = &
                    Tdomain%sFace(nf)%sSimu(0)%Veloc(1:ngll1-2,1:ngll2-2,0:2) + veloc_loc(0,1:ngll1-2,1:ngll2-2,0:2)
                case (5)
                   Tdomain%sFace(nf)%sSimu(0)%Displ(1:ngll1-2,1:ngll2-2,0:2) = &
                    Tdomain%sFace(nf)%sSimu(0)%Displ(1:ngll1-2,1:ngll2-2,0:2) + displ_loc(1:ngll1-2,1:ngll2-2,ngllz-1,0:2)
                   Tdomain%sFace(nf)%sSimu(0)%Veloc(1:ngll1-2,1:ngll2-2,0:2) = &
                    Tdomain%sFace(nf)%sSimu(0)%Veloc(1:ngll1-2,1:ngll2-2,0:2) + veloc_loc(1:ngll1-2,1:ngll2-2,ngllz-1,0:2)
               end select
           endif
       endif
   enddo

   do i = 0,11
       ne = Tdomain%specel(n)%Near_Edges(i)
       if (L_Edge(ne)) then
           L_Edge(ne) = .false.
           if (.not.Tdomain%sedge(ne)%PML) then
               ngll = Tdomain%sEdge(ne)%ngll
               select case (i)
                case (0)
                   Tdomain%sEdge(ne)%sSimu(0)%Displ(1:ngll-2,0:2) = &
                    Tdomain%sEdge(ne)%sSimu(0)%Displ(1:ngll-2,0:2) + displ_loc(1:ngll-2,0,0,0:2)
                   Tdomain%sEdge(ne)%sSimu(0)%Veloc(1:ngll-2,0:2) = &
                    Tdomain%sEdge(ne)%sSimu(0)%Veloc(1:ngll-2,0:2) + veloc_loc(1:ngll-2,0,0,0:2)
                case (1)
                   Tdomain%sEdge(ne)%sSimu(0)%Displ(1:ngll-2,0:2) = &
                    Tdomain%sEdge(ne)%sSimu(0)%Displ(1:ngll-2,0:2) + displ_loc(ngllx-1,1:ngll-2,0,0:2)
                   Tdomain%sEdge(ne)%sSimu(0)%Veloc(1:ngll-2,0:2) = &
                    Tdomain%sEdge(ne)%sSimu(0)%Veloc(1:ngll-2,0:2) + veloc_loc(ngllx-1,1:ngll-2,0,0:2)
                case (2)
                   Tdomain%sEdge(ne)%sSimu(0)%Displ(1:ngll-2,0:2) = &
                    Tdomain%sEdge(ne)%sSimu(0)%Displ(1:ngll-2,0:2) + displ_loc(1:ngll-2,nglly-1,0,0:2)
                   Tdomain%sEdge(ne)%sSimu(0)%Veloc(1:ngll-2,0:2) = &
                    Tdomain%sEdge(ne)%sSimu(0)%Veloc(1:ngll-2,0:2) + veloc_loc(1:ngll-2,nglly-1,0,0:2)
                case (3)
                   Tdomain%sEdge(ne)%sSimu(0)%Displ(1:ngll-2,0:2) = &
                    Tdomain%sEdge(ne)%sSimu(0)%Displ(1:ngll-2,0:2) + displ_loc(0,1:ngll-2,0,0:2)
                   Tdomain%sEdge(ne)%sSimu(0)%Veloc(1:ngll-2,0:2) = &
                    Tdomain%sEdge(ne)%sSimu(0)%Veloc(1:ngll-2,0:2) + veloc_loc(0,1:ngll-2,0,0:2)
                case (4)
                   Tdomain%sEdge(ne)%sSimu(0)%Displ(1:ngll-2,0:2) = &
                    Tdomain%sEdge(ne)%sSimu(0)%Displ(1:ngll-2,0:2) + displ_loc(ngllx-1,0,1:ngll-2,0:2)
                   Tdomain%sEdge(ne)%sSimu(0)%Veloc(1:ngll-2,0:2) = &
                    Tdomain%sEdge(ne)%sSimu(0)%Veloc(1:ngll-2,0:2) + veloc_loc(ngllx-1,0,1:ngll-2,0:2)
                case (5)
                   Tdomain%sEdge(ne)%sSimu(0)%Displ(1:ngll-2,0:2) = &
                    Tdomain%sEdge(ne)%sSimu(0)%Displ(1:ngll-2,0:2) + displ_loc(1:ngll-2,0,ngllz-1,0:2)
                   Tdomain%sEdge(ne)%sSimu(0)%Veloc(1:ngll-2,0:2) = &
                    Tdomain%sEdge(ne)%sSimu(0)%Veloc(1:ngll-2,0:2) + veloc_loc(1:ngll-2,0,ngllz-1,0:2)
                case (6)
                   Tdomain%sEdge(ne)%sSimu(0)%Displ(1:ngll-2,0:2) = &
                    Tdomain%sEdge(ne)%sSimu(0)%Displ(1:ngll-2,0:2) + displ_loc(0,0,1:ngll-2,0:2)
                   Tdomain%sEdge(ne)%sSimu(0)%Veloc(1:ngll-2,0:2) = &
                    Tdomain%sEdge(ne)%sSimu(0)%Veloc(1:ngll-2,0:2) + veloc_loc(0,0,1:ngll-2,0:2)
                case (7)
                   Tdomain%sEdge(ne)%sSimu(0)%Displ(1:ngll-2,0:2) = &
                    Tdomain%sEdge(ne)%sSimu(0)%Displ(1:ngll-2,0:2) + displ_loc(ngllx-1,nglly-1,1:ngll-2,0:2)
                   Tdomain%sEdge(ne)%sSimu(0)%Veloc(1:ngll-2,0:2) = &
                    Tdomain%sEdge(ne)%sSimu(0)%Veloc(1:ngll-2,0:2) + veloc_loc(ngllx-1,nglly-1,1:ngll-2,0:2)
                case (8)
                   Tdomain%sEdge(ne)%sSimu(0)%Displ(1:ngll-2,0:2) = &
                    Tdomain%sEdge(ne)%sSimu(0)%Displ(1:ngll-2,0:2) + displ_loc(ngllx-1,1:ngll-2,ngllz-1,0:2)
                   Tdomain%sEdge(ne)%sSimu(0)%Veloc(1:ngll-2,0:2) = &
                    Tdomain%sEdge(ne)%sSimu(0)%Veloc(1:ngll-2,0:2) + veloc_loc(ngllx-1,1:ngll-2,ngllz-1,0:2)
                case (9)
                   Tdomain%sEdge(ne)%sSimu(0)%Displ(1:ngll-2,0:2) = &
                    Tdomain%sEdge(ne)%sSimu(0)%Displ(1:ngll-2,0:2) + displ_loc(1:ngll-2,nglly-1,ngllz-1,0:2)
                   Tdomain%sEdge(ne)%sSimu(0)%Veloc(1:ngll-2,0:2) = &
                    Tdomain%sEdge(ne)%sSimu(0)%Veloc(1:ngll-2,0:2) + veloc_loc(1:ngll-2,nglly-1,ngllz-1,0:2)
                case (10)
                   Tdomain%sEdge(ne)%sSimu(0)%Displ(1:ngll-2,0:2) = &
                    Tdomain%sEdge(ne)%sSimu(0)%Displ(1:ngll-2,0:2) + displ_loc(0,nglly-1,1:ngll-2,0:2)
                   Tdomain%sEdge(ne)%sSimu(0)%Veloc(1:ngll-2,0:2) = &
                    Tdomain%sEdge(ne)%sSimu(0)%Veloc(1:ngll-2,0:2) + veloc_loc(0,nglly-1,1:ngll-2,0:2)
                case (11)
                   Tdomain%sEdge(ne)%sSimu(0)%Displ(1:ngll-2,0:2) = &
                    Tdomain%sEdge(ne)%sSimu(0)%Displ(1:ngll-2,0:2) + displ_loc(0,1:ngll-2,ngllz-1,0:2)
                   Tdomain%sEdge(ne)%sSimu(0)%Veloc(1:ngll-2,0:2) = &
                    Tdomain%sEdge(ne)%sSimu(0)%Veloc(1:ngll-2,0:2) + veloc_loc(0,1:ngll-2,ngllz-1,0:2)
               end select
           endif
       endif
   enddo

   do i = 0,7
       nv = Tdomain%specel(n)%Near_Vertices(i)
       if (L_Vertex(nv)) then
           L_Vertex(nv) = .false.
           if (.not.Tdomain%svertex(nv)%PML) then
               select case (i)
                case (0)
                   Tdomain%sVertex(nv)%sSimu(0)%Displ(0:2) = &
                    Tdomain%sVertex(nv)%sSimu(0)%Displ(0:2) + displ_loc(0,0,0,0:2)
                   Tdomain%sVertex(nv)%sSimu(0)%Veloc(0:2) = &
                    Tdomain%sVertex(nv)%sSimu(0)%Veloc(0:2) + veloc_loc(0,0,0,0:2)
                case (1)
                   Tdomain%sVertex(nv)%sSimu(0)%Displ(0:2) = &
                    Tdomain%sVertex(nv)%sSimu(0)%Displ(0:2) + displ_loc(ngllx-1,0,0,0:2)
                   Tdomain%sVertex(nv)%sSimu(0)%Veloc(0:2) = &
                    Tdomain%sVertex(nv)%sSimu(0)%Veloc(0:2) + veloc_loc(ngllx-1,0,0,0:2)
                case (2)
                   Tdomain%sVertex(nv)%sSimu(0)%Displ(0:2) = &
                    Tdomain%sVertex(nv)%sSimu(0)%Displ(0:2) + displ_loc(ngllx-1,nglly-1,0,0:2)
                   Tdomain%sVertex(nv)%sSimu(0)%Veloc(0:2) = &
                    Tdomain%sVertex(nv)%sSimu(0)%Veloc(0:2) + veloc_loc(ngllx-1,nglly-1,0,0:2)
                case (3)
                   Tdomain%sVertex(nv)%sSimu(0)%Displ(0:2) = &
                    Tdomain%sVertex(nv)%sSimu(0)%Displ(0:2) + displ_loc(0,nglly-1,0,0:2)
                   Tdomain%sVertex(nv)%sSimu(0)%Veloc(0:2) = &
                    Tdomain%sVertex(nv)%sSimu(0)%Veloc(0:2) + veloc_loc(0,nglly-1,0,0:2)
                case (4)
                   Tdomain%sVertex(nv)%sSimu(0)%Displ(0:2) = &
                    Tdomain%sVertex(nv)%sSimu(0)%Displ(0:2) + displ_loc(0,0,ngllz-1,0:2)
                   Tdomain%sVertex(nv)%sSimu(0)%Veloc(0:2) = &
                    Tdomain%sVertex(nv)%sSimu(0)%Veloc(0:2) + veloc_loc(0,0,ngllz-1,0:2)
                case (5)
                   Tdomain%sVertex(nv)%sSimu(0)%Displ(0:2) = &
                    Tdomain%sVertex(nv)%sSimu(0)%Displ(0:2) + displ_loc(ngllx-1,0,ngllz-1,0:2)
                   Tdomain%sVertex(nv)%sSimu(0)%Veloc(0:2) = &
                    Tdomain%sVertex(nv)%sSimu(0)%Veloc(0:2) + veloc_loc(ngllx-1,0,ngllz-1,0:2)
                case (6)
                   Tdomain%sVertex(nv)%sSimu(0)%Displ(0:2) = &
                    Tdomain%sVertex(nv)%sSimu(0)%Displ(0:2) + displ_loc(ngllx-1,nglly-1,ngllz-1,0:2)
                   Tdomain%sVertex(nv)%sSimu(0)%Veloc(0:2) = &
                    Tdomain%sVertex(nv)%sSimu(0)%Veloc(0:2) + veloc_loc(ngllx-1,nglly-1,ngllz-1,0:2)
                case (7)
                   Tdomain%sVertex(nv)%sSimu(0)%Displ(0:2) = &
                    Tdomain%sVertex(nv)%sSimu(0)%Displ(0:2) + displ_loc(0,nglly-1,ngllz-1,0:2)
                   Tdomain%sVertex(nv)%sSimu(0)%Veloc(0:2) = &
                    Tdomain%sVertex(nv)%sSimu(0)%Veloc(0:2) + veloc_loc(0,nglly-1,ngllz-1,0:2)
               end select
           endif
       endif
   enddo

   deallocate (displ_loc,veloc_loc)
   if (Tdomain%ellipticity) deallocate (Rsph)

enddo
deallocate (L_Face,L_Edge,L_Vertex)


call deallocate_init_cond()

!--------------------------------------------------------------
end subroutine put_init_cond
!--------------------------------------------------------------
