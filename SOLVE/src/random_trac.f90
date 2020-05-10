subroutine random_trac (Tdomain, rg, isrc)


use sdomains
use angles

implicit none

include 'mpif.h'

integer, intent(IN) :: rg, isrc
type (domain), intent(INOUT) :: Tdomain

logical :: pouet, return_sgn = .true.
integer :: nstep,nstep2, ios, n, proc, elem, x,y,z, ngllx,nglly,ngllz, comp, i,j, ne,nv, ngll, &
           envoyeur,receveur, code
integer, parameter :: etiquette = 100
integer, dimension(mpi_status_size) :: statut
integer, dimension(:), allocatable :: L_Edge, L_Vertex
doubleprecision :: lat,long, duree,duree2, wt, f1h,f2h,f3h,f4h, freq, taper_debut,taper_fin
doubleprecision, dimension(0:2) :: Amp, normale
doubleprecision, dimension(:), allocatable :: spectre
doubleprecision, dimension(:,:), allocatable :: Ffilt, tmp
character(len=20) :: filename

doubleprecision :: dt = 1.d0   ! Pas de temps des signaux aleatoires
taper_debut = Tdomain%sSource(isrc)%tau_b   ! Longueur du signal qui sera attenuee au debut
taper_fin = 0.d0   ! Longueur du signal qui sera attenuee a la fin


if (Tdomain%curve .eqv. .false.) then
   write (*,*) "Random tractions can't be used in a cuboid for now."
   stop
endif

! Allocation
duree = Tdomain%sTimeParam%duration
nstep = int(duree/dt) + 1
dt = Tdomain%sTimeParam%duration/(nstep-1)
Tdomain%sTimeParam%dt_rand_F = dt
allocate (Ffilt(0:nstep-1,0:2))
allocate (tmp(0:nstep-1,0:2))

nstep2 = int(2.d0**(int(log(dble(nstep))/log(2.d0))+1))
duree2 = (nstep2-1)*dt
allocate (spectre(1:nstep2))

allocate (L_Edge(0:Tdomain%n_edge-1));   L_Edge = 0
allocate (L_Vertex(0:Tdomain%n_vertex-1));   L_Vertex = 0

! Amplitude spectrale du filtre
f1h = Tdomain%sSource(isrc)%fh(0)
f2h = Tdomain%sSource(isrc)%fh(1)
f3h = Tdomain%sSource(isrc)%fh(2)
f4h = Tdomain%sSource(isrc)%fh(3)
do i = 1,nstep2
   if (i<=nstep2/2) then
      freq = (i-1)/duree2
   else if (i==nstep2/2+1) then
      freq = 1/(2.d0*dt)
   else
      freq = -(nstep2-i+1)/duree2
   endif
   call wtcoef(abs(freq),f1h,f2h,f3h,f4h,wt)
   spectre(i) = wt
enddo

! Generation des signaux aleatoires et Comm intraproc
call random_seed
open (15,file='random_trac.table',iostat=ios)
if (ios/=0)   stop 'CANNOT OPEN random_trac.table'
read(15,*) Amp(0:2)
read(15,*) proc, elem
do while (ios==0)
   read(15,*) lat, long, pouet

   if (rg==proc) then
      Tdomain%specel(elem)%random_t = pouet
      ngllx = Tdomain%specel(elem)%ngllx
      nglly = Tdomain%specel(elem)%nglly
      ngllz = Tdomain%specel(elem)%ngllz
      if (Tdomain%specel(elem)%random_t) then
         allocate (Tdomain%specel(elem)%random_coeff(0:ngllx-1, 0:nglly-1, 0:nstep-1, 0:2))
         z = ngllz-1
         do y = 1,nglly-2
          do x = 1,ngllx-2
             ! Generation du signal aleatoire
             call comp_rand_t (Amp, dt, taper_debut, taper_fin, duree, nstep2, spectre, nstep, Ffilt)
             ! Expression du signal dans le repere cartesien
             n = Tdomain%specel(elem)%Iglobnum(x,y,z)
             call vect_sph2cart (Tdomain%Globcoord(:,n), nstep, Ffilt(:,:), tmp(:,:))
             Tdomain%specel(elem)%random_coeff(x,y,:,:) = tmp(:,:)
          enddo
         enddo
      endif
   endif

   if (rg==proc) then

      ne = Tdomain%specel(elem)%Near_Edges(5)
      if (L_Edge(ne)==0) then
         allocate (Tdomain%sEdge(ne)%random_coeff(1:ngllx-2, 0:nstep-1, 0:2))
         z = ngllz-1; y = 0
         do x = 1,ngllx-2 
            ! Generation du signal aleatoire
            call comp_rand_t (Amp, dt, taper_debut, taper_fin, duree, nstep2, spectre, nstep, Ffilt)
            ! Expression du signal dans le repere cartesien
            n = Tdomain%specel(elem)%Iglobnum(x,y,z)
            call vect_sph2cart (Tdomain%Globcoord(:,n), nstep, Ffilt(:,:), tmp(:,:))
            Tdomain%sEdge(ne)%random_coeff(x,:,:) = tmp(:,:)
            if (Tdomain%specel(elem)%random_t) then   ! L'edge donne a l'element
               Tdomain%specel(elem)%random_coeff(x,y,:,:) = Tdomain%sEdge(ne)%random_coeff(x,:,:)
            endif
         enddo
         L_Edge(ne) = L_Edge(ne) + 1
      else
         z = ngllz-1; y = 0
         do x = 1,ngllx-2
            if (Tdomain%specel(elem)%random_t) then   ! L'edge donne a l'element
               Tdomain%specel(elem)%random_coeff(x,y,:,:) = Tdomain%sEdge(ne)%random_coeff(x,:,:)
            endif
         enddo
         L_Edge(ne) = 0;   deallocate (Tdomain%sEdge(ne)%random_coeff)
      endif

      ne = Tdomain%specel(elem)%Near_Edges(8)
      if (L_Edge(ne)==0) then
         allocate (Tdomain%sEdge(ne)%random_coeff(1:nglly-2, 0:nstep-1, 0:2))
         z = ngllz-1; x = ngllx-1
         do y = 1,nglly-2
            ! Generation du signal aleatoire
            call comp_rand_t (Amp, dt, taper_debut, taper_fin, duree, nstep2, spectre, nstep, Ffilt)
            ! Expression du signal dans le repere cartesien
            n = Tdomain%specel(elem)%Iglobnum(x,y,z)
            call vect_sph2cart (Tdomain%Globcoord(:,n), nstep, Ffilt(:,:), tmp(:,:))
            Tdomain%sEdge(ne)%random_coeff(y,:,:) = tmp(:,:)
            if (Tdomain%specel(elem)%random_t) then   ! L'edge donne a l'element
               Tdomain%specel(elem)%random_coeff(x,y,:,:) = Tdomain%sEdge(ne)%random_coeff(y,:,:)
            endif
         enddo
         L_Edge(ne) = L_Edge(ne) + 1
      else
         z = ngllz-1; x = ngllx-1
         do y = 1,nglly-2
            if (Tdomain%specel(elem)%random_t) then   ! L'edge donne a l'element
               Tdomain%specel(elem)%random_coeff(x,y,:,:) = Tdomain%sEdge(ne)%random_coeff(y,:,:)
            endif
         enddo
         L_Edge(ne) = 0;   deallocate (Tdomain%sEdge(ne)%random_coeff)
      endif

      ne = Tdomain%specel(elem)%Near_Edges(9)
      if (L_Edge(ne)==0) then
         allocate (Tdomain%sEdge(ne)%random_coeff(1:ngllx-2, 0:nstep-1, 0:2))
         z = ngllz-1; y = nglly-1
         do x = 1,ngllx-2
            ! Generation du signal aleatoire
            call comp_rand_t (Amp, dt, taper_debut, taper_fin, duree, nstep2, spectre, nstep, Ffilt)
            ! Expression du signal dans le repere cartesien
            n = Tdomain%specel(elem)%Iglobnum(x,y,z)
            call vect_sph2cart (Tdomain%Globcoord(:,n), nstep, Ffilt(:,:), tmp(:,:))
            Tdomain%sEdge(ne)%random_coeff(x,:,:) = tmp(:,:)
            if (Tdomain%specel(elem)%random_t) then   ! L'edge donne a l'element
               Tdomain%specel(elem)%random_coeff(x,y,:,:) = Tdomain%sEdge(ne)%random_coeff(x,:,:)
            endif
         enddo
         L_Edge(ne) = L_Edge(ne) + 1
      else
         z = ngllz-1; y = nglly-1
         do x = 1,ngllx-2
            if (Tdomain%specel(elem)%random_t) then   ! L'edge donne a l'element
               Tdomain%specel(elem)%random_coeff(x,y,:,:) = Tdomain%sEdge(ne)%random_coeff(x,:,:)
            endif
         enddo
         L_Edge(ne) = 0;   deallocate (Tdomain%sEdge(ne)%random_coeff)
      endif

      ne = Tdomain%specel(elem)%Near_Edges(11)
      if (L_Edge(ne)==0) then
         allocate (Tdomain%sEdge(ne)%random_coeff(1:nglly-2, 0:nstep-1, 0:2))
         z = ngllz-1; x = 0
         do y = 1,nglly-2
            ! Generation du signal aleatoire
            call comp_rand_t (Amp, dt, taper_debut, taper_fin, duree, nstep2, spectre, nstep, Ffilt)
            ! Expression du signal dans le repere cartesien
            n = Tdomain%specel(elem)%Iglobnum(x,y,z)
            call vect_sph2cart (Tdomain%Globcoord(:,n), nstep, Ffilt(:,:), tmp(:,:))
            Tdomain%sEdge(ne)%random_coeff(y,:,:) = tmp(:,:)
            if (Tdomain%specel(elem)%random_t) then   ! L'edge donne a l'element
               Tdomain%specel(elem)%random_coeff(x,y,:,:) = Tdomain%sEdge(ne)%random_coeff(y,:,:)
            endif
         enddo
         L_Edge(ne) = L_Edge(ne) + 1
      else
         z = ngllz-1; x = 0
         do y = 1,nglly-2
            if (Tdomain%specel(elem)%random_t) then   ! L'edge donne a l'element
               Tdomain%specel(elem)%random_coeff(x,y,:,:) = Tdomain%sEdge(ne)%random_coeff(y,:,:)
            endif
         enddo
         L_Edge(ne) = 0;   deallocate (Tdomain%sEdge(ne)%random_coeff)
      endif

      nv = Tdomain%specel(elem)%Near_Vertices(4)
      if (L_Vertex(nv)==0) then
         allocate (Tdomain%sVertex(nv)%random_coeff(0:nstep-1, 0:2))
         z = ngllz-1; y = 0; x = 0
         ! Generation du signal aleatoire
         call comp_rand_t (Amp, dt, taper_debut, taper_fin, duree, nstep2, spectre, nstep, Ffilt)
         ! Expression du signal dans le repere cartesien
         n = Tdomain%specel(elem)%Iglobnum(x,y,z)
         call vect_sph2cart (Tdomain%Globcoord(:,n), nstep, Ffilt(:,:), tmp(:,:))
         Tdomain%sVertex(nv)%random_coeff(:,:) = tmp(:,:)
         if (Tdomain%specel(elem)%random_t) then   ! Le vertex donne a l'element
            Tdomain%specel(elem)%random_coeff(x,y,:,:) = Tdomain%sVertex(nv)%random_coeff(:,:)
         endif
         L_Vertex(nv) = L_Vertex(nv) + 1
      else
         z = ngllz-1; y = 0; x = 0
         if (Tdomain%specel(elem)%random_t) then   ! Le vertex donne a l'element
            Tdomain%specel(elem)%random_coeff(x,y,:,:) = Tdomain%sVertex(nv)%random_coeff(:,:)
         endif
         L_Vertex(nv) = L_Vertex(nv) + 1
         if (L_Vertex(nv)==4) then
            L_Vertex(nv) = 0;   deallocate (Tdomain%sVertex(nv)%random_coeff)
         endif
      endif

      nv = Tdomain%specel(elem)%Near_Vertices(5)
      if (L_Vertex(nv)==0) then
         allocate (Tdomain%sVertex(nv)%random_coeff(0:nstep-1, 0:2))
         z = ngllz-1; y = 0; x = ngllx-1
         ! Generation du signal aleatoire
         call comp_rand_t (Amp, dt, taper_debut, taper_fin, duree, nstep2, spectre, nstep, Ffilt)
         ! Expression du signal dans le repere cartesien
         n = Tdomain%specel(elem)%Iglobnum(x,y,z)
         call vect_sph2cart (Tdomain%Globcoord(:,n), nstep, Ffilt(:,:), tmp(:,:))
         Tdomain%sVertex(nv)%random_coeff(:,:) = tmp(:,:)
         if (Tdomain%specel(elem)%random_t) then   ! Le vertex donne a l'element
            Tdomain%specel(elem)%random_coeff(x,y,:,:) = Tdomain%sVertex(nv)%random_coeff(:,:)
         endif
         L_Vertex(nv) = L_Vertex(nv) + 1
      else
         z = ngllz-1; y = 0; x = ngllx-1
         if (Tdomain%specel(elem)%random_t) then   ! Le vertex donne a l'element
            Tdomain%specel(elem)%random_coeff(x,y,:,:) = Tdomain%sVertex(nv)%random_coeff(:,:)
         endif
         L_Vertex(nv) = L_Vertex(nv) + 1
         if (L_Vertex(nv)==4) then
            L_Vertex(nv) = 0;   deallocate (Tdomain%sVertex(nv)%random_coeff)
         endif
      endif

      nv = Tdomain%specel(elem)%Near_Vertices(6)
      if (L_Vertex(nv)==0) then
         allocate (Tdomain%sVertex(nv)%random_coeff(0:nstep-1, 0:2))
         z = ngllz-1; y = nglly-1; x = ngllx-1
         ! Generation du signal aleatoire
         call comp_rand_t (Amp, dt, taper_debut, taper_fin, duree, nstep2, spectre, nstep, Ffilt)
         ! Expression du signal dans le repere cartesien
         n = Tdomain%specel(elem)%Iglobnum(x,y,z)
         call vect_sph2cart (Tdomain%Globcoord(:,n), nstep, Ffilt(:,:), tmp(:,:))
         Tdomain%sVertex(nv)%random_coeff(:,:) = tmp(:,:)
         if (Tdomain%specel(elem)%random_t) then   ! Le vertex donne a l'element
            Tdomain%specel(elem)%random_coeff(x,y,:,:) = Tdomain%sVertex(nv)%random_coeff(:,:)
         endif
         L_Vertex(nv) = L_Vertex(nv) + 1
      else
         z = ngllz-1; y = nglly-1; x = ngllx-1
         if (Tdomain%specel(elem)%random_t) then   ! Le vertex donne a l'element
            Tdomain%specel(elem)%random_coeff(x,y,:,:) = Tdomain%sVertex(nv)%random_coeff(:,:)
         endif
         L_Vertex(nv) = L_Vertex(nv) + 1
         if (L_Vertex(nv)==4) then
            L_Vertex(nv) = 0;   deallocate (Tdomain%sVertex(nv)%random_coeff)
         endif
      endif

      nv = Tdomain%specel(elem)%Near_Vertices(7)
      if (L_Vertex(nv)==0) then
         allocate (Tdomain%sVertex(nv)%random_coeff(0:nstep-1, 0:2))
         z = ngllz-1; y = nglly-1; x = 0
         ! Generation du signal aleatoire
         call comp_rand_t (Amp, dt, taper_debut, taper_fin, duree, nstep2, spectre, nstep, Ffilt)
         ! Expression du signal dans le repere cartesien
         n = Tdomain%specel(elem)%Iglobnum(x,y,z)
         call vect_sph2cart (Tdomain%Globcoord(:,n), nstep, Ffilt(:,:), tmp(:,:))
         Tdomain%sVertex(nv)%random_coeff(:,:) = tmp(:,:)
         if (Tdomain%specel(elem)%random_t) then   ! Le vertex donne a l'element
            Tdomain%specel(elem)%random_coeff(x,y,:,:) = Tdomain%sVertex(nv)%random_coeff(:,:)
         endif
         L_Vertex(nv) = L_Vertex(nv) + 1
      else
         z = ngllz-1; y = nglly-1; x = 0
         if (Tdomain%specel(elem)%random_t) then   ! Le vertex donne a l'element
            Tdomain%specel(elem)%random_coeff(x,y,:,:) = Tdomain%sVertex(nv)%random_coeff(:,:)
         endif
         L_Vertex(nv) = L_Vertex(nv) + 1
         if (L_Vertex(nv)==4) then
            L_Vertex(nv) = 0;   deallocate (Tdomain%sVertex(nv)%random_coeff)
         endif
         if (return_sgn .and. Tdomain%specel(elem)%random_t) then
            do comp = 0,2
               if (rg<10) then
                write (filename(1:10), '(a7,i1,a1,i1)') "random_", rg, "_", comp
               elseif (rg<100) then
                write (filename(1:11), '(a7,i2,a1,i1)') "random_", rg, "_", comp
               elseif (rg<1000) then
                write (filename(1:12), '(a7,i3,a1,i1)') "random_", rg, "_", comp
               else
                write (filename(1:13), '(a7,i4,a1,i1)') "random_", rg, "_", comp
               endif
               open (16,file=trim(filename))
               do n = 0,nstep-1
                  write(16,*) n*dt, Tdomain%specel(elem)%random_coeff(x,y,n,comp)
               enddo
               close (16)
            enddo
            return_sgn = .false.
         endif
      endif

   endif

   read(15,fmt=*,iostat=ios) proc, elem
enddo
close (15)

! Comm interproc (MPI)
do n = 0,Tdomain%n_proc-1
   ngll = 0
   do i = 0,Tdomain%sComm(n)%nb_edges-1
      ne = Tdomain%sComm(n)%edges(i)
      if (L_Edge(ne)>0) then
         ngll = ngll + Tdomain%sEdge(ne)%ngll-2
      endif
   enddo
   do i = 0,Tdomain%sComm(n)%nb_vertices-1
      nv = Tdomain%sComm(n)%vertices(i)
      if (L_Vertex(nv)>0) then
         ngll = ngll + 1
      endif
   enddo
   Tdomain%sComm(n)%ngllrandt = ngll
   if (ngll>0)   allocate (Tdomain%sComm(n)%RandT(0:ngll-1,0:nstep-1,0:2))
enddo
do n = rg+1,Tdomain%n_proc-1
   ngll = 0
   do i = 0,Tdomain%sComm(n)%nb_edges-1
      ne = Tdomain%sComm(n)%edges(i)
      if (L_Edge(ne)>0) then
         do j = 1,Tdomain%sEdge(ne)%ngll-2
            Tdomain%sComm(n)%RandT(ngll,:,:) = Tdomain%sEdge(ne)%random_coeff(j,:,:)
            ngll = ngll + 1
         enddo
      endif
   enddo
   do i = 0,Tdomain%sComm(n)%nb_vertices-1
      nv = Tdomain%sComm(n)%vertices(i)
      if (L_Vertex(nv)>0) then
         Tdomain%sComm(n)%RandT(ngll,:,:) = Tdomain%sVertex(nv)%random_coeff(:,:)
         ngll = ngll + 1
      endif
   enddo
enddo
do envoyeur = 0,Tdomain%n_proc-2
   do receveur = envoyeur+1,Tdomain%n_proc-1
      if (rg==envoyeur .and. Tdomain%sComm(receveur)%ngllrandt>0) then
         call MPI_SEND (Tdomain%sComm(receveur)%RandT, 3*nstep*Tdomain%sComm(receveur)%ngllrandt, &
                        MPI_DOUBLE_PRECISION, receveur, etiquette, MPI_COMM_WORLD, code)
      endif
      if (rg==receveur .and. Tdomain%sComm(envoyeur)%ngllrandt>0) then
         call MPI_RECV (Tdomain%sComm(envoyeur)%RandT, 3*nstep*Tdomain%sComm(envoyeur)%ngllrandt, &
                        MPI_DOUBLE_PRECISION, envoyeur, etiquette, MPI_COMM_WORLD, statut, code)
      endif
   enddo
enddo
do n = 0,rg-1
   ngll = 0
   do i = 0,Tdomain%sComm(n)%nb_edges-1
      ne = Tdomain%sComm(n)%edges(i)
      if (L_Edge(ne)>0) then
         do j = 1,Tdomain%sEdge(ne)%ngll-2
            Tdomain%sEdge(ne)%random_coeff(j,:,:) = Tdomain%sComm(n)%RandT(ngll,:,:)
            ngll = ngll + 1
         enddo
      endif
   enddo
   do i = 0,Tdomain%sComm(n)%nb_vertices-1
      nv = Tdomain%sComm(n)%vertices(i)
      if (L_Vertex(nv)>0) then
         Tdomain%sVertex(nv)%random_coeff(:,:) = Tdomain%sComm(n)%RandT(ngll,:,:)
         ngll = ngll + 1
      endif
   enddo
enddo
do n = 0,Tdomain%n_proc-1
   if (Tdomain%sComm(n)%ngllrandt>0)   deallocate(Tdomain%sComm(n)%RandT)
enddo

do n = 0,Tdomain%n_elem-1
   if (Tdomain%specel(n)%random_t) then
      ngllx = Tdomain%specel(n)%ngllx
      nglly = Tdomain%specel(n)%nglly
      ngllz = Tdomain%specel(n)%ngllz
      ne = Tdomain%specel(n)%Near_Edges(5)
      if (L_Edge(ne)>0) then
         z = ngllz-1; y = 0
         do x = 1,ngllx-2
            Tdomain%specel(n)%random_coeff(x,y,:,:) = Tdomain%sEdge(ne)%random_coeff(x,:,:)
         enddo
      endif
      ne = Tdomain%specel(n)%Near_Edges(8)
      if (L_Edge(ne)>0) then
         z = ngllz-1; x = ngllx-1
         do y = 1,nglly-2
            Tdomain%specel(n)%random_coeff(x,y,:,:) = Tdomain%sEdge(ne)%random_coeff(y,:,:)
         enddo
      endif
      ne = Tdomain%specel(n)%Near_Edges(9)
      if (L_Edge(ne)>0) then
         z = ngllz-1; y = nglly-1
         do x = 1,ngllx-2
            Tdomain%specel(n)%random_coeff(x,y,:,:) = Tdomain%sEdge(ne)%random_coeff(x,:,:)
         enddo
      endif
      ne = Tdomain%specel(n)%Near_Edges(11)
      if (L_Edge(ne)>0) then
         z = ngllz-1; x = 0
         do y = 1,nglly-2
            Tdomain%specel(n)%random_coeff(x,y,:,:) = Tdomain%sEdge(ne)%random_coeff(y,:,:)
         enddo
      endif
      nv = Tdomain%specel(n)%Near_Vertices(4)
      if (L_Vertex(nv)>0) then
         z = ngllz-1; y = 0; x = 0
         Tdomain%specel(n)%random_coeff(x,y,:,:) = Tdomain%sVertex(nv)%random_coeff(:,:)
      endif
      nv = Tdomain%specel(n)%Near_Vertices(5)
      if (L_Vertex(nv)>0) then
         z = ngllz-1; y = 0; x = ngllx-1
         Tdomain%specel(n)%random_coeff(x,y,:,:) = Tdomain%sVertex(nv)%random_coeff(:,:)
      endif
      nv = Tdomain%specel(n)%Near_Vertices(6)
      if (L_Vertex(nv)>0) then
         z = ngllz-1; y = nglly-1; x = ngllx-1
         Tdomain%specel(n)%random_coeff(x,y,:,:) = Tdomain%sVertex(nv)%random_coeff(:,:)
      endif
      nv = Tdomain%specel(n)%Near_Vertices(7)
      if (L_Vertex(nv)>0) then
         z = ngllz-1; y = nglly-1; x = 0
         Tdomain%specel(n)%random_coeff(x,y,:,:) = Tdomain%sVertex(nv)%random_coeff(:,:)
      endif
      ! Multiplication par le Jacobien surfacique et les poids d'integration
      if (Tdomain%specel(n)%PML) then
         print *,"Tractions on PML elements are not allowed"
         stop
      else 
         z = ngllz-1
         do y = 0,nglly-1
          do x = 0,ngllx-1
             normale(0:2) = Tdomain%specel(n)%InvGrad(x,y,z,0:2,2)
             Tdomain%specel(n)%random_coeff(x,y,:,:) = Tdomain%specel(n)%random_coeff(x,y,:,:) * &
                                                       Tdomain%specel(n)%Jacob(x,y,z) * &
                                                       dsqrt(normale(0)**2 + normale(1)**2 + normale(2)**2) * &
                                                       Tdomain%specel(n)%wgtx(x)*Tdomain%specel(n)%wgty(y)
          enddo
         enddo
      endif
   endif
enddo

! Desallocation
do ne = 0,Tdomain%n_edge-1
   if (L_Edge(ne)>0)   deallocate (Tdomain%sEdge(ne)%random_coeff)
enddo
do nv = 0,Tdomain%n_vertex-1
   if (L_Vertex(nv)>0)   deallocate (Tdomain%sVertex(nv)%random_coeff)
enddo
deallocate (spectre,Ffilt,L_Edge,L_Vertex)


end subroutine random_trac
! ##############################################################################

! ##############################################################################
subroutine comp_rand_t (Amp, dt, taper_debut, taper_fin, duree, nstep2, spectre, nstep, Ffilt)


use fft_util
use angles

implicit none

integer, intent(IN) :: nstep, nstep2
doubleprecision, intent(IN) ::dt, taper_debut, taper_fin, duree
doubleprecision, dimension(0:2), intent(IN) :: Amp
doubleprecision, dimension(1:nstep2), intent(IN) :: spectre
doubleprecision, dimension(0:nstep-1,0:2), intent(INOUT) :: Ffilt

integer :: n, comp
doubleprecision :: hasard, wt, t0,t1,t2,t3, xi,yi, tanphi, epsil
complex*16, dimension(:,:), allocatable :: F


! Definition du zero des reels
epsil = tiny(hasard)

! Allocation
allocate (F(1:nstep2,0:2))

! Tirage des signaux aleatoires pour le GLL considere
do comp = 0,2
    if (Amp(comp)>epsil) then
        do n = 1,nstep2
            call random_number(hasard)
            hasard = hasard - 0.5
            F(n,comp) = cmplx(hasard,0.d0)
        enddo
    else
        F(:,comp) = cmplx(0.d0,0.d0)
    endif
enddo

! Filtrage
do comp = 0,2
    if (Amp(comp)>epsil) then
        call dfour1(F(:,comp),nstep2,-1)
        do n = 1,nstep2
            tanphi = aimag(F(n,comp))/dble(F(n,comp))
            xi = spectre(n)/dsqrt(1+tanphi**2)
            if (dble(F(n,comp))<0) xi = -xi
            yi = xi*tanphi
            F(n,comp) = cmplx(xi,yi)
        enddo
        call dfour1(F(:,comp),nstep2,1)
    endif
enddo

! Taper (on attenue le debut et la fin du signal genere aleatoirement)
t0 = 0; t1 = taper_debut; t2 = duree-taper_fin; t3 = duree
do n = 0,nstep-1
    call wtcoef(n*dt,t0,t1,t2,t3,wt)
    Ffilt(n,:) = Amp(:) * wt * dble(F(n+1,:))
enddo

! Desallocation
deallocate (F)


end subroutine comp_rand_t
