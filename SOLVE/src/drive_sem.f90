program  drive_sem


use sdomains

implicit none

include 'mpif.h'

type(domain), target :: Tdomain

logical :: check_random_t = .true.
integer :: code, rg, nb_procs, ntime, i
doubleprecision :: t
!character*60 :: command, command2


!call suspendinit

call mpi_init(code)
call mpi_comm_rank(mpi_comm_world, rg, code)
call mpi_comm_size(mpi_comm_world, nb_procs, code)

! ##############  Begin of the program  #######################
write (*,*) "Define mesh properties ",rg
call read_input (Tdomain, rg)

write (*,*) "Compute Gauss-Lobatto-Legendre weights and zeroes ",rg
call compute_GLL (Tdomain)

write (*,*) "Define a global numbering for the collocation points ",rg
call global_numbering (Tdomain)

write (*,*) "Compute shape functions within their derivatives ",rg
if (Tdomain%n_nodes==8) then
   call shape8(Tdomain)   ! Linear interpolation
else if (Tdomain%n_nodes==20) then
   write (*,*) "20 control points not yet implemented in the code. Wait for an upgrade"
   stop
else
   call shape27(Tdomain)   ! Quadratic interpolation
endif

call PML_definition (Tdomain) 

call ocean_definition (Tdomain)

call mirror_definition (Tdomain)

write (*,*) "Compute source parameters ",rg
call SourcePosition (Tdomain, rg)
call double_couple (Tdomain, rg)
call pulse (Tdomain, rg)

if (Tdomain%save_trace) then
   write (*,*) "Compute receiver locations ",rg
   call ReceiverPosition (Tdomain, rg)
endif

write (*,*) "Allocate fields ",rg
call allocate_domain (Tdomain, rg)

write (*,*) "Elastic parameters, Mass matrix and Time step ",rg
call define_arrays (Tdomain, rg)

if (Tdomain%n_sls>0) then
   write (*,*) "Compute attenuation features ",rg
   call set_attenuation_param (Tdomain)
endif

do i = 0, Tdomain%n_source-1
   if (Tdomain%sSource(i)%i_type_source==3) then
      write (*,*) "Load initial conditions ",rg
      call put_init_cond (Tdomain)
   endif
   if (Tdomain%sSource(i)%i_type_source==4) then
      if (check_random_t) then
         write (*,*) "Compute random tractions ",rg
         call random_trac (Tdomain, rg, i)
         check_random_t = .false.
      else
         stop 'ONLY ONE SET OF RANDOM SOURCES IS ALLOWED !!!'
      endif
   endif
enddo

if (Tdomain%adjoint) then
   write (*,*) "Compute adjoint source ",rg
   call adj_input (Tdomain)
   call SourceAdjoint (Tdomain, rg)
   call forceAdjoint (Tdomain, rg)
endif

write (*,*) "Entering the time evolution ",rg
Tdomain%sTimeParam%Nsnap = 1
Tdomain%sTimeParam%rtime = 0.d0
do ntime = 0, Tdomain%sTimeParam%ntime-1
   t = mpi_wtime()
   call suspendcheck
   call Newmark (Tdomain, rg, ntime)
!   if (Tdomain%save_snapshot .and. mod(ntime,Tdomain%sTimeParam%dnsnap)==0 .and. ntime/=0) then
   if (Tdomain%save_snapshot .and. &
       ( (Tdomain%t_reversal_mirror<=1 .and. mod(ntime,Tdomain%sTimeParam%dnsnap)==0) .or. &
         (Tdomain%t_reversal_mirror==2 .and. mod(Tdomain%sTimeParam%ntime-ntime-2,Tdomain%sTimeParam%dnsnap)==0) ) ) then
      call MPI_BARRIER(MPI_COMM_WORLD, code)
      if (rg==0) then
         call merge_snapshot_files(nb_procs,Tdomain%sTimeParam%Nsnap,'snapshot_forward_')
         call merge_snapshot_files(nb_procs,Tdomain%sTimeParam%Nsnap,'snapshot_backward_')
         call merge_snapshot_files(nb_procs,Tdomain%sTimeParam%Nsnap,'snapshot_adjoint_')
!         !!! CAUTION: This can make the run unstable !!!
!         write (command,'(a,I3.3,a,I3.3)') "cat snapshot",Tdomain%sTimeParam%Nsnap,"_* > shot",Tdomain%sTimeParam%Nsnap
!!         call system(command)
!!         call system("rm snapshot*")
!         if (Tdomain%adjoint) then
!            write (command2,'(a,I3.3,a,I3.3)') "cat adjsnap_",Tdomain%sTimeParam%Nsnap,"_* > shot_adj_",Tdomain%sTimeParam%Nsnap
!!            call system(command2)
!!            call system("rm adjsnap_*")
!         endif
      endif
      Tdomain%sTimeParam%Nsnap = Tdomain%sTimeParam%Nsnap + 1
      if (Tdomain%sTimeParam%Nsnap==1000) then
         print *,"WARNING : TOO MANY SNAPSHOTS"
         Tdomain%save_snapshot = .false.
      endif
   endif
   Tdomain%sTimeParam%rtime = Tdomain%sTimeParam%rtime + Tdomain%sTimeParam%dt
   t = mpi_wtime() - t
   if (rg==0)   write (*,*) ntime, t
enddo
write (*,*) "Simulation is finished ",rg

!if (Tdomain%t_reversal_mirror==1) then
!   write (*,*) "One more prediction for the time-reversal mirror ",rg
!   call one_more_predic (Tdomain, rg, Tdomain%sTimeParam%ntime)
!endif

write (*,*) "Deallocate fields ",rg
call deallocate_domain (Tdomain, rg)

call mpi_finalize(code)
write (*,*) "END ",rg


end program drive_sem
