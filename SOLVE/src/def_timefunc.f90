subroutine def_timefunc (Tdomain,rg)


use sdomains
use angles
use fft_util

implicit none

type (domain), intent(INOUT) :: Tdomain
integer, intent(IN) :: rg

integer :: i,j, nstep,nstep2, isrc, type_src
doubleprecision :: wt,freq, t0,t1,tmax, dt, f1h,f2h,f3h,f4h, rtime
complex*16 :: dphi
complex*16, dimension(:), allocatable :: spectre,tmp

doubleprecision, parameter :: taper = 0.2   ! a cosine-taper is applied between t=0 and t=taper*t0
                                            ! (t0 being the time-shift defined in input.spec) 


nstep = Tdomain%sTimeParam%ntime
dt = Tdomain%sTimeParam%dt

do isrc = 0,Tdomain%n_source-1
 type_src = Tdomain%sSource(isrc)%i_type_source
 
 if ( (type_src==1 .or. type_src==2) .and. rg==Tdomain%sSource(isrc)%proc ) then
  t0 = Tdomain%sSource(isrc)%tau_b
  allocate (Tdomain%sSource(isrc)%timefunc(0:nstep-1))

  if (Tdomain%sSource(isrc)%i_time_function==3) then
     nstep2 = int(2.d0**(int(log(dble(nstep))/log(2.d0))+1))
     allocate(spectre(0:nstep2-1))
     allocate(tmp(1:nstep2))
     f1h = Tdomain%sSource(isrc)%fh(0)
     f2h = Tdomain%sSource(isrc)%fh(1)
     f3h = Tdomain%sSource(isrc)%fh(2)
     f4h = Tdomain%sSource(isrc)%fh(3)
     spectre(:)=cmplx(0.d0,0.d0)
     do j = 1,nstep2
        if (j<=nstep2/2) then
           freq = (j-1)/(dt*nstep2)
        else if (j==nstep2/2+1) then
           freq = 1/(2.d0*dt)
        else
           freq = -(nstep2-j+1)/(dt*nstep2)
        endif
        dphi = exp(-2.d0*pi*freq*t0*cmplx(0.d0,1.d0))
        call wtcoef(abs(freq),f1h,f2h,f3h,f4h,wt)
        if (j/=0)   spectre(j-1) = wt*dphi
     enddo
     do j = 1,nstep2
        tmp(j) = spectre(j-1)
     enddo
     call dfour1(tmp,nstep2,1)
     do j = 1,nstep2
        spectre(j-1) = tmp(j)
     enddo
     Tdomain%sSource(isrc)%timefunc(0:nstep-1) = real(spectre(0:nstep-1))/nstep2/dt
     deallocate(spectre,tmp)
  else
     do i = 0,nstep-1
        rtime = i*dt
        Tdomain%sSource(isrc)%timefunc(i) = CompSource(Tdomain%sSource(isrc),rtime)
     enddo
  endif

  t1 = taper*t0
  tmax = nstep*dt
  do i = 0,nstep-1
     rtime = i*dt
     call wtcoef(rtime,0.d0,t1,tmax,tmax,wt)
     Tdomain%sSource(isrc)%timefunc(i) = Tdomain%sSource(isrc)%timefunc(i) * wt
  enddo

  if (isrc==0) then
     open (60,file="sgn_src",status="unknown",form="formatted")
     do i = 0,nstep-1
        rtime = i*dt
        write (60,*) rtime, Tdomain%sSource(isrc)%timefunc(i)
     enddo
     close (60)
  endif
     
 endif

enddo


end subroutine def_timefunc
