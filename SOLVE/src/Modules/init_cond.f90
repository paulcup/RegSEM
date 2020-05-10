!--------------------------------------------------------------------
module init_cond 
!--------------------------------------------------------------------
  implicit none
  public:: init_init_cond,deallocate_init_cond,get_init_cond
  private
!
  character(len=30) :: datafile='init_cond_dep'
  integer :: NBSP,NBST,NBD
  doubleprecision :: phi_deb,phi_end,theta_deb,theta_end,depth_bottom,depth_top,dr,dp,dt
  doubleprecision, dimension(:,:,:,:), allocatable :: map_init,map_initv,y2depth,y2depthv
  doubleprecision, dimension(:), allocatable :: theta,phi,depth
!
  logical :: init_done=.false.

!---------------------------------
contains
!---------------------------------
!---------------------------------
  subroutine init_init_cond()
!---------------------------------
    use module_spline
    implicit none
    integer :: i,j,ic,ird,ii
    doubleprecision :: yp1,ypn
    doubleprecision :: a,b
!
    if (init_done) STOP 'init_init_cond has already been called ???'
!
    open(17,file=datafile,status='old',form='unformatted')
    read(17)  NBSP,NBST,NBD,phi_deb,phi_end,theta_deb,theta_end,depth_bottom,depth_top
    allocate(map_init(NBD,NBSP,NBST,3),map_initv(NBD,NBSP,NBST,3))
    do i=1,NBSP
       do ird=1,NBD       
          do j=1,NBST
             do ic=1,3
                read(17) map_init (ird,i,j,ic),map_initv(ird,i,j,ic)
             enddo
          enddo
       enddo
    enddo
    close(17)
    allocate(theta(NBST),phi(NBSP),depth(NBD))
    dt=(theta_end-theta_deb)/(NBST-1)
    do i=1,NBST
       theta(i)=theta_deb+(i-1)*dt
    enddo
    dp=(phi_end-phi_deb)/(NBSP-1)
    do i=1,NBSP
       phi(i)=phi_deb+(i-1)*dp
    enddo
    dr=(depth_top-depth_bottom)/(NBD-1)
    do i=1,NBD
       depth(i)=depth_bottom+(i-1)*dr
    enddo  
    allocate(y2depth(NBD,NBSP,NBST,3),y2depthv(NBD,NBSP,NBST,3))
    do ic=1,3
       do j=1,NBST
          do i=1,NBSP
             yp1=(map_init(2,i,j,ic)  -map_init(1,i,j,ic)    )/dr
             ypn=(map_init(NBD,i,j,ic)-map_init(NBD-1,i,j,ic))/dr
             call spline(depth,map_init(:,i,j,ic),yp1,ypn,y2depth(:,i,j,ic))
!
             yp1=(map_initv(2,i,j,ic)  -map_initv(1,i,j,ic)    )/dr
             ypn=(map_initv(NBD,i,j,ic)-map_initv(NBD-1,i,j,ic))/dr
             call spline(depth,map_initv(:,i,j,ic),yp1,ypn,y2depthv(:,i,j,ic))
          enddo
       enddo
    enddo
          
    init_done=.true.
!---------------------------------
  end subroutine init_init_cond
!---------------------------------
!---------------------------------
  subroutine get_init_cond(r,t,p,U,V)
!---------------------------------
    use module_spline
    implicit none
    doubleprecision, intent(inout) :: r,t,p
    doubleprecision, dimension(3), intent(out):: U,V
!
    logical :: out_of_range
    integer :: i,j,itd,itf,ipd,ipf,ic,it,ip
    integer, parameter :: N=5,N2=2*N,N2p2=N2+2
    doubleprecision :: yp1,ypn
    doubleprecision, dimension(:,:), allocatable :: tmpu,tmpv 
    doubleprecision, dimension(:), allocatable :: tmp2u,tmp2v 
    doubleprecision, dimension(:), allocatable :: y2u,y2v
!
    if (.not.init_done) STOP 'get_int_cond: init_init_cond must be called first!'
!
    out_of_range = .false.
    if (t<theta_deb .or. t > theta_end) then
!       print *,t,theta_deb,theta_end
!       STOP 'get_init_cond: theta is out of range'
       out_of_range = .true.
    endif
    if (p<phi_deb .or. p > phi_end) then
!       print *,p,phi_deb,phi_end
!       STOP 'get_init_cond: phi is out of range'
       out_of_range = .true.
    endif
    if (r<depth_bottom) then
!       print *,r,depth_bottom,depth_top
!       STOP 'get_init_cond: r is out of range'
       out_of_range = .true.
    endif
    if (r > depth_top) then
       r = depth_top
    endif
!
if (out_of_range) then

    do ic=1,3
        U(ic) = 0
        V(ic) = 0
    enddo

else

    allocate(tmp2u(N2p2),tmp2v(N2p2),y2u(N2),y2v(N2),tmpu(N2p2,N2p2),tmpv(N2p2,N2p2))
    ip=locate(phi  ,NBSP,p)-1
    it=locate(theta,NBST,t)-1
!
    itd=it-N-1
    itf=it+N
    ipd=ip-N-1
    ipf=ip+N
    if (itd<1) then
       itf=itf-itd+1
       itd=1
    endif
    if (ipd<1) then
       ipf=ipf-ipd+1
       ipd=1
    endif
    if (itf>NBST) then
       itd=itd-(itf-NBST)
       itf=NBST
    endif
    if (ipf>NBSP) then
       ipd=ipd-(ipf-NBSP)
       ipf=NBSP
    endif
!
    do ic=1,3
       do j=1,N2p2
          it=itd+j-1
          do i=1,N2p2
             ip=ipd+i-1
             tmpu(i,j)=splint(depth,map_init (:,ip,it,ic),y2depth (:,ip,it,ic),r)
             tmpv(i,j)=splint(depth,map_initv(:,ip,it,ic),y2depthv(:,ip,it,ic),r)
          enddo
       enddo
       do j=1,N2p2
          yp1=(tmpu(3,j)   -tmpu(1,j))     /(2.*dp)
          ypn=(tmpu(N2p2,j)-tmpu(N2p2-2,j))/(2.*dp)
          call spline(phi(ipd+1:ipf-1),tmpu(2:N2p2-1,j),yp1,ypn,y2u)
          tmp2u(j)=splint(phi(ipd+1:ipf-1),tmpu(2:N2p2-1,j),y2u,p)
!
          yp1=(tmpv(3,j)   -tmpv(1,j))     /(2.*dp)
          ypn=(tmpv(N2p2,j)-tmpv(N2p2-2,j))/(2.*dp)
          call spline(phi(ipd+1:ipf-1),tmpv(2:N2p2-1,j),yp1,ypn,y2v)
          tmp2v(j)=splint(phi(ipd+1:ipf-1),tmpv(2:N2p2-1,j),y2v,p)
       enddo
       yp1=(tmp2u(3)   -tmp2u(1))     /(2.*dt)
       ypn=(tmp2u(N2p2)-tmp2u(N2p2-2))/(2.*dt)
       call spline(theta(itd+1:itf-1),tmp2u(2:N2p2-1),yp1,ypn,y2u)
       U(ic)=splint(theta(itd+1:itf-1),tmp2u(2:N2p2-1),y2u,t)
       yp1=(tmp2v(3)   -tmp2v(1))     /(2.*dt)
       ypn=(tmp2v(N2p2)-tmp2v(N2p2-2))/(2.*dt)
       call spline(theta(itd+1:itf-1),tmp2v(2:N2p2-1),yp1,ypn,y2v)
       V(ic)=splint(theta(itd+1:itf-1),tmp2v(2:N2p2-1),y2v,t)
    enddo

endif
!---------------------------------
  end subroutine get_init_cond
!---------------------------------
!---------------------------------
  subroutine deallocate_init_cond
!---------------------------------
    implicit none
    deallocate(map_init,map_initv,theta,phi,depth)
!---------------------------------
  end subroutine deallocate_init_cond
!---------------------------------
!-----------------------------------------------------------------
      integer function locate(inter,n,x)
!-----------------------------------------------------------------
         implicit none
         integer, intent(in) :: n
         doubleprecision, dimension(1:n), intent(in) :: inter
         doubleprecision, intent(in) :: x
!
         integer :: i
         i=1
         if (abs (x-inter(1))/max(abs(x),abs(inter(1))) <1.d-10) then
            i=1
         else if ( abs (x-inter(n))/max(abs(x),abs(inter(n)))  <1.d-10) then
            i=n-1
         else
            do while (  x < inter(i) .or. x > inter(i+1) )
               i=i+1
               if (i > n) then
                  print*,'i=',i,'n=',n
                  print*,'x=',x
                  print*,'inter=',inter
                  stop 'locate: failed to locate x in inter!'
               endif
            enddo
         endif
         locate=i+1
!-----------------------------------------------------------------
      end function locate
!-----------------------------------------------------------------

!--------------------------------------------------------------------
end module init_cond
!--------------------------------------------------------------------
