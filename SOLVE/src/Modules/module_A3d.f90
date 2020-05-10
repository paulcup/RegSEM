!-----------------------------------------------------------
module module_A3d !!! for BSL models !!!
!-----------------------------------------------------------

  use def_gparam
  use earth_modele

  implicit none

  public :: init_A3d, A3d_full, &
            init_crust, get_crust

  private

  type :: nodes_modele
  doubleprecision, dimension(:,:), pointer :: param
  end type

  type :: long_modele
  type(nodes_modele), dimension(:), pointer :: longitude
  end type

  type(long_modele), dimension(:), allocatable :: latitude

  doubleprecision :: latmin_vel,latmax_vel,dlat_vel, &
                     longmin_vel,longmax_vel,dlong_vel
  doubleprecision, dimension(:), pointer :: mem_lat,mem_long
  doubleprecision, dimension(:), allocatable :: mdl,kntrad,aknot,oknot, &
                                                aknotreg,oknotreg,disc_rad,aknot2,oknot2
  doubleprecision, dimension(:,:,:), allocatable :: y2
  integer :: npar,ndisc,nknotreg,nvsplreg,nb_disc,surface,NBPARAM,n_lat,n_long
  integer, parameter :: MAXPAR = 4
  integer, dimension(MAXPAR) :: nknotA1,nknotA2
  integer, dimension(:), allocatable :: level,NBF,reg,indexreg, &
                                        vindexreg,levelreg,disc_index,level2
  integer, dimension (:,:), allocatable :: nb_nodes
  logical :: model1D_exist=.false.,modele1D_read=.false.,hknots2_exist=.false., &
             unconformal=.false.,interp=.false.,TI_interp=.false.
  character, dimension(:), allocatable :: parblock
  character(len=52), dimension(:,:), allocatable:: name_rec_tab
  character(len=20), parameter :: A3d_dat            ='A3d.dat'
  character(len=20), parameter :: hknots_dat         ='hknots.dat'
  character(len=20), parameter :: hknots2_dat        ='hknots2.dat'
  character(len=20), parameter :: model1Dfile        ='model1D.dat'
  character(len=20), parameter :: crust_dat          ='Crust.asc'

  type(modele) :: modele1D


  contains

!--------------------------------------------------------
subroutine init_crust
!--------------------------------------------------------

    integer :: i,j,k, openstatus, n

    open(18, file=crust_dat, status="old", iostat=openstatus)
    if(openstatus>0) STOP "***Cannot open the file Crust.asc***"
    read(18,*), latmin_vel, latmax_vel, dlat_vel
    if (mod(latmax_vel-latmin_vel,dlat_vel)/=0) then
        print *,"In the file Crust.asc the latitude step doesn't fit the latitude limits"
        stop
    endif
    read(18,*), longmin_vel, longmax_vel, dlong_vel
    if (longmin_vel<0.d0)   longmin_vel = 360.d0 + longmin_vel
    if (longmax_vel<0.d0)   longmax_vel = 360.d0 + longmax_vel
    if (mod(longmax_vel-longmin_vel,dlong_vel)/=0) then
        print *,"In the file Crust.asc the longitude step doesn't fit the longitude limits"
        stop
    endif
    n_lat = int((latmax_vel-latmin_vel)/dlat_vel) + 1
    n_long = int((longmax_vel-longmin_vel)/dlong_vel) + 1
    allocate (latitude(0:n_lat-1))
    allocate (mem_lat(0:n_lat-1))
    allocate (mem_long(0:n_long-1))
    allocate (nb_nodes(0:n_lat-1,0:n_long-1))
    do i = 0,n_lat-1
        allocate (latitude(i)%longitude(0:n_long-1))
        do j = 0,n_long-1
            read(18,*) mem_lat(i), mem_long(j), nb_nodes(i,j)
            mem_lat(i) = 90.d0 - mem_lat(i)
            if (mem_long(j)<0.d0)   mem_long(j) = 360.d0 + mem_long(j)
            allocate (latitude(i)%longitude(j)%param(1:nb_nodes(i,j),1:9))
            do k = 1,nb_nodes(i,j)
                read(18,*) (latitude(i)%longitude(j)%param(k,n), n=1,9)
            enddo
        enddo
    enddo    
    close(18)

!---------------------------------------------------
end subroutine init_crust
!---------------------------------------------------

!---------------------------------------------------
subroutine init_A3d
!---------------------------------------------------

    integer :: unit1=41,unit2=42,dum1,dum2,i,j,k,n,mdim
    real :: dum3
    character :: trash

    inquire(file=model1Dfile,exist=model1D_exist)
    if(model1D_exist) then
       call read_modele(model1Dfile,modele1D)
       modele1D_read=.true.
    endif

    open(unit1,file=hknots_dat,status='old')
    read(unit1,*) nknotA2(1)
    if(allocated(oknot).or.allocated(aknot).or.allocated(level)) then
       print*,'A3d already initiated'
       return
    endif
    allocate(oknot(nknotA2(1)),aknot(nknotA2(1)),level(nknotA2(1)))
    do i=1,nknotA2(1)
       read(unit1,*) oknot(i),aknot(i),level(i)
    enddo
    close(unit1)

    inquire(file=hknots2_dat,exist=hknots2_exist)
    if(hknots2_exist) then
       open(unit1,file=hknots2_dat,status='old')
       read(unit1,*) nknotA2(2)
       if(allocated(oknot2).or.allocated(aknot2).or.allocated(level2)) then
          print*,'A3d already initiated'
          return
       endif
       allocate(oknot2(nknotA2(2)),aknot2(nknotA2(2)),level2(nknotA2(2)))
       do i=1,nknotA2(2)
          read(unit1,*) oknot2(i),aknot2(i),level2(i)
       enddo
       close(unit1)
    else
       nknotA2(2)=nknotA2(1)
    endif

    !open model file and read in model
    open(unit2,file=A3d_dat,status='old')
    read(unit2,*) npar
    NBPARAM=npar
    if(npar>MAXPAR)   stop 'npar greater than MAXPAR'
    allocate(parblock(npar))
    do i=1,npar
       read(unit2,*) dum2,nknotA1(i),parblock(i)
       if(i>1.and.nknotA1(i)/=nknotA1(1)) then
          stop 'Inconsistent A1 splines between parameters'
       endif
       if(i>2) then
          nknotA2(i) = dum2
          if (nknotA2(i)/=nknotA2(2))   stop 'Param 3 and 4 need the same A2 splines than param 2'
       elseif (dum2/=nknotA2(i)) then
          stop 'Inconsistent hknots.dat and A3d.dat'
       endif
       if(i==2.and.nknotA2(i)/=nknotA2(1)) then
          unconformal=.true.
          if(.not.hknots2_exist)   stop 'unconformal grid requires hknots2.dat'
       endif
    enddo
    read(unit2,*) ndisc
    surface=0
    if (ndisc>0) then
       surface=ndisc
       if(unconformal)   print*,'discontinuities assumed same grid as first par'
       do i=1,surface
          read(unit2,*) dum2, trash
       enddo
    endif
    allocate(kntrad(nknotA1(1)))
    read(unit2,*) (kntrad(i),i=1,nknotA1(1))
    mdim=0
    do i=1,npar
       mdim=mdim+nknotA1(i)*nknotA2(i)
    enddo
    mdim=mdim+ndisc*nknotA2(1)
    allocate(mdl(mdim))
    n=0
    do i=1,npar
       do j=1,nknotA1(i)
          read(unit2,*) (mdl(k+n),k=1,nknotA2(i))
          n=n+nknotA2(i)
       enddo
    enddo
    do i=1,ndisc
       read(unit2,*) (mdl(k+n),k=1,nknotA2(1))
       n=n+nknotA2(1)
    enddo

    if(n/=mdim) stop 'init_A3d dimension error'
    close(unit2)

!---------------------------------------------------
end subroutine init_A3d
!---------------------------------------------------

!---------------------------------------------------
subroutine A3d_full(r,theta,phi,ifanis,vs,vp,rho,Qs,moho,xi,fi,eta,Gc,Gs,Hc,Hs,Bc,Bs)
!---------------------------------------------------
!returns isotropic vs, vp, and rho assuming scaling dlnVs/dlnVp=2 dlnVs/dlnrho=3
!also returns anisotropic parameters xi,fi,eta,Gc,Gs,Hc,Hs,Bc,Bs if ifanis=1

    doubleprecision, intent(in) :: r,theta,phi
    integer, intent(in) :: ifanis,moho
    doubleprecision, intent(out) :: vs,vp,rho,Qs
    doubleprecision, optional, intent(inout) :: xi,fi,eta,Gc,Gs,Hc,Hs,Bc,Bs

    integer :: jump,effnknot,i,j,k
    integer, dimension(:), allocatable :: kindex
    doubleprecision :: lat,lon,del,dr,dv,AA,CC,FF,LL,NN,eta1,adel1,r_
    doubleprecision, dimension(7) :: adel
    doubleprecision, dimension(:), allocatable :: dh

    adel=(/63.4, 31.7, 15.8, 7.9, 3.95, 1.98, 0.99/)

    if(.not.model1D_exist) &
        stop 'no 1D model file'

    if(ifanis==1) then
        if(.not.(present(xi).and.present(fi).and.present(eta))) &
            stop 'A3d_full: ifanis inconsistent'
    endif

    if(r>6371.0) then
        r_=6371.0_DP
    else
        r_=r
    endif
    call get_1Dmodel_TI(AA,CC,FF,LL,NN,rho,Qs,r_*1000.0_DP)
    if (rho<1200.d0)   call get_1Dmodel_TI(AA,CC,FF,LL,NN,rho,Qs,6367999.d0)   ! No water in RegSEM please

    eta1=FF/(AA-2.*LL)
    !Voigt average
    vp=dsqrt((3.*CC+(8.+4.*eta1)*AA+8.*(1.-eta1)*LL)/(15.*rho))
    vs=dsqrt((CC+(1.-2.*eta1)*AA+(6.+4.*eta1)*LL+5.*NN)/(15.*rho))

    if(ifanis==1) then
       eta=eta1
       xi=NN/LL
       fi=CC/AA
    endif

    ! Vs perturbation
    if (r_>kntrad(nknotA1(1)).or.r_<kntrad(1)) then
       dv=0.0
    else
       jump=0
       allocate(dh(nknotA2(1)),kindex(nknotA2(1)))

       lat=90-rad2deg*theta
       lon=rad2deg*phi
       effnknot=0
       do i=1,nknotA2(1)
          del=getdel(lat,lon,aknot(i),oknot(i))
          if(del<=adel(level(i))*2.0) then
             effnknot=effnknot+1
             kindex(effnknot)=i
             dh(effnknot)=spbsp(del,adel(level(i)))
          endif
       enddo

       dv=0.0
       do i=1,nknotA1(1)
          call fspl(i,nknotA1(1),kntrad,r_,dr)
          do j=1,effnknot
             dv=dv+dr*dh(j)*mdl(jump+kindex(j)+nknotA2(1)*(i-1))
          enddo
       enddo
       deallocate(dh,kindex)
    endif

    ! Scaling
    vs=vs+dv*vs
    vp=vp+0.5_DP*dv*vp
    rho=rho+0.33333_DP*dv*rho

    ! Perturbation of other parameters
    if(ifanis==1) then
       if(npar<2) then   ! no other parameters
          dv=0.
       else
          do k=2,npar
             if (r_>kntrad(nknotA1(k)).or.r_<kntrad(1)) then
                dv=0.0
             else
                jump = jump + nknotA1(k-1)*nknotA2(k-1)
                allocate(dh(nknotA2(k)),kindex(nknotA2(k)))

                lat=90-rad2deg*theta
                lon=rad2deg*phi
                effnknot=0
                do i=1,nknotA2(k)
                   if(unconformal) then
                      del=getdel(lat,lon,aknot2(i),oknot2(i))
                      adel1=adel(level2(i))
                   else
                      del=getdel(lat,lon,aknot(i),oknot(i))
                      adel1=adel(level(i))
                   endif
                   if(del<=adel1*2.0) then
                      effnknot=effnknot+1
                      kindex(effnknot)=i
                      dh(effnknot)=spbsp(del,adel1)
                   endif
                enddo

                dv=0.0
                do i=1,nknotA1(k)
                   call fspl(i,nknotA1(k),kntrad,r_,dr)
                   do j=1,effnknot
                      dv=dv+dr*dh(j)*mdl(jump+kindex(j)+nknotA2(k)*(i-1))
                   enddo
                enddo
                deallocate(dh,kindex)
             endif

             if (k==2) then
                xi=xi+dv*xi
                fi=fi-1.5*dv*fi
                eta=eta-2.5*dv*eta
             elseif (k==3) then
                Gc = dv
             elseif (k==4) then
                Gs = dv
             endif
             ! Here we can add a scaling to get Hc, Hs, Bc and Bs
          enddo
       endif
    endif

!---------------------------------------------------
end subroutine A3d_full
!---------------------------------------------------

!-----------------------------------------------------------------
subroutine get_1Dmodel_TI(AA,CC,FF,LL,NN,rho,Qs,rad)
!-----------------------------------------------------------------

        use module_spline

        implicit none

        doubleprecision, intent(in ) :: rad
        doubleprecision, intent(out) :: AA,CC,FF,LL,NN,rho,Qs

        integer :: i,k,ka,kb,nd,si,kas,kbs
        doubleprecision :: yp1,ypn

        if (modele1D_read.and..not.interp) then
           interp=.true.
           nb_disc=0
           do i=1,modele1D%nbcou
              if (abs(modele1D%r(i)-modele1D%r(i+1))<1e-8) &
                 nb_disc=nb_disc+1
           enddo
           allocate(disc_index(0:nb_disc+1),disc_rad(0:nb_disc+1))
           disc_index(0)=0
           disc_index(nb_disc+1)=modele1D%nbcou
           disc_rad(0)  =modele1D%r(1)
           disc_rad(nb_disc+1) =modele1D%r(modele1D%nbcou)
           k=0
           do i=1,modele1D%nbcou-1
              if (abs(modele1D%r(i)-modele1D%r(i+1))<1e-8) then
                 k=k+1
                 disc_index(k)=i
                 disc_rad  (k)=modele1D%r(i)
              endif
           enddo

           if(.not.allocated(y2)) allocate(y2(modele1D%nbcou,nb_disc+1,7))

           do nd=1,nb_disc+1
              ka=disc_index(nd-1)+1
              kb=disc_index(nd)

              yp1=(modele1D%rho(ka+1)-modele1D%rho(ka  ))/(modele1D%r(ka+1)-modele1D%r(ka  ))
              ypn=(modele1D%rho(kb  )-modele1D%rho(kb-1))/(modele1D%r(kb  )-modele1D%r(kb-1))
              call spline(modele1D%r(ka:kb),modele1D%rho(ka:kb),yp1,ypn,y2(ka:kb,nd,1)) 

              yp1=(modele1D%vpv(ka+1)-modele1D%vpv(ka  ))/(modele1D%r(ka+1)-modele1D%r(ka  ))
              ypn=(modele1D%vpv(kb  )-modele1D%vpv(kb-1))/(modele1D%r(kb  )-modele1D%r(kb-1))
              call spline(modele1D%r(ka:kb),modele1D%vpv(ka:kb),yp1,ypn,y2(ka:kb,nd,2))

              yp1=(modele1D%vsv(ka+1)-modele1D%vsv(ka  ))/(modele1D%r(ka+1)-modele1D%r(ka  ))
              ypn=(modele1D%vsv(kb  )-modele1D%vsv(kb-1))/(modele1D%r(kb  )-modele1D%r(kb-1))
              call spline(modele1D%r(ka:kb),modele1D%vsv(ka:kb),yp1,ypn,y2(ka:kb,nd,3)) 

              yp1=(modele1D%qshear(ka+1)-modele1D%qshear(ka  ))/(modele1D%r(ka+1)-modele1D%r(ka  ))
              ypn=(modele1D%qshear(kb  )-modele1D%qshear(kb-1))/(modele1D%r(kb  )-modele1D%r(kb-1))
              call spline(modele1D%r(ka:kb),modele1D%qshear(ka:kb),yp1,ypn,y2(ka:kb,nd,4)) 

              yp1=(modele1D%vph(ka+1)-modele1D%vph(ka  ))/(modele1D%r(ka+1)-modele1D%r(ka  ))
              ypn=(modele1D%vph(kb  )-modele1D%vph(kb-1))/(modele1D%r(kb  )-modele1D%r(kb-1))
              call spline(modele1D%r(ka:kb),modele1D%vph(ka:kb),yp1,ypn,y2(ka:kb,nd,5))

              yp1=(modele1D%vsh(ka+1)-modele1D%vsh(ka  ))/(modele1D%r(ka+1)-modele1D%r(ka  ))
              ypn=(modele1D%vsh(kb  )-modele1D%vsh(kb-1))/(modele1D%r(kb  )-modele1D%r(kb-1))
              call spline(modele1D%r(ka:kb),modele1D%vsv(ka:kb),yp1,ypn,y2(ka:kb,nd,6))

              yp1=(modele1D%eta(ka+1)-modele1D%eta(ka  ))/(modele1D%r(ka+1)-modele1D%r(ka  ))
              ypn=(modele1D%eta(kb  )-modele1D%eta(kb-1))/(modele1D%r(kb  )-modele1D%r(kb-1))
              call spline(modele1D%r(ka:kb),modele1D%eta(ka:kb),yp1,ypn,y2(ka:kb,nd,7))

           enddo
        endif
        nd=locate(disc_rad,nb_disc,rad)
        ka=disc_index(nd-1)+1
        kb=disc_index(nd)
        rho=splint(modele1D%r(ka:kb),modele1D%rho(ka:kb),y2(ka:kb,nd,1) &
             ,rad)
        AA=splint(modele1D%r(ka:kb),modele1D%vph(ka:kb),y2(ka:kb,nd,5) &
             ,rad)**2*rho
        CC=splint(modele1D%r(ka:kb),modele1D%vpv(ka:kb),y2(ka:kb,nd,2) &
             ,rad)**2*rho
        LL=splint(modele1D%r(ka:kb),modele1D%vsv(ka:kb),y2(ka:kb,nd,3) &
             ,rad)**2*rho
        NN=splint(modele1D%r(ka:kb),modele1D%vsh(ka:kb),y2(ka:kb,nd,6) &
             ,rad)**2*rho
        FF=splint(modele1D%r(ka:kb),modele1D%eta(ka:kb),y2(ka:kb,nd,7) &
             ,rad)*(AA-2._DP*LL)
        qs=splint(modele1D%r(ka:kb),modele1D%qshear(ka:kb) &
                   ,y2(ka:kb,nd,4),rad)

!-----------------------------------------------------------------
end subroutine get_1Dmodel_TI
!-----------------------------------------------------------------

!--------------------------------------------------------
subroutine get_crust(r,theta_rad,phi_rad,rho,vpv,vph,vsv,vsh,eta_aniso,Qmu)
!--------------------------------------------------------

implicit none

doubleprecision, intent(IN) :: r,theta_rad,phi_rad
doubleprecision, intent(OUT) :: rho,vpv,vph,vsv,vsh,eta_aniso,Qmu

logical :: special
integer :: i,j, ok
doubleprecision :: depth, Wcolat,Wlong,Wrad, theta,phi, epsil
doubleprecision, dimension(1:9) :: param
doubleprecision, dimension(1:2,1:9) :: pt_intermediaire
doubleprecision, dimension(1:8,1:9) :: pt
doubleprecision, parameter :: pi = 3.141592653, Rterre = 6371.d0

epsil = dlong_vel/1000.d0
theta = rad2deg*theta_rad
phi = rad2deg*phi_rad
if (phi>360.d0)   phi = 359.999999

depth = Rterre-r
if (depth<=0)   depth = 0.001

j = 0
do while (phi - mem_long(j) >= dlong_vel)
    j = j + 1
    if (j==n_long) then
        print *, "THE LONGITUDE RANGE DEFINED IN Crust.asc IS TOO SMALL !!!"
        STOP
    endif
enddo
special = .false.
if (j==0 .and. phi<longmin_vel) then
    if (dabs((longmax_vel+dlong_vel-360.d0)-longmin_vel) > epsil) then
        print *, "THE LONGITUDE RANGE DEFINED IN Crust.asc IS TOO SMALL !!!"
        STOP
    else
        j = n_long-1
        special = .true.
    endif
else if (j==n_long-1) then
    if (dabs((longmax_vel+dlong_vel-360.d0)-longmin_vel) > epsil) then
        print *, "THE LONGITUDE RANGE DEFINED IN Crust.asc IS TOO SMALL !!!"
        STOP
    endif
endif

ok = 0
if (theta > mem_lat(0)) then
    i = 0
    ok = 1
else
    i = 0
    boucle : do while (mem_lat(i)-theta >= dlat_vel)
        i = i + 1
        if (i==n_lat-1) then
            ok = 1
            exit boucle
        endif
    enddo boucle
endif

if (ok==1) then

    if (i==0) then
        if (mem_lat(i)+dlat_vel > 180.d0) then
            call lookfordepths(i,j,pt(1,1:9),pt(2,1:9),depth)
            if (j==n_long-1) then
                call lookfordepths(i,0,pt(3,1:9),pt(4,1:9),depth)
            else
                call lookfordepths(i,j+1,pt(3,1:9),pt(4,1:9),depth)
            endif
        else
            print *, "THE LATITUDE RANGE DEFINED IN Crust.asc IS TOO SMALL !!!"
            STOP
        endif
    else if (i==n_lat-1) then
        if (mem_lat(i)-dlat_vel < 0.d0) then
            call lookfordepths(i,j,pt(1,1:9),pt(2,1:9),depth)
            if (j==n_long-1) then
                call lookfordepths(i,0,pt(3,1:9),pt(4,1:9),depth)
            else
                call lookfordepths(i,j+1,pt(3,1:9),pt(4,1:9),depth)
            endif
        else
            print *, "THE LATITUDE RANGE DEFINED IN Crust.asc IS TOO SMALL !!!"
            STOP
        endif
    endif
    Wlong = abs(phi-mem_long(j)) / dlong_vel
    if (special)   Wlong = dabs(dlong_vel-(longmin_vel-phi)) / dlong_vel
    pt_intermediaire(2,1:9) = Wlong*pt(4,1:9) + (1-Wlong)*pt(2,1:9)
    if (depth>=pt_intermediaire(2,1)) then
        param(1:9) = pt_intermediaire(2,1:9)
    else
        pt_intermediaire(1,1:9) = Wlong*pt(3,1:9) + (1-Wlong)*pt(1,1:9)
        Wrad = abs(depth-pt_intermediaire(1,1)) / abs(pt_intermediaire(2,1)-pt_intermediaire(1,1))
        param(1:9) = Wrad*pt_intermediaire(2,1:9) + (1-Wrad)*pt_intermediaire(1,1:9)
    endif

else

    call lookfordepths(i,j,pt(1,1:9),pt(2,1:9),depth)
    if (j==n_long-1) then
        call lookfordepths(i,0,pt(3,1:9),pt(4,1:9),depth)
    else
        call lookfordepths(i,j+1,pt(3,1:9),pt(4,1:9),depth)
    endif
    call lookfordepths(i+1,j,pt(5,1:9),pt(6,1:9),depth)
    if (j==n_long-1) then
        call lookfordepths(i+1,0,pt(7,1:9),pt(8,1:9),depth)
    else
        call lookfordepths(i+1,j+1,pt(7,1:9),pt(8,1:9),depth)
    endif

    Wcolat = abs(theta-mem_lat(i)) / dlat_vel
    Wlong = abs(phi-mem_long(j)) / dlong_vel
    if (special)   Wlong = dabs(dlong_vel-(longmin_vel-phi)) / dlong_vel
    pt_intermediaire(2,1:9) = Wcolat * (Wlong*pt(8,1:9) + (1-Wlong)*pt(6,1:9)) + &
                              (1-Wcolat) * (Wlong*pt(4,1:9) + (1-Wlong)*pt(2,1:9))
    if (depth>=pt_intermediaire(2,1)) then
        param(1:9) = pt_intermediaire(2,1:9)
    else
        pt_intermediaire(1,1:9) = Wcolat * (Wlong*pt(7,1:9) + (1-Wlong)*pt(5,1:9)) + &
                                  (1-Wcolat) * (Wlong*pt(3,1:9) + (1-Wlong)*pt(1,1:9))
        Wrad = abs(depth-pt_intermediaire(1,1)) / abs(pt_intermediaire(2,1)-pt_intermediaire(1,1))
        param(1:9) = Wrad*pt_intermediaire(2,1:9) + (1-Wrad)*pt_intermediaire(1,1:9)
    endif

endif

if (param(4)<2.4) then
    rho = 2400.d0;   vpv = 4200.d0; vph = 4200.d0;   vsv = 2400.d0; vsh = 2400.d0
    eta_aniso = 1.d0
else
    rho = param(2) * 1000.d0
    vpv = param(3) * 1000.d0
    vph = param(7) * 1000.d0
    vsv = param(4) * 1000.d0
    vsh = param(8) * 1000.d0
    eta_aniso = param(9)
endif
Qmu = param(6)

!------------------------------------------------------
end subroutine get_crust
!------------------------------------------------------

!---------------------------------------------------------------
subroutine lookfordepths (latnum,longnum,ptA,ptB,depth)
!---------------------------------------------------------------

implicit none

integer, intent(IN) :: latnum,longnum
doubleprecision, intent(IN) :: depth
doubleprecision, dimension(1:9), intent(OUT) :: ptA, ptB

integer :: k

k = 1
cherche : do while (latitude(latnum)%longitude(longnum)%param(k,1) < depth)
    k = k + 1
    if (k==nb_nodes(latnum,longnum)) exit cherche
enddo cherche
ptA(:) = latitude(latnum)%longitude(longnum)%param(k-1,:)
ptB(:) = latitude(latnum)%longitude(longnum)%param(k,:)

return

!----------------------------------------------------------------
end subroutine lookfordepths
!----------------------------------------------------------------

!---------------------------------------------------
doubleprecision function getdel(a0,o0,a,o)
!---------------------------------------------------

    use def_gparam

    doubleprecision :: a0,o0,a,o
    doubleprecision :: q0,sq0,cq0,q,sq,cq,ff,sff,cff,arg

    q0=(90.0-a0)*deg2rad;
    sq0=sin(q0);
    cq0=cos(q0);

    q=(90.0-a)*deg2rad;
    sq=sin(q);
    cq=cos(q);

    ff=(o-o0)*deg2rad;
    sff=sin(ff);
    cff=cos(ff);

    arg=cq*cq0+sq*sq0*cff;
    if(arg > 1.) arg= 1.;
    if(arg < -1.) arg=-1.;
    getdel=rad2deg*acos(arg);

!---------------------------------------------------
end function getdel
!---------------------------------------------------

!---------------------------------------------------
doubleprecision function spbsp(hdel,ahdel)
!---------------------------------------------------

    doubleprecision :: hdel,ahdel

    if(hdel < ahdel) then
       spbsp=0.75*(hdel/ahdel)*(hdel/ahdel)*(hdel/ahdel)-1.5*(hdel/ahdel)*(hdel/ahdel)+1.
    else if(hdel <= ahdel*2.) then
       spbsp=0.25*(2.-hdel/ahdel)*(2.-hdel/ahdel)*(2.-hdel/ahdel)
       if(spbsp < 0) spbsp=0.0
    else
       spbsp=0.0
    endif

!---------------------------------------------------
end function spbsp
!---------------------------------------------------

!-----------------------------------------------------------------
integer function locate(inter,n,x)
!-----------------------------------------------------------------

         use def_gparam

         implicit none

         integer, intent(in) :: n
         doubleprecision, dimension(0:n+1), intent(in) :: inter
         doubleprecision, intent(in) :: x

         integer :: i

         i=0
         if (abs (x-inter(0))/max(x,inter(0)) <1.E-10_DP) then
            i=0
         else if ( abs (x-inter(n+1))/max(x,inter(n+1))  <1.E-10_DP) then
            i=n
	 else if (x > inter(n) ) then
	    i=n
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

!-----------------------------------------------------------
end module module_A3d
!-----------------------------------------------------------
