!==========================================
subroutine save_model(Tdomain,rg)
!==========================================
  
  use sdomains
  
  implicit none
  
  type (Domain), intent (INOUT) :: Tdomain
  integer, intent (IN) :: rg
  integer :: nrec, n, i, ngll1, ngll2, ngll3, nxint, nyint, nzint, xmiddle, ymiddle, zmiddle
  integer :: x, y, z, xint, yint, zint, indx, indy, indz, imin, imax, jmin, jmax, kmin, kmax
  real :: xloc,yloc,zloc
  real, dimension(0:2) ::  Vp, Vs
  doubleprecision :: rho, rhox2
  character*100 :: Vpfile, Vsfile

  integer, parameter :: which_pts = 1 ! Defines the output points (see below)

 
  if (Tdomain%save_model) then

     write (Vpfile,"(a,I3.3)") "Vp",rg
     open (12,file=trim(Vpfile),form='UNFORMATTED',access='stream')
     do i = 1,7
        write(12)nrec
     enddo

     write (Vsfile,"(a,I3.3)") "Vs",rg
     open (13,file=trim(Vsfile),form='UNFORMATTED',access='stream')
     do i = 1,7
        write(13)nrec
     enddo
     
     nrec=0

     imin =  999999.
     imax = -999999.
     jmin =  999999.
     jmax = -999999.
     kmin =  999999.
     kmax = -999999.

     do n = 0,Tdomain%n_elem-1

        if (.not. Tdomain%specel(n)%PML) then

           ngll1 = Tdomain%specel(n)%ngllx
           ngll2 = Tdomain%specel(n)%nglly
           ngll3 = Tdomain%specel(n)%ngllz

if (which_pts==1) then ! GLL interieurs
           nxint = ngll1 - 2
           nyint = ngll2 - 2
           nzint = ngll3 - 2
else if (which_pts==2) then ! Corners
           nxint = 2
           nyint = 2
           nzint = 2
else if (which_pts==3) then ! Corners + mid-edges
           if (modulo(ngll1,2)==0 .or. modulo(ngll2,2)==0 .or. modulo(ngll3,2)==0) then
              print *, "ERROR: AN ODD NUMBER OF GLLs IS NEEDED WHEN which_pts=3"
           else
              nxint = 3
              nyint = 3
              nzint = 3
           endif
else
           stop 'ERROR: WRONG INPUT FOR which_pts'
endif


           indz = Tdomain%specel(n)%indz*(nzint-1)
           do zint = 0,nzint-1
              indy = Tdomain%specel(n)%indy*(nyint-1)
              do yint = 0,nyint-1
                 indx = Tdomain%specel(n)%indx*(nxint-1)
                 do xint = 0,nxint-1
                    
                    nrec = nrec+1

                    kmin = min(kmin,indz)
                    kmax = max(kmax,indz)
                    jmin = min(jmin,indy)
                    jmax = max(jmax,indy)
                    imin = min(imin,indx)
                    imax = max(imax,indx)

if (which_pts==1) then ! GLL interieurs
                    x = xint + 1
                    y = yint + 1
                    z = zint + 1
else if (which_pts==2) then ! Corners
                    x = 0;   if (xint==1) x = ngll1-1
                    y = 0;   if (yint==1) y = ngll2-1
                    z = 0;   if (zint==1) z = ngll3-1
else if (which_pts==3) then ! Corners + mid-edges
                    x = 0;   if (xint==1) x = int(ngll1/2);   if (xint==2) x = ngll1-1
                    y = 0;   if (yint==1) y = int(ngll2/2);   if (yint==2) y = ngll2-1
                    z = 0;   if (zint==1) z = int(ngll3/2);   if (zint==2) z = ngll3-1
endif

                    i = Tdomain%specel(n)%Iglobnum(x,y,z)
                    xloc = Tdomain%Globcoord(0,i)
                    yloc = Tdomain%Globcoord(1,i)
                    zloc = Tdomain%Globcoord(2,i)

                    rho = Tdomain%specel(n)%Density(x,y,z)
                    if (Tdomain%aniso) then
                       Vp(0) = sqrt(Tdomain%specel(n)%Cij(0,x,y,z)/rho)
                       Vp(1) = sqrt(Tdomain%specel(n)%Cij(6,x,y,z)/rho)
                       Vp(2) = sqrt(Tdomain%specel(n)%Cij(11,x,y,z)/rho)
                    else
                       Vp(:) = sqrt((Tdomain%specel(n)%Lambda(x,y,z) + &
                                     2.d0*Tdomain%specel(n)%Mu(x,y,z)) &
                                    / rho)
                    endif

                    rhox2 = 2.d0*Tdomain%specel(n)%Density(x,y,z)
                    if (Tdomain%aniso) then
                       Vs(0) = sqrt(Tdomain%specel(n)%Cij(15,x,y,z)/rhox2)
                       Vs(1) = sqrt(Tdomain%specel(n)%Cij(18,x,y,z)/rhox2)
                       Vs(2) = sqrt(Tdomain%specel(n)%Cij(20,x,y,z)/rhox2)
                    else
                       Vs(:) = sqrt(Tdomain%specel(n)%Mu(x,y,z)/rho)
                    endif

                    write (12) indx, indy, indz, xloc, yloc, zloc
                    write (12)Vp(:)

                    write (13) indx, indy, indz, xloc, yloc, zloc
                    write (13)Vs(:)

                    indx = indx+1
                 enddo
                 indy = indy+1
              enddo
              indz = indz+1
           enddo
           
        endif

     enddo

     write(12,pos=1)nrec,imin,imax,jmin,jmax,kmin,kmax
     close (12)

     write(13,pos=1)nrec,imin,imax,jmin,jmax,kmin,kmax
     close (13)

  endif

  return

!============================
end subroutine save_model
!============================
