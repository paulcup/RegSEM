!==========================================
subroutine save_snapshots(Tdomain,ntime,rg)
!==========================================
  
  use sdomains
  use angles
  
  implicit none
  
  type (Domain), intent (INOUT) :: Tdomain
  integer, intent (IN) :: ntime, rg
  integer :: nrec, n, i,j, ngll1, ngll2, ngll3, nsimu, nxint, nyint, nzint, mat
  integer :: x, y, z, xint, yint, zint, isimu, nf, ne, nv, indx, indy, indz, imin, imax, jmin, jmax, kmin, kmax
  double precision :: xi, eta, zeta, xi_int, eta_int, zeta_int, xa, ya, za, a
  real :: xloc,yloc,zloc
  real, allocatable :: displ(:,:,:,:,:),displ_gll(:,:)
  double precision, allocatable :: pol(:,:,:)
  character*100 :: snapfile, adjfile

 
!  if (Tdomain%save_snapshots .and.((Tdomain%t_reversal_mirror<=1.and.mod(ntime,Tdomain%sTimeParam%dnsnap)==0)&
!       &.or.(Tdomain%t_reversal_mirror==2.and.mod(Tdomain%sTimeParam%ntime-ntime-2,Tdomain%sTimeParam%dnsnap)==0)))then
     !if (Tdomain%save_snapshots .and. mod(ntime,Tdomain%sTimeParam%dnsnap)==0 .and. ntime/=0) then
  if (Tdomain%save_snapshot .and. &
      ( (Tdomain%t_reversal_mirror<=1 .and. mod(ntime,Tdomain%sTimeParam%dnsnap)==0) .or. &
        (Tdomain%t_reversal_mirror==2 .and. mod(Tdomain%sTimeParam%ntime-ntime-2,Tdomain%sTimeParam%dnsnap)==0) ) ) then

     if(Tdomain%t_reversal_mirror==2)then
        write (snapfile,"(a,I3.3,a,I3.3)") "snapshot_backward_",Tdomain%sTimeParam%Nsnap,"_",rg
     else
        write (snapfile,"(a,I3.3,a,I3.3)") "snapshot_forward_",Tdomain%sTimeParam%Nsnap,"_",rg
     endif
     open (12,file=trim(snapfile),form='UNFORMATTED',access='stream')
     do i = 1,7
        write(12)nrec
     enddo
     if (Tdomain%adjoint) then
        write (adjfile,"(a,I3.3,a,I3.3)") "snapshot_adjoint_",Tdomain%sTimeParam%Nsnap,"_",rg
        open (22,file=trim(adjfile),form='UNFORMATTED',access='stream')
        do i = 1,7
           write(22)nrec
        enddo
     endif
     
     nrec=0
     
     if(Tdomain%adjoint)then
        nsimu = 1
        allocate(displ_gll(0:1,0:2))
     else
        nsimu = 0
        allocate(displ_gll(0:0,0:2))
     endif

     imin =  999999.
     imax = -999999.
     jmin =  999999.
     jmax = -999999.
     kmin =  999999.
     kmax = -999999.

     !     Elements

     do n = 0,Tdomain%n_elem-1

        if (.not. Tdomain%specel(n)%PML) then

           mat = Tdomain%specel(n)%mat_index
           
           ngll1 = Tdomain%specel(n)%ngllx
           ngll2 = Tdomain%specel(n)%nglly
           ngll3 = Tdomain%specel(n)%ngllz
           
           allocate (displ(0:ngll1-1,0:ngll2-1,0:ngll3-1,0:2,0:nsimu))
           
           displ = 0
           
           ! get displ
           
           do isimu = 0,nsimu
       
           displ(1:ngll1-2,1:ngll2-2,1:ngll3-2,0:2,isimu) = Tdomain%specel(n)%sSimu(isimu)%Displ(:,:,:,:)
           
           ! faces
           
           nf = Tdomain%specel(n)%near_faces(0)
           Displ(1:ngll1-2,1:ngll2-2,0,0:2,isimu) = Tdomain%sFace(nf)%sSimu(isimu)%Displ(1:ngll1-2,1:ngll2-2,0:2) 
           
           nf = Tdomain%specel(n)%near_faces(1)
           Displ(1:ngll1-2,0,1:ngll3-2,0:2,isimu) = Tdomain%sFace(nf)%sSimu(isimu)%Displ(1:ngll1-2,1:ngll3-2,0:2) 
           
           nf = Tdomain%specel(n)%near_faces(2)
           Displ(ngll1-1,1:ngll2-2,1:ngll3-2,0:2,isimu) = Tdomain%sFace(nf)%sSimu(isimu)%Displ(1:ngll2-2,1:ngll3-2,0:2) 
           
           nf = Tdomain%specel(n)%near_faces(3)
           Displ(1:ngll1-2,ngll2-1,1:ngll3-2,0:2,isimu) = Tdomain%sFace(nf)%sSimu(isimu)%Displ(1:ngll1-2,1:ngll3-2,0:2) 
           
           nf = Tdomain%specel(n)%near_faces(4)
           Displ(0,1:ngll2-2,1:ngll3-2,0:2,isimu) = Tdomain%sFace(nf)%sSimu(isimu)%Displ(1:ngll2-2,1:ngll3-2,0:2) 
           
           nf = Tdomain%specel(n)%near_faces(5)
           Displ(1:ngll1-2,1:ngll2-2,ngll3-1,0:2,isimu) = Tdomain%sFace(nf)%sSimu(isimu)%Displ(1:ngll1-2,1:ngll2-2,0:2) 
           
           ! edges
           
           ne = Tdomain%specel(n)%near_edges(0)
           Displ(1:ngll1-2,0,0,0:2,isimu) = Tdomain%sEdge(ne)%sSimu(isimu)%Displ(1:ngll1-2,0:2) 
           
           ne = Tdomain%specel(n)%near_edges(1)
           Displ(ngll1-1,1:ngll2-2,0,0:2,isimu) = Tdomain%sEdge(ne)%sSimu(isimu)%Displ(1:ngll2-2,0:2)
           
           ne = Tdomain%specel(n)%near_edges(2)
           Displ(1:ngll1-2,ngll2-1,0,0:2,isimu) = Tdomain%sEdge(ne)%sSimu(isimu)%Displ(1:ngll1-2,0:2)
           
           ne = Tdomain%specel(n)%near_edges(3)
           Displ(0,1:ngll2-2,0,0:2,isimu) = Tdomain%sEdge(ne)%sSimu(isimu)%Displ(1:ngll2-2,0:2)
           
           ne = Tdomain%specel(n)%near_edges(4)
           Displ(ngll1-1,0,1:ngll3-2,0:2,isimu) = Tdomain%sEdge(ne)%sSimu(isimu)%Displ(1:ngll3-2,0:2)
           
           ne = Tdomain%specel(n)%near_edges(5)
           Displ(1:ngll1-2,0,ngll3-1,0:2,isimu) = Tdomain%sEdge(ne)%sSimu(isimu)%Displ(1:ngll1-2,0:2)
           
           ne = Tdomain%specel(n)%near_edges(6)
           Displ(0,0,1:ngll3-2,0:2,isimu) = Tdomain%sEdge(ne)%sSimu(isimu)%Displ(1:ngll3-2,0:2)
           
           ne = Tdomain%specel(n)%near_edges(7)
           Displ(ngll1-1,ngll2-1,1:ngll3-2,0:2,isimu) = Tdomain%sEdge(ne)%sSimu(isimu)%Displ(1:ngll3-2,0:2)
           
           ne = Tdomain%specel(n)%near_edges(8)
           Displ(ngll1-1,1:ngll2-2,ngll3-1,0:2,isimu) = Tdomain%sEdge(ne)%sSimu(isimu)%Displ(1:ngll2-2,0:2)
           
           ne = Tdomain%specel(n)%near_edges(9)
           Displ(1:ngll1-2,ngll2-1,ngll3-1,0:2,isimu) = Tdomain%sEdge(ne)%sSimu(isimu)%Displ(1:ngll1-2,0:2)
           
           ne = Tdomain%specel(n)%near_edges(10)
           Displ(0,ngll2-1,1:ngll3-2,0:2,isimu) = Tdomain%sEdge(ne)%sSimu(isimu)%Displ(1:ngll3-2,0:2)
           
           ne = Tdomain%specel(n)%near_edges(11)
           Displ(0,1:ngll2-2,ngll3-1,0:2,isimu) = Tdomain%sEdge(ne)%sSimu(isimu)%Displ(1:ngll2-2,0:2)
           
           ! vertices
           
           nv = Tdomain%specel(n)%near_vertices(0)
           Displ(0,0,0,0:2,isimu) = Tdomain%sVertex(nv)%sSimu(isimu)%Displ(0:2) 
           
           nv = Tdomain%specel(n)%near_vertices(1)
           Displ(ngll1-1,0,0,0:2,isimu) = Tdomain%sVertex(nv)%sSimu(isimu)%Displ(0:2)
           
           nv = Tdomain%specel(n)%near_vertices(2)
           Displ(ngll1-1,ngll2-1,0,0:2,isimu) = Tdomain%sVertex(nv)%sSimu(isimu)%Displ(0:2)
           
           nv = Tdomain%specel(n)%near_vertices(3)
           Displ(0,ngll2-1,0,0:2,isimu) = Tdomain%sVertex(nv)%sSimu(isimu)%Displ(0:2)
           
           nv = Tdomain%specel(n)%near_vertices(4)
           Displ(0,0,ngll3-1,0:2,isimu) = Tdomain%sVertex(nv)%sSimu(isimu)%Displ(0:2)
           
           nv = Tdomain%specel(n)%near_vertices(5)
           Displ(ngll1-1,0,ngll3-1,0:2,isimu) = Tdomain%sVertex(nv)%sSimu(isimu)%Displ(0:2)
           
           nv = Tdomain%specel(n)%near_vertices(6)
           Displ(ngll1-1,ngll2-1,ngll3-1,0:2,isimu) = Tdomain%sVertex(nv)%sSimu(isimu)%Displ(0:2)
           
           nv = Tdomain%specel(n)%near_vertices(7)
           Displ(0,ngll2-1,ngll3-1,0:2,isimu) = Tdomain%sVertex(nv)%sSimu(isimu)%Displ(0:2)
       
!	   Displ(:,:,:,0,isimu)=Tdomain%specel(n)%win_mirror(:,:,:)
!	   Displ(:,:,:,1,isimu)=Tdomain%specel(n)%win_mirror(:,:,:)
!	   Displ(:,:,:,2,isimu)=Tdomain%specel(n)%win_mirror(:,:,:)

           enddo

           nxint = Tdomain%nx_interp
           nyint = Tdomain%ny_interp
           nzint = Tdomain%nz_interp

           ! write field values
           
           indz = Tdomain%specel(n)%indz*(nzint-1)
           do zint = 0,nzint-1
              indy = Tdomain%specel(n)%indy*(nyint-1)
              zeta_int = -1+2*real(zint)/real(nzint-1) 
              do yint = 0,nyint-1
                 indx = Tdomain%specel(n)%indx*(nxint-1)
                 eta_int  =  -1+2*real(yint)/real(nyint-1)
                 do xint = 0,nxint-1
                    xi_int   = -1+2*real(xint)/real(nxint-1)
                    
                    nrec = nrec+1

                    kmin = min(kmin,indz)
                    kmax = max(kmax,indz)
                    jmin = min(jmin,indy)
                    jmax = max(jmax,indy)
                    imin = min(imin,indx)
                    imax = max(imax,indx)

                    ! compute coefs for interpolation
                    
                    allocate (pol(0:ngll1-1,0:ngll2-1,0:ngll3-1))
                    pol = 1
                    do x = 0,ngll1-1
                       do y = 0,ngll2-1
                          do z = 0,ngll3-1
                             xi = Tdomain%sSubdomain(mat)%GLLcx(x)
                             eta = Tdomain%sSubdomain(mat)%GLLcy(y)
                             zeta = Tdomain%sSubdomain(mat)%GLLcz(z)
                             do j = 0,ngll1-1
                                if (j/=x) then
                                   xa = Tdomain%sSubdomain(mat)%GLLcx(j)
                                   pol(x,y,z) = pol(x,y,z) * (xi_int-xa)/(xi-xa)
                                endif
                             enddo
                             do j = 0,ngll2-1
                                if (j/=y) then
                                   ya = Tdomain%sSubdomain(mat)%GLLcy(j)
                                   pol(x,y,z) = pol(x,y,z) * (eta_int-ya)/(eta-ya)
                                endif
                             enddo
                             do j = 0,ngll3-1
                                if (j/=z) then
                                   za = Tdomain%sSubdomain(mat)%GLLcz(j)
                                   pol(x,y,z) = pol(x,y,z) * (zeta_int-za)/(zeta-za)
                                endif
                             enddo
                          enddo
                       enddo
                    enddo

                    ! interpolate coordinates
                    
                    xloc = 0
                    yloc = 0
                    zloc = 0

                    displ_gll = 0

                    do x = 0,ngll1-1
                       do y = 0,ngll2-1
                          do z = 0,ngll3-1
                             i = Tdomain%specel(n)%Iglobnum(x,y,z)
                             xloc = xloc + Tdomain%Globcoord(0,i) * pol(x,y,z)
                             yloc = yloc + Tdomain%Globcoord(1,i) * pol(x,y,z)
                             zloc = zloc + Tdomain%Globcoord(2,i) * pol(x,y,z)

                             do isimu = 0,nsimu
                                displ_gll(isimu,:) = displ_gll(isimu,:) + displ(x,y,z,:,isimu) * pol(x,y,z)
                             enddo

                          enddo
                       enddo
                    enddo
                    
                    deallocate(pol)
                    
                    write (12) indx, indy, indz, xloc, yloc, zloc
                    if (Tdomain%adjoint) write (22) indx, indy, indz, xloc, yloc, zloc
                    
                    write (12)displ_gll(0,:)

                    ! save adjoint snapshots
                    
                    if (Tdomain%adjoint) then
                       write (22)displ_gll(1,:)
                    endif
                    indx = indx+1
                 enddo
                 indy = indy+1
              enddo
              indz = indz+1
           enddo
           
           deallocate(displ)
           
        endif
     enddo
     write(12,pos=1)nrec,imin,imax,jmin,jmax,kmin,kmax
     close (12)
     if (Tdomain%adjoint) then
        write(22,pos=1)nrec,imin,imax,jmin,jmax,kmin,kmax
        close (22)
     endif
  endif
  !
  return
  !
!============================
end subroutine save_snapshots
!============================
