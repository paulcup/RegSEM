subroutine deallocate_domain (Tdomain, rg)


use sdomains

implicit none

type(domain), intent (INOUT) :: Tdomain
integer, intent (IN) :: rg

integer :: n, i


deallocate (Tdomain%GlobCoord)
deallocate (Tdomain%Coord_Nodes)

do n = 0,Tdomain%n_elem-1
   deallocate (Tdomain%specel(n)%Density)
   deallocate (Tdomain%specel(n)%MassMat)
   deallocate (Tdomain%specel(n)%IglobNum)
   deallocate (Tdomain%specel(n)%Control_Nodes)
   deallocate (Tdomain%specel(n)%Jacob)
   do i = 0,Tdomain%nb_simu-1
      deallocate (Tdomain%specel(n)%sSimu(i)%Veloc)
      deallocate (Tdomain%specel(n)%sSimu(i)%Accel)
      deallocate (Tdomain%specel(n)%sSimu(i)%Forces)
   enddo
   if (Tdomain%specel(n)%PML) then
      deallocate (Tdomain%specel(n)%Lambda)
      deallocate (Tdomain%specel(n)%Mu)
      deallocate (Tdomain%specel(n)%DumpSx)
      deallocate (Tdomain%specel(n)%DumpSy)
      deallocate (Tdomain%specel(n)%DumpSz)
      deallocate (Tdomain%specel(n)%DumpVx)
      deallocate (Tdomain%specel(n)%DumpVy)
      deallocate (Tdomain%specel(n)%DumpVz)
      deallocate (Tdomain%specel(n)%Acoeff)
      do i = 0,Tdomain%nb_simu-1
          deallocate (Tdomain%specel(n)%sSimu(i)%Diagonal_Stress)
          deallocate (Tdomain%specel(n)%sSimu(i)%Diagonal_Stress1)
          deallocate (Tdomain%specel(n)%sSimu(i)%Diagonal_Stress2)
          deallocate (Tdomain%specel(n)%sSimu(i)%Diagonal_Stress3)
          deallocate (Tdomain%specel(n)%sSimu(i)%Residual_Stress)
          deallocate (Tdomain%specel(n)%sSimu(i)%Residual_Stress1)
          deallocate (Tdomain%specel(n)%sSimu(i)%Residual_Stress2)
          deallocate (Tdomain%specel(n)%sSimu(i)%Veloc1)
          deallocate (Tdomain%specel(n)%sSimu(i)%Veloc2)
          deallocate (Tdomain%specel(n)%sSimu(i)%Veloc3)
          deallocate (Tdomain%specel(n)%sSimu(i)%Forces1)
          deallocate (Tdomain%specel(n)%sSimu(i)%Forces2)
          deallocate (Tdomain%specel(n)%sSimu(i)%Forces3)
      enddo
      if (Tdomain%curve) then
          deallocate (Tdomain%specel(n)%Normales)
          deallocate (Tdomain%specel(n)%Inv_Normales)
          do i = 0,Tdomain%nb_simu-1
             deallocate (Tdomain%specel(n)%sSimu(i)%Residual_Stress3)
          enddo
      endif
   else
      deallocate (Tdomain%specel(n)%wgtx)
      deallocate (Tdomain%specel(n)%wgty)
      deallocate (Tdomain%specel(n)%wgtz)
      do i = 0,Tdomain%nb_simu-1
          deallocate (Tdomain%specel(n)%sSimu(i)%Displ)
          deallocate (Tdomain%specel(n)%sSimu(i)%V0)
      enddo
      if (Tdomain%aniso) then
         deallocate (Tdomain%specel(n)%Cij)
         if (Tdomain%n_sls>0) then
             deallocate (Tdomain%specel(n)%Lambda)
             deallocate (Tdomain%specel(n)%Mu)
         endif
         if (Tdomain%adjoint) then
            deallocate (Tdomain%specel(n)%Kern_rho)
            deallocate (Tdomain%specel(n)%Kern_aniso)
            do i = 0,Tdomain%nb_simu-1
               deallocate (Tdomain%specel(n)%sSimu(i)%save_strain)
            enddo
         endif
      else
         deallocate (Tdomain%specel(n)%Lambda)
         deallocate (Tdomain%specel(n)%Mu)
         if (Tdomain%adjoint) then
            deallocate (Tdomain%specel(n)%Kern_rho)
            deallocate (Tdomain%specel(n)%Kern_lambda)
            deallocate (Tdomain%specel(n)%Kern_mu)
            do i = 0,Tdomain%nb_simu-1
               deallocate (Tdomain%specel(n)%sSimu(i)%save_strain)
            enddo
         endif
      endif
      if (Tdomain%n_sls>0) then
         deallocate (Tdomain%specel(n)%Q)
         deallocate (Tdomain%specel(n)%onemSbeta)
         deallocate (Tdomain%specel(n)%factor_common_3)
         deallocate (Tdomain%specel(n)%alphaval_3)
         deallocate (Tdomain%specel(n)%betaval_3)
         deallocate (Tdomain%specel(n)%gammaval_3)
         do i = 0,Tdomain%nb_simu-1 
            deallocate (Tdomain%specel(n)%sSimu(i)%epsilondev_xx_)
            deallocate (Tdomain%specel(n)%sSimu(i)%epsilondev_yy_)
            deallocate (Tdomain%specel(n)%sSimu(i)%epsilondev_xy_)
            deallocate (Tdomain%specel(n)%sSimu(i)%epsilondev_xz_)
            deallocate (Tdomain%specel(n)%sSimu(i)%epsilondev_yz_)
            deallocate (Tdomain%specel(n)%sSimu(i)%R_xx_)
            deallocate (Tdomain%specel(n)%sSimu(i)%R_yy_)
            deallocate (Tdomain%specel(n)%sSimu(i)%R_xy_)
            deallocate (Tdomain%specel(n)%sSimu(i)%R_xz_)
            deallocate (Tdomain%specel(n)%sSimu(i)%R_yz_)
         enddo
      endif
   endif
   deallocate (Tdomain%specel(n)%sSimu)
   if (Tdomain%specel(n)%random_t)   deallocate (Tdomain%specel(n)%random_coeff)
enddo

do n = 0, Tdomain%n_face-1
   deallocate (Tdomain%sFace(n)%MassMat)
   do i = 0,Tdomain%nb_simu-1
      deallocate (Tdomain%sFace(n)%sSimu(i)%Veloc)
      deallocate (Tdomain%sFace(n)%sSimu(i)%Forces)
      deallocate (Tdomain%sFace(n)%sSimu(i)%Accel)
   enddo
   if (Tdomain%sFace(n)%ocean) then
      deallocate (Tdomain%sFace(n)%M33ocean)
   endif  
   if (Tdomain%sFace(n)%PML) then 
      deallocate (Tdomain%sFace(n)%DumpVx) 
      deallocate (Tdomain%sFace(n)%DumpVy) 
      deallocate (Tdomain%sFace(n)%DumpVz)
      do i = 0,Tdomain%nb_simu-1
         deallocate (Tdomain%sFace(n)%sSimu(i)%Forces1)
         deallocate (Tdomain%sFace(n)%sSimu(i)%Forces2)
         deallocate (Tdomain%sFace(n)%sSimu(i)%Forces3) 
         deallocate (Tdomain%sFace(n)%sSimu(i)%Veloc1)
         deallocate (Tdomain%sFace(n)%sSimu(i)%Veloc2)
         deallocate (Tdomain%sFace(n)%sSimu(i)%Veloc3)
      enddo
   else
      do i = 0,Tdomain%nb_simu-1
         deallocate (Tdomain%sFace(n)%sSimu(i)%Displ)
         deallocate (Tdomain%sFace(n)%sSimu(i)%V0)
      enddo
   endif
   deallocate (Tdomain%sFace(n)%sSimu)
enddo

do n = 0,Tdomain%n_edge-1
   deallocate (Tdomain%sEdge(n)%MassMat)
   do i = 0,Tdomain%nb_simu-1
      deallocate (Tdomain%sEdge(n)%sSimu(i)%Veloc)
      deallocate (Tdomain%sEdge(n)%sSimu(i)%Forces)
      deallocate (Tdomain%sEdge(n)%sSimu(i)%Accel)
   enddo
   if (Tdomain%sEdge(n)%ocean) then
      deallocate (Tdomain%sEdge(n)%M33ocean)
   endif
   if (Tdomain%sEdge(n)%PML) then
      deallocate (Tdomain%sEdge(n)%DumpVx)
      deallocate (Tdomain%sEdge(n)%DumpVy)
      deallocate (Tdomain%sEdge(n)%DumpVz)
      do i = 0,Tdomain%nb_simu-1
         deallocate (Tdomain%sEdge(n)%sSimu(i)%Forces1)
         deallocate (Tdomain%sEdge(n)%sSimu(i)%Forces2)
         deallocate (Tdomain%sEdge(n)%sSimu(i)%Forces3)
         deallocate (Tdomain%sEdge(n)%sSimu(i)%Veloc1)
         deallocate (Tdomain%sEdge(n)%sSimu(i)%Veloc2)
         deallocate (Tdomain%sEdge(n)%sSimu(i)%Veloc3)
      enddo
   else
      do i = 0,Tdomain%nb_simu-1
         deallocate (Tdomain%sEdge(n)%sSimu(i)%Displ)
         deallocate (Tdomain%sEdge(n)%sSimu(i)%V0)
      enddo
   endif
   deallocate (Tdomain%sEdge(n)%sSimu)
enddo

do n = 0,Tdomain%n_vertex-1
   do i = 0,Tdomain%nb_simu-1
      deallocate (Tdomain%sVertex(n)%sSimu(i)%Veloc)
      deallocate (Tdomain%sVertex(n)%sSimu(i)%Accel)
      deallocate (Tdomain%sVertex(n)%sSimu(i)%Forces)
   enddo
   if (Tdomain%sVertex(n)%ocean) then
      deallocate (Tdomain%sVertex(n)%M33ocean)
   endif
   if (Tdomain%sVertex(n)%PML) then
      deallocate (Tdomain%sVertex(n)%DumpVx)
      deallocate (Tdomain%sVertex(n)%DumpVy)       
      deallocate (Tdomain%sVertex(n)%DumpVz)
      do i = 0,Tdomain%nb_simu-1
         deallocate (Tdomain%sVertex(n)%sSimu(i)%Forces1)
         deallocate (Tdomain%sVertex(n)%sSimu(i)%Forces2)
         deallocate (Tdomain%sVertex(n)%sSimu(i)%Forces3)
         deallocate (Tdomain%sVertex(n)%sSimu(i)%Veloc1)
         deallocate (Tdomain%sVertex(n)%sSimu(i)%Veloc2)
         deallocate (Tdomain%sVertex(n)%sSimu(i)%Veloc3)
      enddo
   else
      do i = 0,Tdomain%nb_simu-1
         deallocate (Tdomain%sVertex(n)%sSimu(i)%Displ)
         deallocate (Tdomain%sVertex(n)%sSimu(i)%V0)
      enddo
   endif
   deallocate (Tdomain%sVertex(n)%sSimu)
enddo

do n = 0,Tdomain%n_proc-1
    if (Tdomain%sComm(n)%ngll>0) then
        deallocate (Tdomain%sComm(n)%GiveForces)
        deallocate (Tdomain%sComm(n)%TakeForces)
    endif
    if (Tdomain%sComm(n)%ngllPML>0) then
        deallocate (Tdomain%sComm(n)%GiveForcesPML)
        deallocate (Tdomain%sComm(n)%TakeForcesPML)
    endif
enddo

do n = 0, Tdomain%n_mat-1 
    if (associated(Tdomain%sSubdomain(n)%GLLcz, Tdomain%sSubdomain(n)%GLLcx) .or. &
        associated(Tdomain%sSubdomain(n)%GLLcz, Tdomain%sSubdomain(n)%GLLcy )) then
        nullify (Tdomain%sSubdomain(n)%GLLcz)
        nullify(Tdomain%sSubdomain(n)%GLLwz)
        nullify(Tdomain%sSubdomain(n)%hprimez)
        nullify(Tdomain%sSubdomain(n)%hTprimez)
    else
        deallocate (Tdomain%sSubdomain(n)%GLLcz)
        deallocate (Tdomain%sSubdomain(n)%GLLwz)
        deallocate (Tdomain%sSubdomain(n)%hprimez)
        deallocate (Tdomain%sSubdomain(n)%hTprimez)
    endif
    if (associated(Tdomain%sSubdomain(n)%GLLcy, Tdomain%sSubdomain(n)%GLLcx) ) then
        nullify (Tdomain%sSubdomain(n)%GLLcy)
        nullify(Tdomain%sSubdomain(n)%GLLwy)
        nullify(Tdomain%sSubdomain(n)%hprimey)
        nullify(Tdomain%sSubdomain(n)%hTprimey)
    else
        deallocate (Tdomain%sSubdomain(n)%GLLcy)
        deallocate (Tdomain%sSubdomain(n)%GLLwy)
        deallocate (Tdomain%sSubdomain(n)%hprimey)
        deallocate (Tdomain%sSubdomain(n)%hTprimey)
    endif
    deallocate (Tdomain%sSubdomain(n)%GLLcx)
    deallocate (Tdomain%sSubdomain(n)%GLLwx)
    deallocate (Tdomain%sSubdomain(n)%hprimex)
    deallocate (Tdomain%sSubdomain(n)%hTprimex)
enddo

do n = 0, Tdomain%n_receivers-1
    if (rg == Tdomain%sReceiver(n)%proc) then
        deallocate (Tdomain%sReceiver(n)%StoreTrace)
        deallocate (Tdomain%sReceiver(n)%pol)
        deallocate (Tdomain%sReceiver(n)%coeff)
    endif
enddo

do n = 0, Tdomain%n_source-1
    if (Tdomain%sSource(n)%i_type_source==1 .or. Tdomain%sSource(n)%i_type_source==2) then
     if (rg==Tdomain%sSource(n)%proc) then
         deallocate (Tdomain%sSource(n)%coeff)
         deallocate (Tdomain%sSource(n)%timefunc)
     endif
    endif
enddo

return
end subroutine deallocate_domain
