subroutine allocate_domain (Tdomain, rg)


use sdomains

implicit none

type(domain), intent (INOUT) :: Tdomain
integer, intent (IN) :: rg

integer :: n, nf,ne,nv, i, ngllx,nglly,ngllz, ngll1,ngll2, ngll, ngllPML, ngllocean, n_solid


n_solid = Tdomain%n_sls

do n = 0,Tdomain%n_elem-1
    ngllx = Tdomain%specel(n)%ngllx
    nglly = Tdomain%specel(n)%nglly
    ngllz = Tdomain%specel(n)%ngllz
    allocate (Tdomain%specel(n)%Density (0:ngllx-1, 0:nglly-1, 0:ngllz-1))
    allocate (Tdomain%specel(n)%MassMat (0:ngllx-1, 0:nglly-1, 0:ngllz-1))
    allocate (Tdomain%specel(n)%sSimu(0:Tdomain%nb_simu-1))
    do i = 0,Tdomain%nb_simu-1
        allocate (Tdomain%specel(n)%sSimu(i)%Veloc (1:ngllx-2, 1:nglly-2, 1:ngllz-2, 0:2))
        allocate (Tdomain%specel(n)%sSimu(i)%Accel (1:ngllx-2, 1:nglly-2, 1:ngllz-2, 0:2))
        allocate (Tdomain%specel(n)%sSimu(i)%Forces (0:ngllx-1, 0:nglly-1, 0:ngllz-1, 0:2))
        Tdomain%specel(n)%sSimu(i)%Veloc = 0.d0
        Tdomain%specel(n)%sSimu(i)%Accel = 0.d0
    enddo
    if (Tdomain%specel(n)%PML) then
        allocate (Tdomain%specel(n)%Lambda (0:ngllx-1, 0:nglly-1, 0:ngllz-1))
        allocate (Tdomain%specel(n)%Mu (0:ngllx-1, 0:nglly-1, 0:ngllz-1))
        allocate (Tdomain%specel(n)%DumpSx (0:ngllx-1, 0:nglly-1, 0:ngllz-1, 0:1))
        allocate (Tdomain%specel(n)%DumpSy (0:ngllx-1, 0:nglly-1, 0:ngllz-1, 0:1))
        allocate (Tdomain%specel(n)%DumpSz (0:ngllx-1, 0:nglly-1, 0:ngllz-1, 0:1))
        allocate (Tdomain%specel(n)%DumpMass (0:ngllx-1, 0:nglly-1, 0:ngllz-1, 0:2))
        allocate (Tdomain%specel(n)%DumpVx (1:ngllx-2, 1:nglly-2, 1:ngllz-2, 0:1))
        allocate (Tdomain%specel(n)%DumpVy (1:ngllx-2, 1:nglly-2, 1:ngllz-2, 0:1))
        allocate (Tdomain%specel(n)%DumpVz (1:ngllx-2, 1:nglly-2, 1:ngllz-2, 0:1))
        allocate (Tdomain%specel(n)%Acoeff (0:ngllx-1, 0:nglly-1, 0:ngllz-1, 0:35))
        do i = 0,Tdomain%nb_simu-1
            allocate (Tdomain%specel(n)%sSimu(i)%Diagonal_Stress (0:ngllx-1, 0:nglly-1, 0:ngllz-1, 0:2))
            allocate (Tdomain%specel(n)%sSimu(i)%Diagonal_Stress1 (0:ngllx-1, 0:nglly-1, 0:ngllz-1, 0:2))
            allocate (Tdomain%specel(n)%sSimu(i)%Diagonal_Stress2 (0:ngllx-1, 0:nglly-1, 0:ngllz-1, 0:2))
            allocate (Tdomain%specel(n)%sSimu(i)%Diagonal_Stress3 (0:ngllx-1, 0:nglly-1, 0:ngllz-1, 0:2))
            allocate (Tdomain%specel(n)%sSimu(i)%Residual_Stress (0:ngllx-1, 0:nglly-1, 0:ngllz-1, 0:2))
            allocate (Tdomain%specel(n)%sSimu(i)%Residual_Stress1 (0:ngllx-1, 0:nglly-1, 0:ngllz-1, 0:2))
            allocate (Tdomain%specel(n)%sSimu(i)%Residual_Stress2 (0:ngllx-1, 0:nglly-1, 0:ngllz-1, 0:2))
            allocate (Tdomain%specel(n)%sSimu(i)%Veloc1 (1:ngllx-2, 1:nglly-2,1:ngllz-2, 0:2))
            allocate (Tdomain%specel(n)%sSimu(i)%Veloc2 (1:ngllx-2, 1:nglly-2,1:ngllz-2, 0:2))
            allocate (Tdomain%specel(n)%sSimu(i)%Veloc3 (1:ngllx-2, 1:nglly-2,1:ngllz-2, 0:2))
            allocate (Tdomain%specel(n)%sSimu(i)%Forces1 (0:ngllx-1, 0:nglly-1, 0:ngllz-1, 0:2))
            allocate (Tdomain%specel(n)%sSimu(i)%Forces2 (0:ngllx-1, 0:nglly-1, 0:ngllz-1, 0:2))
            allocate (Tdomain%specel(n)%sSimu(i)%Forces3 (0:ngllx-1, 0:nglly-1, 0:ngllz-1, 0:2))
            Tdomain%specel(n)%sSimu(i)%Diagonal_Stress = 0.d0
            Tdomain%specel(n)%sSimu(i)%Diagonal_Stress1 = 0.d0
            Tdomain%specel(n)%sSimu(i)%Diagonal_Stress2 = 0.d0
            Tdomain%specel(n)%sSimu(i)%Diagonal_Stress3 = 0.d0
            Tdomain%specel(n)%sSimu(i)%Residual_Stress = 0.d0
            Tdomain%specel(n)%sSimu(i)%Residual_Stress1 = 0.d0
            Tdomain%specel(n)%sSimu(i)%Residual_Stress2 = 0.d0
            Tdomain%specel(n)%sSimu(i)%Veloc1 = 0.d0
            Tdomain%specel(n)%sSimu(i)%Veloc2 = 0.d0
            Tdomain%specel(n)%sSimu(i)%Veloc3 = 0.d0
        enddo
        if (Tdomain%curve) then
            allocate (Tdomain%specel(n)%Normales (0:2, 0:2))
            allocate (Tdomain%specel(n)%Inv_Normales (0:2, 0:2))
            do i = 0,Tdomain%nb_simu-1
                allocate (Tdomain%specel(n)%sSimu(i)%Residual_Stress3 (0:ngllx-1, 0:nglly-1, 0:ngllz-1, 0:2))
                Tdomain%specel(n)%sSimu(i)%Residual_Stress3 = 0.d0
            enddo
        endif
    else
        allocate (Tdomain%specel(n)%wgtx (0:ngllx-1))
        allocate (Tdomain%specel(n)%wgty (0:nglly-1))
        allocate (Tdomain%specel(n)%wgtz (0:ngllz-1))
!        allocate (Tdomain%specel(n)%Acoeff (0:ngllx-1, 0:nglly-1, 0:ngllz-1, 0:44))
        do i = 0,Tdomain%nb_simu-1
            allocate (Tdomain%specel(n)%sSimu(i)%Displ (1:ngllx-2, 1:nglly-2, 1:ngllz-2, 0:2))
            allocate (Tdomain%specel(n)%sSimu(i)%V0 (1:ngllx-2, 1:nglly-2, 1:ngllz-2, 0:2))
            Tdomain%specel(n)%sSimu(i)%Displ = 0.d0
            Tdomain%specel(n)%sSimu(i)%V0 = 0.d0
        enddo
        if (Tdomain%aniso) then
            allocate (Tdomain%specel(n)%Cij (0:20, 0:ngllx-1, 0:nglly-1, 0:ngllz-1))
            if (n_solid>0) then
                allocate (Tdomain%specel(n)%Lambda (0:ngllx-1, 0:nglly-1, 0:ngllz-1))
                allocate (Tdomain%specel(n)%Mu (0:ngllx-1, 0:nglly-1, 0:ngllz-1))
            endif
            if (Tdomain%adjoint) then
                allocate (Tdomain%specel(n)%save_TIparam (1:5, 0:ngllx-1, 0:nglly-1, 0:ngllz-1))
                allocate (Tdomain%specel(n)%Kern_rho (0:ngllx-1, 0:nglly-1, 0:ngllz-1))
                allocate (Tdomain%specel(n)%Kern_aniso (1:5, 0:ngllx-1, 0:nglly-1, 0:ngllz-1))
                Tdomain%specel(n)%Kern_rho = 0.d0
                Tdomain%specel(n)%Kern_aniso = 0.d0
                do i = 0,Tdomain%nb_simu-1
                    allocate (Tdomain%specel(n)%sSimu(i)%save_strain (0:ngllx-1, 0:nglly-1, 0:ngllz-1, 0:5))
                enddo
!                allocate (Tdomain%specel(n)%anomaly (0:ngllx-1, 0:nglly-1, 0:ngllz-1))
            endif
        else
            allocate (Tdomain%specel(n)%Lambda (0:ngllx-1, 0:nglly-1, 0:ngllz-1))
            allocate (Tdomain%specel(n)%Mu (0:ngllx-1, 0:nglly-1, 0:ngllz-1))
            if (Tdomain%adjoint) then
                allocate (Tdomain%specel(n)%Kern_rho (0:ngllx-1, 0:nglly-1, 0:ngllz-1))
                allocate (Tdomain%specel(n)%Kern_lambda (0:ngllx-1, 0:nglly-1, 0:ngllz-1))
                allocate (Tdomain%specel(n)%Kern_mu (0:ngllx-1, 0:nglly-1, 0:ngllz-1))
                Tdomain%specel(n)%Kern_rho = 0.d0
                Tdomain%specel(n)%Kern_lambda = 0.d0
                Tdomain%specel(n)%Kern_mu = 0.d0
                do i = 0,Tdomain%nb_simu-1
                    allocate (Tdomain%specel(n)%sSimu(i)%save_strain (0:ngllx-1, 0:nglly-1, 0:ngllz-1, 0:5))
                enddo
!                allocate (Tdomain%specel(n)%anomaly (0:ngllx-1, 0:nglly-1, 0:ngllz-1))
            endif
        endif
        if (n_solid>0) then
            allocate (Tdomain%specel(n)%Q (0:ngllx-1, 0:nglly-1, 0:ngllz-1))
            allocate (Tdomain%specel(n)%onemSbeta (0:ngllx-1, 0:nglly-1, 0:ngllz-1))
            allocate (Tdomain%specel(n)%factor_common_3 (0:n_solid-1, 0:ngllx-1, 0:nglly-1, 0:ngllz-1))
            allocate (Tdomain%specel(n)%alphaval_3 (0:n_solid-1, 0:ngllx-1, 0:nglly-1, 0:ngllz-1))
            allocate (Tdomain%specel(n)%betaval_3 (0:n_solid-1, 0:ngllx-1, 0:nglly-1, 0:ngllz-1))
            allocate (Tdomain%specel(n)%gammaval_3 (0:n_solid-1, 0:ngllx-1, 0:nglly-1, 0:ngllz-1))
            do i = 0,Tdomain%nb_simu-1
                allocate (Tdomain%specel(n)%sSimu(i)%epsilondev_xx_ (0:ngllx-1, 0:nglly-1, 0:ngllz-1))
                allocate (Tdomain%specel(n)%sSimu(i)%epsilondev_yy_ (0:ngllx-1, 0:nglly-1, 0:ngllz-1))
                allocate (Tdomain%specel(n)%sSimu(i)%epsilondev_xy_ (0:ngllx-1, 0:nglly-1, 0:ngllz-1))
                allocate (Tdomain%specel(n)%sSimu(i)%epsilondev_xz_ (0:ngllx-1, 0:nglly-1, 0:ngllz-1))
                allocate (Tdomain%specel(n)%sSimu(i)%epsilondev_yz_ (0:ngllx-1, 0:nglly-1, 0:ngllz-1))
                allocate (Tdomain%specel(n)%sSimu(i)%R_xx_ (0:n_solid-1, 0:ngllx-1, 0:nglly-1, 0:ngllz-1))
                allocate (Tdomain%specel(n)%sSimu(i)%R_yy_ (0:n_solid-1, 0:ngllx-1, 0:nglly-1, 0:ngllz-1))
                allocate (Tdomain%specel(n)%sSimu(i)%R_xy_ (0:n_solid-1, 0:ngllx-1, 0:nglly-1, 0:ngllz-1))
                allocate (Tdomain%specel(n)%sSimu(i)%R_xz_ (0:n_solid-1, 0:ngllx-1, 0:nglly-1, 0:ngllz-1))
                allocate (Tdomain%specel(n)%sSimu(i)%R_yz_ (0:n_solid-1, 0:ngllx-1, 0:nglly-1, 0:ngllz-1))
                Tdomain%specel(n)%sSimu(i)%epsilondev_xx_ = 0.d0
                Tdomain%specel(n)%sSimu(i)%epsilondev_yy_ = 0.d0
                Tdomain%specel(n)%sSimu(i)%epsilondev_xy_ = 0.d0
                Tdomain%specel(n)%sSimu(i)%epsilondev_xz_ = 0.d0
                Tdomain%specel(n)%sSimu(i)%epsilondev_yz_ = 0.d0
                Tdomain%specel(n)%sSimu(i)%R_xx_ = 0.d0
                Tdomain%specel(n)%sSimu(i)%R_yy_ = 0.d0
                Tdomain%specel(n)%sSimu(i)%R_xy_ = 0.d0
                Tdomain%specel(n)%sSimu(i)%R_xz_ = 0.d0
                Tdomain%specel(n)%sSimu(i)%R_yz_ = 0.d0
            enddo
        endif
    endif
enddo

do n = 0,Tdomain%n_face-1
    ngll1 = Tdomain%sFace(n)%ngll1
    ngll2 = Tdomain%sFace(n)%ngll2
    allocate (Tdomain%sFace(n)%MassMat (1:ngll1-2, 1:ngll2-2))
    Tdomain%sFace(n)%MassMat = 0.d0
    allocate (Tdomain%sFace(n)%sSimu(0:Tdomain%nb_simu-1))
    do i = 0,Tdomain%nb_simu-1
        allocate (Tdomain%sFace(n)%sSimu(i)%Veloc (1:ngll1-2, 1:ngll2-2, 0:2))
        allocate (Tdomain%sFace(n)%sSimu(i)%Accel (1:ngll1-2, 1:ngll2-2, 0:2))
        allocate (Tdomain%sFace(n)%sSimu(i)%Forces (1:ngll1-2, 1:ngll2-2, 0:2))
        Tdomain%sFace(n)%sSimu(i)%Veloc = 0.d0
        Tdomain%sFace(n)%sSimu(i)%Accel = 0.d0
        Tdomain%sFace(n)%sSimu(i)%Forces = 0.d0
    enddo
    if (Tdomain%sFace(n)%ocean) then
        allocate (Tdomain%sFace(n)%Mocean (1:ngll1-2, 1:ngll2-2))
        Tdomain%sFace(n)%Mocean = 0.d0
        allocate (Tdomain%sFace(n)%verticale (1:ngll1-2, 1:ngll2-2, 0:2, 0:2))
        allocate (Tdomain%sFace(n)%M33ocean (1:ngll1-2, 1:ngll2-2, 0:2, 0:2))
    endif
    if (Tdomain%sFace(n)%PML) then 
        allocate (Tdomain%sFace(n)%DumpMass (1:ngll1-2, 1:ngll2-2, 0:2))
        allocate (Tdomain%sFace(n)%DumpVx (1:ngll1-2, 1:ngll2-2, 0:1)) 
        allocate (Tdomain%sFace(n)%DumpVy (1:ngll1-2, 1:ngll2-2, 0:1)) 
        allocate (Tdomain%sFace(n)%DumpVz (1:ngll1-2, 1:ngll2-2, 0:1)) 
        Tdomain%sFace(n)%DumpMass = 0.d0
        do i = 0,Tdomain%nb_simu-1
            allocate (Tdomain%sFace(n)%sSimu(i)%Forces1 (1:ngll1-2, 1:ngll2-2, 0:2))
            allocate (Tdomain%sFace(n)%sSimu(i)%Forces2 (1:ngll1-2, 1:ngll2-2, 0:2))
            allocate (Tdomain%sFace(n)%sSimu(i)%Forces3 (1:ngll1-2, 1:ngll2-2, 0:2)) 
            allocate (Tdomain%sFace(n)%sSimu(i)%Veloc1 (1:ngll1-2, 1:ngll2-2, 0:2))
            allocate (Tdomain%sFace(n)%sSimu(i)%Veloc2 (1:ngll1-2, 1:ngll2-2, 0:2))
            allocate (Tdomain%sFace(n)%sSimu(i)%Veloc3 (1:ngll1-2, 1:ngll2-2, 0:2)) 
            Tdomain%sFace(n)%sSimu(i)%Veloc1 = 0.d0
            Tdomain%sFace(n)%sSimu(i)%Veloc2 = 0.d0
            Tdomain%sFace(n)%sSimu(i)%Veloc3 = 0.d0
        enddo
    else
        do i = 0,Tdomain%nb_simu-1
            allocate (Tdomain%sFace(n)%sSimu(i)%Displ (1:ngll1-2, 1:ngll2-2, 0:2))
            allocate (Tdomain%sFace(n)%sSimu(i)%V0 (1:ngll1-2, 1:ngll2-2, 0:2))
            Tdomain%sFace(n)%sSimu(i)%Displ = 0.d0
            Tdomain%sFace(n)%sSimu(i)%V0 = 0.d0
        enddo
    endif
enddo

do n = 0,Tdomain%n_edge-1
    ngll = Tdomain%sEdge(n)%ngll
    allocate (Tdomain%sEdge(n)%MassMat (1:ngll-2))
    Tdomain%sEdge(n)%MassMat = 0.d0
    allocate (Tdomain%sEdge(n)%sSimu(0:Tdomain%nb_simu-1))
    do i = 0,Tdomain%nb_simu-1
        allocate (Tdomain%sEdge(n)%sSimu(i)%Veloc (1:ngll-2, 0:2))
        allocate (Tdomain%sEdge(n)%sSimu(i)%Accel (1:ngll-2, 0:2))
        allocate (Tdomain%sEdge(n)%sSimu(i)%Forces (1:ngll-2, 0:2))
        Tdomain%sEdge(n)%sSimu(i)%Veloc = 0.d0
        Tdomain%sEdge(n)%sSimu(i)%Accel = 0.d0
        Tdomain%sEdge(n)%sSimu(i)%Forces = 0.d0
    enddo
    if (Tdomain%sEdge(n)%ocean) then
        allocate (Tdomain%sEdge(n)%Mocean (1:ngll-2))
        Tdomain%sEdge(n)%Mocean = 0.d0
        allocate (Tdomain%sEdge(n)%verticale (1:ngll-2, 0:2, 0:2))
        allocate (Tdomain%sEdge(n)%M33ocean (1:ngll-2, 0:2, 0:2))
    endif
    if (Tdomain%sEdge(n)%PML) then
        allocate (Tdomain%sEdge(n)%DumpMass (1:ngll-2, 0:2))
        allocate (Tdomain%sEdge(n)%DumpVx (1:ngll-2, 0:1))
        allocate (Tdomain%sEdge(n)%DumpVy (1:ngll-2, 0:1))
        allocate (Tdomain%sEdge(n)%DumpVz (1:ngll-2, 0:1))
        Tdomain%sEdge(n)%DumpMass = 0.d0
        do i = 0,Tdomain%nb_simu-1
            allocate (Tdomain%sEdge(n)%sSimu(i)%Forces1 (1:ngll-2, 0:2))
            allocate (Tdomain%sEdge(n)%sSimu(i)%Forces2 (1:ngll-2, 0:2))
            allocate (Tdomain%sEdge(n)%sSimu(i)%Forces3 (1:ngll-2, 0:2))
            allocate (Tdomain%sEdge(n)%sSimu(i)%Veloc1 (1:ngll-2, 0:2))
            allocate (Tdomain%sEdge(n)%sSimu(i)%Veloc2 (1:ngll-2, 0:2))
            allocate (Tdomain%sEdge(n)%sSimu(i)%Veloc3 (1:ngll-2, 0:2))
            Tdomain%sEdge(n)%sSimu(i)%Veloc1 = 0.d0
            Tdomain%sEdge(n)%sSimu(i)%Veloc2 = 0.d0
            Tdomain%sEdge(n)%sSimu(i)%Veloc3 = 0.d0
        enddo
    else
        do i = 0,Tdomain%nb_simu-1
            allocate (Tdomain%sEdge(n)%sSimu(i)%Displ (1:ngll-2, 0:2))
            allocate (Tdomain%sEdge(n)%sSimu(i)%V0 (1:ngll-2, 0:2))
            Tdomain%sEdge(n)%sSimu(i)%Displ = 0.d0
            Tdomain%sEdge(n)%sSimu(i)%V0 = 0.d0
        enddo
    endif
enddo

do n = 0,Tdomain%n_vertex-1
    Tdomain%sVertex(n)%MassMat = 0.d0
    allocate (Tdomain%sVertex(n)%sSimu(0:Tdomain%nb_simu-1))
    do i = 0,Tdomain%nb_simu-1
        allocate (Tdomain%sVertex(n)%sSimu(i)%Veloc (0:2))
        allocate (Tdomain%sVertex(n)%sSimu(i)%Accel (0:2))
        allocate (Tdomain%sVertex(n)%sSimu(i)%Forces (0:2))
        Tdomain%sVertex(n)%sSimu(i)%Veloc = 0.d0
        Tdomain%sVertex(n)%sSimu(i)%Accel = 0.d0
        Tdomain%sVertex(n)%sSimu(i)%Forces = 0.d0
    enddo
    if (Tdomain%sVertex(n)%ocean) then
        Tdomain%sVertex(n)%Mocean = 0.d0
        allocate (Tdomain%sVertex(n)%verticale (0:2, 0:2))
        allocate (Tdomain%sVertex(n)%M33ocean( 0:2, 0:2))
    endif
    if (Tdomain%sVertex(n)%PML) then
        allocate (Tdomain%sVertex(n)%DumpMass (0:2))
        allocate (Tdomain%sVertex(n)%DumpVx (0:1))
        allocate (Tdomain%sVertex(n)%DumpVy (0:1))
        allocate (Tdomain%sVertex(n)%DumpVz (0:1))
        Tdomain%sVertex(n)%DumpMass = 0.d0
        do i = 0,Tdomain%nb_simu-1
            allocate (Tdomain%sVertex(n)%sSimu(i)%Forces1 (0:2))
            allocate (Tdomain%sVertex(n)%sSimu(i)%Forces2 (0:2))
            allocate (Tdomain%sVertex(n)%sSimu(i)%Forces3 (0:2))
            allocate (Tdomain%sVertex(n)%sSimu(i)%Veloc1 (0:2))
            allocate (Tdomain%sVertex(n)%sSimu(i)%Veloc2 (0:2))
            allocate (Tdomain%sVertex(n)%sSimu(i)%Veloc3 (0:2))
            Tdomain%sVertex(n)%sSimu(i)%Veloc1 = 0.d0
            Tdomain%sVertex(n)%sSimu(i)%Veloc2 = 0.d0
            Tdomain%sVertex(n)%sSimu(i)%Veloc3 = 0.d0
        enddo
    else
        do i = 0,Tdomain%nb_simu-1
            allocate (Tdomain%sVertex(n)%sSimu(i)%Displ (0:2))
            allocate (Tdomain%sVertex(n)%sSimu(i)%V0 (0:2))
            Tdomain%sVertex(n)%sSimu(i)%Displ = 0.d0
            Tdomain%sVertex(n)%sSimu(i)%V0 = 0.d0
        enddo
    endif
enddo

do n = 0,Tdomain%n_proc-1
    ngll = 0
    ngllPML = 0
    ngllocean = 0
    do i = 0,Tdomain%sComm(n)%nb_faces-1
        nf = Tdomain%sComm(n)%faces(i)
        ngll1 = Tdomain%sFace(nf)%ngll1; ngll2 = Tdomain%sFace(nf)%ngll2 
        ngll = ngll + (ngll1-2)*(ngll2-2)
        if (Tdomain%sFace(nf)%PML)   ngllPML = ngllPML + (ngll1-2)*(ngll2-2)
    enddo
    do i = 0,Tdomain%sComm(n)%nb_edges-1
        ne = Tdomain%sComm(n)%edges(i)
        ngll = ngll + Tdomain%sEdge(ne)%ngll-2
        if (Tdomain%sEdge(ne)%PML)   ngllPML = ngllPML + Tdomain%sEdge(ne)%ngll-2
        if (Tdomain%sEdge(ne)%ocean)   ngllocean = ngllocean + Tdomain%sEdge(ne)%ngll-2
    enddo
    do i = 0,Tdomain%sComm(n)%nb_vertices-1
        nv = Tdomain%sComm(n)%vertices(i) 
        ngll = ngll + 1
        if (Tdomain%sVertex(nv)%PML)   ngllPML = ngllPML + 1
        if (Tdomain%sVertex(nv)%ocean)   ngllocean = ngllocean + 1
    enddo
    if (ngll>0) then
        allocate (Tdomain%sComm(n)%Give (0:ngll-1))
        allocate (Tdomain%sComm(n)%Take (0:ngll-1))
        allocate (Tdomain%sComm(n)%GiveForces (0:ngll-1, 0:2))
        allocate (Tdomain%sComm(n)%TakeForces (0:ngll-1, 0:2))
    endif
    if (ngllocean>0) then
        allocate (Tdomain%sComm(n)%Giveocean (0:ngllocean-1))
        allocate (Tdomain%sComm(n)%Takeocean (0:ngllocean-1))
    endif
    if (ngllPML>0) then
        allocate (Tdomain%sComm(n)%GivePML (0:ngllPML-1, 0:2))
        allocate (Tdomain%sComm(n)%TakePML (0:ngllPML-1, 0:2))
        allocate (Tdomain%sComm(n)%GiveForcesPML (0:ngllPML-1, 1:3, 0:2))
        allocate (Tdomain%sComm(n)%TakeForcesPML (0:ngllPML-1, 1:3, 0:2))
    endif
    Tdomain%sComm(n)%ngll = ngll
    Tdomain%sComm(n)%ngllPML = ngllPML
    Tdomain%sComm(n)%ngllocean = ngllocean
enddo


return
end subroutine allocate_domain
