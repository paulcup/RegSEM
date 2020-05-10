module selements


use ssimus

type :: element

type(simu), dimension (:), pointer :: sSimu

integer :: mat_index, cov_index, moho_position, topo_position, mirror_position, ngllx,nglly,ngllz, indx,indy,indz
integer, dimension (:), pointer :: Control_nodes
integer, dimension (0:5) :: Near_Faces
integer, dimension (0:11) :: Near_Edges
integer, dimension (0:7) :: Near_Vertices
integer, dimension (:,:,:), pointer :: Iglobnum, win_mirror

doubleprecision, dimension (:), pointer :: wgtx, wgty, wgtz
doubleprecision, dimension (:,:,:), pointer :: Jacob, Density, Lambda, Mu, Q, onemSbeta, MassMat, &
                                               Kern_rho, Kern_lambda, Kern_mu, anomaly
doubleprecision, dimension (:,:,:,:), pointer :: Cij, Acoeff, factor_common_3, alphaval_3,betaval_3,gammaval_3, &
                                                 Kern_aniso, save_TIparam
doubleprecision, dimension (:,:,:,:,:), pointer :: InvGrad

! PML
logical :: PML
doubleprecision, dimension (:,:), pointer :: Normales, Inv_Normales
doubleprecision, dimension (:,:,:,:), pointer :: DumpSx,DumpSy,DumpSz
doubleprecision, dimension (:,:,:,:), pointer :: DumpVx,DumpVy,DumpVz, DumpMass

! Ocean
logical :: ocean
doubleprecision, dimension (:,:), pointer :: Mocean, hocean

! Random Forces
logical :: random_t = .false.
doubleprecision, dimension(:,:,:,:), pointer :: random_coeff

end type


contains 

! ############################################################
subroutine Prediction_Elem_Veloc (Elem, alpha, bega, gam1, dt, i_simu)


implicit none

type (Element), intent (INOUT) :: Elem
integer, intent (IN) :: i_simu
doubleprecision, intent (IN) :: alpha, bega, gam1, dt

integer :: ngllx, nglly, ngllz

   
ngllx = Elem%ngllx ; nglly =Elem%nglly;  ngllz = Elem%ngllz

Elem%sSimu(i_simu)%Forces(1:ngllx-2,1:nglly-2,1:ngllz-2,0:2) = Elem%sSimu(i_simu)%Displ + &
                                                               dt * Elem%sSimu(i_simu)%Veloc + &
                                                               dt**2 * (0.5 - bega) * Elem%sSimu(i_simu)%Accel
Elem%sSimu(i_simu)%V0 = Elem%sSimu(i_simu)%Veloc
Elem%sSimu(i_simu)%Forces(1:ngllx-2,1:nglly-2,1:ngllz-2,0:2) = &
   alpha * Elem%sSimu(i_simu)%Forces(1:ngllx-2,1:nglly-2,1:ngllz-2,0:2) + (1-alpha) * Elem%sSimu(i_simu)%Displ 


return
end subroutine Prediction_Elem_Veloc

! ###########################################################
subroutine Correction_Elem_Veloc (Elem, i_simu, bega, gam1, dt)


implicit none

type (Element), intent (INOUT) :: Elem
doubleprecision, intent (IN) :: bega, gam1, dt
integer, intent (IN) :: i_simu

integer :: ngllx, nglly, ngllz, i


ngllx = Elem%ngllx ;  nglly = Elem%nglly ;  ngllz = Elem%ngllz
do i = 0,2
    Elem%sSimu(i_simu)%Forces(1:ngllx-2,1:nglly-2,1:ngllz-2,i) = Elem%MassMat(:,:,:) * &
                              Elem%sSimu(i_simu)%Forces(1:ngllx-2,1:nglly-2,1:ngllz-2,i)
enddo
Elem%sSimu(i_simu)%Veloc(:,:,:,:) = Elem%sSimu(i_simu)%V0(:,:,:,:) + &
                                    dt * Elem%sSimu(i_simu)%Forces(1:ngllx-2,1:nglly-2,1:ngllz-2,:)
Elem%sSimu(i_simu)%Accel = Elem%sSimu(i_simu)%Accel + &
                           gam1/dt * (Elem%sSimu(i_simu)%Veloc-Elem%sSimu(i_simu)%V0)
Elem%sSimu(i_simu)%Displ = Elem%sSimu(i_simu)%Displ + &
                           bega*dt * (Elem%sSimu(i_simu)%Veloc+Elem%sSimu(i_simu)%V0)


return
end subroutine Correction_Elem_Veloc

! ###########################################################
subroutine  compute_InternalForces_Elem (Elem, i_simu, hprimex, hTprimex, hprimey, hTprimey, hprimez, hTprimez)


implicit none

type (Element), intent (INOUT) :: Elem
doubleprecision, dimension (0:Elem%ngllx-1, 0:Elem%ngllx-1), intent (IN) :: hprimex, htprimex
doubleprecision, dimension (0:Elem%nglly-1, 0:Elem%nglly-1), intent (IN) :: hprimey, htprimey
doubleprecision, dimension (0:Elem%ngllz-1, 0:Elem%ngllz-1), intent (IN) :: hprimez, hTprimez
integer, intent (IN) :: i_simu

integer :: n_z, m1, m2, m3

doubleprecision, dimension (0:Elem%ngllx-1, 0:Elem%nglly-1) :: s0z
doubleprecision, dimension (0:Elem%nglly-1, 0:Elem%ngllz-1) :: s0x
doubleprecision, dimension (0:Elem%ngllx-1, 0:Elem%nglly-1, 0:Elem%ngllz-1) :: dUx_dxi, dUx_deta, dUx_dzeta, &
                                                                    dUy_dxi, dUy_deta, dUy_dzeta, &
                                                                    dUz_dxi, dUz_deta, dUz_dzeta, &
                                                                    t1, s0, Uxloc, Uyloc, Uzloc


m1 = Elem%ngllx;   m2 = Elem%nglly;   m3 = Elem%ngllz

call DGEMM ('N', 'N', m1, m2*m3, m1, 1.d0, htprimex, m1, Elem%sSimu(i_simu)%Forces(0,0,0,0), m1, 0.d0, dUx_dxi, m1)
do n_z = 0,Elem%ngllz-1
   call DGEMM ('N', 'N', m1, m2, m2, 1.d0, Elem%sSimu(i_simu)%Forces(0,0,n_z,0), m1, hprimey, m2, 0.d0, dUx_deta(0,0,n_z), m1)
enddo
call DGEMM ('N', 'N', m1*m2, m3, m3, 1.d0, Elem%sSimu(i_simu)%Forces(0,0,0,0), m1*m2, hprimez, m3, 0.d0, dUx_dzeta, m1*m2)

call DGEMM ('N', 'N', m1, m2*m3, m1, 1.d0, htprimex, m1, Elem%sSimu(i_simu)%Forces(0,0,0,1), m1, 0.d0, dUy_dxi, m1)
do n_z = 0,Elem%ngllz-1
   call DGEMM ('N', 'N', m1, m2, m2, 1.d0, Elem%sSimu(i_simu)%Forces(0,0,n_z,1), m1, hprimey, m2, 0.d0, dUy_deta(0,0,n_z), m1)
enddo
call DGEMM ('N', 'N', m1*m2, m3, m3, 1.d0, Elem%sSimu(i_simu)%Forces(0,0,0,1), m1*m2, hprimez, m3, 0.d0, dUy_dzeta, m1*m2)

call DGEMM ('N', 'N', m1, m2*m3, m1, 1.d0, htprimex, m1, Elem%sSimu(i_simu)%Forces(0,0,0,2), m1, 0.d0, dUz_dxi, m1)
do n_z = 0,Elem%ngllz-1
   call DGEMM ('N', 'N', m1, m2, m2, 1.d0, Elem%sSimu(i_simu)%Forces(0,0,n_z,2), m1, hprimey, m2, 0.d0, dUz_deta(0,0,n_z), m1)
enddo
call DGEMM ('N', 'N', m1*m2, m3, m3, 1.d0, Elem%sSimu(i_simu)%Forces(0,0,0,2), m1*m2, hprimez, m3, 0.d0, dUz_dzeta, m1*m2)

t1 = Elem%Acoeff(:,:,:,0)*dUx_dxi + Elem%Acoeff(:,:,:,1)*dUx_deta + Elem%Acoeff(:,:,:,2)*dUx_dzeta + &
     Elem%Acoeff(:,:,:,3)*dUy_dxi + Elem%Acoeff(:,:,:,4)*dUy_deta + Elem%Acoeff(:,:,:,5)*dUy_dzeta +  &
     Elem%Acoeff(:,:,:,6)*dUz_dxi + Elem%Acoeff(:,:,:,7)*dUz_deta + Elem%Acoeff(:,:,:,8)*dUz_dzeta 

call DGEMM ('N', 'N', m1, m2*m3, m1, 1.d0, hprimex, m1, t1(0,0,0), m1, 0.d0, Uxloc, m1)

t1 = Elem%Acoeff(:,:,:,1)*dUx_dxi + Elem%Acoeff(:,:,:,9)*dUx_deta + Elem%Acoeff(:,:,:,10)*dUx_dzeta + &
     Elem%Acoeff(:,:,:,11)*dUy_dxi + Elem%Acoeff(:,:,:,12)*dUy_deta + Elem%Acoeff(:,:,:,13)*dUy_dzeta + &
     Elem%Acoeff(:,:,:,14)*dUz_dxi + Elem%Acoeff(:,:,:,15)*dUz_deta + Elem%Acoeff(:,:,:,16)*dUz_dzeta

do n_z = 0,Elem%ngllz-1
     call DGEMM ('N', 'N', m1, m2, m2, 1.d0, t1(0,0,n_z), m1, htprimey, m2, 0.d0, s0(0,0,n_z), m1)
enddo
Uxloc = s0 + Uxloc

t1 = Elem%Acoeff(:,:,:,2)*dUx_dxi + Elem%Acoeff(:,:,:,10)*dUx_deta + Elem%Acoeff(:,:,:,17)*dUx_dzeta + &
     Elem%Acoeff(:,:,:,18)*dUy_dxi + Elem%Acoeff(:,:,:,19)*dUy_deta + Elem%Acoeff(:,:,:,20)*dUy_dzeta +&
     Elem%Acoeff(:,:,:,21)*dUz_dxi + Elem%Acoeff(:,:,:,22)*dUz_deta + Elem%Acoeff(:,:,:,23)*dUz_dzeta

call DGEMM ('N', 'N', m1*m2, m3, m3, 1.d0, t1(0,0,0), m1*m2, htprimez, m3, 0.d0, s0, m1*m2)
Uxloc = s0 + Uxloc

t1 = Elem%Acoeff(:,:,:,3)*dUx_dxi + Elem%Acoeff(:,:,:,11)*dUx_deta + Elem%Acoeff(:,:,:,18)*dUx_dzeta + &
     Elem%Acoeff(:,:,:,24)*dUy_dxi + Elem%Acoeff(:,:,:,25)*dUy_deta + Elem%Acoeff(:,:,:,26)*dUy_dzeta + &
     Elem%Acoeff(:,:,:,27)*dUz_dxi + Elem%Acoeff(:,:,:,28)*dUz_deta + Elem%Acoeff(:,:,:,29)*dUz_dzeta

call DGEMM ('N', 'N', m1, m2*m3, m1, 1.d0, hprimex, m1, t1(0,0,0), m1, 0.d0, Uyloc, m1)

t1 = Elem%Acoeff(:,:,:,4)*dUx_dxi + Elem%Acoeff(:,:,:,12)*dUx_deta + Elem%Acoeff(:,:,:,19)*dUx_dzeta + &
     Elem%Acoeff(:,:,:,25)*dUy_dxi + Elem%Acoeff(:,:,:,30)*dUy_deta + Elem%Acoeff(:,:,:,31)*dUy_dzeta + &
     Elem%Acoeff(:,:,:,32)*dUz_dxi + Elem%Acoeff(:,:,:,33)*dUz_deta + Elem%Acoeff(:,:,:,34)*dUz_dzeta 

do n_z = 0,Elem%ngllz-1
     call DGEMM ('N', 'N', m1, m2, m2, 1.d0, t1(0,0,n_z), m1, htprimey, m2, 0.d0, s0(0,0,n_z), m1)
enddo
Uyloc = s0 + Uyloc

t1 = Elem%Acoeff(:,:,:,5)*dUx_dxi + Elem%Acoeff(:,:,:,13)*dUx_deta + Elem%Acoeff(:,:,:,20)* dUx_dzeta +&
     Elem%Acoeff(:,:,:,26)*dUy_dxi + Elem%Acoeff(:,:,:,31)*dUy_deta + Elem%Acoeff(:,:,:,35)*dUy_dzeta + &
     Elem%Acoeff(:,:,:,36)*dUz_dxi + Elem%Acoeff(:,:,:,37)*dUz_deta + Elem%Acoeff(:,:,:,38)*dUz_dzeta

call DGEMM ('N', 'N', m1*m2, m3, m3, 1.d0, t1(0,0,0), m1*m2, htprimez ,m3, 0.d0, s0, m1*m2)
Uyloc = s0 + Uyloc

t1 = Elem%Acoeff(:,:,:,6)*dUx_dxi + Elem%Acoeff(:,:,:,14)*dUx_deta + Elem%Acoeff(:,:,:,21)*dUx_dzeta + &
     Elem%Acoeff(:,:,:,27)*dUy_dxi + Elem%Acoeff(:,:,:,32)*dUy_deta + Elem%Acoeff(:,:,:,36)*dUy_dzeta + &
     Elem%Acoeff(:,:,:,39)*dUz_dxi + Elem%Acoeff(:,:,:,40)*dUz_deta + Elem%Acoeff(:,:,:,41)*dUz_dzeta

call DGEMM ('N', 'N', m1, m2*m3, m1, 1.d0, hprimex, m1, t1(0,0,0), m1, 0.d0, Uzloc, m1)

t1 = Elem%Acoeff(:,:,:,7)*dUx_dxi + Elem%Acoeff(:,:,:,15)*dUx_deta + Elem%Acoeff(:,:,:,22)*dUx_dzeta + &
     Elem%Acoeff(:,:,:,28)*dUy_dxi + Elem%Acoeff(:,:,:,33)*dUy_deta + Elem%Acoeff(:,:,:,37)*dUy_dzeta + &
     Elem%Acoeff(:,:,:,40)*dUz_dxi + Elem%Acoeff(:,:,:,42)*dUz_deta + Elem%Acoeff(:,:,:,43)*dUz_dzeta
      
do n_z = 0,Elem%ngllz-1
     call DGEMM ('N', 'N', m1, m2, m2, 1.d0, t1(0,0,n_z), m1, htprimey, m2, 0.d0, s0(0,0,n_z), m1)
enddo
Uzloc = s0 + Uzloc
     
t1 = Elem%Acoeff(:,:,:,8)*dUx_dxi + Elem%Acoeff(:,:,:,16)*dUx_deta + Elem%Acoeff(:,:,:,23)*dUx_dzeta + &
     Elem%Acoeff(:,:,:,29)*dUy_dxi + Elem%Acoeff(:,:,:,34)*dUy_deta + Elem%Acoeff(:,:,:,38)*dUy_dzeta + &
     Elem%Acoeff(:,:,:,41)*dUz_dxi + Elem%Acoeff(:,:,:,43)*dUz_deta + Elem%Acoeff(:,:,:,44)*dUz_dzeta        

call DGEMM ('N', 'N', m1*m2, m3, m3, 1.d0, t1(0,0,0), m1*m2, htprimez, m3, 0.d0, s0, m1*m2)
Uzloc = Uzloc + s0

Elem%sSimu(i_simu)%Forces(:,:,:,0) = Uxloc
Elem%sSimu(i_simu)%Forces(:,:,:,1) = Uyloc
Elem%sSimu(i_simu)%Forces(:,:,:,2) = Uzloc


return
end subroutine compute_InternalForces_Elem

! ###########################################################
subroutine Prediction_Elem_PML_Veloc (Elem, bega, dt, hTprimex, hprimey, hprimez, Sim)


implicit none

type (Element), intent (INOUT) :: Elem
type (Simu), intent (INOUT) :: Sim
doubleprecision, intent (IN) :: bega, dt
doubleprecision, dimension (0:Elem%ngllx-1, 0:Elem%ngllx-1), intent (IN) :: hTprimex
doubleprecision, dimension (0:Elem%nglly-1, 0:Elem%nglly-1), intent (IN) :: hprimey
doubleprecision, dimension (0:Elem%ngllz-1, 0:Elem%ngllz-1), intent (IN) :: hprimez

doubleprecision, dimension (0:Elem%ngllx-1, 0:Elem%nglly-1, 0:Elem%ngllz-1) :: dVx_dxi,dVx_deta,dVx_dzeta, & 
                                                                               dVy_dxi,dVy_deta,dVy_dzeta, &
                                                                               dVz_dxi,dVz_deta,dVz_dzeta
integer :: m1,m2,m3, n_z


m1 = Elem%ngllx;  m2 = Elem%nglly;  m3= Elem%ngllz

Sim%Forces(1:m1-2, 1:m2-2, 1:m3-2, 0:2) = Sim%Veloc(:,:,:,:) + dt * (0.5-bega) * Sim%Accel(:,:,:,:)

call DGEMM ( 'N', 'N', m1, m2*m3, m1, 1.d0, htprimex, m1, Sim%Forces(0,0,0,0) ,m1, 0.d0, dVx_dxi, m1 )
do n_z = 0,Elem%ngllz-1
   call DGEMM ( 'N', 'N', m1, m2, m2, 1.d0, Sim%Forces(0,0,n_z,0), m1, hprimey ,m2, 0.d0, dVx_deta(0,0,n_z), m1 )
enddo
call DGEMM ( 'N', 'N', m1*m2, m3, m3, 1.d0, Sim%Forces(0,0,0,0), m1*m2, hprimez ,m3, 0.d0, dVx_dzeta, m1*m2 )

call DGEMM ( 'N', 'N', m1, m2*m3, m1, 1.d0, htprimex, m1, Sim%Forces(0,0,0,1) ,m1, 0.d0, dVy_dxi, m1 )
do n_z = 0,Elem%ngllz-1
   call DGEMM ( 'N', 'N', m1, m2, m2, 1.d0, Sim%Forces(0,0,n_z,1), m1, hprimey ,m2, 0.d0, dVy_deta(0,0,n_z), m1 )
enddo
call DGEMM ( 'N', 'N', m1*m2, m3, m3, 1.d0, Sim%Forces(0,0,0,1), m1*m2, hprimez ,m3, 0.d0, dVy_dzeta, m1*m2 )

call DGEMM ( 'N', 'N', m1, m2*m3, m1, 1.d0, htprimex, m1, Sim%Forces(0,0,0,2) ,m1, 0.d0, dVz_dxi, m1 )
do n_z = 0,Elem%ngllz-1
   call DGEMM ( 'N', 'N', m1, m2, m2, 1.d0, Sim%Forces(0,0,n_z,2), m1, hprimey ,m2, 0.d0, dVz_deta(0,0,n_z), m1)
enddo
call DGEMM ( 'N', 'N', m1*m2, m3, m3, 1.d0, Sim%Forces(0,0,0,2), m1*m2, hprimez ,m3, 0.d0, dVz_dzeta, m1*m2 )


Sim%Diagonal_Stress1 (:,:,:,0) = Elem%DumpSx(:,:,:,0) * Sim%Diagonal_Stress1 (:,:,:,0) + &
   Elem%DumpSx(:,:,:,1) * Dt * (Elem%Acoeff(:,:,:,0) * dVx_dxi + Elem%Acoeff(:,:,:,1) * dVx_deta + Elem%Acoeff(:,:,:,2) * dVx_dzeta)
Sim%Diagonal_Stress2 (:,:,:,0) = Elem%DumpSy(:,:,:,0) * Sim%Diagonal_Stress2 (:,:,:,0) + &
   Elem%DumpSy(:,:,:,1) * Dt * (Elem%Acoeff(:,:,:,3) * dVy_dxi + Elem%Acoeff(:,:,:,4) * dVy_deta + Elem%Acoeff(:,:,:,5) * dVy_dzeta)
Sim%Diagonal_Stress3 (:,:,:,0) = Elem%DumpSz(:,:,:,0) * Sim%Diagonal_Stress3 (:,:,:,0) + &
   Elem%DumpSz(:,:,:,1) * Dt * (Elem%Acoeff(:,:,:,6) * dVz_dxi + Elem%Acoeff(:,:,:,7) * dVz_deta + Elem%Acoeff(:,:,:,8) * dVz_dzeta)

Sim%Diagonal_Stress1 (:,:,:,1) = Elem%DumpSx(:,:,:,0) * Sim%Diagonal_Stress1 (:,:,:,1) + &
   Elem%DumpSx(:,:,:,1) * Dt * (Elem%Acoeff(:,:,:,9) * dVx_dxi + Elem%Acoeff(:,:,:,10) * dVx_deta + Elem%Acoeff(:,:,:,11) * dVx_dzeta)
Sim%Diagonal_Stress2 (:,:,:,1) = Elem%DumpSy(:,:,:,0) * Sim%Diagonal_Stress2 (:,:,:,1) + &
   Elem%DumpSy(:,:,:,1) * Dt * (Elem%Acoeff(:,:,:,12) * dVy_dxi + Elem%Acoeff(:,:,:,13) * dVy_deta + Elem%Acoeff(:,:,:,14) * dVy_dzeta)
Sim%Diagonal_Stress3 (:,:,:,1) = Elem%DumpSz(:,:,:,0) * Sim%Diagonal_Stress3 (:,:,:,1) + &
   Elem%DumpSz(:,:,:,1) * Dt * (Elem%Acoeff(:,:,:,6) * dVz_dxi + Elem%Acoeff(:,:,:,7) * dVz_deta + Elem%Acoeff(:,:,:,8) * dVz_dzeta)

Sim%Diagonal_Stress1 (:,:,:,2) = Elem%DumpSx(:,:,:,0) * Sim%Diagonal_Stress1 (:,:,:,2) + &
   Elem%DumpSx(:,:,:,1) * Dt * (Elem%Acoeff(:,:,:,9) * dVx_dxi + Elem%Acoeff(:,:,:,10) * dVx_deta + Elem%Acoeff(:,:,:,11) * dVx_dzeta)
Sim%Diagonal_Stress2 (:,:,:,2) = Elem%DumpSy(:,:,:,0) * Sim%Diagonal_Stress2 (:,:,:,2) + &
   Elem%DumpSy(:,:,:,1) * Dt * (Elem%Acoeff(:,:,:,3) * dVy_dxi + Elem%Acoeff(:,:,:,4) * dVy_deta + Elem%Acoeff(:,:,:,5) * dVy_dzeta)
Sim%Diagonal_Stress3 (:,:,:,2) = Elem%DumpSz(:,:,:,0) * Sim%Diagonal_Stress3 (:,:,:,2) + &
   Elem%DumpSz(:,:,:,1) * Dt * (Elem%Acoeff(:,:,:,15) * dVz_dxi + Elem%Acoeff(:,:,:,16) * dVz_deta + Elem%Acoeff(:,:,:,17) * dVz_dzeta)

Sim%Diagonal_Stress = Sim%Diagonal_Stress1 + Sim%Diagonal_Stress2 + Sim%Diagonal_Stress3

Sim%Residual_Stress1 (:,:,:,0) = Elem%DumpSx(:,:,:,0) * Sim%Residual_Stress1 (:,:,:,0) + &
   Elem%DumpSx(:,:,:,1) * Dt * (Elem%Acoeff(:,:,:,18) * dVy_dxi + Elem%Acoeff(:,:,:,19) * dVy_deta + Elem%Acoeff(:,:,:,20) * dVy_dzeta)
Sim%Residual_Stress2 (:,:,:,0) = Elem%DumpSy(:,:,:,0) * Sim%Residual_Stress2 (:,:,:,0) + &
   Elem%DumpSy(:,:,:,1) * Dt * (Elem%Acoeff(:,:,:,21) * dVx_dxi + Elem%Acoeff(:,:,:,22) * dVx_deta + Elem%Acoeff(:,:,:,23) * dVx_dzeta)

Sim%Residual_Stress1 (:,:,:,1) = Elem%DumpSx(:,:,:,0) * Sim%Residual_Stress1 (:,:,:,1) + &
   Elem%DumpSx(:,:,:,1) * Dt * (Elem%Acoeff(:,:,:,18) * dVz_dxi + Elem%Acoeff(:,:,:,19) * dVz_deta + Elem%Acoeff(:,:,:,20) * dVz_dzeta)
Sim%Residual_Stress2 (:,:,:,1) = Elem%DumpSz(:,:,:,0) * Sim%Residual_Stress2 (:,:,:,1) + &
   Elem%DumpSz(:,:,:,1) * Dt * (Elem%Acoeff(:,:,:,24) * dVx_dxi + Elem%Acoeff(:,:,:,25) * dVx_deta + Elem%Acoeff(:,:,:,26) * dVx_dzeta)

Sim%Residual_Stress1 (:,:,:,2) = Elem%DumpSy(:,:,:,0) * Sim%Residual_Stress1 (:,:,:,2) + &
   Elem%DumpSy(:,:,:,1) * Dt * (Elem%Acoeff(:,:,:,21) * dVz_dxi + Elem%Acoeff(:,:,:,22) * dVz_deta + Elem%Acoeff(:,:,:,23) * dVz_dzeta)
Sim%Residual_Stress2 (:,:,:,2) = Elem%DumpSz(:,:,:,0) * Sim%Residual_Stress2 (:,:,:,2) + &
   Elem%DumpSz(:,:,:,1) * Dt * (Elem%Acoeff(:,:,:,24) * dVy_dxi + Elem%Acoeff(:,:,:,25) * dVy_deta + Elem%Acoeff(:,:,:,26) * dVy_dzeta)

Sim%Residual_Stress = Sim%Residual_Stress1 + Sim%Residual_Stress2 


return
end subroutine Prediction_Elem_PML_Veloc

! ###########################################################
subroutine Prediction_Elem_PML_Veloc_curve (Elem, bega, dt, hTprimex, hprimey, hprimez, Sim)


implicit none

type (Element), intent (IN) :: Elem
type (Simu), intent (INOUT) :: Sim
doubleprecision, intent (IN) :: bega, dt
doubleprecision, dimension (0:Elem%ngllx-1, 0:Elem%ngllx-1), intent (IN) :: hTprimex
doubleprecision, dimension (0:Elem%nglly-1, 0:Elem%nglly-1), intent (IN) :: hprimey
doubleprecision, dimension (0:Elem%ngllz-1, 0:Elem%ngllz-1), intent (IN) :: hprimez

integer :: m1,m2,m3, n_x,n_y,n_z
doubleprecision, dimension(0:2,0:2) :: inv_norm, tmp, res, norm
doubleprecision, dimension (0:Elem%ngllx-1, 0:Elem%nglly-1, 0:Elem%ngllz-1) :: dVx_dxi,dVx_deta,dVx_dzeta, &
                                                                               dVy_dxi,dVy_deta,dVy_dzeta, &
                                                                               dVz_dxi,dVz_deta,dVz_dzeta
doubleprecision, dimension (0:Elem%ngllx-1, 0:Elem%nglly-1, 0:Elem%ngllz-1, 0:2, 0:2) :: VxMat,VyMat,VzMat, &
                                                                                         VxTerms,VyTerms,VzTerms


! On calcule les derivees de la vitesse par rapport a xi, eta et zeta

m1 = Elem%ngllx; m2 = Elem%nglly; m3= Elem%ngllz

Sim%Forces(1:m1-2, 1:m2-2, 1:m3-2, 0:2) = Sim%Veloc(:,:,:,:) + dt * (0.5-bega) * Sim%Accel(:,:,:,:)

call DGEMM ('N', 'N', m1, m2*m3, m1, 1.d0, htprimex, m1, Sim%Forces(0,0,0,0), m1, 0.d0, dVx_dxi, m1)
do n_z = 0,m3-1
   call DGEMM ('N', 'N', m1, m2, m2, 1.d0, Sim%Forces(0,0,n_z,0), m1, hprimey, m2, 0.d0, dVx_deta(0,0,n_z), m1)
enddo
call DGEMM ('N', 'N', m1*m2, m3, m3, 1.d0, Sim%Forces(0,0,0,0), m1*m2, hprimez, m3, 0.d0, dVx_dzeta, m1*m2)

call DGEMM ('N', 'N', m1, m2*m3, m1, 1.d0, htprimex, m1, Sim%Forces(0,0,0,1), m1, 0.d0, dVy_dxi, m1)
do n_z = 0,m3-1
   call DGEMM ('N', 'N', m1, m2, m2, 1.d0, Sim%Forces(0,0,n_z,1), m1, hprimey, m2, 0.d0, dVy_deta(0,0,n_z), m1)
enddo
call DGEMM ('N', 'N', m1*m2, m3, m3, 1.d0, Sim%Forces(0,0,0,1), m1*m2, hprimez, m3, 0.d0, dVy_dzeta, m1*m2)

call DGEMM ('N', 'N', m1, m2*m3, m1, 1.d0, htprimex, m1, Sim%Forces(0,0,0,2), m1, 0.d0, dVz_dxi, m1)
do n_z = 0,m3-1
   call DGEMM ('N', 'N', m1, m2, m2, 1.d0, Sim%Forces(0,0,n_z,2), m1, hprimey, m2, 0.d0, dVz_deta(0,0,n_z), m1)
enddo
call DGEMM ('N', 'N', m1*m2, m3, m3, 1.d0, Sim%Forces(0,0,0,2), m1*m2, hprimez, m3, 0.d0, dVz_dzeta, m1*m2)


! On cree VxMat, VyMat et VzMat:

VxMat(:,:,:,0,0) = Elem%Acoeff(:,:,:,0) * dVx_dxi + Elem%Acoeff(:,:,:,1) * dVx_deta + Elem%Acoeff(:,:,:,2) * dVx_dzeta
VxMat(:,:,:,0,1) = Elem%Acoeff(:,:,:,12) * dVx_dxi + Elem%Acoeff(:,:,:,13) * dVx_deta + Elem%Acoeff(:,:,:,14) * dVx_dzeta
VxMat(:,:,:,0,2) = Elem%Acoeff(:,:,:,15) * dVx_dxi + Elem%Acoeff(:,:,:,16) * dVx_deta + Elem%Acoeff(:,:,:,17) * dVx_dzeta
VxMat(:,:,:,1,0) = Elem%Acoeff(:,:,:,9) * dVx_dxi + Elem%Acoeff(:,:,:,10) * dVx_deta + Elem%Acoeff(:,:,:,11) * dVx_dzeta
VxMat(:,:,:,1,1) = Elem%Acoeff(:,:,:,3) * dVx_dxi + Elem%Acoeff(:,:,:,4) * dVx_deta + Elem%Acoeff(:,:,:,5) * dVx_dzeta
VxMat(:,:,:,1,2) = Elem%Acoeff(:,:,:,6) * dVx_dxi + Elem%Acoeff(:,:,:,7) * dVx_deta + Elem%Acoeff(:,:,:,8) * dVx_dzeta
VxMat(:,:,:,2,0) = Elem%Acoeff(:,:,:,18) * dVx_dxi + Elem%Acoeff(:,:,:,19) * dVx_deta + Elem%Acoeff(:,:,:,20) * dVx_dzeta
VxMat(:,:,:,2,1) = Elem%Acoeff(:,:,:,21) * dVx_dxi + Elem%Acoeff(:,:,:,22) * dVx_deta + Elem%Acoeff(:,:,:,23) * dVx_dzeta
VxMat(:,:,:,2,2) = Elem%Acoeff(:,:,:,24) * dVx_dxi + Elem%Acoeff(:,:,:,25) * dVx_deta + Elem%Acoeff(:,:,:,26) * dVx_dzeta

VyMat(:,:,:,0,0) = Elem%Acoeff(:,:,:,0) * dVy_dxi + Elem%Acoeff(:,:,:,1) * dVy_deta + Elem%Acoeff(:,:,:,2) * dVy_dzeta
VyMat(:,:,:,0,1) = Elem%Acoeff(:,:,:,12) * dVy_dxi + Elem%Acoeff(:,:,:,13) * dVy_deta + Elem%Acoeff(:,:,:,14) * dVy_dzeta
VyMat(:,:,:,0,2) = Elem%Acoeff(:,:,:,15) * dVy_dxi + Elem%Acoeff(:,:,:,16) * dVy_deta + Elem%Acoeff(:,:,:,17) * dVy_dzeta
VyMat(:,:,:,1,0) = Elem%Acoeff(:,:,:,9) * dVy_dxi + Elem%Acoeff(:,:,:,10) * dVy_deta + Elem%Acoeff(:,:,:,11) * dVy_dzeta
VyMat(:,:,:,1,1) = Elem%Acoeff(:,:,:,3) * dVy_dxi + Elem%Acoeff(:,:,:,4) * dVy_deta + Elem%Acoeff(:,:,:,5) * dVy_dzeta
VyMat(:,:,:,1,2) = Elem%Acoeff(:,:,:,6) * dVy_dxi + Elem%Acoeff(:,:,:,7) * dVy_deta + Elem%Acoeff(:,:,:,8) * dVy_dzeta
VyMat(:,:,:,2,0) = Elem%Acoeff(:,:,:,18) * dVy_dxi + Elem%Acoeff(:,:,:,19) * dVy_deta + Elem%Acoeff(:,:,:,20) * dVy_dzeta
VyMat(:,:,:,2,1) = Elem%Acoeff(:,:,:,21) * dVy_dxi + Elem%Acoeff(:,:,:,22) * dVy_deta + Elem%Acoeff(:,:,:,23) * dVy_dzeta
VyMat(:,:,:,2,2) = Elem%Acoeff(:,:,:,24) * dVy_dxi + Elem%Acoeff(:,:,:,25) * dVy_deta + Elem%Acoeff(:,:,:,26) * dVy_dzeta

VzMat(:,:,:,0,0) = Elem%Acoeff(:,:,:,0) * dVz_dxi + Elem%Acoeff(:,:,:,1) * dVz_deta + Elem%Acoeff(:,:,:,2) * dVz_dzeta
VzMat(:,:,:,0,1) = Elem%Acoeff(:,:,:,12) * dVz_dxi + Elem%Acoeff(:,:,:,13) * dVz_deta + Elem%Acoeff(:,:,:,14) * dVz_dzeta
VzMat(:,:,:,0,2) = Elem%Acoeff(:,:,:,15) * dVz_dxi + Elem%Acoeff(:,:,:,16) * dVz_deta + Elem%Acoeff(:,:,:,17) * dVz_dzeta
VzMat(:,:,:,1,0) = Elem%Acoeff(:,:,:,9) * dVz_dxi + Elem%Acoeff(:,:,:,10) * dVz_deta + Elem%Acoeff(:,:,:,11) * dVz_dzeta
VzMat(:,:,:,1,1) = Elem%Acoeff(:,:,:,3) * dVz_dxi + Elem%Acoeff(:,:,:,4) * dVz_deta + Elem%Acoeff(:,:,:,5) * dVz_dzeta
VzMat(:,:,:,1,2) = Elem%Acoeff(:,:,:,6) * dVz_dxi + Elem%Acoeff(:,:,:,7) * dVz_deta + Elem%Acoeff(:,:,:,8) * dVz_dzeta
VzMat(:,:,:,2,0) = Elem%Acoeff(:,:,:,18) * dVz_dxi + Elem%Acoeff(:,:,:,19) * dVz_deta + Elem%Acoeff(:,:,:,20) * dVz_dzeta
VzMat(:,:,:,2,1) = Elem%Acoeff(:,:,:,21) * dVz_dxi + Elem%Acoeff(:,:,:,22) * dVz_deta + Elem%Acoeff(:,:,:,23) * dVz_dzeta
VzMat(:,:,:,2,2) = Elem%Acoeff(:,:,:,24) * dVz_dxi + Elem%Acoeff(:,:,:,25) * dVz_deta + Elem%Acoeff(:,:,:,26) * dVz_dzeta


! On cree VxTerms, VyTerms et VzTerms:

inv_norm = transpose(Elem%Inv_Normales)
do n_x = 0,m1-1
   do n_y = 0,m2-1
      do n_z = 0,m3-1
         tmp(0:2,0:2) = VxMat(n_x,n_y,n_z,0:2,0:2)
         call DGEMM ('N', 'N', 3, 3, 3, 1.d0, tmp, 3, inv_norm, 3, 0.d0, res, 3)
         VxTerms(n_x,n_y,n_z,0:2,0:2) = res(0:2,0:2)
         
         tmp(0:2,0:2) = VyMat(n_x,n_y,n_z,0:2,0:2)
         call DGEMM ('N', 'N', 3, 3, 3, 1.d0, tmp, 3, inv_norm, 3, 0.d0, res, 3)
         VyTerms(n_x,n_y,n_z,0:2,0:2) = res(0:2,0:2)

         tmp(0:2,0:2) = VzMat(n_x,n_y,n_z,0:2,0:2)
         call DGEMM ('N', 'N', 3, 3, 3, 1.d0, tmp, 3, inv_norm, 3, 0.d0, res, 3)
         VzTerms(n_x,n_y,n_z,0:2,0:2) = res(0:2,0:2)
      enddo
   enddo
enddo


! On calcule Diagonal_Stress:

norm = Elem%Normales

Sim%Diagonal_Stress1(:,:,:,0) = Elem%DumpSx(:,:,:,0) * Sim%Diagonal_Stress1(:,:,:,0) + &
   Elem%DumpSx(:,:,:,1) * Dt * (norm(0,0)*VxTerms(:,:,:,0,0) + norm(1,0)*VyTerms(:,:,:,1,0) + norm(2,0)*VzTerms(:,:,:,1,0))
Sim%Diagonal_Stress1(:,:,:,1) = Elem%DumpSx(:,:,:,0) * Sim%Diagonal_Stress1(:,:,:,1) + &
   Elem%DumpSx(:,:,:,1) * Dt * (norm(0,0)*VxTerms(:,:,:,1,0) + norm(1,0)*VyTerms(:,:,:,0,0) + norm(2,0)*VzTerms(:,:,:,1,0))
Sim%Diagonal_Stress1(:,:,:,2) = Elem%DumpSx(:,:,:,0) * Sim%Diagonal_Stress1(:,:,:,2) + &
   Elem%DumpSx(:,:,:,1) * Dt * (norm(0,0)*VxTerms(:,:,:,1,0) + norm(1,0)*VyTerms(:,:,:,1,0) + norm(2,0)*VzTerms(:,:,:,0,0))

Sim%Diagonal_Stress2(:,:,:,0) = Elem%DumpSy(:,:,:,0) * Sim%Diagonal_Stress2(:,:,:,0) + &
   Elem%DumpSy(:,:,:,1) * Dt * (norm(0,1)*VxTerms(:,:,:,0,1) + norm(1,1)*VyTerms(:,:,:,1,1) + norm(2,1)*VzTerms(:,:,:,1,1))
Sim%Diagonal_Stress2(:,:,:,1) = Elem%DumpSy(:,:,:,0) * Sim%Diagonal_Stress2(:,:,:,1) + &
   Elem%DumpSy(:,:,:,1) * Dt * (norm(0,1)*VxTerms(:,:,:,1,1) + norm(1,1)*VyTerms(:,:,:,0,1) + norm(2,1)*VzTerms(:,:,:,1,1))
Sim%Diagonal_Stress2(:,:,:,2) = Elem%DumpSy(:,:,:,0) * Sim%Diagonal_Stress2(:,:,:,2) + &
   Elem%DumpSy(:,:,:,1) * Dt * (norm(0,1)*VxTerms(:,:,:,1,1) + norm(1,1)*VyTerms(:,:,:,1,1) + norm(2,1)*VzTerms(:,:,:,0,1))

Sim%Diagonal_Stress3(:,:,:,0) = Elem%DumpSz(:,:,:,0) * Sim%Diagonal_Stress3(:,:,:,0) + &
   Elem%DumpSz(:,:,:,1) * Dt * (norm(0,2)*VxTerms(:,:,:,0,2) + norm(1,2)*VyTerms(:,:,:,1,2) + norm(2,2)*VzTerms(:,:,:,1,2))
Sim%Diagonal_Stress3(:,:,:,1) = Elem%DumpSz(:,:,:,0) * Sim%Diagonal_Stress3(:,:,:,1) + &
   Elem%DumpSz(:,:,:,1) * Dt * (norm(0,2)*VxTerms(:,:,:,1,2) + norm(1,2)*VyTerms(:,:,:,0,2) + norm(2,2)*VzTerms(:,:,:,1,2))
Sim%Diagonal_Stress3(:,:,:,2) = Elem%DumpSz(:,:,:,0) * Sim%Diagonal_Stress3(:,:,:,2) + &
   Elem%DumpSz(:,:,:,1) * Dt * (norm(0,2)*VxTerms(:,:,:,1,2) + norm(1,2)*VyTerms(:,:,:,1,2) + norm(2,2)*VzTerms(:,:,:,0,2))

Sim%Diagonal_Stress = Sim%Diagonal_Stress1 + Sim%Diagonal_Stress2 + Sim%Diagonal_Stress3


! On calcule Residual_Stress:

Sim%Residual_Stress1(:,:,:,0) = Elem%DumpSx(:,:,:,0) * Sim%Residual_Stress1(:,:,:,0) + &
   Elem%DumpSx(:,:,:,1) * Dt * (norm(0,0)*VyTerms(:,:,:,2,0) + norm(1,0)*VxTerms(:,:,:,2,0))
Sim%Residual_Stress1(:,:,:,1) = Elem%DumpSx(:,:,:,0) * Sim%Residual_Stress1(:,:,:,1) + &
   Elem%DumpSx(:,:,:,1) * Dt * (norm(0,0)*VzTerms(:,:,:,2,0) + norm(2,0)*VxTerms(:,:,:,2,0))
Sim%Residual_Stress1(:,:,:,2) = Elem%DumpSx(:,:,:,0) * Sim%Residual_Stress1(:,:,:,2) + &
   Elem%DumpSx(:,:,:,1) * Dt * (norm(1,0)*VzTerms(:,:,:,2,0) + norm(2,0)*VyTerms(:,:,:,2,0))

Sim%Residual_Stress2(:,:,:,0) = Elem%DumpSy(:,:,:,0) * Sim%Residual_Stress2(:,:,:,0) + &
   Elem%DumpSy(:,:,:,1) * Dt * (norm(0,1)*VyTerms(:,:,:,2,1) + norm(1,1)*VxTerms(:,:,:,2,1))
Sim%Residual_Stress2(:,:,:,1) = Elem%DumpSy(:,:,:,0) * Sim%Residual_Stress2(:,:,:,1) + &
   Elem%DumpSy(:,:,:,1) * Dt * (norm(0,1)*VzTerms(:,:,:,2,1) + norm(2,1)*VxTerms(:,:,:,2,1))
Sim%Residual_Stress2(:,:,:,2) = Elem%DumpSy(:,:,:,0) * Sim%Residual_Stress2(:,:,:,2) + &
   Elem%DumpSy(:,:,:,1) * Dt * (norm(1,1)*VzTerms(:,:,:,2,1) + norm(2,1)*VyTerms(:,:,:,2,1))

Sim%Residual_Stress3(:,:,:,0) = Elem%DumpSz(:,:,:,0) * Sim%Residual_Stress3(:,:,:,0) + &
   Elem%DumpSz(:,:,:,1) * Dt * (norm(0,2)*VyTerms(:,:,:,2,2) + norm(1,2)*VxTerms(:,:,:,2,2))
Sim%Residual_Stress3(:,:,:,1) = Elem%DumpSz(:,:,:,0) * Sim%Residual_Stress3(:,:,:,1) + &
   Elem%DumpSz(:,:,:,1) * Dt * (norm(0,2)*VzTerms(:,:,:,2,2) + norm(2,2)*VxTerms(:,:,:,2,2))
Sim%Residual_Stress3(:,:,:,2) = Elem%DumpSz(:,:,:,0) * Sim%Residual_Stress3(:,:,:,2) + &
   Elem%DumpSz(:,:,:,1) * Dt * (norm(1,2)*VzTerms(:,:,:,2,2) + norm(2,2)*VyTerms(:,:,:,2,2))

Sim%Residual_Stress = Sim%Residual_Stress1 + Sim%Residual_Stress2 + Sim%Residual_Stress3


return
end subroutine Prediction_Elem_PML_Veloc_curve

! ###########################################################
subroutine Correction_Elem_PML_Veloc (Elem, i_simu, dt)


implicit none

type (Element), intent (INOUT) :: Elem
doubleprecision, intent (IN) :: dt
integer, intent (IN) :: i_simu

integer :: ngllx, nglly, ngllz, i 


ngllx = Elem%ngllx; nglly = Elem%nglly; ngllz = Elem%ngllz

do i = 0,2
   Elem%sSimu(i_simu)%Veloc1(:,:,:,i) = Elem%DumpVx(:,:,:,0) * Elem%sSimu(i_simu)%Veloc1(:,:,:,i) + &
              dt * Elem%DumpVx(:,:,:,1) * Elem%sSimu(i_simu)%Forces1(1:ngllx-2,1:nglly-2,1:ngllz-2,i)
   Elem%sSimu(i_simu)%Veloc2(:,:,:,i) = Elem%DumpVy(:,:,:,0) * Elem%sSimu(i_simu)%Veloc2(:,:,:,i) + &
              dt * Elem%DumpVy(:,:,:,1) * Elem%sSimu(i_simu)%Forces2(1:ngllx-2,1:nglly-2,1:ngllz-2,i)
   Elem%sSimu(i_simu)%Veloc3(:,:,:,i) = Elem%DumpVz(:,:,:,0) * Elem%sSimu(i_simu)%Veloc3(:,:,:,i) + &
              dt * Elem%DumpVz(:,:,:,1) * Elem%sSimu(i_simu)%Forces3(1:ngllx-2,1:nglly-2,1:ngllz-2,i)
enddo

Elem%sSimu(i_simu)%Veloc = Elem%sSimu(i_simu)%Veloc1 + Elem%sSimu(i_simu)%Veloc2 + Elem%sSimu(i_simu)%Veloc3


return
end subroutine Correction_Elem_PML_Veloc

! ###########################################################
subroutine  compute_InternalForces_PML_Elem (Elem, i_simu, hprimex, hTprimey, hTprimez)


implicit none

type (Element), intent (INOUT) :: Elem
doubleprecision, dimension (0:Elem%ngllx-1, 0:Elem%ngllx-1), intent (IN) :: hprimex
doubleprecision, dimension (0:Elem%nglly-1, 0:Elem%nglly-1), intent (IN) :: hTprimey
doubleprecision, dimension (0:Elem%ngllz-1, 0:Elem%ngllz-1), intent (IN) :: hTprimez
integer, intent(IN) :: i_simu

integer :: m1, m2, m3, n_z
doubleprecision, dimension (0:Elem%ngllx-1, 0:Elem%nglly-1, 0:Elem%ngllz-1) :: s0,s1


m1 = Elem%ngllx;  m2 = Elem%nglly;  m3 = Elem%ngllz

s0 = Elem%Acoeff(:,:,:,27) * Elem%sSimu(i_simu)%Diagonal_Stress(:,:,:,0) + &
     Elem%Acoeff(:,:,:,28) * Elem%sSimu(i_simu)%Residual_Stress(:,:,:,0) + & 
     Elem%Acoeff(:,:,:,29) * Elem%sSimu(i_simu)%Residual_Stress(:,:,:,1)
call DGEMM ( 'N', 'N', m1, m2*m3, m1, 1.d0, hprimex, m1, s0(0,0,0), m1, 0.d0, s1, m1 )
Elem%sSimu(i_simu)%Forces1(:,:,:,0) = s1

s0 = Elem%Acoeff(:,:,:,30) * Elem%sSimu(i_simu)%Diagonal_Stress(:,:,:,0) + &
     Elem%Acoeff(:,:,:,31) * Elem%sSimu(i_simu)%Residual_Stress(:,:,:,0) + & 
     Elem%Acoeff(:,:,:,32) * Elem%sSimu(i_simu)%Residual_Stress(:,:,:,1)
do n_z = 0,m3-1
     call DGEMM ( 'N', 'N', m1, m2, m2, 1.d0, s0(0,0,n_z), m1, htprimey ,m2, 0.d0, s1(0,0,n_z), m1 )
enddo
Elem%sSimu(i_simu)%Forces2(:,:,:,0) = s1

s0 = Elem%Acoeff(:,:,:,33) * Elem%sSimu(i_simu)%Diagonal_Stress(:,:,:,0) + &
     Elem%Acoeff(:,:,:,34) * Elem%sSimu(i_simu)%Residual_Stress(:,:,:,0) + & 
     Elem%Acoeff(:,:,:,35) * Elem%sSimu(i_simu)%Residual_Stress(:,:,:,1)
call DGEMM ( 'N', 'N', m1*m2, m3, m3, 1.d0, s0(0,0,0), m1*m2, htprimez ,m3, 0.d0, s1, m1*m2 )
Elem%sSimu(i_simu)%Forces3(:,:,:,0) = s1

s0 = Elem%Acoeff(:,:,:,27) * Elem%sSimu(i_simu)%Residual_Stress(:,:,:,0) + &
     Elem%Acoeff(:,:,:,28) * Elem%sSimu(i_simu)%Diagonal_Stress(:,:,:,1) + & 
     Elem%Acoeff(:,:,:,29) * Elem%sSimu(i_simu)%Residual_Stress(:,:,:,2)
call DGEMM ( 'N', 'N', m1, m2*m3, m1, 1.d0, hprimex, m1, s0(0,0,0) ,m1, 0.d0, s1, m1 )
Elem%sSimu(i_simu)%Forces1(:,:,:,1) = s1

s0 = Elem%Acoeff(:,:,:,30) * Elem%sSimu(i_simu)%Residual_Stress(:,:,:,0) + &
     Elem%Acoeff(:,:,:,31) * Elem%sSimu(i_simu)%Diagonal_Stress(:,:,:,1) + & 
     Elem%Acoeff(:,:,:,32) * Elem%sSimu(i_simu)%Residual_Stress(:,:,:,2)
do n_z = 0,m3-1
     call DGEMM ( 'N', 'N', m1, m2, m2, 1.d0, s0(0,0,n_z), m1, htprimey ,m2, 0.d0, s1(0,0,n_z), m1 )
enddo
Elem%sSimu(i_simu)%Forces2(:,:,:,1) = s1

s0 = Elem%Acoeff(:,:,:,33) * Elem%sSimu(i_simu)%Residual_Stress(:,:,:,0) + &
     Elem%Acoeff(:,:,:,34) * Elem%sSimu(i_simu)%Diagonal_Stress(:,:,:,1) + &
     Elem%Acoeff(:,:,:,35) * Elem%sSimu(i_simu)%Residual_Stress(:,:,:,2)
call DGEMM ( 'N', 'N', m1*m2, m3, m3, 1.d0, s0(0,0,0), m1*m2, htprimez, m3, 0.d0, s1, m1*m2 )
Elem%sSimu(i_simu)%Forces3(:,:,:,1) = s1

s0 = Elem%Acoeff(:,:,:,27) * Elem%sSimu(i_simu)%Residual_Stress(:,:,:,1) + &
     Elem%Acoeff(:,:,:,28) * Elem%sSimu(i_simu)%Residual_Stress(:,:,:,2) + & 
     Elem%Acoeff(:,:,:,29) * Elem%sSimu(i_simu)%Diagonal_Stress(:,:,:,2)
call DGEMM ( 'N', 'N', m1, m2*m3, m1, 1.d0, hprimex, m1, s0(0,0,0) ,m1, 0.d0, s1, m1 )
Elem%sSimu(i_simu)%Forces1(:,:,:,2) = s1

s0 = Elem%Acoeff(:,:,:,30) * Elem%sSimu(i_simu)%Residual_Stress(:,:,:,1) + &
     Elem%Acoeff(:,:,:,31) * Elem%sSimu(i_simu)%Residual_Stress(:,:,:,2) + & 
     Elem%Acoeff(:,:,:,32) * Elem%sSimu(i_simu)%Diagonal_Stress(:,:,:,2)
do n_z = 0,m3-1
    call DGEMM ( 'N', 'N', m1, m2, m2, 1.d0,s0(0,0,n_z), m1, htprimey ,m2, 0.d0, s1(0,0,n_z), m1 )
enddo
Elem%sSimu(i_simu)%Forces2(:,:,:,2) = s1

s0 = Elem%Acoeff(:,:,:,33) * Elem%sSimu(i_simu)%Residual_Stress(:,:,:,1) + &
     Elem%Acoeff(:,:,:,34) * Elem%sSimu(i_simu)%Residual_Stress(:,:,:,2) + & 
     Elem%Acoeff(:,:,:,35) * Elem%sSimu(i_simu)%Diagonal_Stress(:,:,:,2)
call DGEMM ( 'N', 'N', m1*m2, m3, m3, 1.d0, s0(0,0,0), m1*m2, htprimez ,m3, 0.d0, s1, m1*m2 )
Elem%sSimu(i_simu)%Forces3(:,:,:,2) = s1

Elem%sSimu(i_simu)%Forces = Elem%sSimu(i_simu)%Forces1 + &
                            Elem%sSimu(i_simu)%Forces2 + &
                            Elem%sSimu(i_simu)%Forces3


return
end subroutine compute_InternalForces_PML_Elem

! ###########################################################
end module selements
