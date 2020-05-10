subroutine internal_forces (Elem, i_simu, hprimex, htprimex, hprimey, htprimey, hprimez, htprimez, n_solid, aniso, adj, trm)


use sdomains

implicit none

type (Element), intent (INOUT) :: Elem
doubleprecision, dimension (0:Elem%ngllx-1, 0:Elem%ngllx-1), intent (IN) :: hprimex, htprimex
doubleprecision, dimension (0:Elem%nglly-1, 0:Elem%nglly-1), intent (IN) :: hprimey, htprimey
doubleprecision, dimension (0:Elem%ngllz-1, 0:Elem%ngllz-1), intent (IN) :: hprimez, hTprimez
integer, intent(IN) :: n_solid, i_simu, trm
logical, intent(IN) :: aniso, adj

integer :: n_z, m1,m2,m3, i,j,k, mat
doubleprecision :: epsilon_trace_over_3
doubleprecision, dimension (0:Elem%ngllx-1, 0:Elem%nglly-1, 0:Elem%ngllz-1) :: dUx_dxi, dUx_deta, dUx_dzeta, &
                                                                               dUy_dxi, dUy_deta, dUy_dzeta, &
                                                                               dUz_dxi, dUz_deta, dUz_dzeta, &
                                                                               DXX,DXY,DXZ,DYX,DYY,DYZ,DZX,DZY,DZZ, &
                                                                               Fox,Foy,Foz
doubleprecision, dimension(:,:,:), allocatable :: epsilondev_xx_loc, epsilondev_yy_loc, &
                                                  epsilondev_xy_loc, epsilondev_xz_loc, epsilondev_yz_loc


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

do i = 0,m1-1
 do j = 0,m2-1
  do k = 0,m3-1

     dxx(i,j,k) = dUx_dxi(i,j,k)*Elem%InvGrad(i,j,k,0,0) + dUx_deta(i,j,k)*Elem%InvGrad(i,j,k,0,1) + dUx_dzeta(i,j,k)*Elem%InvGrad(i,j,k,0,2)
     dyy(i,j,k) = dUy_dxi(i,j,k)*Elem%InvGrad(i,j,k,1,0) + dUy_deta(i,j,k)*Elem%InvGrad(i,j,k,1,1) + dUy_dzeta(i,j,k)*Elem%InvGrad(i,j,k,1,2)
     dzz(i,j,k) = dUz_dxi(i,j,k)*Elem%InvGrad(i,j,k,2,0) + dUz_deta(i,j,k)*Elem%InvGrad(i,j,k,2,1) + dUz_dzeta(i,j,k)*Elem%InvGrad(i,j,k,2,2)

     dyx(i,j,k) = dUx_dxi(i,j,k)*Elem%InvGrad(i,j,k,1,0) + dUx_deta(i,j,k)*Elem%InvGrad(i,j,k,1,1) + dUx_dzeta(i,j,k)*Elem%InvGrad(i,j,k,1,2)
     dzx(i,j,k) = dUx_dxi(i,j,k)*Elem%InvGrad(i,j,k,2,0) + dUx_deta(i,j,k)*Elem%InvGrad(i,j,k,2,1) + dUx_dzeta(i,j,k)*Elem%InvGrad(i,j,k,2,2)

     dxy(i,j,k) = dUy_dxi(i,j,k)*Elem%InvGrad(i,j,k,0,0) + dUy_deta(i,j,k)*Elem%InvGrad(i,j,k,0,1) + dUy_dzeta(i,j,k)*Elem%InvGrad(i,j,k,0,2)
     dzy(i,j,k) = dUy_dxi(i,j,k)*Elem%InvGrad(i,j,k,2,0) + dUy_deta(i,j,k)*Elem%InvGrad(i,j,k,2,1) + dUy_dzeta(i,j,k)*Elem%InvGrad(i,j,k,2,2)

     dxz(i,j,k) = dUz_dxi(i,j,k)*Elem%InvGrad(i,j,k,0,0) + dUz_deta(i,j,k)*Elem%InvGrad(i,j,k,0,1) + dUz_dzeta(i,j,k)*Elem%InvGrad(i,j,k,0,2)
     dyz(i,j,k) = dUz_dxi(i,j,k)*Elem%InvGrad(i,j,k,1,0) + dUz_deta(i,j,k)*Elem%InvGrad(i,j,k,1,1) + dUz_dzeta(i,j,k)*Elem%InvGrad(i,j,k,1,2)

  enddo
 enddo
enddo

if (adj) then
   do i = 0,m1-1
    do j = 0,m2-1
     do k = 0,m3-1
        Elem%sSimu(i_simu)%save_strain(i,j,k,0) = DXX(i,j,k)
        Elem%sSimu(i_simu)%save_strain(i,j,k,1) = DYY(i,j,k)
        Elem%sSimu(i_simu)%save_strain(i,j,k,2) = DZZ(i,j,k)
        Elem%sSimu(i_simu)%save_strain(i,j,k,3) = 0.5 * (DXY(i,j,k) + DYX(i,j,k))
        Elem%sSimu(i_simu)%save_strain(i,j,k,4) = 0.5 * (DZX(i,j,k) + DXZ(i,j,k))
        Elem%sSimu(i_simu)%save_strain(i,j,k,5) = 0.5 * (DZY(i,j,k) + DYZ(i,j,k))
     enddo
    enddo
   enddo
endif

if (n_solid>0) then
   allocate (epsilondev_xx_loc(0:m1-1,0:m2-1,0:m3-1))
   allocate (epsilondev_yy_loc(0:m1-1,0:m2-1,0:m3-1))
   allocate (epsilondev_xy_loc(0:m1-1,0:m2-1,0:m3-1))
   allocate (epsilondev_xz_loc(0:m1-1,0:m2-1,0:m3-1))
   allocate (epsilondev_yz_loc(0:m1-1,0:m2-1,0:m3-1))
   do i = 0,m1-1
    do j = 0,m2-1
     do k = 0,m3-1
        epsilon_trace_over_3 = 0.3333333333333333d0 * (DXX(i,j,k) + DYY(i,j,k) + DZZ(i,j,k))
        epsilondev_xx_loc(i,j,k) = DXX(i,j,k) - epsilon_trace_over_3
        epsilondev_yy_loc(i,j,k) = DYY(i,j,k) - epsilon_trace_over_3
        epsilondev_xy_loc(i,j,k) = 0.5 * (DXY(i,j,k) + DYX(i,j,k))
        epsilondev_xz_loc(i,j,k) = 0.5 * (DZX(i,j,k) + DXZ(i,j,k))
        epsilondev_yz_loc(i,j,k) = 0.5 * (DZY(i,j,k) + DYZ(i,j,k))
     enddo
    enddo
   enddo
endif

if (aniso) then
   if (n_solid>0) then
      call calcul_forces_aniso_att(Fox,Foy,Foz, &
                                   Elem%Invgrad(:,:,:,0,0), &
                                   Elem%Invgrad(:,:,:,1,0), &
                                   Elem%Invgrad(:,:,:,2,0), &
                                   Elem%Invgrad(:,:,:,0,1), &
                                   Elem%Invgrad(:,:,:,1,1), &
                                   Elem%Invgrad(:,:,:,2,1), &
                                   Elem%Invgrad(:,:,:,0,2), &
                                   Elem%Invgrad(:,:,:,1,2), &
                                   Elem%Invgrad(:,:,:,2,2), &
                                   htprimex, htprimey, htprimez, &
                                   Elem%Jacob, Elem%wgtx, Elem%wgty, Elem%wgtz, &
                                   DXX,DXY,DXZ, &
                                   DYX,DYY,DYZ, &
                                   DZX,DZY,DZZ, &
                                   Elem%Mu, Elem%Lambda, Elem%Cij, &
                                   m1,m2,m3, n_solid, &
                                   Elem%onemSbeta, Elem%sSimu(i_simu)%R_xx_, Elem%sSimu(i_simu)%R_yy_, &
                                   Elem%sSimu(i_simu)%R_xy_, Elem%sSimu(i_simu)%R_xz_, Elem%sSimu(i_simu)%R_yz_)
      call attenuation_update(Elem%sSimu(i_simu)%epsilondev_xx_,Elem%sSimu(i_simu)%epsilondev_yy_, &
                              Elem%sSimu(i_simu)%epsilondev_xy_,Elem%sSimu(i_simu)%epsilondev_xz_,Elem%sSimu(i_simu)%epsilondev_yz_, &
                              epsilondev_xx_loc,epsilondev_yy_loc, &
                              epsilondev_xy_loc,epsilondev_xz_loc,epsilondev_yz_loc, &
                              Elem%sSimu(i_simu)%R_xx_,Elem%sSimu(i_simu)%R_yy_, &
                              Elem%sSimu(i_simu)%R_xy_,Elem%sSimu(i_simu)%R_xz_,Elem%sSimu(i_simu)%R_yz_, &
                              Elem%factor_common_3,Elem%alphaval_3,Elem%betaval_3,Elem%gammaval_3, &
                              Elem%Mu, m1,m2,m3, n_solid, i_simu, trm)
      deallocate(epsilondev_xx_loc,epsilondev_yy_loc,epsilondev_xy_loc,epsilondev_xz_loc,epsilondev_yz_loc)
   else
      call calcul_forces_aniso(Fox,Foy,Foz,  &
                               Elem%Invgrad(:,:,:,0,0), &
                               Elem%Invgrad(:,:,:,1,0), &
                               Elem%Invgrad(:,:,:,2,0), &
                               Elem%Invgrad(:,:,:,0,1), &
                               Elem%Invgrad(:,:,:,1,1), &
                               Elem%Invgrad(:,:,:,2,1), &
                               Elem%Invgrad(:,:,:,0,2), &
                               Elem%Invgrad(:,:,:,1,2), &
                               Elem%Invgrad(:,:,:,2,2), &
                               htprimex, htprimey, htprimez, &
                               Elem%Jacob, Elem%wgtx, Elem%wgty, Elem%wgtz, &
                               DXX,DXY,DXZ, &
                               DYX,DYY,DYZ, &
                               DZX,DZY,DZZ, &
                               Elem%Cij, &
                               m1,m2,m3)
   endif
else
   if (n_solid>0) then
      call calcul_forces_att(Fox,Foy,Foz, &
                             Elem%Invgrad(:,:,:,0,0), &
                             Elem%Invgrad(:,:,:,1,0), &
                             Elem%Invgrad(:,:,:,2,0), &
                             Elem%Invgrad(:,:,:,0,1), &
                             Elem%Invgrad(:,:,:,1,1), &
                             Elem%Invgrad(:,:,:,2,1), &
                             Elem%Invgrad(:,:,:,0,2), &
                             Elem%Invgrad(:,:,:,1,2), &
                             Elem%Invgrad(:,:,:,2,2), &
                             htprimex, htprimey, htprimez, &
                             Elem%Jacob, Elem%wgtx, Elem%wgty, Elem%wgtz, &
                             DXX,DXY,DXZ, &
                             DYX,DYY,DYZ, &
                             DZX,DZY,DZZ, &
                             Elem%Mu, Elem%Lambda, &
                             m1,m2,m3, n_solid, &
                             Elem%onemSbeta, Elem%sSimu(i_simu)%R_xx_, Elem%sSimu(i_simu)%R_yy_, &
                             Elem%sSimu(i_simu)%R_xy_, Elem%sSimu(i_simu)%R_xz_, Elem%sSimu(i_simu)%R_yz_)
      call attenuation_update(Elem%sSimu(i_simu)%epsilondev_xx_,Elem%sSimu(i_simu)%epsilondev_yy_, &
                              Elem%sSimu(i_simu)%epsilondev_xy_,Elem%sSimu(i_simu)%epsilondev_xz_,Elem%sSimu(i_simu)%epsilondev_yz_, &
                              epsilondev_xx_loc,epsilondev_yy_loc, &
                              epsilondev_xy_loc,epsilondev_xz_loc,epsilondev_yz_loc, &
                              Elem%sSimu(i_simu)%R_xx_,Elem%sSimu(i_simu)%R_yy_, &
                              Elem%sSimu(i_simu)%R_xy_,Elem%sSimu(i_simu)%R_xz_,Elem%sSimu(i_simu)%R_yz_, &
                              Elem%factor_common_3,Elem%alphaval_3,Elem%betaval_3,Elem%gammaval_3, &
                              Elem%Mu, m1,m2,m3, n_solid, i_simu, trm)
      deallocate(epsilondev_xx_loc,epsilondev_yy_loc,epsilondev_xy_loc,epsilondev_xz_loc,epsilondev_yz_loc)
   else
      call calcul_forces(Fox,Foy,Foz,  &
                         Elem%Invgrad(:,:,:,0,0), &
                         Elem%Invgrad(:,:,:,1,0), &
                         Elem%Invgrad(:,:,:,2,0), &
                         Elem%Invgrad(:,:,:,0,1), &
                         Elem%Invgrad(:,:,:,1,1), &
                         Elem%Invgrad(:,:,:,2,1), &
                         Elem%Invgrad(:,:,:,0,2), &
                         Elem%Invgrad(:,:,:,1,2), &
                         Elem%Invgrad(:,:,:,2,2), &
                         htprimex, htprimey, htprimez, &
                         Elem%Jacob, Elem%wgtx, Elem%wgty, Elem%wgtz, &
                         DXX,DXY,DXZ, &
                         DYX,DYY,DYZ, &
                         DZX,DZY,DZZ, &
                         Elem%Mu, Elem%Lambda, &
                         m1,m2,m3)
   endif
endif

Elem%sSimu(i_simu)%Forces(:,:,:,0) = -Fox
Elem%sSimu(i_simu)%Forces(:,:,:,1) = -Foy
Elem%sSimu(i_simu)%Forces(:,:,:,2) = -Foz


return
end subroutine internal_forces
!------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------
subroutine calcul_forces_aniso_att(Fox,Foy,Foz, xi1,xi2,xi3,et1,et2,et3,ga1,ga2,ga3, &
                             dx,dy,dz, jac, poidsx,poidsy,poidsz, DXX,DXY,DXZ,DYX,DYY,DYZ,DZX,DZY,DZZ, &
                             mu_,la_, Cij, ngllx,nglly,ngllz, n_solid, onemSbeta_, R_xx_,R_yy_,R_xy_,R_xz_,R_yz_)


use sdomains

implicit none

integer, intent(in) :: ngllx,nglly,ngllz, n_solid
doubleprecision, dimension(0:20, 0:ngllx-1,0:nglly-1,0:ngllz-1), intent(in) :: Cij
doubleprecision, dimension(0:n_solid-1,0:ngllx-1,0:nglly-1,0:ngllz-1), intent(in) :: R_xx_,R_yy_,R_xy_,R_xz_,R_yz_
doubleprecision, dimension(0:ngllx-1,0:nglly-1,0:ngllz-1), intent(in) :: xi1,xi2,xi3, et1,et2,et3, &
                                                                         ga1,ga2,ga3, jac, mu_,la_, &
                                                                         DXX,DXY,DXZ,DYX,DYY,DYZ,DZX,DZY,DZZ, &
                                                                         onemSbeta_
doubleprecision, dimension(0:ngllx-1,0:nglly-1,0:ngllz-1), intent(out) :: Fox,Foz,Foy
doubleprecision, dimension(0:ngllx-1,0:ngllx-1), intent(in) :: dx
doubleprecision, dimension(0:nglly-1,0:nglly-1), intent(in) :: dy
doubleprecision, dimension(0:ngllz-1,0:ngllz-1), intent(in) :: dz
doubleprecision, dimension(0:ngllx-1), intent(in) :: poidsx
doubleprecision, dimension(0:nglly-1), intent(in) :: poidsy
doubleprecision, dimension(0:ngllz-1), intent(in) :: poidsz

doubleprecision :: sxx,sxy,sxz,syy,syz,szz,xdiv,t4,F1,F2,F3
doubleprecision :: t41,t42,t43,t11,t51,t52,t53,t12,t61,t62,t63,t13
doubleprecision :: xt1,xt2,xt3,xt5,xt6,xt7,xt8,xt9,xt10
doubleprecision, parameter :: deuxtiers = 0.666666666666667, &
                              quatretiers = 1.333333333333337, &
                              s2 = 1.414213562373095, &
                              s2o2 = 0.707106781186547, &
                              zero = 0.d0, two = 2.d0
integer :: i,j,k,l, i_sls, jj, rank, ier
doubleprecision, dimension(0:5) :: eij
doubleprecision, dimension(0:ngllx-1,0:nglly-1,0:ngllz-1) :: xmu,xla,xla2mu,kappa
doubleprecision, dimension(0:ngllx-1,0:nglly-1,0:ngllz-1) :: t1,t5,t8,t2,t6,t9
doubleprecision, dimension(0:ngllx-1,0:nglly-1,0:ngllz-1) :: t3,t7,t10
doubleprecision, dimension(0:5, 0:5, 0:ngllx-1,0:nglly-1,0:ngllz-1) :: C


xmu(:,:,:) = mu_(:,:,:)
xmu(:,:,:) = xmu(:,:,:) * onemSbeta_(:,:,:)
!ATTENTION a l'histoire kappa/lambda !!!
kappa(:,:,:) = la_(:,:,:)
xla(:,:,:) = kappa(:,:,:) - deuxtiers * xmu(:,:,:)
xla2mu(:,:,:) = xla(:,:,:) + two * xmu(:,:,:)

!on remplit la partie superieure isotrope
C(:,:,:,:,:) = 0.d0
C(0,0,:,:,:) = xla2mu(:,:,:)
C(1,1,:,:,:) = xla2mu(:,:,:)
C(2,2,:,:,:) = xla2mu(:,:,:)
C(3,3,:,:,:) = two*xmu(:,:,:)
C(4,4,:,:,:) = two*xmu(:,:,:)
C(5,5,:,:,:) = two*xmu(:,:,:)
C(0,1,:,:,:) = xla(:,:,:)
C(0,2,:,:,:) = xla(:,:,:)
C(1,2,:,:,:) = xla(:,:,:)
!on ajoute la partie anisotrope et on symetrise:
k = 0
do i = 0,5
   do j = i,5
      C(i,j,:,:,:) = C(i,j,:,:,:) + Cij(k,:,:,:)
      k = k + 1
   enddo
enddo
do i = 1,5
   do j = 0,i-1
      C(i,j,:,:,:) = C(j,i,:,:,:)
   enddo
enddo

do k = 0,ngllz-1
 do j = 0,nglly-1
  do i = 0,ngllx-1

      eij(0) = DXX(i,j,k)
      eij(1) = DYY(i,j,k)
      eij(2) = DZZ(i,j,k)
      eij(3) = s2o2*(DYZ(i,j,k)+DZY(i,j,k))
      eij(4) = s2o2*(DXZ(i,j,k)+DZX(i,j,k))
      eij(5) = s2o2*(DXY(i,j,k)+DYX(i,j,k))

      sxx = 0.d0
      syy = 0.d0
      szz = 0.d0
      syz = 0.d0
      sxz = 0.d0
      sxy = 0.d0
      do l = 0,5
          sxx = sxx + C(0,l,i,j,k)*eij(l)
          syy = syy + C(1,l,i,j,k)*eij(l)
          szz = szz + C(2,l,i,j,k)*eij(l)
          syz = syz + C(3,l,i,j,k)*eij(l)
          sxz = sxz + C(4,l,i,j,k)*eij(l)
          sxy = sxy + C(5,l,i,j,k)*eij(l)
      enddo
      syz=syz/s2
      sxz=sxz/s2
      sxy=sxy/s2

      do i_sls = 0,n_solid-1
         sxx = sxx - R_xx_(i_sls,i,j,k)
         syy = syy - R_yy_(i_sls,i,j,k)
         ! ici on utilise le fait que la trace est nulle
         szz = szz + R_xx_(i_sls,i,j,k) + R_yy_(i_sls,i,j,k)
         sxy = sxy - R_xy_(i_sls,i,j,k)
         sxz = sxz - R_xz_(i_sls,i,j,k)
         syz = syz - R_yz_(i_sls,i,j,k)
      enddo
!
!=====================
!       FX 
!=====================
!
      xt1 = sxx * xi1(i,j,k)
      xt2 = sxx * et1(i,j,k)
      xt3 = sxx * ga1(i,j,k)

      xt1 = xt1 + sxy * xi2(i,j,k)
      xt2 = xt2 + sxy * et2(i,j,k)
      xt3 = xt3 + sxy * ga2(i,j,k)

      xt1 = xt1 + sxz * xi3(i,j,k)
      xt2 = xt2 + sxz * et3(i,j,k)
      xt3 = xt3 + sxz * ga3(i,j,k)
!
!=====================
!       FY 
!=====================
!
      xt5 = syy * xi2(i,j,k)
      xt6 = syy * et2(i,j,k)
      xt7 = syy * ga2(i,j,k)

      xt5 = xt5 + sxy * xi1(i,j,k)
      xt6 = xt6 + sxy * et1(i,j,k)
      xt7 = xt7 + sxy * ga1(i,j,k)

      xt5 = xt5 + syz * xi3(i,j,k)
      xt6 = xt6 + syz * et3(i,j,k)
      xt7 = xt7 + syz * ga3(i,j,k)
!
!=====================
!       FZ 
!=====================
!
      xt8  = szz * xi3(i,j,k)
      xt9  = szz * et3(i,j,k)
      xt10 = szz * ga3(i,j,k)

      xt8  = xt8  + sxz * xi1(i,j,k)
      xt9  = xt9  + sxz * et1(i,j,k)
      xt10 = xt10 + sxz * ga1(i,j,k)

      xt8  = xt8  + syz * xi2(i,j,k)
      xt9  = xt9  + syz * et2(i,j,k)
      xt10 = xt10 + syz * ga2(i,j,k)

!
!- Multiplication par le Jacobien et le poids d'integration
!
      t4 = jac(i,j,k) * poidsx(i)
      xt1  =  xt1 * t4
      xt5  =  xt5 * t4
      xt8  =  xt8 * t4

      t4 = jac(i,j,k) * poidsy(j)
      xt2  =  xt2 * t4
      xt6  =  xt6 * t4
      xt9  =  xt9 * t4

      t4 = jac(i,j,k) * poidsz(k)
      xt3  =  xt3 * t4
      xt7  =  xt7 * t4
      xt10 = xt10 * t4

      t1(i,j,k) = xt1
      t5(i,j,k) = xt5
      t8(i,j,k) = xt8

      t2(j,i,k) = xt2
      t6(j,i,k) = xt6
      t9(j,i,k) = xt9

      t3(k,i,j) = xt3
      t7(k,i,j) = xt7
      t10(k,i,j) = xt10

  enddo
 enddo
enddo

!
!- Multiplication par la matrice de derivation puis par les poids
!

!=-=-=-=-=-=-=-=-=-=-
do k = 0,ngllz-1
 do j = 0,nglly-1
  do i = 0,ngllx-1
!=-=-=-=-=-=-=-=-=-=-
!
	t11 = poidsy(j) * poidsz(k)
	t12 = poidsx(i) * poidsz(k)
	t13 = poidsx(i) * poidsy(j)
!
	t41 = zero
	t42 = zero
	t43 = zero
	t51 = zero
	t52 = zero
	t53 = zero
	t61 = zero
	t62 = zero
	t63 = zero
!
	do l = 0,ngllx-1
	   t41 = t41 + dx(l,i) * t1(l,j,k)
	   t42 = t42 + dx(l,i) * t5(l,j,k)
	   t43 = t43 + dx(l,i) * t8(l,j,k)
	enddo

	do l = 0,nglly-1
	   t51 = t51 + dy(l,j) * t2(l,i,k)
	   t52 = t52 + dy(l,j) * t6(l,i,k)
	   t53 = t53 + dy(l,j) * t9(l,i,k)
	enddo
! FX
	F1 = t41*t11 + t51*t12
! FY
	F2 = t42*t11 + t52*t12
! FZ
	F3 = t43*t11 + t53*t12
!
	do l = 0,ngllz-1
	   t61 = t61 + dz(l,k) * t3(l,i,j)
	   t62 = t62 + dz(l,k) * t7(l,i,j)
	   t63 = t63 + dz(l,k) * t10(l,i,j)
	enddo
! FX
	F1 = F1 + t61*t13
! FY
	F2 = F2 + t62*t13
! FZ
	F3 = F3 + t63*t13
!
	Fox(i,j,k) = F1
	Foy(i,j,k) = F2
	Foz(i,j,k) = F3

!=-=-=-=-=-=-=-=-=-=-
  enddo
 enddo
enddo
!=-=-=-=-=-=-=-=-=-=-

end subroutine calcul_forces_aniso_att
!------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------
subroutine calcul_forces_aniso(Fox,Foy,Foz, xi1,xi2,xi3,et1,et2,et3,ga1,ga2,ga3, &
                               dx,dy,dz, jac, poidsx,poidsy,poidsz, DXX,DXY,DXZ,DYX,DYY,DYZ,DZX,DZY,DZZ, &
                               Cij, ngllx,nglly,ngllz)


use sdomains

implicit none

integer, intent(in) :: ngllx,nglly,ngllz
doubleprecision, dimension(0:20, 0:ngllx-1,0:nglly-1,0:ngllz-1), intent(in) :: Cij
doubleprecision, dimension(0:ngllx-1,0:nglly-1,0:ngllz-1), intent(in) :: xi1,xi2,xi3, et1,et2,et3, &
                                                                         ga1,ga2,ga3, jac, &
                                                                         DXX,DXY,DXZ,DYX,DYY,DYZ,DZX,DZY,DZZ
doubleprecision, dimension(0:ngllx-1,0:nglly-1,0:ngllz-1), intent(out) :: Fox,Foz,Foy
doubleprecision, dimension(0:ngllx-1,0:ngllx-1), intent(in) :: dx
doubleprecision, dimension(0:nglly-1,0:nglly-1), intent(in) :: dy
doubleprecision, dimension(0:ngllz-1,0:ngllz-1), intent(in) :: dz
doubleprecision, dimension(0:ngllx-1), intent(in) :: poidsx
doubleprecision, dimension(0:nglly-1), intent(in) :: poidsy
doubleprecision, dimension(0:ngllz-1), intent(in) :: poidsz

integer :: i,j,k,l
doubleprecision :: sxx,sxy,sxz,syy,syz,szz,t4,F1,F2,F3
doubleprecision :: t41,t42,t43,t11,t51,t52,t53,t12,t61,t62,t63,t13
doubleprecision :: xt1,xt2,xt3,xt5,xt6,xt7,xt8,xt9,xt10
doubleprecision, parameter :: s2 = 1.414213562373095, &
                              s2o2 = 0.707106781186547, &
                              zero = 0.d0
doubleprecision, dimension(0:5) :: eij
doubleprecision, dimension(0:ngllx-1,0:nglly-1,0:ngllz-1) :: t1,t5,t8,t2,t6,t9
doubleprecision, dimension(0:ngllx-1,0:nglly-1,0:ngllz-1) :: t3,t7,t10
doubleprecision, dimension(0:5, 0:5, 0:ngllx-1,0:nglly-1,0:ngllz-1) :: C


k = 0
do i = 0,5
    do j = i,5
        C(i,j,:,:,:) = Cij(k,:,:,:)
        k = k + 1
    enddo
enddo
do i = 1,5
    do j = 0,i-1
        C(i,j,:,:,:) = C(j,i,:,:,:)
    enddo
enddo

do k = 0,ngllz-1
 do j = 0,nglly-1
  do i = 0,ngllx-1

      eij(0) = DXX(i,j,k)
      eij(1) = DYY(i,j,k)
      eij(2) = DZZ(i,j,k)
      eij(3) = s2o2*(DYZ(i,j,k)+DZY(i,j,k))
      eij(4) = s2o2*(DXZ(i,j,k)+DZX(i,j,k))
      eij(5) = s2o2*(DXY(i,j,k)+DYX(i,j,k))

      sxx = 0.
      syy = 0.
      szz = 0.
      syz = 0.
      sxz = 0.
      sxy = 0.
      do l = 0,5
          sxx = sxx + C(0,l,i,j,k)*eij(l)
          syy = syy + C(1,l,i,j,k)*eij(l)
          szz = szz + C(2,l,i,j,k)*eij(l)
          syz = syz + C(3,l,i,j,k)*eij(l)
          sxz = sxz + C(4,l,i,j,k)*eij(l)
          sxy = sxy + C(5,l,i,j,k)*eij(l)
      enddo
      syz=syz/s2
      sxz=sxz/s2
      sxy=sxy/s2

!
!=====================
!       FX 
!=====================
!
      xt1 = sxx * xi1(i,j,k)
      xt2 = sxx * et1(i,j,k)
      xt3 = sxx * ga1(i,j,k)
!---
      xt1 = xt1 + sxy * xi2(i,j,k)
      xt2 = xt2 + sxy * et2(i,j,k)
      xt3 = xt3 + sxy * ga2(i,j,k)
!---
      xt1 = xt1 + sxz * xi3(i,j,k)
      xt2 = xt2 + sxz * et3(i,j,k)
      xt3 = xt3 + sxz * ga3(i,j,k)
!
!=====================
!       FY 
!=====================
!
      xt5 = syy * xi2(i,j,k)
      xt6 = syy * et2(i,j,k)
      xt7 = syy * ga2(i,j,k)
!---
      xt5 = xt5 + sxy * xi1(i,j,k)
      xt6 = xt6 + sxy * et1(i,j,k)
      xt7 = xt7 + sxy * ga1(i,j,k)
!---
      xt5 = xt5 + syz * xi3(i,j,k)
      xt6 = xt6 + syz * et3(i,j,k)
      xt7 = xt7 + syz * ga3(i,j,k)
!
!=====================
!       FZ 
!=====================
!
      xt8  = szz * xi3(i,j,k)
      xt9  = szz * et3(i,j,k)
      xt10 = szz * ga3(i,j,k)
!---
      xt8  = xt8  + sxz * xi1(i,j,k)
      xt9  = xt9  + sxz * et1(i,j,k)
      xt10 = xt10 + sxz * ga1(i,j,k)
!---
      xt8  = xt8  + syz * xi2(i,j,k)
      xt9  = xt9  + syz * et2(i,j,k)
      xt10 = xt10 + syz * ga2(i,j,k)
!---

!
!- Multiplication par le Jacobien et le poids d'integration
!
      t4 = jac(i,j,k) * poidsx(i)
      xt1  =  xt1 * t4
      xt5  =  xt5 * t4
      xt8  =  xt8 * t4

      t4 = jac(i,j,k) * poidsy(j)
      xt2  =  xt2 * t4
      xt6  =  xt6 * t4
      xt9  =  xt9 * t4

      t4 = jac(i,j,k) * poidsz(k)
      xt3  =  xt3 * t4
      xt7  =  xt7 * t4
      xt10 = xt10 * t4


      t1(i,j,k) = xt1
      t5(i,j,k) = xt5
      t8(i,j,k) = xt8

      t2(j,i,k) = xt2
      t6(j,i,k) = xt6
      t9(j,i,k) = xt9

      t3(k,i,j) = xt3
      t7(k,i,j) = xt7
      t10(k,i,j) = xt10

  enddo
 enddo
enddo

!
!- Multiplication par la matrice de derivation puis par les poids
!

!=-=-=-=-=-=-=-=-=-=-
   do k = 0,ngllz-1
    do j = 0,nglly-1
     do i = 0,ngllx-1
!=-=-=-=-=-=-=-=-=-=-
!
	t11 = poidsy(j) * poidsz(k)
	t12 = poidsx(i) * poidsz(k)
	t13 = poidsx(i) * poidsy(j)
!
	t41 = zero
	t42 = zero
	t43 = zero
	t51 = zero
	t52 = zero
	t53 = zero
	t61 = zero
	t62 = zero
	t63 = zero
!
	do l = 0,ngllx-1
	   t41 = t41 + dx(l,i) * t1(l,j,k)
	   t42 = t42 + dx(l,i) * t5(l,j,k)
	   t43 = t43 + dx(l,i) * t8(l,j,k)
	enddo

	do l = 0,nglly-1
	   t51 = t51 + dy(l,j) * t2(l,i,k)
	   t52 = t52 + dy(l,j) * t6(l,i,k)
	   t53 = t53 + dy(l,j) * t9(l,i,k)
	enddo
! FX
	F1 = t41*t11 + t51*t12
! FY
	F2 = t42*t11 + t52*t12
! FZ
	F3 = t43*t11 + t53*t12
!
!
	do l = 0,ngllz-1
	   t61 = t61 + dz(l,k) * t3(l,i,j)
	   t62 = t62 + dz(l,k) * t7(l,i,j)
	   t63 = t63 + dz(l,k) * t10(l,i,j)
	enddo
! FX
	F1 = F1 + t61*t13
! FY
	F2 = F2 + t62*t13
! FZ
	F3 = F3 + t63*t13
!
	Fox(i,j,k) = F1
	Foy(i,j,k) = F2
	Foz(i,j,k) = F3

!=-=-=-=-=-=-=-=-=-=-
     enddo
    enddo
   enddo
!=-=-=-=-=-=-=-=-=-=-

end subroutine calcul_forces_aniso
!------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------
subroutine calcul_forces_att(Fox,Foy,Foz, xi1,xi2,xi3,et1,et2,et3,ga1,ga2,ga3, &
                             dx,dy,dz, jac, poidsx,poidsy,poidsz, DXX,DXY,DXZ,DYX,DYY,DYZ,DZX,DZY,DZZ, &
                             mu_,la_, ngllx,nglly,ngllz, n_solid, onemSbeta_, R_xx_,R_yy_,R_xy_,R_xz_,R_yz_)


use sdomains

implicit none

integer, intent(in) :: ngllx,nglly,ngllz, n_solid
doubleprecision, dimension(0:n_solid-1,0:ngllx-1,0:nglly-1,0:ngllz-1), intent(in) :: R_xx_,R_yy_,R_xy_,R_xz_,R_yz_
doubleprecision, dimension(0:ngllx-1,0:nglly-1,0:ngllz-1), intent(in) :: xi1,xi2,xi3, et1,et2,et3, &
                                                                         ga1,ga2,ga3, jac, mu_,la_, &
                                                                         DXX,DXY,DXZ,DYX,DYY,DYZ,DZX,DZY,DZZ, &
                                                                         onemSbeta_
doubleprecision, dimension(0:ngllx-1,0:nglly-1,0:ngllz-1), intent(out) :: Fox,Foz,Foy
doubleprecision, dimension(0:ngllx-1,0:ngllx-1), intent(in) :: dx
doubleprecision, dimension(0:nglly-1,0:nglly-1), intent(in) :: dy
doubleprecision, dimension(0:ngllz-1,0:ngllz-1), intent(in) :: dz
doubleprecision, dimension(0:ngllx-1), intent(in) :: poidsx
doubleprecision, dimension(0:nglly-1), intent(in) :: poidsy
doubleprecision, dimension(0:ngllz-1), intent(in) :: poidsz

doubleprecision :: sxx,sxy,sxz,syy,syz,szz,xdiv,t4,F1,F2,F3
doubleprecision :: t41,t42,t43,t11,t51,t52,t53,t12,t61,t62,t63,t13
doubleprecision :: xt1,xt2,xt3,xt5,xt6,xt7,xt8,xt9,xt10
doubleprecision, parameter :: deuxtiers = 0.666666666666667, &
                              quatretiers = 1.333333333333337, &
                              s2 = 1.414213562373095, &
                              s2o2 = 0.707106781186547, &
                              zero = 0.d0, two = 2.d0
integer :: i,j,k,l, i_sls, jj, rank, ier
doubleprecision, dimension(0:ngllx-1,0:nglly-1,0:ngllz-1) :: xmu,xla,xla2mu,kappa
doubleprecision, dimension(0:ngllx-1,0:nglly-1,0:ngllz-1) :: t1,t5,t8,t2,t6,t9
doubleprecision, dimension(0:ngllx-1,0:nglly-1,0:ngllz-1) :: t3,t7,t10


xmu(:,:,:) = mu_(:,:,:)
xmu(:,:,:) = xmu(:,:,:) * onemSbeta_(:,:,:)
!ATTENTION a l'histoire kappa/lambda !!!
kappa(:,:,:) = la_(:,:,:)
xla (:,:,:) = kappa(:,:,:) - deuxtiers * xmu(:,:,:)
xla2mu(:,:,:) = xla(:,:,:) + two * xmu(:,:,:)

do k = 0,ngllz-1
 do j = 0,nglly-1
  do i = 0,ngllx-1

      sxx = xla2mu(i,j,k) * DXX(i,j,k) +		&
	    xla(i,j,k)  * ( DYY(i,j,k) +		&
		     DZZ(i,j,k) )

      sxy = xmu(i,j,k)  *( DXY(i,j,k)  +		&
		    DYX(i,j,k)  )

      sxz = xmu(i,j,k) * ( DXZ(i,j,k)  +		&
		    DZX(i,j,k)  )

      syy = xla2mu(i,j,k) * DYY(i,j,k) +		&
	    xla(i,j,k)  * ( DXX(i,j,k) +		&
		     DZZ(i,j,k) )

      syz = xmu(i,j,k) * ( DYZ(i,j,k)  +		&
		    DZY(i,j,k)  )

      szz = xla2mu(i,j,k) * DZZ(i,j,k) +		&
	    xla(i,j,k)  * ( DXX(i,j,k) +		&
		     DYY(i,j,k) )

      do i_sls = 0,n_solid-1
         sxx = sxx - R_xx_(i_sls,i,j,k)
         syy = syy - R_yy_(i_sls,i,j,k)
         ! ici on utilise le fait que la trace est nulle
         szz = szz + R_xx_(i_sls,i,j,k) + R_yy_(i_sls,i,j,k)
         sxy = sxy - R_xy_(i_sls,i,j,k)
         sxz = sxz - R_xz_(i_sls,i,j,k)
         syz = syz - R_yz_(i_sls,i,j,k)
      enddo
!
!=====================
!       FX 
!=====================
!
      xt1 = sxx * xi1(i,j,k)
      xt2 = sxx * et1(i,j,k)
      xt3 = sxx * ga1(i,j,k)

      xt1 = xt1 + sxy * xi2(i,j,k)
      xt2 = xt2 + sxy * et2(i,j,k)
      xt3 = xt3 + sxy * ga2(i,j,k)

      xt1 = xt1 + sxz * xi3(i,j,k)
      xt2 = xt2 + sxz * et3(i,j,k)
      xt3 = xt3 + sxz * ga3(i,j,k)
!
!=====================
!       FY 
!=====================
!
      xt5 = syy * xi2(i,j,k)
      xt6 = syy * et2(i,j,k)
      xt7 = syy * ga2(i,j,k)

      xt5 = xt5 + sxy * xi1(i,j,k)
      xt6 = xt6 + sxy * et1(i,j,k)
      xt7 = xt7 + sxy * ga1(i,j,k)

      xt5 = xt5 + syz * xi3(i,j,k)
      xt6 = xt6 + syz * et3(i,j,k)
      xt7 = xt7 + syz * ga3(i,j,k)
!
!=====================
!       FZ 
!=====================
!
      xt8  = szz * xi3(i,j,k)
      xt9  = szz * et3(i,j,k)
      xt10 = szz * ga3(i,j,k)

      xt8  = xt8  + sxz * xi1(i,j,k)
      xt9  = xt9  + sxz * et1(i,j,k)
      xt10 = xt10 + sxz * ga1(i,j,k)

      xt8  = xt8  + syz * xi2(i,j,k)
      xt9  = xt9  + syz * et2(i,j,k)
      xt10 = xt10 + syz * ga2(i,j,k)

!
!- Multiplication par le Jacobien et le poids d'integration
!
      t4 = jac(i,j,k) * poidsx(i)
      xt1  =  xt1 * t4
      xt5  =  xt5 * t4
      xt8  =  xt8 * t4

      t4 = jac(i,j,k) * poidsy(j)
      xt2  =  xt2 * t4
      xt6  =  xt6 * t4
      xt9  =  xt9 * t4

      t4 = jac(i,j,k) * poidsz(k)
      xt3  =  xt3 * t4
      xt7  =  xt7 * t4
      xt10 = xt10 * t4


      t1(i,j,k) = xt1
      t5(i,j,k) = xt5
      t8(i,j,k) = xt8

      t2(j,i,k) = xt2
      t6(j,i,k) = xt6
      t9(j,i,k) = xt9

      t3(k,i,j) = xt3
      t7(k,i,j) = xt7
      t10(k,i,j) = xt10

  enddo
 enddo
enddo

!
!- Multiplication par la matrice de derivation puis par les poids
!

!=-=-=-=-=-=-=-=-=-=-
do k = 0,ngllz-1
 do j = 0,nglly-1
  do i = 0,ngllx-1
!=-=-=-=-=-=-=-=-=-=-
!
	t11 = poidsy(j) * poidsz(k)
	t12 = poidsx(i) * poidsz(k)
	t13 = poidsx(i) * poidsy(j)
!
	t41 = zero
	t42 = zero
	t43 = zero
	t51 = zero
	t52 = zero
	t53 = zero
	t61 = zero
	t62 = zero
	t63 = zero
!
	do l = 0,ngllx-1
	   t41 = t41 + dx(l,i) * t1(l,j,k)
	   t42 = t42 + dx(l,i) * t5(l,j,k)
	   t43 = t43 + dx(l,i) * t8(l,j,k)
	enddo

	do l = 0,nglly-1
	   t51 = t51 + dy(l,j) * t2(l,i,k)
	   t52 = t52 + dy(l,j) * t6(l,i,k)
	   t53 = t53 + dy(l,j) * t9(l,i,k)
	enddo
! FX
	F1 = t41*t11 + t51*t12
! FY
	F2 = t42*t11 + t52*t12
! FZ
	F3 = t43*t11 + t53*t12
!
!
	do l = 0,ngllz-1
	   t61 = t61 + dz(l,k) * t3(l,i,j)
	   t62 = t62 + dz(l,k) * t7(l,i,j)
	   t63 = t63 + dz(l,k) * t10(l,i,j)
	enddo
! FX
	F1 = F1 + t61*t13
! FY
	F2 = F2 + t62*t13
! FZ
	F3 = F3 + t63*t13
!
	Fox(i,j,k) = F1
	Foy(i,j,k) = F2
	Foz(i,j,k) = F3

!=-=-=-=-=-=-=-=-=-=-
  enddo
 enddo
enddo
!=-=-=-=-=-=-=-=-=-=-

end subroutine calcul_forces_att
!------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------
subroutine calcul_forces(Fox,Foy,Foz, xi1,xi2,xi3,et1,et2,et3,ga1,ga2,ga3, &
                         dx,dy,dz, jac, poidsx,poidsy,poidsz, DXX,DXY,DXZ,DYX,DYY,DYZ,DZX,DZY,DZZ, &
                         mu_,la_, ngllx,nglly,ngllz)


use sdomains

implicit none

integer, intent(in) :: ngllx,nglly,ngllz
doubleprecision, dimension(0:ngllx-1,0:nglly-1,0:ngllz-1), intent(in) :: xi1,xi2,xi3, et1,et2,et3, &
                                                                         ga1,ga2,ga3, jac, mu_,la_, &
                                                                         DXX,DXY,DXZ,DYX,DYY,DYZ,DZX,DZY,DZZ
doubleprecision, dimension(0:ngllx-1,0:nglly-1,0:ngllz-1), intent(out) :: Fox,Foz,Foy
doubleprecision, dimension(0:ngllx-1,0:ngllx-1), intent(in) :: dx
doubleprecision, dimension(0:nglly-1,0:nglly-1), intent(in) :: dy
doubleprecision, dimension(0:ngllz-1,0:ngllz-1), intent(in) :: dz
doubleprecision, dimension(0:ngllx-1), intent(in) :: poidsx
doubleprecision, dimension(0:nglly-1), intent(in) :: poidsy
doubleprecision, dimension(0:ngllz-1), intent(in) :: poidsz

integer :: i,j,k,l
doubleprecision :: sxx,sxy,sxz,syy,syz,szz,t4,F1,F2,F3
doubleprecision :: t41,t42,t43,t11,t51,t52,t53,t12,t61,t62,t63,t13
doubleprecision :: xt1,xt2,xt3,xt5,xt6,xt7,xt8,xt9,xt10
doubleprecision, parameter :: zero = 0.
doubleprecision, dimension(0:ngllx-1,0:nglly-1,0:ngllz-1) :: xmu,xla,xla2mu
doubleprecision, dimension(0:ngllx-1,0:nglly-1,0:ngllz-1) :: t1,t5,t8,t2,t6,t9
doubleprecision, dimension(0:ngllx-1,0:nglly-1,0:ngllz-1) :: t3,t7,t10


xmu(:,:,:) = mu_(:,:,:)
xla(:,:,:) = la_(:,:,:)
xla2mu(:,:,:) = xla(:,:,:) + 2. * xmu(:,:,:)

do k = 0,ngllz-1
 do j = 0,nglly-1
  do i = 0,ngllx-1

      sxx = xla2mu(i,j,k) * DXX(i,j,k) +		&
	    xla(i,j,k)  * ( DYY(i,j,k) +		&
		     DZZ(i,j,k) )
!---
      sxy = xmu(i,j,k) * ( DXY(i,j,k)  +		&
		    DYX(i,j,k)  )
!---
      sxz = xmu(i,j,k) * ( DXZ(i,j,k)  +		&
		    DZX(i,j,k)  )
!---
      syy = xla2mu(i,j,k) * DYY(i,j,k) +		&
	    xla(i,j,k)  * ( DXX(i,j,k) +		&
		     DZZ(i,j,k) )
!---
      syz = xmu(i,j,k) * ( DYZ(i,j,k)  +		&
		    DZY(i,j,k)  )
!---
      szz = xla2mu(i,j,k) * DZZ(i,j,k) +		&
	    xla(i,j,k)  * ( DXX(i,j,k) +		&
		     DYY(i,j,k) )

!
!=====================
!       FX 
!=====================
!
      xt1 = sxx * xi1(i,j,k)
      xt2 = sxx * et1(i,j,k)
      xt3 = sxx * ga1(i,j,k)
!---
      xt1 = xt1 + sxy * xi2(i,j,k)
      xt2 = xt2 + sxy * et2(i,j,k)
      xt3 = xt3 + sxy * ga2(i,j,k)
!---
      xt1 = xt1 + sxz * xi3(i,j,k)
      xt2 = xt2 + sxz * et3(i,j,k)
      xt3 = xt3 + sxz * ga3(i,j,k)
!
!=====================
!       FY 
!=====================
!
      xt5 = syy * xi2(i,j,k)
      xt6 = syy * et2(i,j,k)
      xt7 = syy * ga2(i,j,k)
!---
      xt5 = xt5 + sxy * xi1(i,j,k)
      xt6 = xt6 + sxy * et1(i,j,k)
      xt7 = xt7 + sxy * ga1(i,j,k)
!---
      xt5 = xt5 + syz * xi3(i,j,k)
      xt6 = xt6 + syz * et3(i,j,k)
      xt7 = xt7 + syz * ga3(i,j,k)
!
!=====================
!       FZ 
!=====================
!
      xt8  = szz * xi3(i,j,k)
      xt9  = szz * et3(i,j,k)
      xt10 = szz * ga3(i,j,k)
!---
      xt8  = xt8  + sxz * xi1(i,j,k)
      xt9  = xt9  + sxz * et1(i,j,k)
      xt10 = xt10 + sxz * ga1(i,j,k)
!---
      xt8  = xt8  + syz * xi2(i,j,k)
      xt9  = xt9  + syz * et2(i,j,k)
      xt10 = xt10 + syz * ga2(i,j,k)
!---

!
!- Multiplication par le Jacobien et le poids d'integration
!
      t4 = jac(i,j,k) * poidsx(i)
      xt1  =  xt1 * t4
      xt5  =  xt5 * t4
      xt8  =  xt8 * t4

      t4 = jac(i,j,k) * poidsy(j)
      xt2  =  xt2 * t4
      xt6  =  xt6 * t4
      xt9  =  xt9 * t4

      t4 = jac(i,j,k) * poidsz(k)
      xt3  =  xt3 * t4
      xt7  =  xt7 * t4
      xt10 = xt10 * t4


      t1(i,j,k) = xt1
      t5(i,j,k) = xt5
      t8(i,j,k) = xt8

      t2(j,i,k) = xt2
      t6(j,i,k) = xt6
      t9(j,i,k) = xt9

      t3(k,i,j) = xt3
      t7(k,i,j) = xt7
      t10(k,i,j) = xt10

  enddo
 enddo
enddo

!
!- Multiplication par la matrice de derivation puis par les poids
!

!=-=-=-=-=-=-=-=-=-=-
   do k = 0,ngllz-1
    do j = 0,nglly-1
     do i = 0,ngllx-1
!=-=-=-=-=-=-=-=-=-=-
!
	t11 = poidsy(j) * poidsz(k)
	t12 = poidsx(i) * poidsz(k)
	t13 = poidsx(i) * poidsy(j)
!
	t41 = zero
	t42 = zero
	t43 = zero
	t51 = zero
	t52 = zero
	t53 = zero
	t61 = zero
	t62 = zero
	t63 = zero
!
	do l = 0,ngllx-1
	   t41 = t41 + dx(l,i) * t1(l,j,k)
	   t42 = t42 + dx(l,i) * t5(l,j,k)
	   t43 = t43 + dx(l,i) * t8(l,j,k)
	enddo

	do l = 0,nglly-1
	   t51 = t51 + dy(l,j) * t2(l,i,k)
	   t52 = t52 + dy(l,j) * t6(l,i,k)
	   t53 = t53 + dy(l,j) * t9(l,i,k)
	enddo
! FX
	F1 = t41*t11 + t51*t12
! FY
	F2 = t42*t11 + t52*t12
! FZ
	F3 = t43*t11 + t53*t12
!
!
	do l = 0,ngllz-1
	   t61 = t61 + dz(l,k) * t3(l,i,j)
	   t62 = t62 + dz(l,k) * t7(l,i,j)
	   t63 = t63 + dz(l,k) * t10(l,i,j)
	enddo
! FX
	F1 = F1 + t61*t13
! FY
	F2 = F2 + t62*t13
! FZ
	F3 = F3 + t63*t13
!
	Fox(i,j,k) = F1
	Foy(i,j,k) = F2
	Foz(i,j,k) = F3

!=-=-=-=-=-=-=-=-=-=-
     enddo
    enddo
   enddo
!=-=-=-=-=-=-=-=-=-=-

end subroutine calcul_forces
