subroutine SourceAdjoint (Tdomain,rg)


use sdomains
use angles
use module_ellipticity

implicit none

include 'mpif.h'

type (domain), intent(inout) :: Tdomain
integer, intent(in) :: rg

integer :: n_src, near_node, i,j,k, n_around_elem, x,y,z, ngllx,nglly,ngllz, code, mat, num, n, i_count
integer , dimension (0:20) :: el_around_node
doubleprecision :: xs,ys,zs, d,dmin, R,colat,long, ct,st,cp,sp, epsil, xa,ya,za, dist_xi,dist_eta,dist_zeta, xi,eta,zeta, f
doubleprecision, dimension(0:2) :: xi_search,eta_search,zeta_search, centre
doubleprecision, dimension(:), allocatable :: distance
doubleprecision, dimension(0:2,0:2) :: tRot, LocInvGrad
doubleprecision, dimension(:,:), allocatable :: coord
logical :: src_hors_domaine


allocate (distance(0:Tdomain%n_proc-1))


if (Tdomain%ellipticity)   call init_ellipticity()


do n_src = 0, Tdomain%n_source_adj-1


 ! On calcule les coordonnees des src dans le chunk de reference
 if (Tdomain%curve) then
    ! Passage de latitude a colatitude dans le chunk reel
    Tdomain%sAdjoint(n_src)%realcolat = 90.d0 - Tdomain%sAdjoint(n_src)%realcolat
    if (Tdomain%sAdjoint(n_src)%reallong<0.d0) Tdomain%sAdjoint(n_src)%reallong = 360.d0 + Tdomain%sAdjoint(n_src)%reallong
    colat = Tdomain%sAdjoint(n_src)%realcolat
    long = Tdomain%sAdjoint(n_src)%reallong
    ! Calcul de la distance au centre
    if (Tdomain%ellipticity) then
        R = Rterre * (1.d0 - get_ellipticity(Rterre)*(3.d0*dcos(pi*colat/180.d0)**2.d0 - 1.d0)/3.d0)
    else
        R = Rterre
    endif
    R = R - Tdomain%sAdjoint(n_src)%depth
    ! Passage de spherique a cartesien dans le chunk reel
    ct = dcos(pi*colat/180.d0)
    st = dsin(pi*colat/180.d0)
    cp = dcos(pi*long/180.d0)
    sp = dsin(pi*long/180.d0)
    xa = R * st * cp
    ya = R * st * sp
    za = R * ct
    ! Passage dans le chunk de reference
    tRot = transpose(Tdomain%rot)
    xs = tRot(0,0)*xa + tRot(0,1)*ya + tRot(0,2)*za
    ys = tRot(1,0)*xa + tRot(1,1)*ya + tRot(1,2)*za
    zs = tRot(2,0)*xa + tRot(2,1)*ya + tRot(2,2)*za
    ! Passage de cartesien a spherique dans le chunk de reference
    call cart2sph(xs,ys,zs,R,colat,long)
    Tdomain%sAdjoint(n_src)%refcolat = colat
    Tdomain%sAdjoint(n_src)%reflong = long
 else
    xs = Tdomain%sAdjoint(n_src)%realcolat
    ys = Tdomain%sAdjoint(n_src)%reallong
    zs = Tdomain%sAdjoint(n_src)%depth
 endif

 ! On trouve le noeud le plus proche de la src
 dmin = huge(dmin)
 do i = 0,Tdomain%n_glob_nodes-1
    d = dsqrt ((Tdomain%Coord_nodes(0,i)-xs)**2 + (Tdomain%Coord_nodes(1,i)-ys)**2 + (Tdomain%Coord_nodes(2,i)-zs)**2)
    if (d <= dmin) then
       dmin = d
       near_node = i
    endif
 enddo

 ! On trouve l'element dans lequel est la src ainsi que le GLL interne le plus proche
 n_around_elem = 0
 do i = 0,Tdomain%n_elem-1
    do j = 0,Tdomain%n_nodes-1
       if (Tdomain%specel(i)%Control_nodes(j)==near_node) then
          el_around_node(n_around_elem) = i
          n_around_elem = n_around_elem + 1
       endif
    enddo
 enddo
 dmin = huge(dmin)
 do i = 0,n_around_elem-1
    ngllx = Tdomain%specel(el_around_node(i))%ngllx
    nglly = Tdomain%specel(el_around_node(i))%nglly
    ngllz = Tdomain%specel(el_around_node(i))%ngllz
    do x = 1,ngllx-2
     do y = 1,nglly-2
      do z = 1,ngllz-2
         j = Tdomain%specel(el_around_node(i))%Iglobnum(x,y,z)
         d = dsqrt ((Tdomain%Globcoord(0,j)-xs)**2 + (Tdomain%Globcoord(1,j)-ys)**2 + (Tdomain%Globcoord(2,j)-zs)**2)
         if (d <= dmin) then
            dmin = d
            Tdomain%sAdjoint(n_src)%elem = el_around_node(i)
            Tdomain%sAdjoint(n_src)%gll(0) = x
            Tdomain%sAdjoint(n_src)%gll(1) = y
            Tdomain%sAdjoint(n_src)%gll(2) = z
         endif
      enddo
     enddo
    enddo
 enddo

 ! On trouve le processeur qui contient la src
 call mpi_allgather(dmin,1,mpi_double_precision,distance,1,mpi_double_precision,mpi_comm_world,code)
 do i = 0,Tdomain%n_proc-1
    if (distance(i) <= dmin) then
       dmin = distance(i)
       Tdomain%sAdjoint(n_src)%proc = i
    endif
 enddo


 if (rg==Tdomain%sAdjoint(n_src)%proc) then

    ! Dichotomie pour passer de S(x,y,z) a S(xi,eta,zeta)
    n = Tdomain%sAdjoint(n_src)%elem
    mat = Tdomain%specel(n)%mat_index
    src_hors_domaine = .false.
    if (Tdomain%sSubDomain(mat)%material_type == "P") then
        print *,"ADJOINT SOURCE ", n_src, " OUT OF THE DOMAIN"
        src_hors_domaine = .true.
    endif
    x = Tdomain%sAdjoint(n_src)%gll(0)
    y = Tdomain%sAdjoint(n_src)%gll(1)
    z = Tdomain%sAdjoint(n_src)%gll(2)
    dist_xi = Tdomain%sSubdomain(mat)%GLLcx(x+1) - Tdomain%sSubdomain(mat)%GLLcx(x-1)
    dist_eta = Tdomain%sSubdomain(mat)%GLLcy(y+1) - Tdomain%sSubdomain(mat)%GLLcy(y-1)
    dist_zeta = Tdomain%sSubdomain(mat)%GLLcz(z+1) - Tdomain%sSubdomain(mat)%GLLcz(z-1)
    do i = 0,2
       xi_search(i) = Tdomain%sSubdomain(mat)%GLLcx(x-1) + (i+1)*dist_xi/4
       eta_search(i) = Tdomain%sSubdomain(mat)%GLLcy(y-1) + (i+1)*dist_eta/4
       zeta_search(i) = Tdomain%sSubdomain(mat)%GLLcz(z-1) + (i+1)*dist_zeta/4
    enddo
    allocate (coord(0:Tdomain%n_nodes-1,0:2))
    do i = 0,Tdomain%n_nodes-1
        j = Tdomain%specel(n)%Control_Nodes(i)
        coord(i,0:2) = Tdomain%Coord_Nodes(0:2,j)
    enddo
    dmin = max ( distance3d(coord(0,0:2),coord(6,0:2)), distance3d(coord(1,0:2),coord(7,0:2)), &
                 distance3d(coord(2,0:2),coord(4,0:2)), distance3d(coord(3,0:2),coord(5,0:2)) )
    epsil = dmin/10000.
    num = 0
    dicho : do while (dmin > epsil)
       do i = 0,2
        do j = 0,2
         do k = 0,2
            ! calcul des coords dans l'espace physique
            xi = xi_search(i);   eta = eta_search(j);   zeta = zeta_search(k)
            if (Tdomain%n_nodes==27) then
               xa = 0;   ya = 0;   za = 0;
               do i_count = 0,Tdomain%n_nodes-1
                   f = Comp_shapefunc(i_count,xi,eta,zeta)
                   xa = xa + coord(i_count,0)*f
                   ya = ya + coord(i_count,1)*f
                   za = za + coord(i_count,2)*f
               enddo
            else if (Tdomain%n_nodes==8) then
               xa = 0.125 * (coord(0,0)*(1-xi)*(1-eta)*(1-zeta) + coord(1,0)*(1+xi)*(1-eta)*(1-zeta) + &
                             coord(2,0)*(1+xi)*(1+eta)*(1-zeta) + coord(3,0)*(1-xi)*(1+eta)*(1-zeta) + &
                             coord(4,0)*(1-xi)*(1-eta)*(1+zeta) + coord(5,0)*(1+xi)*(1-eta)*(1+zeta) + &
                             coord(6,0)*(1+xi)*(1+eta)*(1+zeta) + coord(7,0)*(1-xi)*(1+eta)*(1+zeta))
               ya = 0.125 * (coord(0,1)*(1-xi)*(1-eta)*(1-zeta) + coord(1,1)*(1+xi)*(1-eta)*(1-zeta) + &
                             coord(2,1)*(1+xi)*(1+eta)*(1-zeta) + coord(3,1)*(1-xi)*(1+eta)*(1-zeta) + &
                             coord(4,1)*(1-xi)*(1-eta)*(1+zeta) + coord(5,1)*(1+xi)*(1-eta)*(1+zeta) + &
                             coord(6,1)*(1+xi)*(1+eta)*(1+zeta) + coord(7,1)*(1-xi)*(1+eta)*(1+zeta))
               za = 0.125 * (coord(0,2)*(1-xi)*(1-eta)*(1-zeta) + coord(1,2)*(1+xi)*(1-eta)*(1-zeta) + &
                             coord(2,2)*(1+xi)*(1+eta)*(1-zeta) + coord(3,2)*(1-xi)*(1+eta)*(1-zeta) + &
                             coord(4,2)*(1-xi)*(1-eta)*(1+zeta) + coord(5,2)*(1+xi)*(1-eta)*(1+zeta) + &
                             coord(6,2)*(1+xi)*(1+eta)*(1+zeta) + coord(7,2)*(1-xi)*(1+eta)*(1+zeta))
            endif
            ! calcul de la distance a la source
            d = dsqrt ((xa-xs)**2 + (ya-ys)**2 + (za-zs)**2)
            if (d < dmin) then
               dmin = d
               centre(0) = xi;   centre(1) = eta;   centre(2) = zeta
            endif
         enddo
        enddo
       enddo
       dist_xi = dist_xi/2.
       dist_eta = dist_eta/2.
       dist_zeta = dist_zeta/2.
       do i = 0,2
          xi_search(i) = centre(0) + (i-1)*dist_xi/4.
          eta_search(i) = centre(1) + (i-1)*dist_eta/4.
          zeta_search(i) = centre(2) + (i-1)*dist_zeta/4.
       enddo
       num = num + 1
       if (num>20) exit dicho
    enddo dicho
    Tdomain%sAdjoint(n_src)%refcoord = centre

    ! Calcul de InvGrad en S(xi,eta,zeta)
    xi = centre(0);   eta = centre(1);   zeta = centre(2)
    if (Tdomain%n_nodes==27) then
       LocInvGrad = 0.
       do i = 0,Tdomain%n_nodes-1
          do j = 0,2
              f = Comp_derivshapefunc(i,xi,eta,zeta,j)
              LocInvGrad(j,0) = LocInvGrad(j,0) + coord(i,0)*f
              LocInvGrad(j,1) = LocInvGrad(j,1) + coord(i,1)*f
              LocInvGrad(j,2) = LocInvGrad(j,2) + coord(i,2)*f
          enddo
       enddo
    else if (Tdomain%n_nodes==8) then
       LocInvGrad(0,0) = 0.125 * ((coord(1,0)-coord(0,0))*(1-eta)*(1-zeta) + (coord(2,0)-coord(3,0))*(1+eta)*(1-zeta) + &
                                  (coord(5,0)-coord(4,0))*(1-eta)*(1+zeta) + (coord(6,0)-coord(7,0))*(1+eta)*(1+zeta))
       LocInvGrad(1,0) = 0.125 * ((coord(3,0)-coord(0,0))*(1-xi)*(1-zeta) + (coord(2,0)-coord(1,0))*(1+xi)*(1-zeta) + &
                                  (coord(7,0)-coord(4,0))*(1-xi)*(1+zeta) + (coord(6,0)-coord(5,0))*(1+xi)*(1+zeta))
       LocInvGrad(2,0) = 0.125 * ((coord(4,0)-coord(0,0))*(1-xi)*(1-eta) + (coord(5,0)-coord(1,0))*(1+xi)*(1-eta) + &
                                  (coord(7,0)-coord(3,0))*(1-xi)*(1+eta) + (coord(6,0)-coord(2,0))*(1+xi)*(1+eta))
       LocInvGrad(0,1) = 0.125 * ((coord(1,1)-coord(0,1))*(1-eta)*(1-zeta) + (coord(2,1)-coord(3,1))*(1+eta)*(1-zeta) + &
                                  (coord(5,1)-coord(4,1))*(1-eta)*(1+zeta) + (coord(6,1)-coord(7,1))*(1+eta)*(1+zeta))
       LocInvGrad(1,1) = 0.125 * ((coord(3,1)-coord(0,1))*(1-xi)*(1-zeta) + (coord(2,1)-coord(1,1))*(1+xi)*(1-zeta) + &
                                  (coord(7,1)-coord(4,1))*(1-xi)*(1+zeta) + (coord(6,1)-coord(5,1))*(1+xi)*(1+zeta))
       LocInvGrad(2,1) = 0.125 * ((coord(4,1)-coord(0,1))*(1-xi)*(1-eta) + (coord(5,1)-coord(1,1))*(1+xi)*(1-eta) + &
                                  (coord(7,1)-coord(3,1))*(1-xi)*(1+eta) + (coord(6,1)-coord(2,1))*(1+xi)*(1+eta))
       LocInvGrad(0,2) = 0.125 * ((coord(1,2)-coord(0,2))*(1-eta)*(1-zeta) + (coord(2,2)-coord(3,2))*(1+eta)*(1-zeta) + &
                                  (coord(5,2)-coord(4,2))*(1-eta)*(1+zeta) + (coord(6,2)-coord(7,2))*(1+eta)*(1+zeta))
       LocInvGrad(1,2) = 0.125 * ((coord(3,2)-coord(0,2))*(1-xi)*(1-zeta) + (coord(2,2)-coord(1,2))*(1+xi)*(1-zeta) + &
                                  (coord(7,2)-coord(4,2))*(1-xi)*(1+zeta) + (coord(6,2)-coord(5,2))*(1+xi)*(1+zeta))
       LocInvGrad(2,2) = 0.125 * ((coord(4,2)-coord(0,2))*(1-xi)*(1-eta) + (coord(5,2)-coord(1,2))*(1+xi)*(1-eta) + &
                                  (coord(7,2)-coord(3,2))*(1-xi)*(1+eta) + (coord(6,2)-coord(2,2))*(1+xi)*(1+eta))
    endif
    call invert_3d (LocInvGrad, R)
    Tdomain%sAdjoint(n_src)%InvGrad(0:2,0:2) = LocInvGrad(0:2,0:2)
    deallocate (coord)

 else

    deallocate (Tdomain%sAdjoint(n_src)%timefunc)

 endif


 call mpi_bcast(src_hors_domaine,1,mpi_logical,Tdomain%sAdjoint(n_src)%proc,mpi_comm_world,code)
 if (src_hors_domaine .eqv. .true.)   stop


enddo


deallocate (distance)


return
end subroutine SourceAdjoint
