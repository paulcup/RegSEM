subroutine ReceiverPosition (Tdomain,rg)


use sdomains
use angles 
use module_ellipticity

implicit none

include 'mpif.h'

type (domain), intent (inout) :: Tdomain
integer, intent(in) :: rg

integer :: n_rcp,n_src, near_node, n_around_elem, i,j,k, x,y,z, ngllx,nglly,ngllz, code, mat, num, n, a, i_count
integer , dimension (0:20) :: el_around_node
doubleprecision :: xs,ys,zs, d,dmin, R,colat,long, ct,st,cp,sp, epsil, coeff, xa,ya,za, &
                   dist_xi,dist_eta,dist_zeta, xi,eta,zeta, f, cosdelta,sindelta
doubleprecision, dimension(0:2) :: xi_search,eta_search,zeta_search, centre
doubleprecision, dimension(:), allocatable :: distance
doubleprecision, dimension(0:2,0:2) :: tRot
doubleprecision, dimension(:,:), allocatable :: coord
logical :: rcp_hors_domaine


allocate (distance(0:Tdomain%n_proc-1))


do n_rcp = 0, Tdomain%n_receivers-1


 ! On calcule les coordonnees du rcp dans le chunk de reference
 if (Tdomain%curve) then
    ! Passage de latitude a colatitude dans le chunk reel
    Tdomain%sReceiver(n_rcp)%realcolat = 90.d0 - Tdomain%sReceiver(n_rcp)%realcolat
    if (Tdomain%sReceiver(n_rcp)%reallong<0.d0) Tdomain%sReceiver(n_rcp)%reallong = 360.d0 + Tdomain%sReceiver(n_rcp)%reallong
    colat = Tdomain%sReceiver(n_rcp)%realcolat
    long = Tdomain%sReceiver(n_rcp)%reallong
    ! Calcul de la distance au centre
    if (Tdomain%ellipticity) then
        R = Rterre * (1.d0 - get_ellipticity(Rterre)*(3.d0*dcos(pi*colat/180.d0)**2 - 1.d0)/3.d0)
    else
        R = Rterre
    endif
    R = R + Tdomain%sReceiver(n_rcp)%elevation
    ! Matrice de passage de spherique (r,theta,phi) a cartesien (x,y,z) dans le chunk reel.
    ! theta pointe vers le nord donc c'est aussi la matrice qui permet d'obtenir les coordonnees Z/N/E dans le chunk reel.
    ct = dcos(pi*colat/180.d0)
    st = dsin(pi*colat/180.d0)
    cp = dcos(pi*long/180.d0)
    sp = dsin(pi*long/180.d0)
    xa = R * st * cp
    ya = R * st * sp
    za = R * ct
    Tdomain%sReceiver(n_rcp)%Pass(0,0) = st*cp; Tdomain%sReceiver(n_rcp)%Pass(0,1) = st*sp; Tdomain%sReceiver(n_rcp)%Pass(0,2) = ct
    Tdomain%sReceiver(n_rcp)%Pass(1,0) =-ct*cp; Tdomain%sReceiver(n_rcp)%Pass(1,1) =-ct*sp; Tdomain%sReceiver(n_rcp)%Pass(1,2) = st
    Tdomain%sReceiver(n_rcp)%Pass(2,0) = -sp  ; Tdomain%sReceiver(n_rcp)%Pass(2,1) = cp   ; Tdomain%sReceiver(n_rcp)%Pass(2,2) = 0.d0
    ! Passage dans le chunk de reference
    tRot = transpose(Tdomain%rot)
    xs = tRot(0,0)*xa + tRot(0,1)*ya + tRot(0,2)*za
    ys = tRot(1,0)*xa + tRot(1,1)*ya + tRot(1,2)*za
    zs = tRot(2,0)*xa + tRot(2,1)*ya + tRot(2,2)*za
    ! Passage de cartesien a spherique dans le chunk de reference
    call cart2sph(xs,ys,zs,R,colat,long)
    Tdomain%sReceiver(n_rcp)%refcolat = colat
    Tdomain%sReceiver(n_rcp)%reflong = long
    ! Matrice de passage du systeme cartesien au systeme Z/N/E dans le chunk de ref
    ct = dcos(colat)
    st = dsin(colat)
    cp = dcos(long)
    sp = dsin(long)
    Tdomain%sReceiver(n_rcp)%Passage(0,0) = st*cp; Tdomain%sReceiver(n_rcp)%Passage(0,1) = st*sp; Tdomain%sReceiver(n_rcp)%Passage(0,2) = ct
    Tdomain%sReceiver(n_rcp)%Passage(1,0) =-ct*cp; Tdomain%sReceiver(n_rcp)%Passage(1,1) =-ct*sp; Tdomain%sReceiver(n_rcp)%Passage(1,2) = st
    Tdomain%sReceiver(n_rcp)%Passage(2,0) = -sp  ; Tdomain%sReceiver(n_rcp)%Passage(2,1) = cp   ; Tdomain%sReceiver(n_rcp)%Passage(2,2) = 0.d0
    if (Tdomain%comp_rot) then   ! Calcul de l'angle de passage du systeme Z/N/E au systeme Z/R/T
       call arc_gd_cercle(long,Tdomain%sSource(0)%reflong,colat,Tdomain%sSource(0)%refcolat,cosdelta,sindelta)
       Tdomain%sReceiver(n_rcp)%cosgamma = (dcos(Tdomain%sSource(0)%refcolat) - cosdelta*ct) / (sindelta*st)
       if (dabs(Tdomain%sReceiver(n_rcp)%cosgamma)<1.d0) then
          Tdomain%sReceiver(n_rcp)%singamma = dsqrt(1-Tdomain%sReceiver(n_rcp)%cosgamma**2)
       else
          Tdomain%sReceiver(n_rcp)%singamma = 0.d0
       endif
    endif
 else
    xs = Tdomain%sReceiver(n_rcp)%realcolat
    ys = Tdomain%sReceiver(n_rcp)%reallong
    zs = Tdomain%sReceiver(n_rcp)%elevation
 endif

 ! On trouve le noeud le plus proche du rcp
 dmin = huge(dmin)
 do i = 0,Tdomain%n_glob_nodes-1
    d = dsqrt ((Tdomain%Coord_nodes(0,i)-xs)**2 + (Tdomain%Coord_nodes(1,i)-ys)**2 + (Tdomain%Coord_nodes(2,i)-zs)**2)
    if (d <= dmin) then
       dmin = d
       near_node = i
    endif
 enddo

 ! On trouve l'element dans lequel est le rcp ainsi que le GLL interne le plus proche
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
             Tdomain%sReceiver(n_rcp)%elem = el_around_node(i)
             Tdomain%sReceiver(n_rcp)%gll(0) = x
             Tdomain%sReceiver(n_rcp)%gll(1) = y
             Tdomain%sReceiver(n_rcp)%gll(2) = z
         endif
      enddo
     enddo
    enddo
 enddo

 ! On trouve le processeur qui contient le rcp
 call mpi_allgather(dmin,1,mpi_double_precision,distance,1,mpi_double_precision,mpi_comm_world,code)
 do i = 0,Tdomain%n_proc-1
    if (distance(i) <= dmin) then
       dmin = distance(i)
       Tdomain%sReceiver(n_rcp)%proc = i
    endif
 enddo


 if (rg==Tdomain%sReceiver(n_rcp)%proc) then

    ! Dichotomie pour passer de S(x,y,z) a S(xi,eta,zeta)
    n = Tdomain%sReceiver(n_rcp)%elem
    mat = Tdomain%specel(n)%mat_index
    rcp_hors_domaine = .false.
    if (Tdomain%sSubDomain(mat)%material_type == "P") then
        print *,"RECEIVER ", n_rcp, " OUT OF THE DOMAIN"
        rcp_hors_domaine = .true.
    endif
    x = Tdomain%sReceiver(n_rcp)%gll(0)
    y = Tdomain%sReceiver(n_rcp)%gll(1)
    z = Tdomain%sReceiver(n_rcp)%gll(2)
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
    epsil = dmin/10000.d0
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
       dist_xi = dist_xi/2.d0
       dist_eta = dist_eta/2.d0
       dist_zeta = dist_zeta/2.d0
       do a = 0,2
          xi_search(a) = centre(0) + (a-1)*dist_xi/4.d0
          eta_search(a) = centre(1) + (a-1)*dist_eta/4.d0
          zeta_search(a) = centre(2) + (a-1)*dist_zeta/4.d0
       enddo
       num = num + 1
       if (num>20) exit dicho
    enddo dicho
    deallocate (coord)

    ! Calcul en S(xi,eta,zeta) des polynomes de Lagrange associes aux GLLs
    ngllx = Tdomain%specel(n)%ngllx
    nglly = Tdomain%specel(n)%nglly
    ngllz = Tdomain%specel(n)%ngllz
    allocate (Tdomain%sReceiver(n_rcp)%coeff(0:ngllx-1,0:nglly-1,0:ngllz-1,0:2))
    allocate (Tdomain%sReceiver(n_rcp)%pol(0:ngllx-1,0:nglly-1,0:ngllz-1))
    Tdomain%sReceiver(n_rcp)%pol = 1
    do x = 0,ngllx-1
     do y = 0,nglly-1
      do z = 0,ngllz-1
         xi = Tdomain%sSubdomain(mat)%GLLcx(x)
         eta = Tdomain%sSubdomain(mat)%GLLcy(y)
         zeta = Tdomain%sSubdomain(mat)%GLLcz(z)
         do a = 0,ngllx-1
            if (a/=x) then
               xa = Tdomain%sSubdomain(mat)%GLLcx(a)
               Tdomain%sReceiver(n_rcp)%pol(x,y,z) = Tdomain%sReceiver(n_rcp)%pol(x,y,z) * (centre(0)-xa)/(xi-xa)
            endif
         enddo
         do a = 0,nglly-1
            if (a/=y) then
               ya = Tdomain%sSubdomain(mat)%GLLcy(a)
               Tdomain%sReceiver(n_rcp)%pol(x,y,z) = Tdomain%sReceiver(n_rcp)%pol(x,y,z) * (centre(1)-ya)/(eta-ya)
            endif
         enddo
         do a = 0,ngllz-1
            if (a/=z) then
               za = Tdomain%sSubdomain(mat)%GLLcz(a)
               Tdomain%sReceiver(n_rcp)%pol(x,y,z) = Tdomain%sReceiver(n_rcp)%pol(x,y,z) * (centre(2)-za)/(zeta-za)
            endif
         enddo
      enddo
     enddo
    enddo

 endif


 call mpi_bcast(rcp_hors_domaine,1,mpi_logical,Tdomain%sReceiver(n_rcp)%proc,mpi_comm_world,code)
 if (rcp_hors_domaine .eqv. .true.) then
    stop
 endif


enddo


deallocate (distance)


end subroutine ReceiverPosition
