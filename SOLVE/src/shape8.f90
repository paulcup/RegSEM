subroutine shape8(Tdomain)


use sdomains

implicit none

type(domain),target, intent (INOUT) :: Tdomain

integer :: n, i_aus, ngllx,nglly,ngllz, mat, i,j,k, ipoint, n_face
doubleprecision :: x0,y0,z0,x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,x5,y5,z5,x6,y6,z6,x7,y7,z7
doubleprecision :: xi,eta,zeta, xp,yp,zp, Jac
doubleprecision, dimension (0:2,0:2) :: LocInvGrad

!---------------------------------------------------------------------------------------------------------------
! Shape functions are derived from "The finite element displayed"
! by Dhatt, G. and Touzot, G.
! John Wiley and Sons, 1984
! --------------------------------------------------------------------------------------------------------------

allocate (Tdomain%GlobCoord(0:2,0:Tdomain%n_glob_points-1))

do n = 0,Tdomain%n_elem - 1
   i_aus = Tdomain%specel(n)%Control_Nodes(0);  x0 = Tdomain%Coord_Nodes(0,i_aus);  
               y0 = Tdomain%Coord_Nodes(1,i_aus);             z0 = Tdomain%Coord_Nodes(2,i_aus)
   i_aus = Tdomain%specel(n)%Control_Nodes(1);  x1 = Tdomain%Coord_Nodes(0,i_aus);  
               y1 = Tdomain%Coord_Nodes(1,i_aus);             z1 = Tdomain%Coord_Nodes(2,i_aus)
   i_aus = Tdomain%specel(n)%Control_Nodes(2);  x2 = Tdomain%Coord_Nodes(0,i_aus);  
               y2 = Tdomain%Coord_Nodes(1,i_aus);             z2 = Tdomain%Coord_Nodes(2,i_aus)
   i_aus = Tdomain%specel(n)%Control_Nodes(3);  x3 = Tdomain%Coord_Nodes(0,i_aus);  
               y3 = Tdomain%Coord_Nodes(1,i_aus);             z3 = Tdomain%Coord_Nodes(2,i_aus)
   i_aus = Tdomain%specel(n)%Control_Nodes(4);  x4 = Tdomain%Coord_Nodes(0,i_aus);  
               y4 = Tdomain%Coord_Nodes(1,i_aus);             z4 = Tdomain%Coord_Nodes(2,i_aus)
   i_aus = Tdomain%specel(n)%Control_Nodes(5);  x5 = Tdomain%Coord_Nodes(0,i_aus);  
               y5 = Tdomain%Coord_Nodes(1,i_aus);             z5 = Tdomain%Coord_Nodes(2,i_aus)
   i_aus = Tdomain%specel(n)%Control_Nodes(6);  x6 = Tdomain%Coord_Nodes(0,i_aus);  
               y6 = Tdomain%Coord_Nodes(1,i_aus);             z6 = Tdomain%Coord_Nodes(2,i_aus)
   i_aus = Tdomain%specel(n)%Control_Nodes(7);  x7 = Tdomain%Coord_Nodes(0,i_aus);  
               y7 = Tdomain%Coord_Nodes(1,i_aus);             z7 = Tdomain%Coord_Nodes(2,i_aus)

   ngllx = Tdomain%specel(n)%ngllx
   nglly = Tdomain%specel(n)%nglly
   ngllz = Tdomain%specel(n)%ngllz
   mat = Tdomain%specel(n)%mat_index

   allocate (Tdomain%specel(n)%Jacob(0:ngllx-1,0:nglly-1,0:ngllz-1) ) 
   allocate (Tdomain%specel(n)%InvGrad(0:ngllx-1,0:nglly-1,0:ngllz-1,0:2,0:2) ) 

   do k = 0,ngllz - 1
       zeta =  Tdomain%sSubdomain(mat)%GLLcz (k)
       do j = 0,nglly - 1
           eta =  Tdomain%sSubdomain(mat)%GLLcy (j)    
           do i = 0,ngllx - 1
               xi = Tdomain%sSubdomain(mat)%GLLcx (i)

               ipoint = Tdomain%specel(n)%Iglobnum(i,j,k)
               xp = 0.125 * (x0*(1-xi)*(1-eta)*(1-zeta) + x1*(1+xi)*(1-eta)*(1-zeta) + x2*(1+xi)*(1+eta)*(1-zeta) + &
                    x3*(1-xi)*(1+eta)*(1-zeta) + x4*(1-xi)*(1-eta)*(1+zeta) + x5*(1+xi)*(1-eta)*(1+zeta) + & 
                    x6*(1+xi)*(1+eta)*(1+zeta) + x7*(1-xi)*(1+eta)*(1+zeta)) 

               yp = 0.125 * (y0*(1-xi)*(1-eta)*(1-zeta) + y1*(1+xi)*(1-eta)*(1-zeta) + y2*(1+xi)*(1+eta)*(1-zeta) + &
                    y3*(1-xi)*(1+eta)*(1-zeta) + y4*(1-xi)*(1-eta)*(1+zeta) + y5*(1+xi)*(1-eta)*(1+zeta) + & 
                    y6*(1+xi)*(1+eta)*(1+zeta) + y7*(1-xi)*(1+eta)*(1+zeta)) 

               zp = 0.125 * (z0*(1-xi)*(1-eta)*(1-zeta) + z1*(1+xi)*(1-eta)*(1-zeta) + z2*(1+xi)*(1+eta)*(1-zeta) + &
                    z3*(1-xi)*(1+eta)*(1-zeta) + z4*(1-xi)*(1-eta)*(1+zeta) + z5*(1+xi)*(1-eta)*(1+zeta) + & 
                    z6*(1+xi)*(1+eta)*(1+zeta) + z7*(1-xi)*(1+eta)*(1+zeta)) 

               Tdomain%GlobCoord (0,ipoint) = xp
               Tdomain%GlobCoord (1,ipoint) = yp
               Tdomain%GlobCoord (2,ipoint) = zp

               !!! Computation of the derivative matrix, dx_(jj)/dxi_(ii) !!!
               LocInvGrad(0,0) = 0.125 * ((x1-x0)*(1-eta)*(1-zeta) + (x2-x3)*(1+eta)*(1-zeta) + &
                                 (x5-x4)*(1-eta)*(1+zeta) + (x6-x7)*(1+eta)*(1+zeta))
               LocInvGrad(1,0) = 0.125 * ((x3-x0)*(1-xi)*(1-zeta) + (x2-x1)*(1+xi)*(1-zeta) + &
                                 (x7-x4)*(1-xi)*(1+zeta) + (x6-x5)*(1+xi)*(1+zeta))
               LocInvGrad(2,0) = 0.125 * ((x4-x0)*(1-xi)*(1-eta) + (x5-x1)*(1+xi)*(1-eta) + &
                                 (x7-x3)*(1-xi)*(1+eta) + (x6-x2)*(1+xi)*(1+eta))

               LocInvGrad(0,1) = 0.125 * ((y1-y0)*(1-eta)*(1-zeta) + (y2-y3)*(1+eta)*(1-zeta) + &
                                 (y5-y4)*(1-eta)*(1+zeta) + (y6-y7)*(1+eta)*(1+zeta))
               LocInvGrad(1,1) = 0.125 * ((y3-y0)*(1-xi)*(1-zeta) + (y2-y1)*(1+xi)*(1-zeta) + &
                                 (y7-y4)*(1-xi)*(1+zeta) + (y6-y5)*(1+xi)*(1+zeta))
               LocInvGrad(2,1) = 0.125 * ((y4-y0)*(1-xi)*(1-eta) + (y5-y1)*(1+xi)*(1-eta) + &
                                 (y7-y3)*(1-xi)*(1+eta) + (y6-y2)*(1+xi)*(1+eta))

               LocInvGrad(0,2) = 0.125 * ((z1-z0)*(1-eta)*(1-zeta) + (z2-z3)*(1+eta)*(1-zeta) + &
                                 (z5-z4)*(1-eta)*(1+zeta) + (z6-z7)*(1+eta)*(1+zeta))
               LocInvGrad(1,2) = 0.125 * ((z3-z0)*(1-xi)*(1-zeta) + (z2-z1)*(1+xi)*(1-zeta) + &
                                 (z7-z4)*(1-xi)*(1+zeta) + (z6-z5)*(1+xi)*(1+zeta))
               LocInvGrad(2,2) = 0.125 * ((z4-z0)*(1-xi)*(1-eta) + (z5-z1)*(1+xi)*(1-eta) + &
                                 (z7-z3)*(1-xi)*(1+eta) + (z6-z2)*(1+xi)*(1+eta))

               call invert_3d (LocInvGrad, Jac)
               Tdomain%specel(n)%Jacob(i,j,k) = Jac
               Tdomain%specel(n)%InvGrad (i,j,k, 0:2,0:2) = LocInvGrad(0:2,0:2)
           enddo
       enddo
   enddo

enddo

return
end subroutine shape8
