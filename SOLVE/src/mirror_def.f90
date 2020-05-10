subroutine mirror_definition (Tdomain)


use sdomains
use read_mirror

implicit none

type (Domain), intent (INOUT) :: Tdomain

integer :: n, i, i1, i2, j1, j2, k1, k2, x, y, z, x1, x2, y1, y2, z1, z2, indx, indy, indz, ngllx, nglly, ngllz


i1 = Tdomain%i1_mirror
i2 = Tdomain%n_elem_x-Tdomain%i2_mirror-1
j1 = Tdomain%j1_mirror
j2 = Tdomain%n_elem_y-Tdomain%j2_mirror-1
k1 = Tdomain%k1_mirror
k2 = Tdomain%n_elem_z-Tdomain%k2_mirror-1

do n = 0,Tdomain%n_elem-1

   indx = Tdomain%specel(n)%indx
   indy = Tdomain%specel(n)%indy
   indz = Tdomain%specel(n)%indz
   
   ngllx = Tdomain%specel(n)%ngllx
   nglly = Tdomain%specel(n)%nglly
   ngllz = Tdomain%specel(n)%ngllz
   
   allocate (Tdomain%specel(n)%win_mirror(0:ngllx-1, 0:nglly-1, 0:ngllz-1))

   Tdomain%specel(n)%win_mirror = 0

   if(indx>=i1 .and. indx<=i2 .and. indy>=j1 .and. indy<=j2 .and. indz>=k1 .and. indz<=k2)then

   x1 = 0
   x2 = ngllx-1
   y1 = 0
   y2 = nglly-1 
   z1 = 0
   z2 = ngllz-1

   if(indx==i1)x1=x2
   if(indx==i2)x2=x1
   if(indy==j1)y1=y2
   if(indy==j2)y2=y1
   if(indz==k1)z1=z2
   if(indz==k2)z2=z1

   do z = z1,z2
      do y = y1,y2
         do x = x1,x2
            Tdomain%specel(n)%win_mirror(x,y,z)=1
         enddo
      enddo
   enddo

   endif


enddo

do n = 0,Tdomain%n_elem-1

ngllx = Tdomain%specel(n)%ngllx
nglly = Tdomain%specel(n)%nglly
ngllz = Tdomain%specel(n)%ngllz

! For a mirror of arbitrary shape, feed the code with subroutine get_mirror (it must be placed in the /SOLVE/src/Modules/read_mirror.f90)

! do z = 0,ngllz-1
!    do y = 0,nglly-1
!       do x = 0,ngllx-1
!           i = Tdomain%specel(n)%Iglobnum(x,y,z)
!          call get_mirror(Tdomain%Globcoord(0,i),Tdomain%Globcoord(1,i),Tdomain%Globcoord(2,i),&
!                          Tdomain%specel(n)%win_mirror(x,y,z))
!       enddo
!    enddo
! enddo

enddo

do n = 0,Tdomain%n_elem-1

 ngllx = Tdomain%specel(n)%ngllx	
 nglly = Tdomain%specel(n)%nglly
 ngllz = Tdomain%specel(n)%ngllz

 i = 0

 do z = 0,ngllz-1
    do y = 0,nglly-1
       do x = 0,ngllx-1
          if(Tdomain%specel(n)%win_mirror(x,y,z)==1)then
             i = i+1
          endif
       enddo
    enddo
 enddo

 if(i==0)then
    Tdomain%specel(n)%mirror_position = 0
 elseif(i==ngllx*nglly*ngllz)then
    Tdomain%specel(n)%mirror_position = 2
 else
    Tdomain%specel(n)%mirror_position = 1
 endif

enddo

i = 0

do n = 0,Tdomain%n_elem-1
   
   if(Tdomain%specel(n)%mirror_position == 1)then

      ngllx = Tdomain%specel(n)%ngllx	
      nglly = Tdomain%specel(n)%nglly
      ngllz = Tdomain%specel(n)%ngllz
      
      do z = 0,ngllz-1
         do y = 0,nglly-1
            do x = 0,ngllx-1
               if(Tdomain%specel(n)%win_mirror(x,y,z)==1) i = i+1
            enddo
         enddo
      enddo

   endif

enddo

Tdomain%recl_mirror = 3*i


return
end subroutine mirror_definition

