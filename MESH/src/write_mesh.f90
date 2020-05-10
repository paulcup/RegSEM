subroutine write_mesh (nb_pts,coord_pts,numface)


implicit none

integer, intent(IN) :: nb_pts, numface
doubleprecision, dimension(0:nb_pts-1,0:2), intent(INOUT) :: coord_pts


coord_pts(:,:) = coord_pts(:,:)/1000.d0

if (nb_pts==8) then

   if (numface==0) then
      write (30,*) coord_pts(0,:)
      write (30,*) coord_pts(1,:)
      write (30,*) coord_pts(2,:)
      write (30,*) coord_pts(3,:)
      write (30,'(a1)') ">"
   else if (numface==1) then
      write (31,*) coord_pts(0,:)
      write (31,*) coord_pts(1,:)
      write (31,*) coord_pts(5,:)
      write (31,*) coord_pts(4,:)
      write (31,'(a1)') ">"
   else if (numface==2) then
      write (32,*) coord_pts(1,:)
      write (32,*) coord_pts(2,:)
      write (32,*) coord_pts(6,:)
      write (32,*) coord_pts(5,:)
      write (32,'(a1)') ">"
   else if (numface==3) then
      write (33,*) coord_pts(2,:)
      write (33,*) coord_pts(3,:)
      write (33,*) coord_pts(7,:)
      write (33,*) coord_pts(6,:)
      write (33,'(a1)') ">"
   else if (numface==4) then
      write (34,*) coord_pts(0,:)
      write (34,*) coord_pts(3,:)
      write (34,*) coord_pts(7,:)
      write (34,*) coord_pts(4,:)
      write (34,'(a1)') ">"
   else if (numface==5) then
      write (35,*) coord_pts(4,:)
      write (35,*) coord_pts(5,:)
      write (35,*) coord_pts(6,:)
      write (35,*) coord_pts(7,:)
      write (35,'(a1)') ">"
   endif

else

   if (numface==0) then
      write (30,*) coord_pts(0,:)
      write (30,*) coord_pts(8,:)
      write (30,*) coord_pts(1,:)
      write (30,*) coord_pts(9,:)
      write (30,*) coord_pts(2,:)
      write (30,*) coord_pts(10,:)
      write (30,*) coord_pts(3,:)
      write (30,*) coord_pts(11,:)
      write (30,'(a1)') ">"
   else if (numface==1) then
      write (31,*) coord_pts(0,:)
      write (31,*) coord_pts(8,:)
      write (31,*) coord_pts(1,:)
      write (31,*) coord_pts(13,:)
      write (31,*) coord_pts(5,:)
      write (31,*) coord_pts(16,:)
      write (31,*) coord_pts(4,:)
      write (31,*) coord_pts(12,:)
      write (31,'(a1)') ">"
   else if (numface==2) then
      write (32,*) coord_pts(1,:)
      write (32,*) coord_pts(9,:)
      write (32,*) coord_pts(2,:)
      write (32,*) coord_pts(14,:)
      write (32,*) coord_pts(6,:)
      write (32,*) coord_pts(17,:)
      write (32,*) coord_pts(5,:)
      write (32,*) coord_pts(13,:)
      write (32,'(a1)') ">"
   else if (numface==3) then
      write (33,*) coord_pts(2,:)
      write (33,*) coord_pts(10,:)
      write (33,*) coord_pts(3,:)
      write (33,*) coord_pts(15,:)
      write (33,*) coord_pts(7,:)
      write (33,*) coord_pts(18,:)
      write (33,*) coord_pts(6,:)
      write (33,*) coord_pts(14,:)
      write (33,'(a1)') ">"
   else if (numface==4) then
      write (34,*) coord_pts(0,:)
      write (34,*) coord_pts(11,:)
      write (34,*) coord_pts(3,:)
      write (34,*) coord_pts(15,:)
      write (34,*) coord_pts(7,:)
      write (34,*) coord_pts(19,:)
      write (34,*) coord_pts(4,:)
      write (34,*) coord_pts(12,:)
      write (34,'(a1)') ">"
   else if (numface==5) then
      write (35,*) coord_pts(4,:)
      write (35,*) coord_pts(16,:)
      write (35,*) coord_pts(5,:)
      write (35,*) coord_pts(17,:)
      write (35,*) coord_pts(6,:)
      write (35,*) coord_pts(18,:)
      write (35,*) coord_pts(7,:)
      write (35,*) coord_pts(19,:)
      write (35,'(a1)') ">"
   endif

endif


end subroutine write_mesh
