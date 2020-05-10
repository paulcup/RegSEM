!=====================================================
SUBROUTINE MERGE_SNAPSHOT_FILES(NPROC,NSNAP,BASE_NAME)
!=====================================================
!
IMPLICIT NONE
!
INTEGER :: IPROC,NPROC,NSNAP,N
INTEGER :: I,J,K,L
INTEGER :: I1,I2,J1,J2,K1,K2,X1,X2,Y1,Y2,Z1,Z2
REAL, DIMENSION(:,:,:), ALLOCATABLE :: X,Y,Z,UX,UY,UZ
CHARACTER :: FILENAME*100
CHARACTER(LEN=*) :: BASE_NAME
LOGICAL EXISTS
!
!==========================
!  MERGE SNAPSHOTS FILES  !
!==========================
!
WRITE(FILENAME,"(a,I3.3,a,I3.3)")TRIM(BASE_NAME),NSNAP,"_",0
INQUIRE(FILE=TRIM(FILENAME),EXIST=EXISTS)
IF(EXISTS)THEN
   !
   ! READ ALL FILE HEADER AND DETERMINE GRID DIMENSIONS
   !
   DO IPROC = 0,NPROC-1
      WRITE(FILENAME,"(A,I3.3,A,I3.3)")TRIM(BASE_NAME),NSNAP,"_",IPROC
      OPEN(26,FILE=TRIM(FILENAME),STATUS='OLD',FORM='UNFORMATTED',ACCESS='STREAM')
      READ(26)N,I1,I2,J1,J2,K1,K2
      IF(IPROC==1)THEN
         X1=I1 ; X2=I2 ; Y1=J1 ; Y2=J2 ; Z1=K1 ; Z2=K2 
      ELSE
         X1 = MIN(X1,I1)
         X2 = MAX(X2,I2)
         Y1 = MIN(Y1,J1)
         Y2 = MAX(Y2,J2)
         Z1 = MIN(Z1,K1)
         Z2 = MAX(Z2,K2)
      ENDIF
      CLOSE(26)
   ENDDO
   !
   ALLOCATE(X(x1:x2,y1:y2,z1:z2))
   ALLOCATE(Y(x1:x2,y1:y2,z1:z2))
   ALLOCATE(Z(x1:x2,y1:y2,z1:z2))
   !
   ALLOCATE(UX(x1:x2,y1:y2,z1:z2))
   ALLOCATE(UY(x1:x2,y1:y2,z1:z2))
   ALLOCATE(UZ(x1:x2,y1:y2,z1:z2))
   !
   ! MERGE FILES, READ/WRITE DATA
   !
   DO IPROC = 0,NPROC-1
      WRITE(FILENAME,"(A,I3.3,A,I3.3)")TRIM(BASE_NAME),NSNAP,"_",IPROC
      OPEN(26,FILE=TRIM(FILENAME),STATUS='OLD',FORM='UNFORMATTED',ACCESS='STREAM')
      READ(26)N,I1,I2,J1,J2,K1,K2
      !      
      DO L = 1,N
         !
         READ(26)I        , J        , K
         READ(26)X(I,J,K) , Y(I,J,K) , Z(I,J,K)
         READ(26)UX(I,J,K), UY(I,J,K), UZ(I,J,K)
         !         
      ENDDO
      !      
      CLOSE(26, STATUS='DELETE')
      !
   ENDDO
   !
   ! METERS TO KILOMETERS
   ! 
   x = X/1000.D0
   y = Y/1000.D0
   z = Z/1000.D0
   !
   WRITE(FILENAME,'(A,I3.3,A)')TRIM(BASE_NAME),NSNAP,'.vts' 
   CALL WRITE_VTS_FILE(X1,X2,Y1,Y2,Z1,Z2,X,Y,Z,UX,UY,UZ,TRIM(FILENAME),27)
   !
   DEALLOCATE(X,Y,Z,UX,UY,UZ)
   !
ENDIF
!
RETURN
!
!==================================
END SUBROUTINE MERGE_SNAPSHOT_FILES
!==================================
