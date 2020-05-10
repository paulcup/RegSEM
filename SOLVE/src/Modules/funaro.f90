!====================
  MODULE FUNARO
!====================
!
 public :: def_xgll, def_wgll, def_heta, def_deriv,def_derivx 
 private
!
 CONTAINS
!

!
!----
!


        SUBROUTINE INLEGL (N, ET, VN, QN, X, QX) 
  !*********************************************************************  
  !   COMPUTES THE VALUE AT A GIVEN POINT OF A POLYNOMIAL INDIVIDUATED    
  !   BY THE VALUES ATTAINED AT THE LEGENDRE GAUSS-LOBATTO NODES          
  !   N  = THE DEGREE OF THE POLYNOMIAL                                   
  !   ET = VECTOR OF THE NODES, ET(I), I=0,N                              
  !   VN = VALUES OF THE LEGENDRE POLYNOMIAL AT THE NODES, VN(I), I=0,N   
  !   QN = VALUES OF THE POLYNOMIAL AT THE NODES, QN(I), I=0,N            
  !   X  = THE POINT IN WHICH THE COMPUTATION IS PERFORMED                
  !   QX = VALUE OF THE POLYNOMIAL IN X                                   
  !*********************************************************************  
        IMPLICIT DOUBLEPRECISION (A - H, O - Z) 
        DIMENSION ET (0: * ), VN (0: * ), QN (0: * ) 
        IF (N.EQ.0) RETURN 
                                                                          
        EPS = 1.D-14 
        CALL VALEPO (N, X, Y, DY, D2Y) 
        DN = DFLOAT (N) 
        C = 1.D0 / (DN * (DN + 1.D0) ) 
        QX = QN (0) * C * DY * (X - 1.D0) / VN (0) 
        QX = QX + QN (N) * C * DY * (X + 1.D0) / VN (N) 
        IF (N.EQ.1) RETURN 
                                                                          
        DO 1 J = 1, N - 1 
           ED = X - ET (J) 
           IF (DABS (ED) .LT.EPS) THEN 
              QX = QN (J) 
              RETURN 
           ELSE 
              QX = QX + QN (J) * C * DY * (X * X - 1.D0) / (VN (J)        &
              * ED)                                                       
           ENDIF 
      1 END DO 
                                                                          
        RETURN 
        END SUBROUTINE INLEGL                         
!
!----
!

        SUBROUTINE INLEGL2 (N, ET, VN, QN, X, QXX) 
  !*********************************************************************  
  !   COMPUTES THE VALUE AT A GIVEN POINT OF A POLYNOMIAL INDIVIDUATED    
  !   BY THE VALUES ATTAINED AT THE LEGENDRE GAUSS-LOBATTO NODES          
  !   N  = THE DEGREE OF THE POLYNOMIAL                                   
  !   ET = VECTOR OF THE NODES, ET(I), I=0,N                              
  !   VN = VALUES OF THE LEGENDRE POLYNOMIAL AT THE NODES, VN(I), I=0,N   
  !   QN = VALUES OF THE POLYNOMIAL AT THE NODES, QN(I), I=0,N            
  !   X  = THE POINT IN WHICH THE COMPUTATION IS PERFORMED                
  !   QXX = FIRST DERIVATIVE VALUE OF THE POLYNOMIAL IN X                                   
  !*********************************************************************  
        IMPLICIT DOUBLEPRECISION (A - H, O - Z) 
        DIMENSION ET (0: * ), VN (0: * ), QN (0: * )
        DOUBLEPRECISION, DIMENSION(0: N)  ::  DQN 
        IF (N.EQ.0) RETURN 
                                                                          
        EPS = 1.D-14 
        CALL VALEPO (N, X, Y, DY, D2Y) 
        DN = DFLOAT (N) 
        C = 1.D0 / (DN * (DN + 1.D0) ) 
        QXX = QN (0) * C/ VN (0) *(DY+D2Y*(X-1.D0))
        QXX = QXX + QN (N) * C/VN(N)*(DY+D2Y*(X + 1.D0))
        IF (N.EQ.1) RETURN 
        call DELEGL (N, ET, VN, QN, DQN)                                                                          
        DO 1 J = 1, N - 1 
           ED = X - ET (J) 
           IF (DABS (ED) .LT.EPS) THEN 
              QXX = DQN (J) 
              RETURN 
           ELSE 
              QXX = QXX + QN (J) * C/VN(J)* (           &
                    (D2Y*(X*X-1.D0)+   DY*2.d0*X)/ED    &
                     -DY*(X*X-1.D0)/ED**2   )
 
           ENDIF 
      1 END DO 
                                                                          
        RETURN 
        END SUBROUTINE INLEGL2
!
!----
!

        SUBROUTINE VAJAPO (N, A, B, X, Y, DY, D2Y) 
  !***********************************************************            
  !   COMPUTES THE VALUE OF THE JACOBI POLYNOMIAL OF DEGREE N             
  !   AND ITS FIRST AND SECOND DERIVATIVES AT A GIVEN POINT               
  !   N  = DEGREE OF THE POLYNOMIAL                                       
  !   A  = PARAMETER >-1                                                  
  !   B  = PARAMETER >-1                                                  
  !   X  = POINT IN WHICH THE COMPUTATION IS PERFORMED                    
  !   Y  = VALUE OF THE POLYNOMIAL IN X                                   
  !   DY = VALUE OF THE FIRST DERIVATIVE IN X                             
  !   D2Y= VALUE OF THE SECOND DERIVATIVE IN X                            
  !***********************************************************            
        IMPLICIT DOUBLEPRECISION (A - H, O - Z) 
                                                                          
        Y = 1.D0 
        DY = 0.D0 
        D2Y = 0.D0 
        IF (N.EQ.0) RETURN 
                                                                          
        AB = A + B 
        Y = .5D0 * (AB + 2.D0) * X + .5D0 * (A - B) 
        DY = .5D0 * (AB + 2.D0) 
        D2Y = 0.D0 
        IF (N.EQ.1) RETURN 
                                                                          
        YP = 1.D0 
        DYP = 0.D0 
        D2YP = 0.D0 
        DO 1 I = 2, N 
           DI = DFLOAT (I) 
           C0 = 2.D0 * DI + AB 
           C1 = 2.D0 * DI * (DI + AB) * (C0 - 2.D0) 
           C2 = (C0 - 1.D0) * (C0 - 2.D0) * C0 
           C3 = (C0 - 1.D0) * (A - B) * AB 
           C4 = 2.D0 * (DI + A - 1.D0) * C0 * (DI + B - 1.D0) 
           YM = Y 
           Y = ( (C2 * X + C3) * Y - C4 * YP) / C1 
           YP = YM 
           DYM = DY 
           DY = ( (C2 * X + C3) * DY - C4 * DYP + C2 * YP) / C1 
           DYP = DYM 
           D2YM = D2Y 
           D2Y = ( (C2 * X + C3) * D2Y - C4 * D2YP + 2.D0 * C2 * DYP)     &
           / C1                                                           
           D2YP = D2YM 
      1 END DO 
                                                                          
        RETURN 
        END SUBROUTINE VAJAPO                         
!
!----
!
        SUBROUTINE VALEPO (N, X, Y, DY, D2Y) 
  !*************************************************************          
  !   COMPUTES THE VALUE OF THE LEGENDRE POLYNOMIAL OF DEGREE N           
  !   AND ITS FIRST AND SECOND DERIVATIVES AT A GIVEN POINT               
  !   N  = DEGREE OF THE POLYNOMIAL                                       
  !   X  = POINT IN WHICH THE COMPUTATION IS PERFORMED                    
  !   Y  = VALUE OF THE POLYNOMIAL IN X                                   
  !   DY = VALUE OF THE FIRST DERIVATIVE IN X                             
  !   D2Y= VALUE OF THE SECOND DERIVATIVE IN X                            
  !*************************************************************          
        IMPLICIT DOUBLEPRECISION (A - H, O - Z) 
                                                                          
        Y = 1.D0 
        DY = 0.D0 
        D2Y = 0.D0 
        IF (N.EQ.0) RETURN 
                                                                          
        Y = X 
        DY = 1.D0 
        D2Y = 0.D0 
        IF (N.EQ.1) RETURN 
                                                                          
        YP = 1.D0 
        DYP = 0.D0 
        D2YP = 0.D0 
        DO 1 I = 2, N 
           C1 = DFLOAT (I) 
           C2 = 2.D0 * C1 - 1.D0 
           C4 = C1 - 1.D0 
           YM = Y 
           Y = (C2 * X * Y - C4 * YP) / C1 
           YP = YM 
           DYM = DY 
           DY = (C2 * X * DY - C4 * DYP + C2 * YP) / C1 
           DYP = DYM 
           D2YM = D2Y 
           D2Y = (C2 * X * D2Y - C4 * D2YP + 2.D0 * C2 * DYP) / C1 
           D2YP = D2YM 
      1 END DO 
                                                                          
        RETURN 
        END SUBROUTINE VALEPO                         
!
!----
!
        SUBROUTINE WELEGL (N, ET, VN, WT) 
  !********************************************************************** 
  !   COMPUTES THE WEIGHTS RELATIVE TO THE LEGENDRE GAUSS-LOBATTO FORMULA 
  !   N  = ORDER OF THE FORMULA                                           
  !   ET = JACOBI GAUSS-LOBATTO NODES, ET(I), I=0,N                       
  !   VN = VALUES OF THE LEGENDRE POLYNOMIAL AT THE NODES, VN(I), I=0,N   
  !   WT = VECTOR OF THE WEIGHTS, WT(I), I=0,N                            
  !********************************************************************** 
        IMPLICIT DOUBLEPRECISION (A - H, O - Z) 
        DIMENSION ET (0: * ), VN (0: * ), WT (0: * ) 
        IF (N.EQ.0) RETURN 
                                                                          
        N2 = (N - 1) / 2 
        DN = DFLOAT (N) 
        C = 2.D0 / (DN * (DN + 1.D0) ) 
        DO 1 I = 0, N2 
           X = ET (I) 
           Y = VN (I) 
           WTX = C / (Y * Y) 
           WT (I) = WTX 
           WT (N - I) = WTX 
      1 END DO 
                                                                          
        IF (N - 1.EQ.2 * N2) RETURN 
        X = 0.D0 
        Y = VN (N2 + 1) 
        WT (N2 + 1) = C / (Y * Y) 
                                                                          
        RETURN 
        END SUBROUTINE WELEGL                         
!
!----
!
        SUBROUTINE ZEJAGA (N, A, B, CS, DZ) 
  !**************************************************************         
  !   COMPUTES THE ZEROES OF THE JACOBI POLYNOMIAL OF DEGREE N            
  !   N  = THE NUMBER OF ZEROES                                           
  !   A  = PARAMETER BETWEEN -1/2 AND 1/2                                 
  !   B  = PARAMETER BETWEEN -1/2 AND 1/2                                 
  !   CS = VECTOR OF THE ZEROES, CS(I), I=1,N                             
  !   DZ = VECTOR OF THE DERIVATIVES AT THE ZEROES, DZ(I), I=1,N          
  !**************************************************************         
        IMPLICIT DOUBLEPRECISION (A - H, O - Z) 
        DIMENSION CS (1), DZ (1) 
        IF (N.EQ.0) RETURN 
                                                                          
        AB = A + B 
        CS (1) = (B - A) / (AB + 2.D0) 
        DZ (1) = .5D0 * AB + 1.D0 
        IF (N.EQ.1) RETURN 
                                                                          
        EPS = 1.E-17 
        PH = 1.57079632679489661923D0 
        C = PH / (2.D0 * DFLOAT (N) + AB + 1.D0) 
        DO 1 I = 1, N 
           DI = DFLOAT (I) 
           CSX = - DCOS (C * (4.D0 * DI + AB - 1.D0) ) 
           DO 2 IT = 1, 8 
              CALL VAJAPO (N, A, B, CSX, Y, DY, D2Y) 
              IF (DABS (Y) .LT.EPS) GOTO 3 
              CSX = CSX - Y / DY 
      2    END DO 
      3    IF (DABS (CSX) .LT.EPS) CSX = 0.D0 
           CS (I) = CSX 
           DZ (I) = DY 
      1 END DO 
                                                                          
        RETURN 
        END SUBROUTINE ZEJAGA                         
!
!----
!
        SUBROUTINE ZELEGL (N, ET, VN) 
  !********************************************************************   
  !   COMPUTES THE NODES RELATIVE TO THE LEGENDRE GAUSS-LOBATTO FORMULA   
  !   N  = ORDER OF THE FORMULA                                           
  !   ET = VECTOR OF THE NODES, ET(I), I=0,N                              
  !   VN = VALUES OF THE LEGENDRE POLYNOMIAL AT THE NODES, VN(I), I=0,N   
  !********************************************************************   
        IMPLICIT DOUBLEPRECISION (A - H, O - Z) 
        DIMENSION ET (0: * ), VN (0: * ) 
        IF (N.EQ.0) RETURN 
                                                                          
        N2 = (N - 1) / 2 
        SN = DFLOAT (2 * N - 4 * N2 - 3) 
        ET (0) = - 1.D0 
        ET (N) = 1.D0 
        VN (0) = SN 
        VN (N) = 1.D0 
        IF (N.EQ.1) RETURN 
                                                                          
        ET (N2 + 1) = 0.D0 
        X = 0.D0 
        CALL VALEPO (N, X, Y, DY, D2Y) 
        VN (N2 + 1) = Y 
        IF (N.EQ.2) RETURN 
                                                                          
        PI = 3.14159265358979323846D0 
        C = PI / DFLOAT (N) 
        DO 1 I = 1, N2 
           ETX = DCOS (C * DFLOAT (I) ) 
           DO 2 IT = 1, 8 
              CALL VALEPO (N, ETX, Y, DY, D2Y) 
              ETX = ETX - DY / D2Y 
      2    END DO 
           ET (I) = - ETX 
           ET (N - I) = ETX 
           VN (I) = Y * SN 
           VN (N - I) = Y 
      1 END DO 
                                                                          
        RETURN 
        END SUBROUTINE ZELEGL                         
!
!----
!
        SUBROUTINE DELEGL (N, ET, VN, QN, DQN) 
  !***********************************************************************
  !  COMPUTES THE DERIVATIVE OF A POLYNOMIAL AT THE LEGENDRE GAUSS-LOBATTO
  !  NODES FROM THE VALUES OF THE POLYNOMIAL ATTAINED AT THE SAME POINTS  
  !   N   = THE DEGREE OF THE POLYNOMIAL                                  
  !   ET  = VECTOR OF THE NODES, ET(I), I=0,N                             
  !   VN  = VALUES OF THE LEGENDRE POLYNOMIAL AT THE NODES, VN(I), I=0,N  
  !   QN  = VALUES OF THE POLYNOMIAL AT THE NODES, QN(I), I=0,N           
  !   DQN = DERIVATIVES OF THE POLYNOMIAL AT THE NODES, DQZ(I), I=0,N     
  !***********************************************************************
        IMPLICIT DOUBLEPRECISION (A - H, O - Z) 
        DIMENSION ET (0: * ), VN (0: * ), QN (0: * ), DQN (0: * ) 
        DQN (0) = 0.D0 
        IF (N.EQ.0) RETURN 
                                                                          
        DO 1 I = 0, N 
           SU = 0.D0 
           VI = VN (I) 
           EI = ET (I) 
           DO 2 J = 0, N 
              IF (I.EQ.J) GOTO 2 
              VJ = VN (J) 
              EJ = ET (J) 
              SU = SU + QN (J) / (VJ * (EI - EJ) ) 
      2    END DO 
           DQN (I) = VI * SU 
      1 END DO 
                                                                          
        DN = DFLOAT (N) 
        C = .25D0 * DN * (DN + 1.D0) 
        DQN (0) = DQN (0) - C * QN (0) 
        DQN (N) = DQN (N) + C * QN (N) 
                                                                          
        RETURN 
        END SUBROUTINE DELEGL                         
!
!----
!
  subroutine def_xgll(xgll,n)
!----------------------------
	doubleprecision, intent(out),dimension(n) 	:: xgll
	doubleprecision, dimension(0:n-1) 		:: ET, VN
	integer						:: n
!
	call ZELEGL(n-1,ET,VN)
	do i = 1,n
		xgll(i) = ET(i-1)
	enddo
!
 end subroutine def_xgll
!
!----
!
  subroutine def_wgll(wgll,n)
!----------------------------
	doubleprecision, intent(out),dimension(n) 	:: wgll
	doubleprecision, dimension(0:n-1) 		:: ET, VN, WT
	integer						:: n
!
	call ZELEGL(n-1,ET,VN)
	call WELEGL(n-1,ET,VN,WT)
	do i = 1,n
		wgll(i) = WT(i-1)
	enddo
!
 end subroutine def_wgll
!
!----
!
  subroutine def_heta(eta,delta,n,heta)
!--------------------------------------
	doubleprecision, intent(in), dimension(n) 	:: eta, delta
	doubleprecision, intent(out), dimension(n) 	:: heta	
	doubleprecision, dimension(0:n-1) 		:: ET, VN, QN
	integer						:: n
!
	call ZELEGL( n-1,ET(:),VN(:) )
	QN(0:n-1) = delta(1:n)
	do i = 1,n
		call INLEGL( n-1,ET(:),VN(:),QN(:),eta(i),heta(i) )
	enddo
!
 end subroutine def_heta
!
!----
!
  subroutine def_deriv(deriv,n)
!------------------------------
	doubleprecision, intent(out), dimension(n,n)	:: deriv
	doubleprecision, dimension(0:n-1) 		:: ET, VN, QN, DQN
	integer						:: i,n
!
	call ZELEGL( n-1,ET(:),VN(:) )
	do i = 0,n-1
		QN(0:n-1) = 0.d0
		QN(i) = 1.d0
		call DELEGL (n-1, ET(0:n-1), VN(0:n-1), QN(0:n-1), DQN(0:n-1))
		deriv(i+1,1:n) = DQN(0:n-1)
	enddo
!
  end subroutine def_deriv

  subroutine def_derivx(x,h,deriv,n)
!------------------------------
	doubleprecision, intent(in)                     :: x
	doubleprecision, intent(out), dimension(n)	:: h,deriv
	doubleprecision, dimension(0:n-1) 		:: ET, VN,QN
        doubleprecision                                 :: Q,DQ
	integer						:: i,n
!
	call ZELEGL( n-1,ET(:),VN(:) )
	do i = 0,n-1
		QN(0:n-1) = 0.d0
		QN(i) = 1.d0
		call INLEGL  (n-1, ET(0:n-1), VN(0:n-1), QN(0:n-1),X,  Q)
		call INLEGL2 (n-1, ET(0:n-1), VN(0:n-1), QN(0:n-1),X, DQ)
                h    (i+1) =  Q
		deriv(i+1) = DQ
	enddo
!
  end subroutine def_derivx
!
!----
!

!========================
  END MODULE FUNARO
!========================
