!#  File name : realaxi-th.f90	!# Last modified : 10 Feb 2016
!#  Fortran program for Gross-Pitaevskii equation in three-dimensional 
!#  axially-symmetric trap by real time propagation (Fortran 90/95 Version)
!#
!# BEC-GP-OMP-FOR programs are developed by:
!#
!# Luis E. Young-S.
!# (Departamento de Ciencias Basicas, Universidad Santo Tomas, Tunja, Boyaca, Colombia)
!#
!# Paulsamy Muruganandam
!# (Department of Physics, Bharathidasan University, Tiruchirappalli, Tamil Nadu, India)
!#
!# Sadhan K. Adhikari
!# (Instituto de Fisica Teorica, UNESP - Sao Paulo State University, Brazil)
!#
!# Vladimir Loncar, Dusan Vudragovic, Antun Balaz
!# (Scientific Computing Laboratory, Center for the Study of Complex Systems,
!# Institute of Physics Belgrade, Serbia)
!#
!#
!# Public use and modification of these codes are allowed provided that the following
!# papers are cited:
!# [1] L. E. Young-S. et al., Comput. Phys. Commun. 220 (2017) 503.
!# [2] P. Muruganandam and S. K. Adhikari, Comput. Phys. Commun. 180 (2009) 1888.
!# [3] D. Vudragovic et al., Comput. Phys. Commun. 183 (2012) 2021.
!# [4] R. Kishor Kumar et al., Comput. Phys. Commun. 195 (2015) 117.
!# [5] B. Sataric et al., Comput. Phys. Commun. 200 (2016) 411.
!# [6] V. Loncar et al., Comput. Phys. Commun. 200 (2016) 406.
!# [7] L. E. Young-S. et al., Comput. Phys. Commun. 204 (2016) 209.
!# [8] V. Loncar et al., Comput. Phys. Commun. 209 (2016) 190.
!#
!# The authors would be grateful for all information and/or comments
!# regarding the use of the programs.
!
!# To compile :
!# (1) Intel Fortran Compiler (newer versions)
! ifort -O3 -qopenmp -w -mcmodel medium -shared-intel -V
!#  (note that some older versions require -openmp instead of -qopenmp)
!!
!# (2) GNU Fortran (gfortran)
! gfortran -O3  -fopenmp  -w  

MODULE COMM_DATA
! NX, NY : Number of space mesh points (X and Y)  
  INTEGER, PARAMETER :: NX = 130, NXX = NX-1 ! No. of X steps
  INTEGER, PARAMETER :: NY = 130, NYY = NY-1, NY2 = NY/2 ! No. of Y steps
! NSTP : Number of iterations to introduce the nonlinearity. 
! NSTP=0 reads the wave function, /=0 calculates the wave function.
! NPAS : Number of subsequent iterations with fixed nonlinearity.
! NRUN : Number of final time iterations with fixed nonlinearity. 
  INTEGER, PARAMETER :: NSTP = 100000, NPAS = 1000, NRUN = 20000
  INTEGER, PARAMETER :: NUMBER_OF_THREADS = 0  ! sets the number of CPU cores to be used
     !     NUMBER_OF_THREADS = 0 deactivates the command and uses all available CPU cores   
END MODULE COMM_DATA

MODULE GPE_DATA
  USE COMM_DATA, ONLY : NX, NY
  REAL (8), PARAMETER :: PI = 3.14159265358979D0
  REAL (8), PARAMETER :: AHO = 1.D-6     		 ! Unit of length (l = 1 MICRON)            
  REAL (8), PARAMETER :: Bohr_a0 =  5.2917720859D-11/AHO ! Bohr radius (scaled with AHO)
  COMPLEX (8), PARAMETER :: CI = (0.0D0, 1.0D0)		 ! Complex i
!  
  REAL (8), PARAMETER :: DX = 0.1D0, DY = 0.1D0, DT = 0.001D0	! DX, DY : Space step and DT : Time step
  INTEGER, PARAMETER :: NATOMS = 400			 ! Number of Atoms
  REAL (8), PARAMETER :: AS = 70.71602D0*Bohr_a0	 ! Scattering length (in units of Bohr_a0)
  REAL (8), PARAMETER :: ANU = 1.0D0, LAMBDA = 4.0D0	 ! ANU, LAMBDA : Parameteres of Trap  
!   
  REAL (8), PARAMETER :: G0 = 4.D0*PI*AS*NATOMS		 ! G0 : 3-dimensional nonlinearity  
  REAL (8), PARAMETER :: GPAR = 0.5D0			 ! Change for dynamics
! 
! OPTION and XOP decides which equation to be solved. 
! OPTION = 1 Solves -psi_xx-psi_yy+V(x,y)psi+G0|psi|^2 psi=i psi_t
! OPTION = 2 Solves [-psi_xx-psi_yy+V(x,y)psi]/2+G0|psi|^2 psi=i psi_t 
  INTEGER, PARAMETER :: OPTION = 2
  REAL (8) :: G, XOP, TPI
! X(0:N), Y(0:N) : Space mesh, V(0:N) : Potential, and 
! CP(0:N) : Wave function (COMPLEX)
 REAL (8), DIMENSION(:), ALLOCATABLE :: X, X2, Y, Y2
  REAL (8), DIMENSION(:, :), ALLOCATABLE :: V ,R2,P,P2
  COMPLEX (8), DIMENSION(:, :), ALLOCATABLE :: CP
!
 ! REAL (8), DIMENSION(0:NX, 0:NY) :: V
 ! REAL (8), DIMENSION(0:NX) :: X, X2
 ! REAL (8), DIMENSION(0:NY) :: Y, Y2
 ! COMPLEX (8), DIMENSION(0:NX, 0:NY) :: CP
END MODULE GPE_DATA

MODULE CN_DATA
  USE COMM_DATA, ONLY : NX, NY
  COMPLEX (8), DIMENSION(:), ALLOCATABLE :: CALA, CGAA, CALB, CGAB, CBE
  COMPLEX (8), DIMENSION(0:NX) :: CALX, CGAX, CAIPX
  COMPLEX (8), DIMENSION(1:NX) :: CPX
  COMPLEX (8), DIMENSION(0:NY) :: CALY, CGAY
  COMPLEX (8) :: CAIPMX, CAIPMY
END MODULE CN_DATA 

PROGRAM GROSS_PITAEVSKII_SSCN_AXI
  USE OMP_LIB
  USE COMM_DATA
  USE GPE_DATA
  IMPLICIT NONE
! Subroutine INITIALIZE() used to initialize the space mesh X(I), 
! potential V(I) and the initial wave function. Subroutine COEF() used to
! generate the coefficients for the Crank-Nicholson Scheme. The routine
! CALCNU() performs time progation for the non-derivative part and LU() performs
! time propagation of derivative part. NORM() calculates the norm and 
! normalizes the wave function, CHEM() and RAD() are used to calculate the 
! chemical potential, energy and the rms radius, respectively. The function 
! DIFF() used to calculate the space derivatives of the wave function used 
! in CHEM() and SIMP() does the integration by Simpson's rule.
!------------------------ INTERFACE BLOCKS -----------------------
  INTERFACE 
    SUBROUTINE ALLOCATE_VARIABLES()
    END SUBROUTINE ALLOCATE_VARIABLES
  END INTERFACE
!
  INTERFACE 
    SUBROUTINE FREE_VARIABLES()
    END SUBROUTINE FREE_VARIABLES
  END INTERFACE
 
 INTERFACE 
    SUBROUTINE INITIALIZE()
      IMPLICIT NONE
    END SUBROUTINE INITIALIZE
  END INTERFACE
!
  INTERFACE 
    SUBROUTINE COEF()
      IMPLICIT NONE
    END SUBROUTINE COEF
  END INTERFACE
!
  INTERFACE 
    SUBROUTINE CALCNU()
      IMPLICIT NONE 
    END SUBROUTINE CALCNU
  END INTERFACE
!
  INTERFACE 
    SUBROUTINE LUX()
      IMPLICIT NONE 
    END SUBROUTINE LUX
  END INTERFACE
!
  INTERFACE 
    SUBROUTINE LUY()
      IMPLICIT NONE 
    END SUBROUTINE LUY
  END INTERFACE
!
  INTERFACE
    SUBROUTINE NORM(CP, ZNORM)
      IMPLICIT NONE
      COMPLEX (8), DIMENSION(0:,0:), INTENT(INOUT) :: CP
      REAL (8), INTENT(OUT) :: ZNORM
    END  SUBROUTINE NORM
  END INTERFACE
!
  INTERFACE
    SUBROUTINE RAD(CP, RHORMS, ZRMS)
      IMPLICIT NONE
      COMPLEX (8), DIMENSION(0:,0:), INTENT(INOUT) :: CP
      REAL (8), INTENT(OUT) :: RHORMS, ZRMS
    END  SUBROUTINE RAD
  END INTERFACE
!

  INTERFACE
    SUBROUTINE CHEM(CP, MU, EN)
      IMPLICIT NONE
      COMPLEX (8), DIMENSION(0:,0:), INTENT(IN) :: CP
      REAL (8), INTENT(OUT) :: MU, EN
    END  SUBROUTINE CHEM
  END INTERFACE
!-----------------------------------------------------------
  INTERFACE
    PURE FUNCTION SIMP(F, DX) RESULT(VALUE)
      IMPLICIT NONE
      REAL (8), DIMENSION(0:), INTENT(IN) :: F
      REAL (8), INTENT(IN) :: DX
      REAL (8) :: VALUE
    END FUNCTION SIMP
  END INTERFACE
!-----------------------------------------------------------
!------------------------ END INTERFACE BLOCKS -------------------
  INTEGER :: I, J, K,NO_OF_THREADS
  REAL (8) :: GSTP, ZNORM, RHORMS, ZRMS, MU, EN, T, T1, T2,TMP
 REAL (8), DIMENSION(:,:), ALLOCATABLE :: CP2
 
 INTEGER :: CLCK_COUNTS_BEG, CLCK_COUNTS_END, CLCK_RATE

!  
!  
  CALL SYSTEM_CLOCK ( CLCK_COUNTS_BEG, CLCK_RATE )
  CALL CPU_TIME (T1)
! 

    IF(NUMBER_OF_THREADS.EQ.0) GO TO 9000
  
  CALL OMP_SET_NUM_THREADS(NUMBER_OF_THREADS)

  9000 CONTINUE

!  print *, omp_get_max_threads()
  !$OMP PARALLEL
    !$OMP MASTER
      NO_OF_THREADS = OMP_GET_NUM_THREADS()
    !$OMP END MASTER
  !$OMP END PARALLEL
!

 CALL ALLOCATE_VARIABLES()
  ALLOCATE(CP2(0:NX, 0:NY))

  SELECT CASE (OPTION)
    CASE (1)
      XOP = 1.0D0
    CASE (2)
      XOP = 2.0D0
    CASE (:0,3:)
      PRINT *, 'ERROR: Wrong option', OPTION
      STOP
  END SELECT
!
! INITIALIZE() initializes the starting normalized wave function 'CP'.
  CALL INITIALIZE()
! COEF() defines the coefficients of Crank-Nicholson Scheme.  
  CALL COEF()
!   
  OPEN(7, FILE = 'realaxi-th-out.txt')
  OPEN(4, FILE = 'realaxi-th-rms.txt')
 Write(7,*)  'Number of Threads =', NO_OF_THREADS
!     
  WRITE(7,900) OPTION 
  WRITE(4,900) OPTION
  WRITE(7,*) 
  WRITE(4,*)
  WRITE(7,901) NATOMS, AHO
  WRITE(7,902) AS/Bohr_a0 
  WRITE(7,903) G0  
  WRITE(7,904) ANU, LAMBDA
  WRITE(7,*)
  WRITE(7,905) NX, NY
  WRITE(7,906) DX, DY 
  WRITE(7,907) NSTP, NPAS, NRUN  
  WRITE(7,908) DT 
  WRITE(7,*) 
  WRITE(7, 1014) GPAR 
  WRITE(7,*)
  WRITE(7, 1001)
  WRITE(7, 1002)
  WRITE(7, 1001)  
  WRITE(4, 1001)
  WRITE(4, 1012)
  WRITE(4, 1001)  
  900 FORMAT(' Real time propagation axially-symmetric trap,   OPTION = ',I3 )
  901 FORMAT('  Number of Atoms N =',I10,', Unit of length AHO =',F12.8,' m')
  902 FORMAT('  Scattering length a = ',F11.6,'*a0 ')
  903 FORMAT('  Nonlinearity G_3D =',F11.6 )
  904 FORMAT('  Parameters of trap: ANU =',F7.2, ', LAMBDA =',F7.2)
  905 FORMAT(' # Space Stp: NX = ', I8, ', NY = ', I8)
  906 FORMAT('  Space Step: DX = ', F10.6, ', DY = ', F10.6)
  907 FORMAT(' # Time Stp : NSTP = ',I9,', NPAS = ',I9,', NRUN = ',I9)
  908 FORMAT('   Time Step:   DT = ', F10.6)
  1014 FORMAT(' * Change for dynamics: GPAR = ',F11.3, ' *')
  1001 FORMAT (19X,'-----------------------------------------------------')
  1002 FORMAT (20X, 'Norm', 6X, 'Chem', 8X, 'Ener/N', 6X, '<r>', 4X, '|Psi(0)|^2')
  1012 FORMAT ('Values of rms size:', 8X, '<r>', 13X, '<rho>', 13X, '<z>')
!   
!   OPEN(1, FILE = 'realaxi-th-den-init.txt')
!  	DO I = 0, NX, 10 ! Writes initial density
!     	 DO J = 0, NY, 10
!        	WRITE (1,999) X(I), Y(J), ABS(CP(I,J))**2
!     	 END DO
!     	 WRITE (1,*)
!	END DO
!   CLOSE(1)
! 
!   CP2 = ABS(CP) * ABS(CP)
! 
!   OPEN(11, FILE = 'realaxi-th-den-init1d_rho.txt')
!     DO I = 0, NX 
! 	 TMP = SIMP(CP2(I,:), DY)
!          WRITE(11, 99) X(I), TMP
!     END DO
!   CLOSE(11)
! 
!   OPEN(12, FILE = 'realaxi-th-den-init1d_z.txt') 
!     DO I = 0, NY 
! 	 TMP = 2.D0*PI* SIMP(CP2(:,I)*X(:), DX)
!          WRITE(12, 99) Y(I), TMP
!     END DO
!   CLOSE(12)
! 
  IF (NSTP /= 0) THEN  
    GSTP = XOP*G0/NSTP
    G = 0.D0
    CALL NORM(CP, ZNORM)	! NORM() calculates and restore normalization
    CALL CHEM(CP, MU, EN)	! CHEM() calculates the chemical potential MU and energy EN.
    CALL RAD(CP, RHORMS, ZRMS)	! RAD() calculates the r.m.s radius RMS
    WRITE (7, 1003) ZNORM, MU/XOP, EN/XOP, SQRT(RHORMS*RHORMS + ZRMS*ZRMS), ABS(CP(0,NY2+1))**2
    WRITE (4, 1013) SQRT(RHORMS*RHORMS + ZRMS*ZRMS), RHORMS, ZRMS
! 
    DO K = 1, NSTP		! NSTP iterations to introduce the nonlinearity
      G = G + GSTP
      CALL CALCNU()	! CALCNU() performs time propagation with non-derivative parts.
      CALL LUX()		! LU() performs the time iteration with space derivative alone.
      CALL LUY()
    END DO
    CALL NORM(CP, ZNORM)
    CALL CHEM(CP, MU, EN)
    CALL RAD(CP, RHORMS, ZRMS)
    WRITE (7, 1004) ZNORM, MU/XOP, EN/XOP, SQRT(RHORMS*RHORMS + ZRMS*ZRMS), ABS(CP(0,NY2+1))**2
    WRITE (4, 1017) SQRT(RHORMS*RHORMS + ZRMS*ZRMS), RHORMS, ZRMS    
    1004 FORMAT('After NSTP iter.:', F8.4, 2F12.4, 2F11.3)
    1017 FORMAT (2X, 'After NSTP iter.:', F14.5, 2F16.5)
!     
!   OPEN(2, FILE = 'realaxi-th-den-nstp.txt')
!       DO I = 0, NX, 10 ! Writes intermediate density
!     	 DO J = 0, NY, 10
!         WRITE (2,999) X(I), Y(J), ABS(CP(I,J))**2
!     	 END DO
!     	 WRITE (2,*)
!  	END DO
!   CLOSE(2)
! 
!   CP2 = ABS(CP) * ABS(CP)
! 
!   OPEN(21, FILE = 'realaxi-th-den-nstp1d_rho.txt')
!     DO I = 0, NX 
! 	 TMP = SIMP(CP2(I,:), DY)
!          WRITE(21, 99) X(I), TMP
!     END DO
!   CLOSE(21)
! 
!   OPEN(22, FILE = 'realaxi-th-den-nstp1d_z.txt') 
!     DO I = 0, NY 
! 	 TMP = 2.D0*PI* SIMP(CP2(:,I)*X(:), DX)
!          WRITE(22, 99) Y(I), TMP
!     END DO
!   CLOSE(22)
! 
  ELSE
    G = XOP*G0
    CALL NORM(CP, ZNORM)
    CALL CHEM(CP, MU, EN)
    CALL RAD(CP, RHORMS, ZRMS)    
    WRITE (7, 1003) ZNORM, MU/XOP, EN/XOP, SQRT(RHORMS*RHORMS + ZRMS*ZRMS), ABS(CP(0,NY2+1))**2
    WRITE (4, 1013) SQRT(RHORMS*RHORMS + ZRMS*ZRMS), RHORMS, ZRMS    
  END IF
! 
  T = 0.D0  
  IF((NPAS /= 0).OR.(NRUN /= 0)) OPEN(8, FILE = 'realaxi-th-dyna.txt')
!   
  DO K = 1, NPAS ! NPAS Iterations transient
     T = T + DT
     CALL CALCNU()
     CALL LUX()
     CALL LUY()
     IF (MOD(K,100).EQ.0) THEN          
        CALL RAD(CP, RHORMS, ZRMS)
        WRITE(8,1999) T*XOP, RHORMS, ZRMS
     END IF      
  END DO
  IF(NPAS /= 0)THEN   
    CALL NORM(CP, ZNORM)
    CALL CHEM(CP, MU, EN)
    CALL RAD(CP, RHORMS, ZRMS)
    WRITE (7, 1005) ZNORM, MU/XOP, EN/XOP, SQRT(RHORMS*RHORMS + ZRMS*ZRMS), ABS(CP(0,NY2+1))**2
    WRITE (4, 1015) SQRT(RHORMS*RHORMS + ZRMS*ZRMS), RHORMS, ZRMS     
    1005 FORMAT('After NPAS iter.:',F8.4, 2F12.4, 2F11.3)
    1015 FORMAT (2X, 'After NPAS iter.:', F14.5, 2F16.5)    
  END IF   
!   
  OPEN(3, FILE = 'realaxi-th-den.txt')
    DO I = 0, NX 
     DO J = 0, NY
        WRITE (3,999) X(I), Y(J), ABS(CP(I,J))**2
     END DO
     WRITE (3,*)
    END DO
  CLOSE(3)

  CP2 = ABS(CP) * ABS(CP)

  OPEN(31, FILE = 'realaxi-th-den1d_rho.txt')
    DO I = 0, NX 
	 TMP = SIMP(CP2(I,:), DY)
         WRITE(31, 99) X(I), TMP
    END DO
  CLOSE(31)

  OPEN(32, FILE = 'realaxi-th-den1d_z.txt') 
    DO I = 0, NY 
	 TMP = 2.D0*PI* SIMP(CP2(:,I)*X(:), DX)
         WRITE(32, 99) Y(I), TMP
    END DO
  CLOSE(32)
! 

!
!  The following line defines a problem which is studied in the time evolution.
  G =GPAR*G          
!   
  DO K = 1, NRUN  ! NRUN iterations to study nonlinear dynamics
     T = T + DT
     CALL CALCNU()
     CALL LUX()
     CALL LUY()
     IF (MOD(K,100).EQ.0) THEN          
        CALL RAD(CP, RHORMS, ZRMS)
        WRITE(8,1999) T*XOP, RHORMS, ZRMS
     END IF       
  END DO
  CLOSE(8)
  IF(NRUN /= 0)THEN 
    CALL NORM(CP, ZNORM)
    CALL CHEM(CP, MU, EN)
    CALL RAD(CP, RHORMS, ZRMS)
    WRITE (7, 1006) ZNORM, MU/XOP, EN/XOP, SQRT(RHORMS*RHORMS + ZRMS*ZRMS), ABS(CP(0,NY2+1))**2
    WRITE (4, 1016) SQRT(RHORMS*RHORMS + ZRMS*ZRMS), RHORMS, ZRMS     
    1006 FORMAT('After NRUN iter.:',F8.4, 2F12.4, 2F11.3)
    1016 FORMAT (2X, 'After NRUN iter.:', F14.5, 2F16.5) 
!     
!   OPEN(40, FILE = 'realaxi-th-den-nrun.txt')
!       DO I = 0, NX, 10 ! Writes intermediate density
!     	 DO J = 0, NY, 10
!         WRITE (40,999) X(I), Y(J), ABS(CP(I,J))**2
!     	 END DO
!     	 WRITE (40,*)
!  	END DO
!   CLOSE(40)
! 
!   CP2 = ABS(CP) * ABS(CP)
! 
!   OPEN(41, FILE = 'realaxi-th-den-nrun1d_rho.txt')
!     DO I = 0, NX 
! 	 TMP = SIMP(CP2(I,:), DY)
!          WRITE(41, 99) X(I), TMP
!     END DO
!   CLOSE(41)
! 
!   OPEN(42, FILE = 'realaxi-th-den-nrun1d_z.txt') 
!     DO I = 0, NY 
! 	 TMP = 2.D0*PI* SIMP(CP2(:,I)*X(:), DX)
!          WRITE(42, 99) Y(I), TMP
!     END DO
!   CLOSE(42)
! 
  END IF    

  CALL NORM(CP, ZNORM)

  99 FORMAT(F10.3, E16.8)

 CALL FREE_VARIABLES()
  IF (ALLOCATED(CP2)) DEALLOCATE(CP2)

! 
  CALL SYSTEM_CLOCK (CLCK_COUNTS_END, CLCK_RATE)
  CALL CPU_TIME(T2)
  WRITE (7, 1001) 
  WRITE (4, 1001) 
  WRITE (7,*)
  WRITE (7,'(A,I7,A)') ' Clock Time: ', (CLCK_COUNTS_END - CLCK_COUNTS_BEG)/INT (CLCK_RATE,8), ' seconds'
  WRITE (7,'(A,I7,A)') '   CPU Time: ', INT(T2-T1), ' seconds' 
!
  CLOSE (7)
  CLOSE (4)
  999 FORMAT (2F12.3, F16.8)
  1003 FORMAT ('Initial : ', 4X, F11.4, 2F12.4, 2F11.3)
  1013 FORMAT (11X, 'Initial:', F14.5, 2F16.5)
  1999 FORMAT(F12.6, 2F16.8)
END PROGRAM GROSS_PITAEVSKII_SSCN_AXI

SUBROUTINE ALLOCATE_VARIABLES()
  USE COMM_DATA, ONLY : NX, NY
  USE GPE_DATA, ONLY : X, Y, X2, Y2, V, CP, R2, P, P2
  USE CN_DATA, ONLY : CALA, CGAA, CALB, CGAB, CBE
 
  IMPLICIT NONE
!
  
 
  ALLOCATE(CBE(0:MAXVAL([NX,NY])-1))
!
  ALLOCATE(X(0:NX), X2(0:NX))   
  ALLOCATE(Y(0:NY), Y2(0:NY))   
  ALLOCATE(V(0:NX, 0:NY))
  ALLOCATE(CP(0:NX, 0:NY))
  ALLOCATE(R2(0:NX, 0:NY))
  ALLOCATE(P(0:NX, 0:NY))
  ALLOCATE(P2(0:NX, 0:NY))
  ALLOCATE(CALA(0:NX))   
  ALLOCATE(CGAA(0:NX))
  ALLOCATE(CALB(0:NY))   
  ALLOCATE(CGAB(0:NY))
!
END SUBROUTINE ALLOCATE_VARIABLES

SUBROUTINE FREE_VARIABLES()
 
  USE GPE_DATA, ONLY : X, Y, X2, Y2, V, CP, R2, P, P2
  USE CN_DATA, ONLY : CALA, CGAA, CALB, CGAB, CBE
  IMPLICIT NONE
  
  
  IF (ALLOCATED(CBE)) DEALLOCATE(CBE)
  IF (ALLOCATED(X)) DEALLOCATE(X)   
  IF (ALLOCATED(X2)) DEALLOCATE(X2)   
  IF (ALLOCATED(Y)) DEALLOCATE(Y)   
  IF (ALLOCATED(Y2)) DEALLOCATE(Y2)   
  IF (ALLOCATED(V)) DEALLOCATE(V)
  IF (ALLOCATED(CP)) DEALLOCATE(CP)
  IF (ALLOCATED(R2)) DEALLOCATE(R2)
  IF (ALLOCATED(P)) DEALLOCATE(P)
  IF (ALLOCATED(P2)) DEALLOCATE(P2)
  IF (ALLOCATED(CALA)) DEALLOCATE(CALA)
  IF (ALLOCATED(CGAA)) DEALLOCATE(CGAA)
  IF (ALLOCATED(CALB)) DEALLOCATE(CALB)
  IF (ALLOCATED(CGAB)) DEALLOCATE(CGAB)
END SUBROUTINE FREE_VARIABLES
!



SUBROUTINE INITIALIZE()
!   Routine that initizlizes the constant and variables.
!   Calculates the potential term V and the initial wave function CP
  USE COMM_DATA, ONLY : NX, NY, NY2,NSTP
  USE GPE_DATA
  IMPLICIT NONE
  REAL (8), DIMENSION(0:NX,0:NY):: TMP2D
  REAL (8) :: PI3, FAC, ANU2, LAMBDA2,T1,T2,TMP
  INTEGER :: I,J
!   
  ANU2 = ANU * ANU
  LAMBDA2 = LAMBDA * LAMBDA
  PI3 = PI**3
  TPI = 2.0D0*PI
  FAC = SQRT(SQRT(LAMBDA*ANU2/(PI3)))
  FORALL (I=0:NX) X(I) = I*DX
  FORALL (I=0:NY) Y(I) = (I-NY2)*DY
  X2 = X*X
  Y2 = Y*Y
!   
  IF (NSTP == 0) THEN
    WRITE(*,'(a)') "Run the program using the input file to read. e.g.: ./real3d-th < imag3d-den.txt"
    DO I = 0, NX; DO J = 0, NY 
       READ(*,*) T1,T2,tmp  !READ OUTPUT IMAGTIME PROGRAM WITH SAME PARAMETERS 
       CP(I,J) = SQRT(tmp) 
    END DO; END DO 
  ELSE


  FORALL(I=0:NX) V(I,:) = XOP * (ANU2*X2(I) + LAMBDA2*Y2)/2.0D0
  FORALL(I=0:NX) TMP2D(I,:) = (ANU*X2(I) + LAMBDA*Y2)/2.0D0
  CP = FAC*EXP(-TMP2D) 

   END IF
END SUBROUTINE INITIALIZE
!
SUBROUTINE COEF()
!   Routine that initizlizes the constant and variables.
!   Calculates the potential term V and the initial wave function CP
  USE COMM_DATA, ONLY : NX, NY, NXX, NYY
  USE GPE_DATA, ONLY : DX, DY, DT, CI, X
  USE CN_DATA
  IMPLICIT NONE
  INTEGER :: J
  REAL (8) :: DX2, DY2
  COMPLEX (8) :: CDT, CAI0, CAIM
! 
  CDT = DT/CI
  DX2 = DX*DX   ! Generating the Coefficients for the
  DY2 = DY*DY   ! C-N method (to solve the spatial part)           
  CAIPMX = CDT/(2.0D0*DX2)
  CAI0 = 1.D0 - CDT/DX2
  CALX(1) = 1.0D0
  CGAX(1) = -1.0D0/(CAI0+CAIPMX-CDT/(4.D0*DX*X(1)))
!
  DO J = 1, NXX
     CPX(J) = CDT/(4.0D0*DX*X(J))
     CAIPX(J) = CAIPMX-CPX(J)
     CAIM = CAIPMX+CPX(J)
     CALX(J+1)=CGAX(J)*CAIM
     CGAX(J+1)=-1.0D0/(CAI0+(CAIPMX-CDT/(4.D0*DX*X(J+1)))*CALX(J+1))
  END DO
!
  CAIPMY = CDT/(2.0D0*DY2)
  CAI0 = 1.D0 - CDT/DY2
  CALY(NYY) = 0.0D0
  CGAY(NYY) = -1.0D0/CAI0
!
  DO J = NYY, 1, -1
     CALY(J-1) = CGAY(J)*CAIPMY
     CGAY(J-1) = -1.0D0/(CAI0+CAIPMY*CALY(J-1))
  END DO
END SUBROUTINE COEF

SUBROUTINE CALCNU() ! Exact solution
  ! Solves the partial differential equation with the potential and 
  ! the nonlinear term.
  USE COMM_DATA, ONLY : NX, NY, NXX, NYY
  USE GPE_DATA, ONLY : DT, CI, V, G, CP, P, P2
  IMPLICIT NONE
  INTEGER :: I, J 
  REAL (8), DIMENSION(0:NX, 0:NY) :: TMP
!
  !$OMP PARALLEL DO PRIVATE(I, J)
    DO J = 0, NY; DO I = 0, NX 
      P(I,J) = ABS(CP(I,J))
      P2(I,J) = P(I,J) * P(I,J)
      TMP(I,J) = DT * ( V(I,J) + G * P2(I,J) )
      CP(I,J) = CP(I,J) * EXP(-CI*TMP(I,J))
    END DO; END DO
  !$OMP END PARALLEL DO
END SUBROUTINE CALCNU

SUBROUTINE LUX()
! Solves the partial differential equation only with the X-space
! derivative term using the Crank-Nicholson method
 USE OMP_LIB
  USE COMM_DATA, ONLY : NX, NY, NXX
  USE GPE_DATA, ONLY : CP
  USE CN_DATA, ONLY : CALX, CGAX, CAIPX, CAIPMX, CPX
  IMPLICIT NONE 
  INTEGER :: I, J
  COMPLEX (8), DIMENSION (0:NX) :: CBE
  COMPLEX (8) :: CXX
! 
  !$OMP PARALLEL DO PRIVATE(J, CXX, CBE)
  DO J = 0, NY
     CBE(1) = 0.0D0
     DO I = 1, NXX
        CXX = CP(I,J)-CAIPMX*(CP(I+1,J)-2.0D0*CP(I,J)+CP(I-1,J)) &
                - CPX(I)*(CP(I+1,J)-CP(I-1,J))
        CBE(I+1) = CGAX(I)*(CAIPX(I)*CBE(I)-CXX)
     END DO
     CP(NX,J) = 0.0D0
     DO I = NX, 1, -1
        CP(I-1,J) = CALX(I)*CP(I,J) + CBE(I)
     END DO
  END DO
 !$OMP END PARALLEL DO
END SUBROUTINE LUX
 
 

SUBROUTINE LUY()
!  Solves the partial differential equation only with the Y-space
!  derivative term using the Crank-Nicholson method
 USE OMP_LIB
  USE COMM_DATA, ONLY : NX, NY, NYY
  USE GPE_DATA, ONLY : CP
  USE CN_DATA, ONLY : CALY, CGAY, CAIPMY
  IMPLICIT NONE 
  INTEGER :: I, J
  COMPLEX (8), DIMENSION (0:NY) :: CBE
  COMPLEX (8) :: CYY
! 
  !$OMP PARALLEL DO PRIVATE(J, CYY, CBE)
  DO I = 0, NX
     CBE(NYY) = CP(I, NY)
     DO J = NYY, 1, -1
        CYY = CP(I,J)-CAIPMY*(CP(I,J+1)-2.0D0*CP(I,J)+CP(I,J-1))
        CBE(J-1) = CGAY(J)*(CAIPMY*CBE(J)-CYY)
     END DO
     DO J = 0, NY-1
        CP(I,J+1) = CALY(J)*CP(I,J)+CBE(J)
     END DO
  END DO
 !$OMP END PARALLEL DO
END SUBROUTINE LUY

SUBROUTINE NORM(CP, ZNORM)
!  Calculates the normalization of the wave function and sets it to
!  unity.
  USE COMM_DATA, ONLY : NX, NY
  USE GPE_DATA, ONLY : DX, DY, X, TPI
  IMPLICIT NONE
  COMPLEX (8), DIMENSION(0:,0:), INTENT(INOUT) :: CP
  REAL (8), INTENT(OUT) :: ZNORM
!-----------------------------------------------------------
  INTERFACE
    PURE FUNCTION SIMP(F, DX) RESULT(VALUE)
      IMPLICIT NONE
      REAL (8), DIMENSION(0:), INTENT(IN) :: F
      REAL (8), INTENT(IN) :: DX
      REAL (8) :: VALUE
    END FUNCTION SIMP
  END INTERFACE
!-----------------------------------------------------------
  INTEGER :: I
  REAL (8), DIMENSION(0:NX,0:NY) :: P, TMP2D
  REAL (8), DIMENSION(0:NX) :: TMP1D
!   
  P = ABS(CP)
  TMP2D = P*P
  FORALL (I = 0:NX) TMP1D(I) = X(I)*SIMP(TMP2D(I,:), DY)
  ZNORM = SQRT(TPI*SIMP(TMP1D, DX))
  CP = CP/ZNORM
END SUBROUTINE NORM

SUBROUTINE RAD(CP, RHORMS, ZRMS)
 ! Calculates the root mean square size RMS
  USE COMM_DATA, ONLY : NX, NY
  USE GPE_DATA, ONLY : DX, DY, X, X2, Y2, TPI
  IMPLICIT NONE
  COMPLEX (8), DIMENSION(0:,0:), INTENT(INOUT) :: CP
  REAL (8), INTENT(OUT) :: RHORMS, ZRMS
!-----------------------------------------------------------
  INTERFACE
    PURE FUNCTION SIMP(F, DX) RESULT(VALUE)
      IMPLICIT NONE
      REAL (8), DIMENSION(0:), INTENT(IN) :: F
      REAL (8), INTENT(IN) :: DX
      REAL (8) :: VALUE
    END FUNCTION SIMP
  END INTERFACE
!-----------------------------------------------------------
  INTEGER :: I
  REAL (8), DIMENSION(0:NX,0:NY) :: P, TMP2D, TMPY1, TMPY2
  REAL (8), DIMENSION(0:NX) :: TMP1D
!   
  P = ABS(CP)
  TMP2D = P*P
  FORALL (I = 0:NX) TMPY1(I,:) = X2(I)*TMP2D(I,:)
  FORALL (I = 0:NY) TMPY2(:,I) = Y2(I)*TMP2D(:,I)
  FORALL (I = 0:NX) TMP1D(I) = X(I)*SIMP(TMPY1(I,:), DY)
  RHORMS = SQRT(TPI*SIMP(TMP1D, DX))
  FORALL (I = 0:NX) TMP1D(I) = X(I)*SIMP(TMPY2(I,:), DY)
  ZRMS = SQRT(TPI*SIMP(TMP1D, DX))
END SUBROUTINE RAD

SUBROUTINE CHEM(CP, MU, EN)
  ! Calculates the chemical potential MU and energy EN.  CP is the wave
  ! function, V is the potential and G is the nonlinearity.
  USE COMM_DATA, ONLY : NX, NY
  USE GPE_DATA, ONLY : DX, DY, X, V, G, TPI
  IMPLICIT NONE
  COMPLEX (8), DIMENSION(0:,0:), INTENT(IN) :: CP
  REAL (8), INTENT(OUT) :: MU, EN
!-----------------------------------------------------------
  INTERFACE
    PURE FUNCTION DIFF(P, DX) RESULT (DP)
      IMPLICIT NONE
      REAL (8), DIMENSION(0:), INTENT(IN) :: P
      REAL (8), INTENT(IN) :: DX
      REAL (8), DIMENSION(0:SIZE(P)-1) :: DP
    END FUNCTION DIFF
  END INTERFACE
!-----------------------------------------------------------
  INTERFACE
    PURE FUNCTION SIMP(F, DX) RESULT(VALUE)
      IMPLICIT NONE
      REAL (8), DIMENSION(0:), INTENT(IN) :: F
      REAL (8), INTENT(IN) :: DX
      REAL (8) :: VALUE
    END FUNCTION SIMP
  END INTERFACE
!-----------------------------------------------------------
  INTEGER :: I
  REAL (8), DIMENSION(0:NX,0:NY) :: DPX, DPY, P, P2, GP2, DP2, TMP2D
  REAL (8), DIMENSION(0:NX) :: TMP1D
!   
  P = ABS(CP)
!   FORALL (I=0:NY) DPX(:,I) = DIFF(P(:,I), DX)
!   FORALL (I=0:NX) DPY(I,:) = DIFF(P(I,:), DY)
  DO I = 0, NX  
    DPY(I,0:NY) = DIFF(P(I,0:NY), DY)
  END DO 
  DO I = 0,NY
    DPX(0:NX,I) = DIFF(P(0:NX,I), DX)
  END DO
! 
  P2 = P*P
  GP2 = G*P2
  DP2 = DPX*DPX + DPY*DPY
  TMP2D = (V + GP2)*P2 + DP2
  FORALL (I = 0:NX) TMP1D(I) = X(I)*SIMP(TMP2D(I,:), DY)
  MU = TPI*SIMP(TMP1D, DX)
  TMP2D = (V + GP2/2.0D0)*P2 + DP2
  FORALL (I = 0:NX) TMP1D(I) = X(I)*SIMP(TMP2D(I,:), DY)
  EN = TPI*SIMP(TMP1D, DX)
END SUBROUTINE CHEM

PURE FUNCTION SIMP(F, DX) RESULT (VALUE)
! Does the spatial integration with Simpson's rule.
! N refer to the number of integration points, DX space step, and
! F is the function to be integrated.
  IMPLICIT NONE
  REAL (8), DIMENSION(0:), INTENT(IN) :: F
  REAL (8), INTENT(IN) :: DX
  REAL (8) :: VALUE
  REAL (8) :: F1, F2
  INTEGER :: I, N
!   
  N = SIZE(F) - 1
  F1 = F(1) + F(N-1) ! N EVEN
  F2 = F(2) 
  DO I = 3, N-3, 2
     F1 = F1 + F(I)
     F2 = F2 + F(I+1)
  END DO
  VALUE = DX*(F(0) + 4.0D0*F1 + 2.0D0*F2 + F(N))/3.0D0
END FUNCTION SIMP

PURE FUNCTION DIFF(P,DX) RESULT (DP)
  ! Computes the first derivative DP of P using
  ! Richardsonextrapolation formula. The derivative at the  
  ! boundaries are assumed to be zero
  IMPLICIT NONE
  REAL (8), DIMENSION(0:), INTENT(IN) :: P
  REAL (8), INTENT(IN) :: DX
  REAL (8), DIMENSION(0:SIZE(P)-1) :: DP
  INTEGER :: I, N
!   
  N = SIZE(P) - 1
  DP(0) = 0.0D0
  DP(1) = (P(2) - P(0))/(2.0D0*DX)
  FORALL(I=2:N-2)
    DP(I) = (P(I-2)-8.0D0*P(I-1)+8.0D0*P(I+1)-P(I+2))/(12.0D0*DX)
  END FORALL
  DP(N-1) = (P(N) - P(N-2))/(2.0D0*DX)
  DP(N) = 0.0D0
END FUNCTION DIFF

!# File name : realaxi-th.f90
