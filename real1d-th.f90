!#  File name : real1d-th.f90	!# Last modified : 26 April 2017
!#  Fortran program for Gross-Pitaevskii equation in one-dimensional 
!#  harmonic trap by real time propagation (Fortran 90/95 Version)
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
! N : Number of space mesh points
  INTEGER, PARAMETER :: N = 3000, N2 = N/2, NX = N-1
! NSTP : Number of iterations to introduce the nonlinearity. 
! NSTP=0 reads the wave function, /=0 calculates the wave function.
! NPAS : Number of subsequent iterations with fixed nonlinearity.
! NRUN : Number of final time iterations with fixed nonlinearity.
! NUMBER_OF_THREADS : Number of processing units to run program.
  INTEGER, PARAMETER :: NSTP = 1000000, NPAS = 100, NRUN = 40000
  INTEGER, PARAMETER :: NUMBER_OF_THREADS = 0  ! sets the number of CPU cores to be used
     !     NUMBER_OF_THREADS =0 deactivates the command and uses all available CPU cores   
  REAL (8), PARAMETER :: PI = 3.14159265358979D0
END MODULE COMM_DATA

MODULE GPE_DATA
  USE COMM_DATA, ONLY : N, PI
  REAL (8), PARAMETER :: AHO = 1.D-6
  REAL (8), PARAMETER :: Bohr_a0 =  5.2917720859D-11/AHO
  COMPLEX (8), PARAMETER :: CI = (0.D0,1.D0)
  REAL (8), PARAMETER :: DX = 0.01D0, DT = 0.0001D0
  INTEGER, PARAMETER  :: NATOMS = 2000
  REAL (8), PARAMETER :: AS = 74.103248D0*Bohr_a0
  REAL (8), PARAMETER :: GAMMA = 1.D0
  REAL (8), PARAMETER :: DRHO = 0.5D0, DRHO2 = DRHO * DRHO
  REAL (8), PARAMETER :: G_1D = 4.D0*PI*AS*NATOMS/(2.D0*PI*DRHO2)
  REAL (8), PARAMETER :: G_3D = G_1D*2.D0*PI*DRHO2
  REAL (8), PARAMETER :: GPAR = 0.5D0
  INTEGER, PARAMETER :: OPTION = 2
  REAL (8), DIMENSION(:), ALLOCATABLE :: X, X2, V
  COMPLEX (8), DIMENSION(:), ALLOCATABLE :: CP
  REAL (8) :: XOP, G
  INTEGER :: NO_OF_THREADS 
END MODULE GPE_DATA

MODULE CN_DATA
  USE COMM_DATA, ONLY : N
  COMPLEX (8), DIMENSION(:), ALLOCATABLE :: CAL, CGA, CBE
  COMPLEX (8) :: CAIPM
END MODULE CN_DATA

PROGRAM GROSS_PITAEVSKII_SSCN_1D
  USE COMM_DATA
  USE GPE_DATA
 ! USE UTIL, ONLY : NORM, RAD
  USE OMP_LIB
  IMPLICIT NONE
! Subroutine INITIALIZE() used to initialize the space mesh X(I) and the initial wave function. 
! Subroutine CALCULATE_TRAP() used to initialize the harmonic oscillator potential V. 
! Subroutine COEF() used to generate the coefficients for the Crank-Nicholson Scheme. 
! The routine CALCNU() performs time progation for the non-derivative part and LU() performs
! time propagation of derivative part. NORM() calculates the norm and 
! normalizes the wave function, CHEM() and RAD() are used to calculate the 
! chemical potential, energy and the rms radius, respectively. The function 
! DIFF() used to calculate the space derivatives of the wave function used 
! in CHEM() and SIMP() does the integration by Simpson's rule.
!------------------------ Interface Blocks -----------------------
 INTERFACE 
    SUBROUTINE INITIALIZE()
      IMPLICIT NONE
    END SUBROUTINE INITIALIZE
  END INTERFACE
!------------------------
  INTERFACE 
    SUBROUTINE CALCULATE_TRAP()
    END SUBROUTINE CALCULATE_TRAP
  END INTERFACE
!------------------------
 INTERFACE 
    SUBROUTINE COEF()
      IMPLICIT NONE
    END SUBROUTINE COEF
  END INTERFACE
!------------------------
  INTERFACE 
    SUBROUTINE CALCNU(CP, DT)
      IMPLICIT NONE
      COMPLEX (8), DIMENSION(0:), INTENT(INOUT) :: CP
      REAL (8), INTENT(IN) :: DT
    END SUBROUTINE CALCNU
  END INTERFACE
!------------------------
  INTERFACE 
    SUBROUTINE LU(CP)
      IMPLICIT NONE

      COMPLEX (8), DIMENSION(0:), INTENT(INOUT) :: CP
    END SUBROUTINE LU
  END INTERFACE
!------------------------
  INTERFACE
    SUBROUTINE NORM(CP, ZNORM)
      IMPLICIT NONE
      COMPLEX (8), DIMENSION(0:), INTENT(INOUT) :: CP
      REAL (8), INTENT(OUT) :: ZNORM
    END  SUBROUTINE NORM
  END INTERFACE
!------------------------
  INTERFACE
    SUBROUTINE RAD(CP2, RMS)
      IMPLICIT NONE
      REAL (8), DIMENSION(0:), INTENT(IN) :: CP2
      REAL (8), INTENT(OUT) :: RMS
    END  SUBROUTINE RAD
  END INTERFACE
!------------------------
  INTERFACE
    SUBROUTINE CHEM(CP, MU, EN)
      IMPLICIT NONE
      COMPLEX (8), DIMENSION(0:), INTENT(IN) :: CP
      REAL (8), INTENT(OUT) :: MU, EN
    END  SUBROUTINE CHEM
  END INTERFACE 
!------------------------
  INTERFACE
    SUBROUTINE WRITE_DENSITY(FUNIT, U2)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: FUNIT
      REAL (8), DIMENSION(0:), INTENT(IN) :: U2
    END SUBROUTINE WRITE_DENSITY
  END INTERFACE
!------------------------ INTERFACE BLOCKS -----------------------
  INTEGER :: K
  REAL (8) :: ZNORM, MU, RMS, EN, T, T1, T2
  REAL (8), DIMENSION(:), ALLOCATABLE :: CP2
  REAL (8) :: GSTP  
  INTEGER :: CLCK_COUNTS_BEG, CLCK_COUNTS_END, CLCK_RATE
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

  SELECT CASE (OPTION)
    CASE (1)
      XOP = 1.D0
    CASE (2)
      XOP = 2.D0
    CASE (:0,3:)
      PRINT *, 'ERROR: Wrong OPTION', OPTION
      STOP
  END SELECT
!

!
  CALL ALLOCATE_VARIABLES()
  ALLOCATE(CP2(0:NX))
!  
! INITIALIZE() initializes the starting normalized wave function 'CP'.
  CALL INITIALIZE()
! CALCULATE_TRAP() initializes the harmonic oscillator potential V.
  CALL CALCULATE_TRAP()
! COEF() defines the coefficients of Crank-Nicholson Scheme.
  CALL COEF()
!
  OPEN(7, FILE = 'real1d-th-out.txt')
 
  Write(7,*)  'Number of Threads =', NO_OF_THREADS
! 
 WRITE(7,900) OPTION 

  WRITE(7,*) 
  WRITE(7,901) NATOMS, AHO
  WRITE(7,902) AS/Bohr_a0 
  WRITE(7,910) G_3D
  WRITE(7,903) G_1D
  WRITE(7,904) GAMMA
  WRITE(7,909) DRHO
  WRITE(7,*)
  WRITE(7,905) N
  WRITE(7,906) DX
  WRITE(7,907) NSTP, NPAS, NRUN  
  WRITE(7,908) DT 
  WRITE(7,*)
  WRITE(7, 1012) GPAR 
  WRITE(7,*)
  WRITE(7,1001)
  WRITE(7,1002)
  WRITE(7,1001)  
  900 FORMAT(' Real time propagation 1d,   OPTION = ',I3 )
  901 FORMAT('  Number of Atoms N =',I10,', Unit of length AHO =',F12.8,' m')
  902 FORMAT('  Scattering length a = ',F9.2,'*a0')
  910 FORMAT('  Nonlinearity G_3D =',F16.7)
  903 FORMAT('  Nonlinearity G_1D =',F16.7)
  904 FORMAT('  Parameter of trap: GAMMA =',F7.2)
  909 FORMAT('  Transverse trap parameter =',F7.3)
  905 FORMAT(' # Space Stp: NX = ', I8)
  906 FORMAT('  Space Step: DX = ', F10.6)
  907 FORMAT(' # Time Stp : NSTP = ',I9,', NPAS = ',I9,', NRUN = ',I9)
  908 FORMAT('   Time Step:   DT = ', F10.6)
  1012 FORMAT(' * Change for dynamics: GPAR = ',F11.3, ' *')
  1001 FORMAT (19X,'-----------------------------------------------------')
  1002 FORMAT (20X, 'Norm', 6X, 'Chem', 8X, 'Ener/N', 6X, '<x>', 6X, '|Psi(0)|^2')   
  1003 FORMAT ('Initial : ', 4X, F11.4, 2F12.4, 2F11.4)
!
!   OPEN(1, FILE = 'real1d-th-den-init.txt')
!    CALL WRITE_DENSITY(1, CP2)
!   CLOSE(1)
!
! NORM() calculates and restore normalization.
! CHEM() calculates the chemical potential MU and energy EN.
! RAD() calculates the r.m.s radius RMS.

  IF (NSTP /= 0) THEN
    GSTP =  XOP * G_1D / DFLOAT(NSTP)    
    G = 0.D0
    CALL NORM(CP,  ZNORM)
    CALL CHEM(CP, MU, EN)
    CP2 = ABS(CP)**2
    CALL RAD(CP2, RMS)
    WRITE (7, 1003) ZNORM, MU/XOP, EN/XOP, RMS, CP2(N2) 
    DO K = 1, NSTP ! NSTP iterations to introduce the nonlinearity
      G = G + GSTP
      CALL CALCNU(CP, DT)
      CALL LU(CP)
    END DO
    CALL NORM(CP,  ZNORM)
    CALL CHEM(CP, MU, EN)
    CP2 = ABS(CP)**2
    CALL RAD(CP2, RMS)
    WRITE (7, 1005) ZNORM, MU/XOP, EN/XOP, RMS, CP2(N2)
    1005 FORMAT('After NSTP iter.:',F8.4, 2F12.4, 2F11.4)
!
!    OPEN(2, FILE = 'real1d-th-den-nstp.txt')
!     CALL WRITE_DENSITY(2, CP2)
!    CLOSE(2)
!
  ELSE
    G = XOP * G_1D
    CALL NORM(CP,  ZNORM)
    CALL CHEM(CP, MU, EN)
    CP2 = ABS(CP)**2
    CALL RAD(CP2, RMS)
    WRITE (7, 1003) ZNORM, MU/XOP, EN/XOP, RMS, CP2(N2) 
  END IF
!
  T = 0.D0  
  IF((NPAS /= 0).OR.(NRUN /= 0)) OPEN(8, FILE = 'real1d-th-dyna.txt')
!
  DO K = 1, NPAS ! NPAS iterations transient
     T = T + DT
    CALL CALCNU(CP, DT)
    CALL LU(CP)
      IF (MOD(K, 100) == 0) THEN
        CP2 = ABS(CP)**2
        CALL RAD(CP2,  RMS)
        WRITE(8,1006) T*XOP, RMS
      END IF
  END DO
  IF(NPAS /= 0)THEN    
    CALL NORM(CP, ZNORM)
    CALL CHEM(CP, MU, EN)
    CP2 = ABS(CP)**2
    CALL RAD(CP2,  RMS)
    WRITE (7, 1008) ZNORM, MU/XOP, EN/XOP, RMS, CP2(N2)
    1008 FORMAT('After NPAS iter.:',F8.4, 2F12.4, 2F11.4)
!
    OPEN(3, FILE = 'real1d-th-den.txt')
     CALL WRITE_DENSITY(3, CP2)
    CLOSE(3)
!
  END IF
!
!  The following line defines a problem which is studied in the time evolution.
  G = GPAR*G
! 
  DO K = 1, NRUN ! NRUN iterations to check convergence
    T = T + DT
    CALL CALCNU(CP, DT)
    CALL LU(CP)
      IF (MOD(K, 50) == 0) THEN
        CP2 = ABS(CP)**2
        CALL RAD(CP2,  RMS)
        WRITE(8,1006) T*XOP, RMS
      END IF
  END DO
  CLOSE (8)
  IF(NRUN /= 0)THEN   
    CALL NORM(CP,  ZNORM) 
    CALL CHEM(CP, MU, EN)
    CP2 = ABS(CP)**2 
    CALL RAD(CP2,  RMS)
    WRITE (7, 1006) ZNORM, MU/XOP, EN/XOP, RMS, CP2(N2)
    1006 FORMAT('After NRUN iter.:',F8.4, 2F12.4, 2F11.4)
!
!    OPEN(4, FILE = 'real1d-th-den-nrun.txt')
!      CALL WRITE_DENSITY(4, CP2)
!    CLOSE(4)
!  
  END IF
!
  CALL SYSTEM_CLOCK (CLCK_COUNTS_END, CLCK_RATE)
  CALL CPU_TIME(T2)
  WRITE (7, 1001)  
  WRITE (7,*)
  WRITE (7,'(A,I7,A)') ' Clock Time:', (CLCK_COUNTS_END - CLCK_COUNTS_BEG)/INT (CLCK_RATE,8), ' seconds'
  WRITE (7,'(A,I7,A)') '   CPU Time:', INT(T2-T1), ' seconds' 
  CLOSE (7)
!
END PROGRAM GROSS_PITAEVSKII_SSCN_1D

!
SUBROUTINE ALLOCATE_VARIABLES()
  USE COMM_DATA, ONLY : N, NX
  USE GPE_DATA, ONLY : X, X2, V, CP
  USE CN_DATA, ONLY : CAL, CGA, CBE
!  USE UTIL, ONLY : TMPX
  IMPLICIT NONE
!
!  ALLOCATE(TMPX(0:N))
  ALLOCATE(CBE(0:N-1))
!
  ALLOCATE(X(0:N))   
  ALLOCATE(X2(0:N))   
  ALLOCATE(V(0:N))
  ALLOCATE(CP(0:NX))
  ALLOCATE(CAL(0:NX))   
  ALLOCATE(CGA(0:NX))
!
END SUBROUTINE ALLOCATE_VARIABLES

SUBROUTINE FREE_VARIABLES()
  USE GPE_DATA, ONLY : X, X2, V, CP
  USE CN_DATA, ONLY : CAL, CGA, CBE
!  USE UTIL, ONLY : TMPX
  IMPLICIT NONE
!  IF (ALLOCATED(TMPX)) DEALLOCATE(TMPX)
  IF (ALLOCATED(CBE)) DEALLOCATE(CBE)
!
  IF (ALLOCATED(X)) DEALLOCATE(X)   
  IF (ALLOCATED(X2)) DEALLOCATE(X2)   
  IF (ALLOCATED(V)) DEALLOCATE(V)
  IF (ALLOCATED(CP)) DEALLOCATE(CP)
  IF (ALLOCATED(CAL)) DEALLOCATE(CAL)
  IF (ALLOCATED(CGA)) DEALLOCATE(CGA)
!
END SUBROUTINE FREE_VARIABLES


SUBROUTINE INITIALIZE()
! Routine that initializes the constant and variables
! Calculates the initial wave function CP
  USE COMM_DATA, ONLY : N, NX, N2, PI,NSTP 
  USE GPE_DATA, ONLY : GAMMA, DX, X, X2, CP 
  IMPLICIT NONE
  REAL (8) :: PI4,TX,TMP
  INTEGER :: I
!
 FORALL (I=0:N) X(I) = (I-N2)*DX
  X2 = X*X
  PI4 = SQRT(SQRT(PI/GAMMA))
!         
  IF (NSTP == 0) THEN
    WRITE(*,'(a)') "Run the program using the input file to read. e.g.: ./real1d < imag1d-den.txt"
     DO I = 0,N
        READ (*,*)TX, TMP	  !READ OUTPUT of IMAGTIME PROGRAM WITH SAME PARAMETERS 
        CP(I) = SQRT(TMP)
     END DO
  ELSE 
   CP = EXP(-X2*GAMMA/2.D0)/PI4
  END IF
END SUBROUTINE INITIALIZE

SUBROUTINE CALCULATE_TRAP()
! Routine that  initializes the harmonic oscillator potential V.
  USE COMM_DATA, ONLY :  NX
  USE GPE_DATA, ONLY : XOP, X2, V, GAMMA
  IMPLICIT NONE
  REAL (8) :: GAMMA2
  INTEGER :: I
  GAMMA2 = GAMMA * GAMMA
  !$OMP PARALLEL DO PRIVATE(I)
  DO I = 0, NX
    V(I) = XOP * GAMMA2 * X2(I) / 2.0D0
  END DO
  !$OMP END PARALLEL DO
!
END SUBROUTINE CALCULATE_TRAP

SUBROUTINE COEF() 
! Calculates the coefficients needed in subroutine lu.
  USE COMM_DATA, ONLY :  NX
  USE GPE_DATA, ONLY : DX, DT, CI
  USE CN_DATA
  IMPLICIT NONE
  INTEGER :: J
  REAL (8) :: DX2
  COMPLEX (8) :: CAI0, CDT
!
  CDT = DT/CI
  DX2 = DX*DX
  CAIPM = CDT/2.D0/DX2
  CAI0 = 1.D0-CDT/DX2
  CAL(NX) = 0.D0
  CGA(NX) = -1.D0/CAI0
  DO J = NX, 1, -1
     CAL(J-1) = CGA(J)*CAIPM
     CGA(J-1) = -1.D0/(CAI0+CAIPM*CAL(J-1))    
  END DO
END SUBROUTINE COEF

SUBROUTINE CALCNU(CP, DT) ! Exact solution
! Solves the partial differential equation with the potential and the
! nonlinear term.
  USE COMM_DATA, ONLY : N
  USE GPE_DATA, ONLY : V, G, CI
  IMPLICIT NONE
!-------------------------------------------------
  COMPLEX (8), DIMENSION(0:), INTENT(INOUT) :: CP
  REAL (8), INTENT(IN) :: DT
  REAL (8) :: P2, TMP
  INTEGER :: I
!
  !$OMP PARALLEL DO PRIVATE(I, P2, TMP)
  DO I = 0, N
    P2 = ABS(CP(I)) * ABS(CP(I))
    TMP = DT * (V(I) + G * P2) 
    CP(I) = CP(I) * EXP(-CI*TMP)
  END DO
  !$OMP END PARALLEL DO
END SUBROUTINE CALCNU

SUBROUTINE LU(CP)
! Solves the partial differential equation only with the space
! derivative term using the Crank-Nicholson method
  USE COMM_DATA, ONLY : N, NX
  USE CN_DATA
  IMPLICIT NONE
  COMPLEX (8), DIMENSION(0:), INTENT(INOUT) :: CP
  COMPLEX (8),DIMENSION(0:N)  :: CXX
  INTEGER :: I
!
  CBE(NX) = CP(N)
   
  DO I = NX, 1, -1
     CXX(I) = CP(I)-(CP(I+1)-2.D0*CP(I)+CP(I-1))*CAIPM
     CBE(I-1) = CGA(I)*(CAIPM*CBE(I) - CXX(I))
  END DO
  
!-----------------------------
! Boundary condition periodic:
!  DO I = 0, NX
!     CP(I+1) = CAL(I)*CP(I) + CBE(I)
!  END DO
!  CP(0) = CP(N)
!  CP(1) = CP(NX)
!-----------------------------
! Boundary condition reflecting:
   CP(0) = 0.D0
 
   DO I = 0, NX-1
      CP(I+1) = CAL(I) * CP(I) + CBE(I)
   END DO
 
   CP(N) = 0.D0
!-----------------------------
END SUBROUTINE LU

SUBROUTINE CHEM(CP, MU, EN)
! Calculates the chemical potential MU and energy EN.  CP is the wave
! function, V is the potential and G is the nonlinearity. 
  USE COMM_DATA, ONLY : N,NX
  USE GPE_DATA, ONLY : DX, V, G 
  IMPLICIT NONE
  COMPLEX (8), DIMENSION(0:), INTENT(IN) :: CP
  REAL (8), INTENT(OUT) :: MU, EN
!----------------------------------------------------------------------
  INTERFACE
    PURE FUNCTION DIFF(P, DX) RESULT (DP)
      IMPLICIT NONE
      REAL (8), DIMENSION(0:), INTENT(IN) :: P
      REAL (8), INTENT(IN) :: DX
      REAL (8), DIMENSION(0:SIZE(P)-1) :: DP
    END FUNCTION DIFF
  END INTERFACE
!----------------------------------------------------------------------
   INTERFACE
    PURE FUNCTION SIMP(F, DX) RESULT(VALUE)
      IMPLICIT NONE
      REAL (8), DIMENSION(0:), INTENT(IN) :: F
      REAL (8), INTENT(IN) :: DX
      REAL (8) :: VALUE
    END FUNCTION SIMP
  END INTERFACE
!----------------------------------------------------------------------
  REAL (8), DIMENSION(0:SIZE(CP)-1) :: TMP1D, P2, DP2, GP2
  REAL (8), DIMENSION(0:SIZE(CP)-1) :: P, DCP, EMP1D
!
  P = ABS(CP)
  DCP = DIFF(P, DX)
  DP2 = DCP * DCP
  P2 = ABS(CP) * ABS(CP)
  GP2 = G * P2
!   
      TMP1D = (V + GP2 ) * P2 + DP2
      EMP1D = (V + GP2/2.D0) * P2 + DP2    
!  
  MU = SIMP(TMP1D, DX)
  EN = SIMP(EMP1D, DX)
END SUBROUTINE CHEM

SUBROUTINE NORM(CP, ZNORM)
! Calculates the normalization of the wave function and sets it to unity.
  USE COMM_DATA, ONLY : N
  USE GPE_DATA, ONLY : DX
  IMPLICIT NONE
  COMPLEX (8), DIMENSION(0:), INTENT(INOUT) :: CP
  REAL (8), INTENT(OUT) :: ZNORM
  !----------------------------------------------------------------------
  INTERFACE
    PURE FUNCTION SIMP(F, DX) RESULT (VALUE)
      IMPLICIT NONE
      REAL (8), DIMENSION(0:), INTENT(IN)  :: F
      REAL (8), INTENT(IN) :: DX
      REAL (8) :: VALUE
    END FUNCTION SIMP
  END INTERFACE
  !----------------------------------------------------------------------
  REAL (8), DIMENSION(0:N) :: P, P2
!
  P = ABS(CP)
  P2 = P*P
  ZNORM = SQRT(SIMP(P2, DX))
!  CP = CP/ZNORM
END SUBROUTINE NORM

SUBROUTINE RAD(CP2, RMS)
! Calculates the root mean square size RMS.
  USE COMM_DATA, ONLY : N
  USE GPE_DATA, ONLY : DX, X2
  IMPLICIT NONE
  REAL (8), DIMENSION(0:), INTENT(IN) :: CP2
  REAL (8), INTENT(OUT) :: RMS
!----------------------------------------------------------------------
   INTERFACE
    FUNCTION SIMP(F, DX) RESULT (VALUE)
      IMPLICIT NONE
      REAL (8), DIMENSION(0:), INTENT(IN) :: F
      REAL (8), INTENT(IN) :: DX
      REAL (8) :: VALUE
    END FUNCTION SIMP
  END INTERFACE
!----------------------------------------------------------------------
 REAL (8), DIMENSION(0:SIZE(CP2)-1) :: P2
!
  P2 = X2 * CP2
  RMS = SQRT(SIMP(P2, DX))
END SUBROUTINE RAD


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
  VALUE = DX*(F(0) + 4.D0*F1 + 2.D0*F2 + F(N))/3.D0
END FUNCTION SIMP

PURE FUNCTION DIFF(P,DX)  RESULT (DP)
! Computes the first derivative DP of P using
! Richardson extrapolation formula. The derivative at the  
! boundaries are assumed to be zero
  IMPLICIT NONE
  REAL (8), DIMENSION(0:), INTENT(IN) :: P
  REAL (8), INTENT(IN) :: DX
  REAL (8), DIMENSION(0:SIZE(P)-1) :: DP
  INTEGER :: I, N
!
  N = SIZE(P) - 1
  DP(0) = 0.D0
  DP(1) = (P(2) - P(0))/(2.D0*DX)
  FORALL(I=2:N-2)
    DP(I) = (P(I-2)-8.D0*P(I-1)+8.D0*P(I+1)-P(I+2))/(12.D0*DX)
  END FORALL
  DP(N-1) = (P(N) - P(N-2))/(2.D0*DX)
  DP(N) = 0.D0
END FUNCTION DIFF

SUBROUTINE WRITE_DENSITY(FUNIT, U2)
  USE COMM_DATA, ONLY : N
  USE GPE_DATA, ONLY : X
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: FUNIT
  REAL (8), DIMENSION(0:), INTENT(IN) :: U2
  INTEGER :: I
!
  DO I = 0, N
    WRITE(FUNIT, 999) X(I), U2(I)
  END DO
  999 FORMAT(F12.4, E17.7E3)
END SUBROUTINE WRITE_DENSITY

!# File name : real1d-th.f90
