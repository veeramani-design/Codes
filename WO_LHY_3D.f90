MODULE COMM_DATA
  ! NX, NY, NZ : Number of space mesh points (X, Y and Z)  
    INTEGER, PARAMETER :: NX = 240, NXX = NX-1, NX2 = NX/2 ! NX=240
    INTEGER, PARAMETER :: NY = 200, NYY = NY-1, NY2 = NY/2 ! NY=200
    INTEGER, PARAMETER :: NZ = 160, NZZ = NZ-1, NZ2 = NZ/2 ! NZ=160
  ! NSTP : Number of iterations to introduce the nonlinearity.
  ! NPAS : Number of subsequent iterations with fixed nonlinearity.
  ! NRUN : Number of final time iterations with fixed nonlinearity.
    INTEGER, PARAMETER :: NSTP = 1, NPAS = 5000, NRUN = 500
    INTEGER, PARAMETER :: NNN = NX * NY * NZ
    INTEGER, PARAMETER :: NUMBER_OF_THREADS = 0 ! sets the number of CPU cores to be used
       !     NUMBER_OF_THREADS =0 deactivates the command and uses all available CPU cores   
    REAL (8), PARAMETER :: PI = 3.14159265358979D0
  END MODULE COMM_DATA
  !===============================================================================
  MODULE GPE_DATA
    USE COMM_DATA, ONLY : PI
    REAL (8), PARAMETER :: AHO = 1.0D-6                    ! Unit of length (= 1 MICRON)            
    REAL (8), PARAMETER :: Bohr_a0 =  5.2917720859D-11/AHO ! Bohr radius (scaled with AHO)
  !  
    REAL (8), PARAMETER :: DX = 0.005D0, DY = DX, DZ = DX   ! DX, DY, DZ : SPACE STEPS dx=0.005
    REAL (8), PARAMETER :: DT = 0.0004D0                   ! DT : TIME STEP  dt=0.0004
    INTEGER, PARAMETER  :: NATOMS = 1000                   ! Number of Atoms
    REAL (8), PARAMETER :: AS = 67.530979D0*Bohr_a0        ! Scattering length (in units of Bohr_a0)
    REAL (8), PARAMETER :: GAM = 1.D0, ANU = 1.414214D0, LAM = 2.D0  ! GAMMA, NU, LAMBDA : Parameteres of Trap in x, y and z directions
  !   
    REAL (8), PARAMETER :: G0 = 4.D0*PI*AS*NATOMS ! 44.907D0   ! Three-dimensional nonlinearity 
    REAL (8), PARAMETER :: LHY = 0  ! LHY correction term = 128*PI**0.5*AS**2.5*NATOMS**2.5/3 
  !
  ! OPTION decides which equation to be solved.
  ! OPTION=1 Solves -psi_xx-psi_yy-psi_zz+V(x,y,z)psi+G0|psi|^2 psi =i psi_t
  ! OPTION=2 Solves [-psi_xx-psi_yy-psi_zz+V(x,y,z)psi]/2+G0|psi|^2 psi =i psi_t
    INTEGER, PARAMETER :: OPTION = 2 
  ! X(0:NX), Y(0:NY), Z(0:NZ) : Space mesh, V(0:NX,0:NY,0:NZ) : Potential, 
  ! CP(0:NX,0:NY,0:NZ) : Wave function  
    REAL (8), DIMENSION(:), ALLOCATABLE :: X, X2, Y, Y2, Z, Z2
    REAL (8), DIMENSION(:,:,:), ALLOCATABLE :: V, R2, CP
    REAL (8), DIMENSION(:), ALLOCATABLE :: TMPX, TMPY, TMPZ
    REAL (8) :: G, XOP
  END MODULE GPE_DATA
  !===============================================================================
  MODULE CN_DATA
    REAL (8), DIMENSION(:), ALLOCATABLE :: CALA, CGAA, CALB, CGAB, CALC, CGAC, CBE
    REAL (8) :: CT0X, CT0Y, CT0Z
    REAL (8) :: CA0, CB0, CC0, CA0R, CB0R, CC0R
  END MODULE CN_DATA 
  !===============================================================================
  PROGRAM GROSS_PITAEVSKII_SSCN_3D
    USE COMM_DATA, ONLY : NX, NY, NZ, NX2, NY2, NZ2, NPAS, NSTP, NRUN, NUMBER_OF_THREADS
    USE GPE_DATA
    USE OMP_LIB
    IMPLICIT NONE
  ! Subroutine INTIALIZE() used to initialize the space mesh X(I), 
  ! potential V(I) and the initial wave function. Subroutine COEF() used to
  ! generate the coefficients for the Crank-Nicholson Scheme. The routine
  ! NU() performs time progation for the non-derivative part and LU() performs
  ! time propagation of derivative part. NORM() calculates the norm and 
  ! normalizes the wave function, CHEM() and RAD() are used to calculate the 
  ! chemical potential, energy and the rms radius, respectively. The function 
  ! DIFF() used to calculate the space derivatives of the wave function used 
  ! in CHEM() and SIMP() does the integration by Simpson's rule.
  !------------------------ interface blocks -------------------------------------
   INTERFACE 
      SUBROUTINE ALLOCATE_VARIABLES()
      END SUBROUTINE ALLOCATE_VARIABLES
    END INTERFACE
  !-------------------------------------------------------------------------------
    INTERFACE 
      SUBROUTINE FREE_VARIABLES()
      END SUBROUTINE FREE_VARIABLES
    END INTERFACE
  !-------------------------------------------------------------------------------
    INTERFACE 
      SUBROUTINE INITIALIZE()
      END SUBROUTINE INITIALIZE
    END INTERFACE
  !-------------------------------------------------------------------------------
    INTERFACE 
      SUBROUTINE CALCULATE_TRAP()
      END SUBROUTINE CALCULATE_TRAP
    END INTERFACE
  !-------------------------------------------------------------------------------
    INTERFACE 
      SUBROUTINE COEF()
      END SUBROUTINE COEF
    END INTERFACE
  !-------------------------------------------------------------------------------
    INTERFACE 
      SUBROUTINE CALCNU(CP, DT)
        REAL (8), DIMENSION(0:, 0:, 0:), INTENT(INOUT) :: CP
        REAL (8), INTENT(IN) :: DT
      END SUBROUTINE CALCNU
    END INTERFACE
  !-------------------------------------------------------------------------------
    INTERFACE 
      SUBROUTINE LUX(CP)
        REAL (8), DIMENSION(0:, 0:, 0:), INTENT(INOUT) :: CP
      END SUBROUTINE LUX
    END INTERFACE
  !-------------------------------------------------------------------------------
    INTERFACE 
      SUBROUTINE LUY(CP)
        REAL (8), DIMENSION(0:, 0:, 0:), INTENT(INOUT) :: CP
      END SUBROUTINE LUY
    END INTERFACE
  !-------------------------------------------------------------------------------
    INTERFACE 
      SUBROUTINE LUZ(CP)
        REAL (8), DIMENSION(0:, 0:, 0:), INTENT(INOUT) :: CP
      END SUBROUTINE LUZ
    END INTERFACE
  !-------------------------------------------------------------------------------
    INTERFACE 
      SUBROUTINE NORM(P, DX, DY, DZ, ZNORM)
        REAL (8), DIMENSION(0:,0:,0:), INTENT(INOUT) :: P
        REAL (8), INTENT(IN) :: DX, DY, DZ
        REAL (8), INTENT(OUT) :: ZNORM
      END SUBROUTINE NORM
    END INTERFACE
  !-------------------------------------------------------------------------------
    INTERFACE 
      SUBROUTINE RAD(P, R2, X2, Y2, Z2, DX, DY, DZ, RMS)
        REAL (8), DIMENSION(0:,0:,0:), INTENT(IN) :: P, R2
        REAL (8), DIMENSION(0:), INTENT(IN) :: X2, Y2, Z2
        REAL (8), INTENT(IN) :: DX, DY, DZ
        REAL (8), DIMENSION(:), INTENT(OUT) :: RMS
      END SUBROUTINE RAD
    END INTERFACE
  !-------------------------------------------------------------------------------
    INTERFACE
      SUBROUTINE CHEM(CP, MU, EN)
        REAL (8), DIMENSION(0:, 0:, 0:), INTENT(IN) :: CP
        REAL (8), INTENT(OUT) :: MU, EN
      END SUBROUTINE CHEM
    END INTERFACE
  !-------------------------------------------------------------------------------
    INTERFACE
      SUBROUTINE WRITE_2DXY(FUNIT, U2)
        INTEGER, INTENT(IN) :: FUNIT
        REAL (8), DIMENSION(0:,0:,0:), INTENT(IN) :: U2
      END  SUBROUTINE WRITE_2DXY
    END INTERFACE
  !-------------------------------------------------------------------------------
    INTERFACE
      SUBROUTINE WRITE_2DXZ(FUNIT, U2)
        INTEGER, INTENT(IN) :: FUNIT
        REAL (8), DIMENSION(0:,0:,0:), INTENT(IN) :: U2
      END  SUBROUTINE WRITE_2DXZ
    END INTERFACE
  !-------------------------------------------------------------------------------
    INTERFACE
      SUBROUTINE DEN2DYZ(FUNIT, U2, DE2DYZ)
        INTEGER, INTENT(IN) :: FUNIT
        REAL (8), DIMENSION(0:,0:,0:), INTENT(IN) :: U2
        REAL (8), DIMENSION(0:,0:), INTENT(OUT) :: DE2DYZ
      END  SUBROUTINE DEN2DYZ
    END INTERFACE
  !-------------------------------------------------------------------------------
    INTERFACE
      SUBROUTINE DEN1DY(FUNIT, U2, DE1DY)
        INTEGER, INTENT(IN) :: FUNIT
        REAL (8), DIMENSION(0:,0:,0:), INTENT(IN) :: U2
        REAL (8), DIMENSION(0:), INTENT(OUT) :: DE1DY
      END  SUBROUTINE DEN1DY
    END INTERFACE
  !-------------------------------------------------------------------------------
    INTERFACE
      SUBROUTINE WRITE_3D(FUNIT, U2)
        INTEGER, INTENT(IN) :: FUNIT
        REAL (8), DIMENSION(0:,0:,0:), INTENT(IN) :: U2
      END  SUBROUTINE WRITE_3D
    END INTERFACE
  !-------------------------------------------------------------------------------
    INTERFACE
      SUBROUTINE DEN1DZ(FUNIT, U2, DE1DZ)
        INTEGER, INTENT(IN) :: FUNIT
        REAL (8), DIMENSION(0:,0:,0:), INTENT(IN) :: U2
        REAL (8), DIMENSION(0:), INTENT(OUT) :: DE1DZ
      END  SUBROUTINE DEN1DZ
    END INTERFACE
  !-------------------------------------------------------------------------------
    INTERFACE
      SUBROUTINE DEN1DX(FUNIT, U2, DE1DX)
        INTEGER, INTENT(IN) :: FUNIT
        REAL (8), DIMENSION(0:,0:,0:), INTENT(IN) :: U2
        REAL (8), DIMENSION(0:), INTENT(OUT) :: DE1DX
      END  SUBROUTINE DEN1DX
    END INTERFACE
  !-------------------------------------------------------------------------------
    INTERFACE
      SUBROUTINE DEN2DXZ(FUNIT, U2, DE2DXZ)
        INTEGER, INTENT(IN) :: FUNIT
        REAL (8), DIMENSION(0:,0:,0:), INTENT(IN) :: U2
        REAL (8), DIMENSION(0:,0:), INTENT(OUT) :: DE2DXZ
      END  SUBROUTINE DEN2DXZ
    END INTERFACE
  !-------------------------------------------------------------------------------
    INTERFACE
      SUBROUTINE DEN2DXY(FUNIT, U2, DE2DXY)
        INTEGER, INTENT(IN) :: FUNIT
        REAL (8), DIMENSION(0:,0:,0:), INTENT(IN) :: U2
        REAL (8), DIMENSION(0:,0:), INTENT(OUT) :: DE2DXY
      END  SUBROUTINE DEN2DXY
    END INTERFACE
  !-------------------- END INTERFACE BLOCKS -------------------------------------
    REAL (8), DIMENSION(:,:,:), ALLOCATABLE :: CP2
    REAL (8), DIMENSION(:,:), ALLOCATABLE :: DEN2XY, DEN2XZ, DEN2YZ
    REAL (8), DIMENSION(:), ALLOCATABLE :: DEN1X, DEN1Y, DEN1Z
    REAL (8), DIMENSION(4) :: RMS
    REAL (8) :: ZNORM, MU, EN, GSTP, T1, T2 
    INTEGER :: I, J, K,NO_OF_THREADS
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
  
   !SELECT CASE (OPTION)
    !    CASE (1)
     !   XOP = 1.0D0
      !! XOP = 2.0D0
      !CASE (:0,3:)
       ! PRINT *, 'ERROR: Wrong option', OPTION
        !STOP
    !END SELECT
  
  !
    CALL ALLOCATE_VARIABLES()
    ALLOCATE(CP2(0:NX, 0:NY, 0:NZ))
    ALLOCATE(DEN2XY(0:NX,0:NY), DEN2XZ(0:NX,0:NZ), DEN2YZ(0:NY,0:NZ))
    ALLOCATE(DEN1X(0:NX), DEN1Y(0:NY), DEN1Z(0:NZ))
  !
      
  !
    PRINT *, "Writing to file..."  
    
    OPEN(7, FILE = 'imag3d-th-out.txt')
    OPEN(4, FILE = 'imag3d-th-rms.txt')
    WRITE(7,*) "This is a test line."
    WRITE(7,*) 'Number of Threads =', NO_OF_THREADS
    WRITE(7,900) OPTION
    WRITE(4,900) OPTION
    WRITE(7,*) 
    WRITE(4,*) 
    WRITE(7,901) NATOMS, AHO
    WRITE(7,902) AS/Bohr_a0 
    WRITE(7,903) G0 
    WRITE(7,904) GAM, ANU, LAM
    WRITE(7,*)
    WRITE(7,905) NX, NY, NZ
    WRITE(7,906) DX, DY, DZ
    WRITE(7,907) NSTP, NPAS, NRUN  
    WRITE(7,908) DT 
    WRITE(7,*)
    900 FORMAT(' Imaginary time propagation 3D, OPTION =',I3)
    901 FORMAT('  Number of Atoms N =',I10,', Unit of length AHO =',F12.8,' m')
    902 FORMAT('  Scattering length a = ',F9.2,'*a0')
    903 FORMAT('  Nonlinearity G_3D =',F16.7 )
    904 FORMAT('  Parameters of trap: GAMMA =',F7.4, ', NU =',F7.4, ', LAMBDA =',F7.4)
    905 FORMAT(' # Space Stp: NX = ', I8, ', NY = ', I8, ', NZ = ', I8)
    906 FORMAT('  Space Step: DX = ', F10.6, ', DY = ', F10.6, ', DZ = ', F10.6)
    907 FORMAT(' # Time Stp : NSTP = ',I9,', NPAS = ',I9,', NRUN = ',I9)
    908 FORMAT('   Time Step:   DT = ', F10.6)
  !
   
  ! INITIALIZE() initializes the starting normalized wave function 'CP'.
    CALL INITIALIZE()
  ! CALCULATE_TRAP() initializes the harmonic oscillator potential V.
    CALL CALCULATE_TRAP()
  ! COEF() defines the coefficients of Crank-Nicholson Scheme.
    CALL COEF()
  ! NORM() calculates norm and restores normalization.
    CALL NORM(CP, DX, DY, DZ, ZNORM)
  ! NORM() calculates norm and restores normalization
    CALL RAD(CP, R2, X2, Y2, Z2, DX, DY, DZ, RMS)
  !  CHEM() calculates the chemical potential MU and energy EN.
    CALL CHEM(CP, MU, EN)
  !
   
    !$OMP PARALLEL DO PRIVATE(I, J, K)
    DO I = 0, NX; DO J = 0, NY;  DO K = 0, NZ
       CP2(I,J,K) = CP(I,J,K) * CP(I,J,K)
    END DO; END DO; END DO
    !$OMP END PARALLEL DO
  !
    WRITE (7, 1001)
    WRITE (7, 1002)
    WRITE (7, 1001)
    WRITE (7, 1003) ZNORM, MU/XOP, EN/XOP, RMS(1), CP2(NX2, NY2, NZ2) 
    WRITE (4, 1001)
    WRITE (4, 1012)
    WRITE (4, 1001)
    WRITE (4, 1013) RMS(1:4) 
    1001 FORMAT (19X,'-----------------------------------------------------')
    1002 FORMAT (20X, 'Norm', 6X, 'Chem', 8X, 'Ener/N', 5X, '<r>', 5X, '|Psi(0,0,0)|^2')
    1003 FORMAT ('Initial : ', 4X, F11.4, 2F12.4, 2F11.4)
    1012 FORMAT ('Values of rms size:', 5X, '<r>', 10X, '<x>', 10X, '<y>', 10X, '<z>')
    1013 FORMAT (11X, 'Initial:', F12.5, 3F13.5)
  ! 
  !    OPEN(31, FILE = 'imag3d-th-den-init.txt')
  !    CALL WRITE_3D(31, CP2)
  !    CLOSE(31)
  !  
  !    OPEN(11, FILE = 'imag3d-th-den-init1d_z.txt')
  !    OPEN(12, FILE = 'imag3d-th-den-init1d_x.txt')
  !    OPEN(13, FILE = 'imag3d-th-den-init1d_y.txt')
  !    CALL DEN1DZ(11, CP2, DEN1Z)
  !    CALL DEN1DX(12, CP2, DEN1X)
  !    CALL DEN1DY(13, CP2, DEN1Y)
  !    CLOSE(11)
  !    CLOSE(12)
  !    CLOSE(13)
  ! 
    !OPEN(21, FILE = 'imag3d-th-den-init2d_xz.txt')
    !OPEN(22, FILE = 'imag3d-th-den-init2d_xy.txt')
    !OPEN(23, FILE = 'imag3d-th-den-init2d_yz.txt')
    !CALL DEN2DXZ(21, CP2, DEN2XZ)
    !CALL DEN2DXY(22, CP2, DEN2XY)
    !CALL DEN2DYZ(23, CP2, DEN2YZ)
    !CLOSE(21)
    !CLOSE(22)
    !CLOSE(23)
  
    !OPEN(32, FILE = 'imag3d-th-den-init3d_x0z.txt')
    !OPEN(33, FILE = 'imag3d-th-den-init3d_xy0.txt')
    !CALL WRITE_2DXZ(32, CP2)
    !CALL WRITE_2DXY(33, CP2)
    !CLOSE(32)
    !CLOSE(33)
  !
    IF (NSTP /= 0) THEN
      GSTP =  XOP * G0 / DFLOAT(NSTP) ! GSTP = Nonlinearity step size 
      DO K = 1, NSTP ! NSTP iterations to introduce the nonlinearity  NSTP=1
        G = G + GSTP
        CALL CALCNU(CP, DT)
        CALL LUX(CP)
        CALL LUY(CP)
        CALL LUZ(CP)
       CALL NORM(CP, DX, DY, DZ, ZNORM)
      END DO
      CALL RAD(CP, R2, X2, Y2, Z2, DX, DY, DZ, RMS)
       CALL CHEM(CP, MU, EN)
  !
      !$OMP PARALLEL DO PRIVATE(I, J, K)
      DO K = 0, NZ; DO J = 0, NY;  DO I = 0, NX
        CP2(I,J,K) = CP(I,J,K) * CP(I,J,K)
      END DO; END DO; END DO
      !$OMP END PARALLEL DO
  
      WRITE (7, 1005) ZNORM, MU/XOP, EN/XOP, RMS(1), CP2(NX2, NY2, NZ2)
      WRITE (4, 1015) RMS(1:4)
      1005 FORMAT('After NSTP iter.:', F8.4, 2F12.4, 2F11.4)
      1015 FORMAT (2X, 'After NSTP iter.:', F12.4, 3F13.4)
  !
  !   OPEN(31, FILE = 'imag3d-th-den-nstp.txt')
  !   CALL WRITE_3D(31, CP2)
  !   CLOSE(31)
  ! 
  !   OPEN(11, FILE = 'imag3d-th-den-nstp1d_z.txt')
  !   OPEN(12, FILE = 'imag3d-th-den-nstp1d_x.txt')
  !   OPEN(13, FILE = 'imag3d-th-den-nstp1d_y.txt')
  !   CALL DEN1DZ(11, CP2, DEN1Z)
  !   CALL DEN1DX(12, CP2, DEN1X)
  !   CALL DEN1DY(13, CP2, DEN1Y)
  !   CLOSE(11)
  !   CLOSE(12)
  !   CLOSE(13)
  ! 
    !OPEN(21, FILE = 'imag3d-th-den-nstp2d_xz.txt')
    !OPEN(22, FILE = 'imag3d-th-den-nstp2d_xy.txt')
    !OPEN(23, FILE = 'imag3d-th-den-nstp2d_yz.txt')
    !CALL DEN2DXZ(21, CP2, DEN2XZ)
    !CALL DEN2DXY(22, CP2, DEN2XY)
    !CALL DEN2DYZ(23, CP2, DEN2YZ)
    !CLOSE(21)
    !CLOSE(22)
    !CLOSE(23)
  
    !OPEN(32, FILE = 'imag3d-th-den-nstp3d_x0z.txt')
    !OPEN(33, FILE = 'imag3d-th-den-nstp3d_xy0.txt')
    !CALL WRITE_2DXZ(32, CP2)
    !CALL WRITE_2DXY(33, CP2)
    !CLOSE(32)
    !CLOSE(33)
  !
    ELSE
       G = XOP * G0    ! when NSTP=0
    END IF
  !
    DO K = 1, NPAS ! NPAS iterations transient (running with fixed non linearity)
       CALL CALCNU(CP, DT)
       CALL LUX(CP)
       CALL LUY(CP)
       CALL LUZ(CP)
       CALL NORM(CP, DX, DY, DZ, ZNORM)
    END DO
    IF(NPAS /= 0) THEN
      CALL RAD(CP, R2, X2, Y2, Z2, DX, DY, DZ, RMS)
       CALL CHEM(CP, MU, EN)
  !
      !$OMP PARALLEL DO PRIVATE(I, J, K)
      DO K = 0, NZ; DO J = 0, NY;  DO I = 0, NX
        CP2(I,J,K) = CP(I,J,K) * CP(I,J,K)
      END DO; END DO; END DO
      !$OMP END PARALLEL DO
      WRITE (7, 1006) ZNORM, MU/XOP, EN/XOP, RMS(1), CP2(NX2, NY2, NZ2)
      WRITE (4, 1016) RMS(1:4)
      1006 FORMAT('After NPAS iter.:',F8.4, 2F12.4, 2F11.4)
      1016 FORMAT (2X, 'After NPAS iter.:', F12.4, 3F13.4)
  !
  !    OPEN(31, FILE = 'imag3d-th-den-npas.txt')
  !    CALL WRITE_3D(31, CP2)
  !    CLOSE(31)
  ! 
  !    OPEN(11, FILE = 'imag3d-th-den-npas1d_z.txt')
  !    OPEN(12, FILE = 'imag3d-th-den-npas1d_x.txt')
  !    OPEN(13, FILE = 'imag3d-th-den-npas1d_y.txt')
  !    CALL DEN1DZ(11, CP2, DEN1Z)
  !    CALL DEN1DX(12, CP2, DEN1X)
  !    CALL DEN1DY(13, CP2, DEN1Y)
  !    CLOSE(11)
  !    CLOSE(12)
  !    CLOSE(13)
  ! 
     !OPEN(21, FILE = 'imag3d-th-den-npas2d_xz.txt')
     !OPEN(22, FILE = 'imag3d-th-den-npas2d_xy.txt')
     !OPEN(23, FILE = 'imag3d-th-den-npas2d_yz.txt')
     !CALL DEN2DXZ(21, CP2, DEN2XZ)
     !CALL DEN2DXY(22, CP2, DEN2XY)
     !CALL DEN2DYZ(23, CP2, DEN2YZ)
     !CLOSE(21)
     !CLOSE(22)
     !CLOSE(23)
  
     !OPEN(32, FILE = 'imag3d-th-den-npas3d_x0z.txt')
     !OPEN(33, FILE = 'imag3d-th-den-npas3d_xy0.txt')
     !CALL WRITE_2DXZ(32, CP2)
     !CALL WRITE_2DXY(33, CP2)
     !CLOSE(32)
     !CLOSE(33)
    END IF  
  !
    DO K = 1, NRUN  ! final time iterations with fixed non linearity 
       CALL CALCNU(CP, DT)
       CALL LUX(CP)
       CALL LUY(CP)
       CALL LUZ(CP)
       CALL NORM(CP, DX, DY, DZ, ZNORM)
    END DO
    IF(NRUN /= 0)THEN    
      CALL RAD(CP, R2, X2, Y2, Z2, DX, DY, DZ, RMS)
       CALL CHEM(CP, MU, EN)
  !
      !$OMP PARALLEL DO PRIVATE(I, J, K)
      DO K = 0, NZ; DO J = 0, NY;  DO I = 0, NX
        CP2(I,J,K) = CP(I,J,K) * CP(I,J,K)
      END DO; END DO; END DO
      !$OMP END PARALLEL DO
      WRITE (7, 1007) ZNORM, MU/XOP, EN/XOP, RMS(1), CP2(NX2, NY2, NZ2)
      WRITE (4, 1017) RMS(1:4)
    END IF 
    1007 FORMAT('After NRUN iter.:',F8.4, 2F12.4, 2F11.4)
    1017 FORMAT (2X, 'After NRUN iter.:', F12.4, 3F13.4)
  ! 
    OPEN(31, FILE = 'imag3d-th-den.txt')
    CALL WRITE_3D(31, CP2)
    CLOSE(31)
  
    OPEN(11, FILE = 'imag3d-th-den1d_z.txt')
    OPEN(12, FILE = 'imag3d-th-den1d_x.txt')
    OPEN(13, FILE = 'imag3d-th-den1d_y.txt')
    CALL DEN1DZ(11, CP2, DEN1Z)
    CALL DEN1DX(12, CP2, DEN1X)
    CALL DEN1DY(13, CP2, DEN1Y)
    CLOSE(11)
    CLOSE(12)
    CLOSE(13)
  !
    !OPEN(21, FILE = 'imag3d-th-den2d_xz.txt')
    !OPEN(22, FILE = 'imag3d-th-den2d_xy.txt')
    !OPEN(23, FILE = 'imag3d-th-den2d_yz.txt')
    !CALL DEN2DXZ(21, CP2, DEN2XZ)
    !CALL DEN2DXY(22, CP2, DEN2XY)
    !CALL DEN2DYZ(23, CP2, DEN2YZ)
    !CLOSE(21)
    !CLOSE(22)
    !CLOSE(23)
  
    !OPEN(32, FILE = 'imag3d-th-den3d_x0z.txt')
    !OPEN(33, FILE = 'imag3d-th-den3d_xy0.txt')
    !CALL WRITE_2DXZ(32, CP2)
    !CALL WRITE_2DXY(33, CP2)
    !CLOSE(32)
    !CLOSE(33)
  !
    CALL FREE_VARIABLES()
    IF (ALLOCATED(CP2))    DEALLOCATE(CP2)
    IF (ALLOCATED(DEN2XY)) DEALLOCATE(DEN2XY)
    IF (ALLOCATED(DEN2XZ)) DEALLOCATE(DEN2XZ)
    IF (ALLOCATED(DEN2YZ)) DEALLOCATE(DEN2YZ)
    IF (ALLOCATED(DEN1X))  DEALLOCATE(DEN1X)
    IF (ALLOCATED(DEN1Y))  DEALLOCATE(DEN1Y)
    IF (ALLOCATED(DEN1Z))  DEALLOCATE(DEN1Z)
  !  
    CALL SYSTEM_CLOCK (CLCK_COUNTS_END, CLCK_RATE)
    CALL CPU_TIME(T2)
    WRITE (7, 1001)  
    WRITE (4, 1001)
    CLOSE (4)
    WRITE (7,*)
    WRITE (7,'(A,I7,A)') ' Clock Time: ', (CLCK_COUNTS_END - CLCK_COUNTS_BEG) / &
       INT(CLCK_RATE,8), ' seconds'
    WRITE (7,'(A,I7,A)') '   CPU Time: ', INT(T2-T1), ' seconds' 
    CLOSE (7)
  END PROGRAM GROSS_PITAEVSKII_SSCN_3D
  !===============================================================================
  SUBROUTINE ALLOCATE_VARIABLES()
    USE COMM_DATA, ONLY : NX, NY, NZ
    USE GPE_DATA, ONLY : X, Y, Z, X2, Y2, Z2, V, CP, R2, TMPX, TMPY, TMPZ
    USE CN_DATA, ONLY : CALA, CGAA, CALB, CGAB, CALC, CGAC, CBE
    IMPLICIT NONE
    ALLOCATE(TMPX(0:NX))
    ALLOCATE(TMPY(0:NY))
    ALLOCATE(TMPZ(0:NZ))
    ALLOCATE(CBE(0:MAXVAL([NX,NY,NZ])-1))
    ALLOCATE(X(0:NX))   
    ALLOCATE(X2(0:NX))   
    ALLOCATE(Y(0:NY))   
    ALLOCATE(Y2(0:NY))   
    ALLOCATE(Z(0:NZ))   
    ALLOCATE(Z2(0:NZ))
    ALLOCATE(V(0:NX, 0:NY, 0:NZ))
    ALLOCATE(CP(0:NX, 0:NY, 0:NZ))
    ALLOCATE(R2(0:NX, 0:NY, 0:NZ))
    ALLOCATE(CALA(0:NX))   
    ALLOCATE(CGAA(0:NX))
    ALLOCATE(CALB(0:NY))   
    ALLOCATE(CGAB(0:NY))
    ALLOCATE(CALC(0:NZ))   
    ALLOCATE(CGAC(0:NZ))
  END SUBROUTINE ALLOCATE_VARIABLES
  !===============================================================================
  SUBROUTINE FREE_VARIABLES()
    USE GPE_DATA, ONLY : X, Y, Z, X2, Y2, Z2, V, CP, R2, TMPX, TMPY, TMPZ
    USE CN_DATA, ONLY : CALA, CGAA, CALB, CGAB, CALC, CGAC, CBE
    IMPLICIT NONE
    IF (ALLOCATED(TMPX)) DEALLOCATE(TMPX)
    IF (ALLOCATED(TMPY)) DEALLOCATE(TMPY)
    IF (ALLOCATED(TMPZ)) DEALLOCATE(TMPZ)
    IF (ALLOCATED(CBE)) DEALLOCATE(CBE)
    IF (ALLOCATED(X)) DEALLOCATE(X)   
    IF (ALLOCATED(X2)) DEALLOCATE(X2)   
    IF (ALLOCATED(Y)) DEALLOCATE(Y)   
    IF (ALLOCATED(Y2)) DEALLOCATE(Y2)   
    IF (ALLOCATED(Z)) DEALLOCATE(Z)   
    IF (ALLOCATED(Z2)) DEALLOCATE(Z2)
    IF (ALLOCATED(V)) DEALLOCATE(V)
    IF (ALLOCATED(CP)) DEALLOCATE(CP)
    IF (ALLOCATED(R2)) DEALLOCATE(R2)
    IF (ALLOCATED(CALA)) DEALLOCATE(CALA)
    IF (ALLOCATED(CGAA)) DEALLOCATE(CGAA)
    IF (ALLOCATED(CALB)) DEALLOCATE(CALB)
    IF (ALLOCATED(CGAB)) DEALLOCATE(CGAB)
    IF (ALLOCATED(CALC)) DEALLOCATE(CALC)
    IF (ALLOCATED(CGAC)) DEALLOCATE(CGAC)
  END SUBROUTINE FREE_VARIABLES
  !===============================================================================
  SUBROUTINE INITIALIZE()
  !   Routine that initizlizes the constant and variables.
  !   Calculates the initial wave function CP
    USE COMM_DATA, ONLY : NX, NY, NZ, NX2, NY2, NZ2, PI
    USE GPE_DATA, ONLY : DX, DY, DZ, ANU, LAM, GAM, X, Y, Z, X2, Y2, Z2, R2, CP
    IMPLICIT NONE
    REAL (8) ::  SP, PI34
    REAL (8), DIMENSION(0:NX,0:NY,0:NZ) :: TMP
    INTEGER :: I, J, K
  
    PI34 = SQRT(PI*SQRT(PI))
    SP = PI34 /SQRT(SQRT(ANU * LAM * GAM))
  
    FORALL (I=0:NX) X(I) = (I-NX2)*DX
    FORALL (J=0:NY) Y(J) = (J-NY2)*DY
    FORALL (K=0:NZ) Z(K) = (K-NZ2)*DZ
    X2 = X*X
    Y2 = Y*Y
    Z2 = Z*Z
  !
    !$OMP PARALLEL DO PRIVATE( I, J, K)
    DO I = 0, NX; DO J = 0, NY; DO K = 0, NZ
       TMP(I,J,K) =  (GAM*X2(I) + ANU * Y2(J) + LAM * Z2(K))
       R2(I,J,K) = X2(I) + Y2(J) + Z2(K)
       CP(I,J,K) = EXP(-TMP(I,J,K)/2.0D0)/SP  ! doubt 
    END DO; END DO; END DO
    !$OMP END PARALLEL DO
  END SUBROUTINE INITIALIZE
  !===============================================================================
  SUBROUTINE CALCULATE_TRAP()
    USE COMM_DATA, ONLY : NX, NY, NZ
    USE GPE_DATA, ONLY : XOP, GAM, ANU, LAM, V, X2, Y2, Z2
    IMPLICIT NONE
    INTEGER :: I, J, K
    REAL (8) :: ANU2, LAM2, GAM2
    GAM2 = GAM * GAM
    ANU2 = ANU * ANU
    LAM2 = LAM * LAM 
    !$OMP PARALLEL DO PRIVATE(I, J, K)
    DO I = 0, NX; DO J = 0, NY; DO K = 0, NZ
       V(I,J,K) = XOP *  (X2(I)*GAM2 + ANU2*Y2(J) + LAM2*Z2(K))/2.0D0
    END DO; END DO; END DO
    !$OMP END PARALLEL DO
  END SUBROUTINE CALCULATE_TRAP!===============================================================================
  SUBROUTINE COEF()
  ! Calculates the coefficients needed in subroutines LUX, LUY, LUZ.      
    USE COMM_DATA, ONLY : NXX, NYY, NZZ
    USE GPE_DATA, ONLY : DX, DY, DZ, DT
    USE CN_DATA
    IMPLICIT NONE
    INTEGER :: I, J, K
    REAL (8) :: DX2, DY2, DZ2
    REAL (8) :: DXX, DYY, DZZ
  !
    DX2 = DX*DX ! generating the coefficients for the
    DY2 = DY*DY ! c-n method (to solve the spatial part)
    DZ2 = DZ*DZ
  !     
    DXX = 1.D0/DX2
    DYY = 1.D0/DY2
    DZZ = 1.D0/DZ2
  !
    CA0 = 1.D0 + DT*DXX
    CA0R = 1.D0 - DT*DXX
    CB0 = 1.D0 + DT*DYY
    CB0R = 1.D0 - DT*DYY
    CC0 = 1.D0 + DT*DZZ
    CC0R = 1.D0 - DT*DZZ
  !
    CT0X = -DT*DXX/2.D0
    CALA(NXX) = 0.D0
    CGAA(NXX) = -1.D0/CA0
    DO I = NXX, 1, -1
       CALA(I-1) = CT0X*CGAA(I)
       CGAA(I-1) = -1.D0/(CA0+CT0X*CALA(I-1))
    END DO
  !
    CT0Y = -DT*DYY/2.D0
    CALB(NYY) = 0.D0
    CGAB(NYY) = -1.D0/CB0
    DO J = NYY, 1, -1
       CALB(J-1) = CT0Y*CGAB(J)
       CGAB(J-1) = -1.D0/(CB0+CT0Y*CALB(J-1))
    END DO
  !
    CT0Z = -DT*DZZ/2.D0
    CALC(NZZ) = 0.D0
    CGAC(NZZ) = -1.D0/CC0
    DO K = NZZ, 1, -1
       CALC(K-1) = CT0Z*CGAC(K)
       CGAC(K-1) = -1.D0/(CC0+CT0Z*CALC(K-1))
    END DO
  END SUBROUTINE COEF
  !===============================================================================
  SUBROUTINE CALCNU(CP, DT) ! Exact solution
  !  Solves the partial differential equation with the potential and the
  !  nonlinear term and the LHY correction term.      
    USE COMM_DATA, ONLY : NX, NY, NZ
    USE GPE_DATA, ONLY : V, G , LHY
    IMPLICIT NONE
    REAL (8), DIMENSION(0:,0:,0:), INTENT(INOUT) :: CP
    REAL (8), INTENT(IN) :: DT
  ! 
   REAL (8), DIMENSION(0:NX, 0:NY, 0:NZ) :: P2, TMP, P3
    INTEGER :: I, J, K
  !
   !$OMP PARALLEL DO PRIVATE(I,J,K)
    DO K = 0, NZ; DO J = 0, NY; DO I = 0, NX
      P2(I,J,K)  = CP(I,J,K) * CP(I,J,K)
      P3(I,J,K)  = CP(I,J,K) * CP(I,J,K) * CP(I,J,K)
      TMP(I,J,K) = V(I,J,K) + G * P2(I,J,K) + LHY * P3(I,J,K)
    END DO; END DO; END DO
    !$OMP END PARALLEL DO
   !$OMP PARALLEL DO PRIVATE(I,J,K)
    DO K = 0, NZ; DO J = 0, NY; DO I = 0, NX
      P2(I,J,K)  =   EXP(-DT * TMP(I,J,K))
       CP(I,J,K) = CP(I,J,K)*  P2(I,J,K)
    END DO; END DO; END DO
    !$OMP END PARALLEL DO
  END SUBROUTINE CALCNU
  !===============================================================================
  SUBROUTINE LUX(CP)
  !  Solves the partial differential equation only with the X-space
  !  derivative term using the Crank-Nicholson method 
    USE COMM_DATA
    USE CN_DATA, ONLY : CA0R, CT0X, CALA, CGAA, CBE
    IMPLICIT NONE
    REAL (8), DIMENSION(0:,0:,0:), INTENT(INOUT) :: CP
    REAL (8), DIMENSION(0:NX,0:NY,0:NZ) :: TMP
    INTEGER :: I, J, K
  !
    !$OMP PARALLEL DO PRIVATE(I,J, K,CBE, TMP)
    DO J = 0, NY; DO K = 0, NZ
      CBE(NXX) = CP(NXX,J,K) ! boundary conditions
      DO I = NXX, 1, -1
        TMP(I,J,K) = -CT0X * CP(I+1,J,K) + CA0R * CP(I,J,K) - CT0X * CP(I-1,J,K)
        CBE(I-1) = CGAA(I) * (CT0X*CBE(I)-TMP(I,J,K))
      END DO
  !-----------------------------
  ! Boundary condition periodic:
  !    DO I = 0, NXX
  !      CP(I+1,J,K) = CALA(I)*CP(I,J,K) + CBE(I,J,K)
  !    END DO
  !    CP(0,J,K) = CP(NX,J,K)
  !    CP(1,J,K) = CP(NXX,J,K)
  !-----------------------------
  ! Boundary condition reflecting:
      CP(0,J,K) = 0.0D0  ! what is J and K
      DO I = 0, NXX
        CP(I+1,J,K) = CALA(I)*CP(I,J,K) + CBE(I)
      END DO
      CP(NX,J,K) = 0.0D0
    END DO; END DO
    !$OMP END PARALLEL DO
  END SUBROUTINE LUX
  !===============================================================================
  SUBROUTINE LUY(CP)
  !  Solves the partial differential equation only with the Y-space
  !  derivative term using the Crank-Nicholson method
    USE COMM_DATA
    USE CN_DATA, ONLY : CB0R, CT0Y, CALB, CGAB, CBE
    IMPLICIT NONE
    REAL (8), DIMENSION(0:,0:,0:), INTENT(INOUT) :: CP
     REAL (8), DIMENSION(0:NX,0:NY,0:NZ) :: TMP
    INTEGER :: I, J, K
  !
    !$OMP PARALLEL DO PRIVATE(I, K,J, CBE,TMP) 
    DO K = 0, NZ; DO I = 0, NX 
      CBE(NYY) = CP(I,NYY,K)
      DO J = NYY, 1, -1
         TMP(I,J,K) = -CT0Y*CP(I,J+1,K)+CB0R*CP(I,J,K)-CT0Y*CP(I,J-1,K)
          CBE(J-1) = CGAB(J)*(CT0Y*CBE(J)-TMP(I,J,K))
      END DO
  !-----------------------------
  ! Boundary condition periodic:
  !    DO J = 0, NYY
  !      CP(I,J+1,K) = CALB(J)*CP(I,J,K) + CBE(I,J,K)
  !    END DO
  !    CP(I,0,K) = CP(I,NY,K)
  !    CP(I,1,K) = CP(I,NYY,K)
  !-----------------------------
  ! Boundary condition reflecting:
      CP(I,0,K) = 0.0D0
      DO J = 0, NYY
        CP(I,J+1,K) = CALB(J)*CP(I,J,K) + CBE(J)
      END DO
      CP(I,NY,K) = 0.0D0
    END DO; END DO
    !$OMP END PARALLEL DO
  END SUBROUTINE LUY
  !===============================================================================
  SUBROUTINE LUZ(CP)
  !  Solves the partial differential equation only with the Y-space
  !  derivative term using the Crank-Nicholson method
    USE COMM_DATA
    USE CN_DATA, ONLY : CC0R, CT0Z, CALC, CGAC, CBE
    IMPLICIT NONE
    REAL (8), DIMENSION(0:,0:,0:), INTENT(INOUT) :: CP
    REAL (8), DIMENSION(0:NX,0:NY,0:NZ) :: TMP
    INTEGER :: I, J, K
  
    !$OMP PARALLEL DO PRIVATE(I, J,  K,CBE,TMP)
    DO J = 0, NY; DO I = 0, NX 
      CBE(NZZ) = CP(I,J,NZZ)
      DO K = NZZ, 1, -1
         TMP(I,J,K) = -CT0Z*CP(I,J,K+1)+CC0R*CP(I,J,K)-CT0Z*CP(I,J,K-1)
         CBE(K-1) = CGAC(K)*(CT0Z*CBE(K)-TMP(I,J,K))
      END DO
  !-----------------------------
  ! Boundary condition periodic:
  !    DO K = 0, NZZ
  !      CP(I,J,K+1) = CALC(K)*CP(I,J,K) + CBE(K)
  !    END DO
  !    CP(I,J,0) = CP(I,J,NZ)
  !    CP(I,J,1) = CP(I,J,NZZ)
  !-----------------------------
  ! Boundary condition reflecting:
      CP(I,J,0) = 0.0D0
      DO K = 0, NZZ
        CP(I,J,K+1) = CALC(K)*CP(I,J,K) + CBE(K)
      END DO
      CP(I,J,NZ) = 0.0D0
    END DO; END DO
    !$OMP END PARALLEL DO
  END SUBROUTINE LUZ
  !===============================================================================
  SUBROUTINE NORM(P, DX, DY, DZ, ZNORM)
    USE COMM_DATA
    IMPLICIT NONE
    REAL (8), DIMENSION(0:,0:,0:), INTENT(INOUT) :: P
    REAL (8), INTENT(IN) :: DX, DY, DZ
    REAL (8), INTENT(OUT) :: ZNORM 
    REAL (8), DIMENSION(0:NX, 0:NY, 0:NZ) :: P2
  !
     INTERFACE
       FUNCTION INTEGRATE(U, DX, DY, DZ) RESULT(RES)
         REAL (8), DIMENSION(0:, 0:, 0:), INTENT(IN) :: U
         REAL (8), INTENT (IN) :: DX, DY, DZ
         REAL (8) :: RES
       END FUNCTION INTEGRATE
     END INTERFACE
  ! 
    INTEGER :: I, J, K  !, NX, NY, NZ
    REAL (8) :: SQ_NORM
   
    !$OMP PARALLEL DO PRIVATE(I, J, K)
    DO K = 0, NZ; DO J = 0, NY; DO I = 0, NX
        P2(I,J,K) = P(I,J,K) * P(I,J,K)
    END DO; END DO; END DO
    !$OMP END PARALLEL DO
  
    ZNORM = INTEGRATE(P2, DX, DY, DZ)
    SQ_NORM = SQRT(ZNORM)
  
    !$OMP PARALLEL DO PRIVATE(I, J, K)
    DO K = 0, NZ; DO J = 0, NY;  DO I = 0, NX
        P(I,J,K) = P(I,J,K)/SQ_NORM
    END DO; END DO; END DO
    !$OMP END PARALLEL DO
  END SUBROUTINE NORM
  !===============================================================================
  SUBROUTINE RAD(P, R2, X2, Y2, Z2, DX, DY, DZ, RMS)
    USE COMM_DATA
    IMPLICIT NONE
    REAL (8), DIMENSION(0:,0:,0:), INTENT(IN) :: P, R2
    REAL (8), DIMENSION(0:), INTENT(IN) :: X2, Y2, Z2
    REAL (8), INTENT(IN) :: DX, DY, DZ
    REAL (8), DIMENSION(:), INTENT(OUT) :: RMS
  !
     INTERFACE
       FUNCTION INTEGRATE(U, DX, DY, DZ) RESULT(RES)
         REAL (8), DIMENSION(0:, 0:, 0:), INTENT(IN) :: U
         REAL (8), INTENT (IN) :: DX, DY, DZ
         REAL (8) :: RES
       END FUNCTION INTEGRATE
     END INTERFACE
  !
  !  REAL (8), DIMENSION(0:SIZE(P,1)-1, 0:SIZE(P,2)-1, 0:SIZE(P,3)-1) :: P2, TMP3D
    REAL (8), DIMENSION(0:NX, 0:NY, 0:NZ) :: P2, TMP3D
    INTEGER :: I, J, K !, NX, NY, NZ
  ! 
    !$OMP PARALLEL DO PRIVATE(I, J, K)
    DO K = 0, NZ; DO J = 0, NY; DO I = 0, NX
       P2(I,J,K) = P(I,J,K) * P(I,J,K)
       TMP3D(I,J,K) = R2(I,J,K) * P2(I,J,K)
    END DO; END DO; END DO
    !$OMP END PARALLEL DO
    RMS(1) = SQRT(INTEGRATE(TMP3D, DX, DY, DZ))
  !
    !$OMP PARALLEL DO PRIVATE(I, J, K)
    DO K = 0, NZ; DO J = 0, NY; DO I = 0, NX
       TMP3D(I,J,K) = X2(I)*P2(I,J,K)
    END DO; END DO; END DO
    !$OMP END PARALLEL DO
    RMS(2) = SQRT(INTEGRATE(TMP3D, DX, DY, DZ))
    
    !$OMP PARALLEL DO PRIVATE(I, J, K)
    DO K = 0, NZ; DO J = 0, NY; DO I = 0, NX
       TMP3D(I,J,K) = Y2(J)*P2(I,J,K)
    END DO; END DO; END DO
    !$OMP END PARALLEL DO
    RMS(3) = SQRT(INTEGRATE(TMP3D, DX, DY, DZ))
  
    !$OMP PARALLEL DO PRIVATE(I, J, K)
    DO K = 0, NZ; DO J = 0, NY; DO I = 0, NX
       TMP3D(I,J,K) = Z2(K)*P2(I,J,K)
    END DO; END DO; END DO
    !$OMP END PARALLEL DO
    RMS(4) = SQRT(INTEGRATE(TMP3D, DX, DY, DZ))
  END SUBROUTINE RAD
  !===============================================================================
  SUBROUTINE CHEM(CP, MU, EN)
  !  Calculates the chemical potential MU and energy EN.  CP is the wave
  !  function, V is the potential and G is the nonlinearity.
     USE COMM_DATA, ONLY : NX, NY, NZ
     USE GPE_DATA, ONLY : DX, DY, DZ, V, G, LHY
     IMPLICIT NONE
     REAL (8), DIMENSION(0:, 0:, 0:), INTENT(IN) :: CP
     REAL (8), INTENT(OUT) :: MU, EN
  !
     INTERFACE
       PURE FUNCTION DIFF(P, DX) RESULT (DP)
         REAL (8), DIMENSION(0:), INTENT(IN) :: P
         REAL (8), INTENT(IN) :: DX
         REAL (8), DIMENSION(0:SIZE(P)-1) :: DP
       END FUNCTION DIFF
     END INTERFACE
  !
     INTERFACE
       FUNCTION INTEGRATE(U, DX, DY, DZ) RESULT(RES)
         REAL (8), DIMENSION(0:, 0:, 0:), INTENT(IN) :: U
         REAL (8), INTENT (IN) :: DX, DY, DZ
         REAL (8) :: RES
       END FUNCTION INTEGRATE
     END INTERFACE
  !
     INTEGER :: I, J, K 
    REAL (8), DIMENSION(0:NX, 0:NY, 0:NZ) :: P2, DPX, DPY, DPZ, P3
     REAL (8), DIMENSION(0:NX, 0:NY, 0:NZ) :: GP2, DP2, TMP3D, EMP3D, LHY2
  !
     !$OMP PARALLEL
     !$OMP DO PRIVATE(I,J)
     DO J = 0, NY; DO I = 0, NX
        DPZ(I,J,:) = DIFF(CP(I,J,:), DZ)
     END DO; END DO
     !$OMP END DO
  
     !$OMP DO PRIVATE(I,K)
     DO K = 0, NZ; DO I = 0, NX
        DPY(I,:,K) = DIFF(CP(I,:,K), DY)
     END DO; END DO
     !$OMP END DO
  
     !$OMP DO PRIVATE(J,K)
     DO J = 0, NY; DO K = 0, NZ
        DPX(:,J,K) = DIFF(CP(:,J,K), DX)
     END DO; END DO
     !$OMP END DO
     !$OMP END PARALLEL
  !
     !$OMP PARALLEL 
     !$OMP DO PRIVATE(I,J,K)
     DO K = 0, NZ; DO J = 0, NY; DO I = 0, NX
       P2(I,J,K) = CP(I,J,K) * CP(I,J,K)
       P3(I,J,K)  = CP(I,J,K) * CP(I,J,K) * CP(I,J,K)
       GP2(I,J,K) = G * P2(I,J,K)
       LHY2(I,J,K) = LHY * P3(I,J,K)
       DP2(I,J,K) = DPX(I,J,K)*DPX(I,J,K) + DPY(I,J,K)*DPY(I,J,K) + DPZ(I,J,K)*DPZ(I,J,K)
       TMP3D(I,J,K) = (V(I,J,K) + GP2(I,J,K) + LHY2(I,J,K)/4.0D0) * P2(I,J,K) + DP2(I,J,K)
       EMP3D(I,J,K) = (V(I,J,K) + GP2(I,J,K)/2.0D0 + LHY2(I,J,K)/5.0D0) * P2(I,J,K) + DP2(I,J,K)
     END DO; END DO; END DO
     !$OMP END DO
     !$OMP END PARALLEL
  !
     MU = INTEGRATE(TMP3D, DX, DY, DZ)
     EN = INTEGRATE(EMP3D, DX, DY, DZ)
  END SUBROUTINE CHEM
  !===============================================================================
  PURE FUNCTION DIFF(P,DX) RESULT (DP)
  ! Computes the first derivative DP of P using
  ! Richardsonextrapolation formula. The derivative at the  
  ! boundaries are assumed to be zero
    IMPLICIT NONE
    REAL (8), DIMENSION(0:), INTENT(IN) :: P
    REAL (8), INTENT(IN) :: DX
    REAL (8), DIMENSION(0:SIZE(P)-1) :: DP
    INTEGER :: I, N
    N = SIZE(P) - 1
    DP(0) = 0.0D0
    DP(1) = (P(2) - P(0))/(2.0D0*DX)
    FORALL(I=2:N-2)
      DP(I) = (P(I-2)-8.0D0*P(I-1)+8.0D0*P(I+1)-P(I+2))/(12.0D0*DX)
    END FORALL
    DP(N-1) = (P(N) - P(N-2))/(2.0D0*DX)
    DP(N) = 0.0D0
  END FUNCTION DIFF
  !===============================================================================
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
  
    N = SIZE(F) - 1
  
    F1 = F(1) + F(N-1) ! N EVEN
    F2 = F(2) 
    DO I = 3, N-3, 2
       F1 = F1 + F(I)
       F2 = F2 + F(I+1)
    END DO
    VALUE = DX*(F(0) + 4.0D0*F1 + 2.0D0*F2 + F(N))/3.0D0
  END FUNCTION SIMP
  !===============================================================================
  FUNCTION INTEGRATE(U, DX, DY, DZ) RESULT(RES)
    USE GPE_DATA, ONLY : TMPX, TMPY, TMPZ
    IMPLICIT NONE
    REAL (8), DIMENSION(0:, 0:, 0:), INTENT(IN) :: U
    REAL (8), INTENT (IN) :: DX, DY, DZ
    REAL (8) :: RES
  !
    INTERFACE
      PURE FUNCTION SIMP(F, DX) RESULT (VALUE)
        REAL (8), DIMENSION(0:), INTENT(IN) :: F
        REAL (8), INTENT(IN) :: DX
        REAL (8) :: VALUE
      END FUNCTION SIMP
    END INTERFACE
  !
    INTEGER :: I, J, K, NX, NY, NZ
  !
    NX = SIZE(U,1)-1
    NY = SIZE(U,2)-1
    NZ = SIZE(U,3)-1
    !$OMP PARALLEL DO PRIVATE(I, J, TMPX, TMPY)
    DO K = 0, NZ
      DO J = 0, NY
        DO I = 0, NX      
          TMPX(I) = U(I,J,K)
        END DO
        TMPY(J) = SIMP(TMPX, DX)
      END DO
      TMPZ(K) = SIMP(TMPY, DY)
    END DO
    !$OMP END PARALLEL DO
    RES = SIMP(TMPZ, DZ)
  END FUNCTION INTEGRATE
  !===============================================================================
  SUBROUTINE WRITE_3D(FUNIT, U2)
    USE COMM_DATA, ONLY : NX, NY, NZ
    INTEGER, INTENT(IN) :: FUNIT
    REAL (8), DIMENSION(0:,0:,0:), INTENT(IN) :: U2
    INTEGER :: I, J, K
  !
    DO I = 0, NX; DO J = 0, NY; DO K = 0, NZ
      WRITE(FUNIT, 1000) U2(I,J,K)
    END DO; END DO; END DO
    1000 FORMAT(E17.6E3)  
  END SUBROUTINE WRITE_3D
  !===============================================================================
  SUBROUTINE WRITE_2DXY(FUNIT, U2)
    USE COMM_DATA, ONLY : NX, NY, NZ2
    USE GPE_DATA, ONLY : X, Y
    INTEGER, INTENT(IN) :: FUNIT
    REAL (8), DIMENSION(0:,0:,0:), INTENT(IN) :: U2
    INTEGER :: I, J, K
  !
    K = NZ2
    DO I = 0, NX
      DO J = 0, NY
        WRITE(FUNIT, 1000) X(I), Y(J), U2(I,J,K)   
      END DO
      WRITE(FUNIT, *)
    END DO
    1000 FORMAT(2F10.2, E17.6E3)
  END SUBROUTINE WRITE_2DXY
  !===============================================================================
  SUBROUTINE WRITE_2DXZ(FUNIT, U2)
    USE COMM_DATA, ONLY : NX, NZ, NY2
    USE GPE_DATA, ONLY : X, Z
    INTEGER, INTENT(IN) :: FUNIT
    REAL (8), DIMENSION(0:,0:,0:), INTENT(IN) :: U2
    INTEGER :: I, J, K
  !
    J = NY2
    DO I = 0, NX
      DO K = 0, NZ
        WRITE(FUNIT, 1000) X(I), Z(K), U2(I,J,K)
      END DO
      WRITE(FUNIT, *)
    END DO
    1000 FORMAT(2F10.2, E17.6E3)
  END SUBROUTINE WRITE_2DXZ
  !===============================================================================
  SUBROUTINE DEN2DXY(FUNIT, U2, DE2DXY)
    USE COMM_DATA, ONLY : NX, NY  
    USE GPE_DATA, ONLY : X,  DZ, Y
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: FUNIT
    REAL (8), DIMENSION(0:, 0:, 0:), INTENT(IN) :: U2
    REAL (8), DIMENSION(0:, 0:), INTENT(OUT) :: DE2DXY
  !
    INTERFACE
      PURE FUNCTION SIMP(F, DX) RESULT (VALUE)
        REAL (8), DIMENSION(0:), INTENT(IN) :: F
        REAL (8), INTENT(IN) :: DX
        REAL (8) :: VALUE
      END FUNCTION SIMP
    END INTERFACE
  !
    INTEGER :: I, J
  !
    FORALL (I = 0:NX, J = 0:NY) DE2DXY(I,J) = SIMP(U2(I,J,0:), DZ)
    DO I = 0, NX
     DO J = 0,NY
       WRITE(FUNIT, 1000) X(I),Y(J), DE2DXY(I,J)
     END DO
      WRITE(FUNIT, 1000)
    END DO 
    1000 FORMAT(2F10.2, E17.6E3)
  END SUBROUTINE DEN2DXY
  !===============================================================================
  SUBROUTINE DEN2DXZ(FUNIT, U2, DE2DXZ)
    USE COMM_DATA, ONLY : NX, NZ  
    USE GPE_DATA, ONLY : X, DY, Z
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: FUNIT
    REAL (8), DIMENSION(0:, 0:, 0:), INTENT(IN) :: U2
    REAL (8), DIMENSION(0:, 0:), INTENT(OUT) :: DE2DXZ
  !
    INTERFACE
      PURE FUNCTION SIMP(F, DX) RESULT (VALUE)
        REAL (8), DIMENSION(0:), INTENT(IN) :: F
        REAL (8), INTENT(IN) :: DX
        REAL (8) :: VALUE
      END FUNCTION SIMP
    END INTERFACE
  !
    INTEGER :: I, K
  !	 
    FORALL (I = 0:NX, K = 0:NZ) DE2DXZ(I,K) = SIMP(U2(I,0:,K), DY)
    DO I = 0, NX
      DO K = 0,NZ
        WRITE(FUNIT, 1000) X(I),Z(K), DE2DXZ(I,K)
      END DO
       WRITE(FUNIT, 1000)
    END DO
    1000 FORMAT(2F10.2, E17.6E3) 
  END SUBROUTINE DEN2DXZ
  !===============================================================================
  SUBROUTINE DEN2DYZ(FUNIT, U2, DE2DYZ)
    USE COMM_DATA, ONLY : NY, NZ  
    USE GPE_DATA, ONLY : Y, DX, Z
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: FUNIT
    REAL (8), DIMENSION(0:, 0:, 0:), INTENT(IN) :: U2
    REAL (8), DIMENSION(0:, 0:), INTENT(OUT) :: DE2DYZ
  !
    INTERFACE
      PURE FUNCTION SIMP(F, DX) RESULT (VALUE)
        REAL (8), DIMENSION(0:), INTENT(IN) :: F
        REAL (8), INTENT(IN) :: DX
        REAL (8) :: VALUE
      END FUNCTION SIMP
    END INTERFACE
  !
    INTEGER :: J, K
  !
    FORALL (J = 0:NY, K = 0:NZ) DE2DYZ(J,K) = SIMP(U2(0:,J,K), DX)
    DO J = 0, NY
    DO K = 0,NZ
          WRITE(FUNIT, 1000) Y(J),Z(K), DE2DYZ(J,K)
    END DO
        WRITE(FUNIT, 1000)
    END DO
    1000 FORMAT(2F10.2, E17.6E3) 
  END SUBROUTINE DEN2DYZ
  !===============================================================================
  SUBROUTINE DEN1DX(FUNIT, U2, DE1DX)
    USE COMM_DATA, ONLY : NX, NZ
    USE GPE_DATA, ONLY : X, DY, DZ
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: FUNIT
    REAL (8), DIMENSION(0:, 0:, 0:), INTENT(IN) :: U2
    REAL (8), DIMENSION(0:), INTENT(OUT) :: DE1DX
    REAL (8), DIMENSION(0:NX, 0:NZ) :: TMP2D
  !
    INTERFACE
      PURE FUNCTION SIMP(F, DX) RESULT (VALUE)
        REAL (8), DIMENSION(0:), INTENT(IN) :: F
        REAL (8), INTENT(IN) :: DX
        REAL (8) :: VALUE
      END FUNCTION SIMP
    END INTERFACE
  !
    INTEGER :: I, K
  !
    FORALL (I = 0:NX, K = 0:NZ) TMP2D(I,K) = SIMP(U2(I,0:,K), DY)
    DO I = 0, NX 
      DE1DX(I) = SIMP(TMP2D(I,0:), DZ)
      WRITE(FUNIT, 1001) X(I), DE1DX(I)
    END DO
    1001 FORMAT(F10.2, E17.6E3)
  END SUBROUTINE DEN1DX
  !===============================================================================
  SUBROUTINE DEN1DY(FUNIT, U2, DE1DY)
    USE COMM_DATA, ONLY : NY, NZ
    USE GPE_DATA, ONLY : Y, DX, DZ  
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: FUNIT
    REAL (8), DIMENSION(0:, 0:, 0:), INTENT(IN) :: U2
    REAL (8), DIMENSION(0:), INTENT(OUT) :: DE1DY
    REAL (8), DIMENSION(0:NY, 0:NZ) :: TMP2D
  !
    INTERFACE
      PURE FUNCTION SIMP(F, DX) RESULT (VALUE)
        REAL (8), DIMENSION(0:), INTENT(IN) :: F
        REAL (8), INTENT(IN) :: DX
        REAL (8) :: VALUE
      END FUNCTION SIMP
    END INTERFACE
  !
    INTEGER :: J, K
  !
    FORALL (J = 0:NY, K = 0:NZ) TMP2D(J,K) = SIMP(U2(0:,J,K), DX)
    DO J = 0, NY 
      DE1DY(J) = SIMP(TMP2D(J,0:), DZ)
      WRITE(FUNIT, 1001) Y(J), DE1DY(J)
    END DO
    1001 FORMAT(F10.2, E17.6E3) 
  END SUBROUTINE DEN1DY
  !===============================================================================
  SUBROUTINE DEN1DZ(FUNIT, U2, DE1DZ)
    USE COMM_DATA, ONLY : NY, NZ
    USE GPE_DATA, ONLY : Z, DX, DY
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: FUNIT
    REAL (8), DIMENSION(0:, 0:, 0:), INTENT(IN) :: U2
    REAL (8), DIMENSION(0:), INTENT(OUT) :: DE1DZ
    REAL (8), DIMENSION(0:NY, 0:NZ) :: TMP2D
  !
    INTERFACE
      PURE FUNCTION SIMP(F, DX) RESULT (VALUE)
        REAL (8), DIMENSION(0:), INTENT(IN) :: F
        REAL (8), INTENT(IN) :: DX
        REAL (8) :: VALUE
      END FUNCTION SIMP
    END INTERFACE
  !
    INTEGER :: J, K
  !
    FORALL (J = 0:NY, K = 0:NZ) TMP2D(J,K) = SIMP(U2(0:,J,K), DX)
    DO K = 0, NZ 
      DE1DZ(K) = SIMP(TMP2D(0:,K), DY)
      WRITE(FUNIT, 1001) Z(K), DE1DZ(K)
    END DO
    1001 FORMAT(F10.2, E17.6E3)
  END SUBROUTINE DEN1DZ
  !===============================================================================
  