!# File name : bec-gp-rot-2d-th.f90
!# Last modified : 20 JANUARY 2019
!# Fortran program for Gross-Pitaevskii equation in two-dimensional ROTATING
!# anisotropic trap by imaginary time propagation (Fortran 90/95 Version)
!#
!# BEC-GP-ROT-OMP programs are developed by:
!#
!# R. Kishor Kumar
!# (Instituto de Fisica, Universidade de Sao Paulo, Sao Paulo, Brazil)
!#
!# Vladimir Loncar, Antun Balaz
!# (Scientific Computing Laboratory, Center for the Study of Complex Systems, Institute of Physics Belgrade, Serbia)
!#
!# Paulsamy Muruganandam
!# (Department of Physics, Bharathidasan University, Tiruchirappalli, Tamil Nadu, India)
!#
!# Sadhan K. Adhikari
!# (Instituto de Fisica Teorica, UNESP - Sao Paulo State University, Brazil)
!#
!# Public use and modification of these codes are allowed provided that the following papers are cited:
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
!# To compile:
!# (1) Intel Fortran Compiler
! ifort -O3 -openmp -w -mcmodel medium -shared-intel
!
!# (2) GNU Fortran (gfortran)
! gfortran -O3  -fopenmp  -w
!
!# (3) PGI Fortran (pgfortran)
! pgfortran -O3 -fast -mp=allcores

MODULE COMM_DATA
! NX, NY : Number of space mesh points (X and Y)
  INTEGER, PARAMETER :: NX = 256, NXX = NX-1, NX2 = NX/2 ! No. of X steps
  INTEGER, PARAMETER :: NY = 256, NYY = NY-1, NY2 = NY/2 ! No. of Y steps
! NSTP : Number of iterations to introduce the nonlinearity.
! NSTP=0 reads the wave function, /=0 calculates the wave function.
! NPAS : Number of subsequent iterations with fixed nonlinearity.
! NRUN : Number of final time iterations with fixed nonlinearity.
! For real time set NSTP=0, and read imag time wave function from im2d-fin-wf.txt
  INTEGER, PARAMETER :: NSTP = 1, NPAS = 2000000, NRUN = 1, ITER = NPAS/10
  INTEGER, PARAMETER :: NUMBER_OF_THREADS = 0 ! sets the number of CPU cores to be used
! NUMBER_OF_THREADS = 0 deactivates the command and uses all available CPU cores
  REAL (8), PARAMETER :: PI = 3.14159265358979D0, SQR_2PI = 2.506628274631D0
  INTEGER :: NO_OF_THREADS
END MODULE COMM_DATA

MODULE GPE_DATA
  USE COMM_DATA, ONLY : NX, NY, SQR_2PI, PI
  REAL (8), PARAMETER :: AHO = 1.D-6 ! Unit of length (l= 1 MICRON)
  REAL (8), PARAMETER :: Bohr_a0 =  5.2917720859D-11/AHO ! Bohr radius (scaled with AHO)
  COMPLEX (8), PARAMETER :: CI = (0.0D0,1.0D0)
  REAL (8), PARAMETER :: DX = 0.05D0, DY = 0.05D0 ! DX, DY : Space step
  INTEGER, PARAMETER  :: NATOMS = 10000 ! Number of Atoms
  REAL (8), PARAMETER :: AS = 3.769458264D0*Bohr_a0 ! Scattering length (in units of Bohr_a0)
  REAL (8), PARAMETER :: GAMMA = 1.D0, NU = 1.D0 ! GAMMA and NU : Parameteres of Trap
  REAL (8), PARAMETER :: D_Z = 0.1D0 ! D_Z : Axial Gaussian Width = l/SQRT(LAMBDA)

  REAL (8), PARAMETER :: G_2D =  4.D0 * PI * AS * NATOMS / (SQR_2PI*D_Z) ! G_2D : Nonlinearity in the 2D GP equation
  REAL (8), PARAMETER :: G_3D = G_2D * SQR_2PI * D_Z ! G_3D : Three-dimensional nonlinearity

  REAL (8), PARAMETER :: OMEGA = 0.8D0 ! Angular Frequency

! OPTION  decides which equation to be solved.
! OPTION=1 Solves -psi_xx-psi_yy+V(x,y)psi+G_2D|psi|^2 psi =i psi_t
! OPTION=2 Solves [-psi_xx-psi_yy+V(x,y)psi]/2+G_2D|psi|^2 psi =i psi_t
  INTEGER, PARAMETER :: OPTION = 2

! INTEGER, PARAMETER :: RANDOM=0 ! No random phase in the initial function
 INTEGER, PARAMETER :: RANDOM=1 ! Random phase included in the initial function

! OPTION_re_im switches between real and imaginary time propagation
  INTEGER, PARAMETER :: OPTION_re_im=1 !  Imaginary time propagation
  REAL (8), PARAMETER ::  DT = 0.00025D0 ! DT : Time step for imag time
!  INTEGER, PARAMETER ::  OPTION_re_im=2   !  Real time propagation
!  REAL (8), PARAMETER :: DT = 0.0001D0 ! DT : Time step for real time

!  SELECT INITIAL FUNCTION 
!  INTEGER, PARAMETER :: FUNCTION=0 !  GAUSSIAN function with no vortex
   INTEGER, PARAMETER :: FUNCTION=1 !  1VORTEX function with 1 vortex

! X(0:NX), Y(0:NY): Space mesh, V(0:NX,0:NY) : Potential, CP(0:NX,0:NY): Wave function
  REAL (8), DIMENSION(:), ALLOCATABLE :: X, X2, Y, Y2
  REAL (8), DIMENSION(:, :), ALLOCATABLE :: V,  P, P2
  COMPLEX (8), DIMENSION(:, :), ALLOCATABLE :: CP
  COMPLEX (8) :: CIJ
  REAL (8) :: G, XOP, SQ2PI_DZ,gstp
END MODULE GPE_DATA

MODULE CN_DATA
  USE COMM_DATA, ONLY : NX, NY
  COMPLEX (8), DIMENSION(0:NX) :: CBP, CBM
  COMPLEX (8), DIMENSION(0:NY) :: CAP, CAM
  COMPLEX (8), DIMENSION(:,:), ALLOCATABLE :: CALA, CALB, CGAA, CGAB
  COMPLEX (8) :: CT0X, CT0Y, CT0, CTMPX, CTMPY
  COMPLEX (8) :: CA0, CB0, CA0R, CB0R, C0
END MODULE CN_DATA

PROGRAM GROSS_PITAEVSKII_SSCN_2D
  USE OMP_LIB
  USE COMM_DATA, ONLY : NX, NY, NX2, NY2, NSTP, NPAS, NRUN, NUMBER_OF_THREADS,ITER
  USE GPE_DATA

  IMPLICIT NONE
! Subroutine INTIIALIZE() used to initialize the space mesh, and the initial
! wave function. Subroutine CALCULATE_TRAP() initializezpotential V(I),
! while subroutine COEF() is used to generate the coefficients for the
! Crank-Nicholson Scheme. The routine NU() performs time progation for the
! non-derivative part and LU() performs time propagation of derivative part.
! NORM() calculates the norm and normalizes the wave function, CHEM() and
! RAD() are used to calculate the chemical potential, energy and the RMS
! radius, respectively. The functions DIFF(), used to calculate the space
! derivatives of the wave function used in CHEM(). The function SIMP() does
! the integration by Simpson's rule.

!------------------------ INTERFACE BLOCKS -----------------------
  INTERFACE
    SUBROUTINE ALLOCATE_VARIABLES()
    END SUBROUTINE ALLOCATE_VARIABLES

    SUBROUTINE FREE_VARIABLES()
    END SUBROUTINE FREE_VARIABLES

    SUBROUTINE INITIALIZE(CP)
      USE COMM_DATA, ONLY : NX, NY
      IMPLICIT NONE
      COMPLEX (8), DIMENSION(0:NX, 0:NY), INTENT(INOUT) :: CP
    END SUBROUTINE INITIALIZE

    SUBROUTINE CALCULATE_TRAP()
    END SUBROUTINE CALCULATE_TRAP

    SUBROUTINE COEF()
      IMPLICIT NONE
    END SUBROUTINE COEF

    SUBROUTINE CALCNU() 
      IMPLICIT NONE
    END SUBROUTINE CALCNU

    SUBROUTINE LUX()
      IMPLICIT NONE
    END SUBROUTINE LUX

    SUBROUTINE LUY()
      IMPLICIT NONE
    END SUBROUTINE LUY

    SUBROUTINE CHEM(CP, MU, EN)
      USE COMM_DATA, ONLY : NX, NY
      IMPLICIT NONE
      COMPLEX (8), DIMENSION(0:NX, 0:NY), INTENT(IN) :: CP
      REAL (8), INTENT(OUT) :: MU, EN
    END SUBROUTINE CHEM

    SUBROUTINE WRITE_DENSITY(FUNIT, U2)
      USE COMM_DATA, ONLY : NX, NY
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: FUNIT
      REAL (8), DIMENSION(0:NX, 0:NY), INTENT(IN) :: U2
    END SUBROUTINE WRITE_DENSITY
  
    SUBROUTINE WRITE_WAVE_FUNCTION(FUNIT, CP)
      USE COMM_DATA, ONLY : NX, NY
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: FUNIT 
      COMPLEX (8), DIMENSION(0:NX,0:NY), INTENT(IN) :: CP
    END SUBROUTINE WRITE_WAVE_FUNCTION

    SUBROUTINE NORM(CP, ZNORM) 
      USE COMM_DATA, ONLY : NX, NY
      IMPLICIT NONE
      COMPLEX (8), DIMENSION(0:NX, 0:NY), INTENT(INOUT) :: CP
      REAL (8), INTENT(OUT) :: ZNORM
    END SUBROUTINE NORM

    SUBROUTINE NORMCHEM(CPR, ZNORM)
      USE COMM_DATA, ONLY : NX, NY
      IMPLICIT NONE
      REAL (8), DIMENSION(0:NX, 0:NY), INTENT(INOUT) :: CPR
      REAL (8), INTENT(OUT) :: ZNORM
    END SUBROUTINE NORMCHEM

    SUBROUTINE RAD(CP2, RMS)
      USE COMM_DATA, ONLY : NX, NY
      IMPLICIT NONE 
      REAL (8), DIMENSION(0:NX, 0:NY), INTENT(IN) :: CP2
      REAL (8), DIMENSION(:), INTENT(OUT) :: RMS
    END SUBROUTINE RAD

    SUBROUTINE DENSITY_1DX(FUNIT, U2, DE1DX)
      USE COMM_DATA, ONLY : NX, NY
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: FUNIT
      REAL (8), DIMENSION(0:NX, 0:NY), INTENT(IN) :: U2 
      REAL (8), DIMENSION(0:NX), INTENT(OUT) :: DE1DX
    END SUBROUTINE DENSITY_1DX

    SUBROUTINE DENSITY_1DY(FUNIT, U2, DE1DY)
      USE COMM_DATA, ONLY : NX, NY
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: FUNIT
      REAL (8), DIMENSION(0:NX, 0:NY), INTENT(IN) :: U2
      REAL (8), DIMENSION(0:NY), INTENT(OUT) :: DE1DY
    END SUBROUTINE DENSITY_1DY
  END INTERFACE
!------------------------ END INTERFACE BLOCKS -------------------

  INTEGER :: I, J, K, NO_OF_THREADS, STEP
  REAL (8) :: ZNORM, T1, T2
  REAL (8) :: MU, EN
  REAL (8), DIMENSION(:,:), ALLOCATABLE :: CP2
  REAL (8), DIMENSION(3) :: RMS
  
  REAL (8), DIMENSION(0:NX) :: DEN1X
  REAL (8), DIMENSION(0:NY) :: DEN1Y
  CHARACTER(256) :: ITER_FILENAME
  REAL(KIND(0.0D0)) :: START_D, END_D

  INTEGER :: CLCK_COUNTS_BEG, CLCK_COUNTS_END, CLCK_RATE, COUNT_MAX

  CALL SYSTEM_CLOCK ( CLCK_COUNTS_BEG, CLCK_RATE )
  CALL CPU_TIME (T1)

  IF (NUMBER_OF_THREADS /= 0) THEN
    CALL OMP_SET_NUM_THREADS(NUMBER_OF_THREADS)
  END IF

  !$OMP PARALLEL
  !$OMP MASTER
  NO_OF_THREADS = OMP_GET_NUM_THREADS()
  !$OMP END MASTER
  !$OMP END PARALLEL

  SELECT CASE (OPTION)
    CASE (1)
      XOP = 1.0D0
    CASE (2) 
      XOP = 2.0D0 
    CASE (:0,3:)
      PRINT *, 'ERROR: Wrong option', OPTION
      STOP
  END SELECT

  CALL ALLOCATE_VARIABLES()
  ALLOCATE(CP2(0:NX, 0:NY))

  CALL INITIALIZE(CP)

  SELECT CASE(OPTION_re_im)
    CASE (1)
      CIJ = CMPLX(1.0D0, 0.0D0)
      OPEN(7, FILE = 'im2d-out.txt')
      OPEN(4, FILE = 'im2d-rms.txt')
    CASE (2)
      CIJ = CMPLX(0.0D0, 1.0D0)
      OPEN(7, FILE = 're2d-out.txt') 
      OPEN(4, FILE = 're2d-rms.txt')
  END SELECT

  SELECT CASE(OPTION_re_im)
    CASE (1)
      WRITE(7,900) OPTION, NO_OF_THREADS
      WRITE(4,900) OPTION, NO_OF_THREADS
    CASE (2)  
      WRITE(7,800) OPTION, NO_OF_THREADS
      WRITE(4,800) OPTION, NO_OF_THREADS
  END SELECT

  WRITE(7,*)
  WRITE(4,*)
  WRITE(7,901) NATOMS, AHO
  WRITE(7,902) AS/Bohr_a0
  WRITE(7,903) G_3D, G_2D
  WRITE(7,904) GAMMA, NU
  WRITE(7,911) OMEGA
  WRITE(7,909) D_Z
  WRITE(7,*)

  WRITE(7,905) NX + 1, NY + 1
  WRITE(7,906) DX, DY
  WRITE(7,907) NSTP, NPAS, NRUN
  WRITE(7,908) DT 
  WRITE(7,*)

  900 FORMAT('Imaginary time propagation 2d,   OPTION = ',I0,', NUM_THREADS = ',I0)
  800 FORMAT('Real time propagation 2d,   OPTION = ',I0,', NUM_THREADS = ',I0)

  901 FORMAT('Number of Atoms N = ',I0,', Unit of length AHO = ',F10.8,' m')
  902 FORMAT('Scattering length a = ',F0.6,'*a0')
  903 FORMAT('Nonlinearity G_3D = ',F0.6, ', G_2D = ',F0.6) 
  904 FORMAT('Parameters of trap: GAMMA = ',F0.2, ', NU = ',F0.2)
  911 FORMAT('Parameters of rotation: ANG VEL = ',F6.4)
  909 FORMAT('Axial trap parameter = ',F0.2)
  905 FORMAT('Space step: NX = ', I0, ', NY = ', I0)
  906 FORMAT('            DX = ', F8.6, ', DY = ', F8.6)
  907 FORMAT('Time step:  NSTP = ',I0,', NPAS = ',I0,', NRUN = ',I0)
  908 FORMAT('            DT = ', F8.6)

  CALL CALCULATE_TRAP() ! CALCULATE_TRAP() initializes the harmonic oscillator potential V.
  CALL COEF() ! COEF() defines the coefficients of Crank-Nicholson Scheme.

  CP2 = CP * CONJG(CP)
  CALL RAD(CP2, RMS) ! RAD() calculates the r.m.s radius RMS
  WRITE (7, 1001)
  WRITE (7, 1002)
  WRITE (7, 1001)

  WRITE (4, 1001)
  WRITE (4, 1012)
  WRITE (4, 1001)
  WRITE (4, 1013) RMS(1:3)
  1001 FORMAT (12X,'-----------------------------------------------------')
  1002 FORMAT (14X,'Iter',5x, 'Norm', 7X, 'Chem', 7X, 'Ener/N', 6X, '<rho>', 3X) 
  1003 FORMAT ('Initial: ', 6X, F13.4, 2F12.5, 2F11.5, F11.4)
  1012 FORMAT ('RMS size:', 4x, 'Iter', 8x, '<x>', 13X, '<y>', 13X, '<rho>') 
  1013 FORMAT ('Initial:',9X, F14.5, 2F16.5)

! IF (OPTION_re_im.EQ.2) THEN
!   OPEN(21, FILE = 're2d-initial-den.txt')
! ELSE
!   OPEN  (21, FILE = 'im2d-initial-den.txt')
! END IF

! CALL WRITE_DENSITY(21, CP2)
! CLOSE(21)

  998 FORMAT(2E20.8)
  999 FORMAT (2F12.6, F16.8)
  IF (NSTP /= 0) THEN
    GSTP =  XOP * G_2D / DFLOAT(NSTP) ! G is divided by Sqrt[2 pi] * LZ
    G = 0.0D0

    CALL NORM(CP, ZNORM) ! NORM() calculates norm and restores normalization
    CALL CHEM(CP, MU, EN) ! CHEM() calculates the chemical potential MU and energy EN.
    CP2 = CP * CONJG(CP)
    CALL RAD(CP2, RMS) ! RAD() calculates the r.m.s radius RMS
    WRITE (7, 1003) ZNORM, REAL(MU)/XOP, REAL(EN)/XOP, RMS(3)

    DO K = 1, NSTP  ! NSTP iterations to introduce the nonlinearity
      G = G + GSTP
      CALL CALCNU() ! CALCNU() performs time propagation with non-derivative parts. 
      CALL LUX()    ! LUX() performs the time iteration with space (x) derivative alone.
      CALL LUY()    ! LUY() performs the time iteration with space (y) derivative alone.
      IF (OPTION_re_im.EQ.1) CALL NORM(CP, ZNORM)
    END DO
    IF (OPTION_re_im.EQ.2) CALL NORM(CP, ZNORM)
    CALL CHEM(CP, MU, EN)
    CP2 = CP * CONJG(CP)
    CALL RAD(CP2, RMS)
    WRITE (7, 1005) ZNORM, REAL(MU)/XOP, REAL(EN)/XOP, RMS(3)
    WRITE (4, 1015) RMS(1:3)
    1005 FORMAT('NSTP iter.:',9x, F8.4, 2F12.5, 2F11.5, F11.4)
    1015 FORMAT('NSTP iter.:'6x, F14.5, 2F16.5)

  ! OPEN(22, FILE = 'den-nstp.txt')
  ! CALL WRITE_DENSITY(22, CP2)
  ! CLOSE(22)

  ! OPEN(11, FILE = 'den-nstp1d_x.txt')
  ! OPEN(12, FILE = 'den-nstp1d_y.txt')
  ! CALL DENSITY_1DX(11, CP2, DEN1X)
  ! CALL DENSITY_1DY(12, CP2, DEN1Y) 
  ! CLOSE(11)
  ! CLOSE(12)
  ELSE
    G = XOP * G_2D
  END IF 

  STEP = 1
  DO K = 1, NPAS ! NPAS iterations transient
    CALL CALCNU()
    CALL LUX()
    CALL LUY()
    IF (OPTION_re_im.EQ.1) CALL NORM(CP, ZNORM)

    IF (MOD(K, ITER) == 0) THEN
      IF (OPTION_RE_IM.EQ.1) WRITE(ITER_FILENAME,'(A,I0,A)') 'im2d-den-', STEP, '.txt'
      IF (OPTION_RE_IM.EQ.2) WRITE(ITER_FILENAME,'(A,I0,A)') 're2d-den-', STEP, '.txt'

      OPEN(101, FILE = ITER_FILENAME)
      CP2 = CP * CONJG(CP)
      CALL WRITE_DENSITY(101, CP2)   ! WRITES DENSITY FOR EVERY ITER
      CLOSE(101)

      CALL CHEM(CP, MU, EN)
      CALL RAD(CP2,  RMS)
      IF (OPTION_re_im.EQ.2) CALL NORM(CP, ZNORM) 
      WRITE (7, 1007)  STEP,ZNORM, REAL(MU)/XOP, REAL(EN)/XOP, RMS(3)
      WRITE (4, 1017) STEP, RMS(1:3)
      STEP = STEP + 1

      IF (OPTION_re_im.EQ.1) OPEN(13, FILE = 'im2d-fin-wf.txt')
      IF (OPTION_re_im.EQ.2) OPEN(13, FILE = 're2d-fin-wf.txt')

      CALL WRITE_WAVE_FUNCTION(13, CP)
      CLOSE(13)
    END IF
  END DO

  1007 FORMAT('NPAS iter.:',I6,3x,F8.4, 2F12.5, 2F11.5, F11.4) 
  1017 FORMAT('NPAS iter.:',I6, F14.5, 2F16.5)

! OPEN(25, FILE = 'den-npas.txt')
! CALL WRITE_DENSITY(25, CP2)
! CLOSE(25)

! OPEN(17, FILE = 'den-npas1d_x.txt')
! OPEN(18, FILE = 'den-npas1d_y.txt')
! CALL DENSITY_1DX(17, CP2, DEN1X)
! CALL DENSITY_1DY(18, CP2, DEN1Y)
! CLOSE(17)
! CLOSE(18) 

  IF (NRUN /= 0) THEN
    DO K = 1, NRUN ! NRUN iterations to check convergence
      CALL CALCNU()
      CALL LUX()
      CALL LUY() 
      IF (OPTION_re_im.EQ.1) CALL NORM(CP, ZNORM)
    END DO

    CALL CHEM(CP, MU, EN)
    CP2 = CP * CONJG(CP)
    CALL RAD(CP2,  RMS)
    IF (OPTION_re_im.EQ.2)  CALL NORM(CP, ZNORM)
    WRITE (7, 1006) ZNORM, REAL(MU)/XOP, REAL(EN)/XOP, RMS(3)
    WRITE (4, 1019) RMS(1:3)
    1006 FORMAT('NRUN iter.:',9x,F8.4, 2F12.5, 2F11.5,   F11.4)
    1019 FORMAT('NRUN iter.:'6x, F14.5, 2F16.5)
  END IF

!  IF (OPTION_re_im.EQ.1) OPEN(13, FILE = 'im2d-fin-den2d.txt')
!  IF (OPTION_re_im.EQ.2) OPEN(13, FILE = 're2d-fin-den2d.txt')
!  CALL WRITE_DENSITY(13, CP2)
!  CLOSE(13)

  IF (OPTION_re_im.EQ.1) OPEN(13, FILE = 'im2d-fin-wf.txt')
  IF (OPTION_re_im.EQ.2) OPEN(13, FILE = 're2d-fin-wf.txt')
  
  CALL WRITE_WAVE_FUNCTION(13, CP)
  CLOSE(13)

  IF (OPTION_re_im.EQ.1) THEN
    OPEN(15, FILE = 'im2d-den1d_x.txt')
    OPEN(16, FILE = 'im2d-den1d_y.txt')
  ELSE
    OPEN(15, FILE = 're2d-den1d_x.txt')
    OPEN(16, FILE = 're2d-den1d_y.txt')
  END IF

  CALL DENSITY_1DX(15, CP2, DEN1X)
  CALL DENSITY_1DY(16, CP2, DEN1Y)

  CLOSE(15)
  CLOSE(16) 

  CALL FREE_VARIABLES()
  IF (ALLOCATED(CP2)) DEALLOCATE(CP2)

  CALL SYSTEM_CLOCK (CLCK_COUNTS_END, CLCK_RATE,COUNT_MAX)
  CALL CPU_TIME(T2)
  END_D = DBLE(CLCK_COUNTS_END)
  START_D = DBLE(CLCK_COUNTS_BEG)
  IF (END_D < START_D) END_D = END_D + DBLE(COUNT_MAX)
  WRITE (7, 1001)
  WRITE (4, 1001)
  CLOSE (4)
  WRITE (7,*)
  WRITE (7,'(A,I7,A)') ' Clock Time: ', INT(( END_D  - START_D )/DBLE(CLCK_RATE)), ' seconds'
  WRITE (7,'(A,I7,A)') '   CPU Time: ', INT(T2-T1), ' seconds'  
  CLOSE (7)
END PROGRAM GROSS_PITAEVSKII_SSCN_2D

SUBROUTINE ALLOCATE_VARIABLES()
  USE COMM_DATA, ONLY : NX, NY, NO_OF_THREADS
  USE GPE_DATA, ONLY : X, Y, X2, Y2, V, CP, P,P2
  USE CN_DATA, ONLY : CALA, CGAA, CALB, CGAB!, CBE
  IMPLICIT NONE

  ALLOCATE(X(0:NX))
  ALLOCATE(X2(0:NX))
  ALLOCATE(Y(0:NY))
  ALLOCATE(Y2(0:NY))
  ALLOCATE(P(0:NX, 0:NY))
  ALLOCATE(P2(0:NX, 0:NY))
  ALLOCATE(V(0:NX, 0:NY))
  ALLOCATE(CP(0:NX, 0:NY))
  ALLOCATE(CALA(0:NX,0:NY))
  ALLOCATE(CGAA(0:NX,0:NY))
  ALLOCATE(CALB(0:NX,0:NY))
  ALLOCATE(CGAB(0:NX,0:NY))
END SUBROUTINE ALLOCATE_VARIABLES

SUBROUTINE FREE_VARIABLES()
  USE GPE_DATA, ONLY : X, Y, X2, Y2, V, CP,  P,P2 
  USE CN_DATA, ONLY : CALA, CGAA, CALB, CGAB!, CBE
  IMPLICIT NONE

  IF (ALLOCATED(X)) DEALLOCATE(X)
  IF (ALLOCATED(X2)) DEALLOCATE(X2)
  IF (ALLOCATED(Y)) DEALLOCATE(Y)
  IF (ALLOCATED(Y2)) DEALLOCATE(Y2)
  IF (ALLOCATED(V)) DEALLOCATE(V)
  IF (ALLOCATED(CP)) DEALLOCATE(CP)
  IF (ALLOCATED(P)) DEALLOCATE(P)
  IF (ALLOCATED(P2)) DEALLOCATE(P2)
  IF (ALLOCATED(CALA)) DEALLOCATE(CALA)
  IF (ALLOCATED(CGAA)) DEALLOCATE(CGAA)
  IF (ALLOCATED(CALB)) DEALLOCATE(CALB)
  IF (ALLOCATED(CGAB)) DEALLOCATE(CGAB)
END SUBROUTINE FREE_VARIABLES

SUBROUTINE INITIALIZE(CP)
  ! Routine that initializes the constant and variables.
  ! Calculates the potential term V and the initial wave function CP
  USE COMM_DATA, ONLY : NX, NY, NX2, NY2, PI ,NSTP
  USE GPE_DATA, ONLY : GAMMA, NU, DX, DY, X, X2, Y, Y2,CI,OPTION_RE_IM,RANDOM,FUNCTION
  IMPLICIT NONE

  COMPLEX (8), DIMENSION(0:NX, 0:NY), INTENT(INOUT) :: CP
  REAL (8) :: QR_AL, TMP, TX, TY 
  INTEGER :: I, J, K
   INTEGER, DIMENSION (100) :: seed   
  real (8), dimension((1+NX)*(1+NY)) :: RANDNUM
 
  QR_AL = SQRT(SQRT(GAMMA * NU)/PI)
  FORALL (I=0:NX) X(I) = (I-NX2)*DX
  FORALL (J=0:NY) Y(J) = (J-NY2)*DY
  X2 = X*X
  Y2 = Y*Y


   seed = 13
  IF(RANDOM.EQ.1)  CALL RANDOM_SEED(PUT=seed)
  
    IF(RANDOM.EQ.1)   call random_number(RANDNUM)
 


  IF (NSTP == 0) THEN
    IF (OPTION_RE_IM.EQ.1) WRITE(*,'(a)') "Run the program using the input file to read. e.g.: ./imag2d < im2d-fin-wf.txt"
    IF (OPTION_RE_IM.EQ.2) WRITE(*,'(a)') "Run the program using the input file to read. e.g.: ./real2d < im2d-fin-wf.txt"

    DO J = 0,NY
      DO I = 0,NX
        READ (*, 998) TX, TY
        CP(I,J) = cmplx(TX,TY)
      END DO
      READ(*,*)
    END DO
  ELSE

 
    !$OMP PARALLEL DO PRIVATE(I,J,TMP)
    DO J = 0, NY; DO I = 0, NX

     K=(I+1)*(J+1)
 

      TMP = (GAMMA * X2(I) +  NU * Y2(J))/2.0D0
 IF(FUNCTION.EQ.1)      CP(I,J) = QR_AL * EXP(-TMP) * (X(I)+CI*Y(J))   * EXP (2.d0*PI*CI*RANDNUM(K))
 IF(FUNCTION.EQ.0)      CP(I,J) = QR_AL * EXP(-TMP)    * EXP (2.d0*PI*CI*RANDNUM(K))

    END DO; END DO
    !$OMP END PARALLEL DO 

  END IF
  998 FORMAT(ES16.9, 1X, ES16.9)
END SUBROUTINE INITIALIZE

SUBROUTINE CALCULATE_TRAP()
  ! Calculates the harmonic oscillator potential term V.
  USE COMM_DATA, ONLY : NX, NY
  USE GPE_DATA, ONLY : XOP, V, X2, Y2, GAMMA, NU
  IMPLICIT NONE
  INTEGER :: I, J
  REAL (8) :: GAMMA2, NU2
  GAMMA2 = GAMMA * GAMMA
  NU2 = NU * NU
  !$OMP PARALLEL DO PRIVATE(I,J)
  DO I = 0, NX; DO J = 0, NY
  V(I,J) = XOP * (GAMMA2 * X2(I) + NU2 * Y2(J)) / 2.0D0
  END DO; END DO
 !$OMP END PARALLEL DO
END SUBROUTINE CALCULATE_TRAP

SUBROUTINE COEF()
  ! Calculates the coefficients needed in subroutine LUX and LUY.
  USE COMM_DATA, ONLY : NXX, NYY
  USE GPE_DATA, ONLY : DX, DY, DT,  X, Y, XOP, OMEGA,OPTION_re_im,CI,CIJ
  USE CN_DATA
  IMPLICIT NONE
  INTEGER :: I, J
  REAL (8) :: DX2, DY2
  REAL (8) :: DXX, DYY
  COMPLEX (8) :: CDT

  DX2 = DX*DX   ! Generating the Coefficients for the
  DY2 = DY*DY   ! C-N method (to solve the spatial part)
  DXX = 1.0D0/DX2
  DYY = 1.0D0/DY2
  CDT = CIJ*DT
  CA0 = 1.0D0 + CDT*DXX
  CA0R = 1.0D0 - CDT*DXX
  CB0 = 1.0D0 + CDT*DYY
  CB0R = 1.0D0 - CDT*DYY

  C0 =CIJ* CI*DT*XOP*OMEGA/2.0D0

  CTMPX = C0/(2.0D0*DX)
  CT0 = cDT*DXX/2.0D0
  CAM = -(CT0 - CTMPX*Y)
  CAP = -(CT0 + CTMPX*Y)

  DO J = 0, NY
    CALA(NXX, J) = 0.0D0
    CGAA(NXX, J) = -1.0D0/CA0
    DO I = NXX, 1, -1
      CALA(I-1,J) = CAM(J)*CGAA(I,J)
      CGAA(I-1,J) = -1.0D0/(CA0+CAP(J)*CALA(I-1,J))
    END DO
  END DO

  CTMPY = C0/(2.0D0*DY)
  CT0 = cDT*DYY/2.0D0
  CBM = -(CT0 + CTMPY*X)
  CBP = -(CT0 - CTMPY*X)

  DO I = 0, NX
    CALB(I, NYY) = 0.0D0
    CGAB(I, NYY) = -1.0D0/CB0
    DO J = NYY, 1, -1
      CALB(I,J-1) = CBM(I)*CGAB(I,J)
      CGAB(I,J-1) = -1.0D0/(CB0+CBP(I)*CALB(I,J-1))
    END DO
  END DO
END SUBROUTINE COEF

SUBROUTINE CALCNU() ! Exact solution
  ! Solves the partial differential equation with the potential and
  ! the nonlinear term.
  USE COMM_DATA, ONLY : NX, NY
  USE GPE_DATA, ONLY : OPTION,  DT, V, G, CP, P2,  OPTION_re_im,CIJ
  IMPLICIT NONE
  REAL (8), DIMENSION(0:NX,0:NY) ::  TMP
  INTEGER :: I, J

  !$OMP PARALLEL DO PRIVATE(I, J,TMP)
  DO J = 0, NY; DO I = 0, NX
    P2(I,J)  = CP(I,J) * CONJG(CP(I,J))
    TMP(I,J) = V(I,J) + G * P2(I,J)
    CP(I,J)  = CP(I,J) * EXP(-CIJ*DT*TMP(I,J))
  END DO; END DO
  !$OMP END PARALLEL DO
END SUBROUTINE CALCNU

SUBROUTINE LUX()
  ! Solves the partial differential equation only with the X-space
  ! derivative term using the Crank-Nicholson method
  USE OMP_LIB
  USE COMM_DATA, ONLY : NX, NY, NXX
  USE GPE_DATA, ONLY : CP, option_re_im
  USE CN_DATA, ONLY : CAP, CAM, CA0, CA0R, CALA, CGAA
  IMPLICIT NONE
  COMPLEX (8), DIMENSION (0:NX,0:NY) :: CBE
  COMPLEX (8) :: CXX
  INTEGER :: I, J

  !$OMP PARALLEL DO PRIVATE(I, J, CXX)
  DO J = 0, NY
    IF (option_re_im.eq.1) CBE(NXX,J) = 0.d0
    IF (option_re_im.eq.2) CBE(NXX,J) = CP(NX,J)

    DO I = NXX, 1, -1
      CXX = -CAP(J) * CP(I+1,J) + CA0R * CP(I,J) - CAM(J) * CP(I-1,J)
      CBE(I-1,J) = CGAA(I,J) * (CAP(J) * CBE(I,J) - CXX)
    END DO
    CP(0,J) = 0.0D0
    DO I = 0, NXX
      CP(I+1,J) = CALA(I,J)*CP(I,J)+CBE(I,J)
    END DO
    CP(NX,J) = 0.0D0
  END DO
  !$OMP END PARALLEL DO
END SUBROUTINE LUX

SUBROUTINE LUY()
  ! Solves the partial differential equation only with the Y-space
  ! derivative term using the Crank-Nicholson method
  USE OMP_LIB
  USE COMM_DATA, ONLY : NX, NY, NYY
  USE CN_DATA, ONLY : CBP, CBM, CB0, CB0R, CALB, CGAB
  USE GPE_DATA, ONLY : CP, option_re_im
  IMPLICIT NONE
  COMPLEX (8), DIMENSION(0:NX,0:NY) :: CBE
  COMPLEX (8) :: CYY
  INTEGER :: I, J

  !$OMP PARALLEL DO PRIVATE(I, J, CYY)
  DO I = 0, NX
    IF (option_re_im.eq.1) CBE(I,NYY) = 0.d0
    IF (option_re_im.eq.2) CBE(I,NYY) =   CP(I,NY)

    DO J = NYY, 1, -1
      CYY = -CBP(I) * CP(I,J+1) + CB0R * CP(I,J) - CBM(I) * CP(I,J-1)
      CBE(I,J-1) = CGAB(I,J) * (CBP(I) * CBE(I,J) - CYY)
    END DO
    CP(I,0) = 0.0D0
    DO J = 0, NYY
      CP(I,J+1) = CALB(I,J) * CP(I,J) + CBE(I,J)
    END DO
    CP(I,NY) = 0.0D0
  END DO
  !$OMP END PARALLEL DO
END SUBROUTINE LUY

SUBROUTINE NORM(CP, ZNORM)
  ! Calculates the normalization of the wave function and sets it to unity.
  USE OMP_LIB
  USE GPE_DATA, ONLY : DX, DY, option_re_im
  USE COMM_DATA, ONLY : NX, NY
  IMPLICIT NONE
  COMPLEX (8), DIMENSION(0:NX, 0:NY), INTENT(INOUT) :: CP
  REAL (8), INTENT(OUT) :: ZNORM

  INTERFACE
    FUNCTION INTEGRATE(P2, DX, DY) RESULT(RES)
      USE COMM_DATA, ONLY : NX, NY
      IMPLICIT NONE
      REAL (8), DIMENSION(0:NX, 0:NY), INTENT(IN) :: P2
      REAL (8), INTENT(IN) :: DX, DY
      REAL (8) :: RES
    END FUNCTION INTEGRATE
  END INTERFACE

  INTEGER :: I,J
  REAL (8), DIMENSION(0:NX,0:NY) :: TMP2D

  !$OMP PARALLEL DO PRIVATE(I, J)
  DO J = 0, NY; DO I = 0, NX
    TMP2D(I,J) = CP(I,J) * CONJG(CP(I,J))
  END DO; END DO
  !$OMP END PARALLEL DO

  ZNORM = SQRT(INTEGRATE(TMP2D, DX, DY))

  IF (option_re_im.EQ.1) THEN
    !$OMP PARALLEL DO PRIVATE(I, J)
    DO J=0,NY; DO I=0,NX
      CP(I,J)=CP(I,J)/ZNORM
    END DO; END DO
    !$OMP END PARALLEL DO
  END IF
END SUBROUTINE NORM

SUBROUTINE NORMCHEM(CPR, ZNORM)
  ! Calculates the normalization of the wave function and sets it to unity.
  USE OMP_LIB
  USE GPE_DATA, ONLY : DX, DY
  USE COMM_DATA, ONLY : NX, NY
  IMPLICIT NONE
  REAL (8), DIMENSION(0:NX, 0:NY), INTENT(INOUT) :: CPR
  REAL (8), INTENT(OUT) :: ZNORM

  INTERFACE
    FUNCTION INTEGRATE(P2, DX, DY) RESULT(RES)
      USE COMM_DATA, ONLY : NX, NY
      IMPLICIT NONE
      REAL (8), DIMENSION(0:NX, 0:NY), INTENT(IN) :: P2
      REAL (8), INTENT(IN) :: DX, DY
      REAL (8) :: RES
    END FUNCTION INTEGRATE
  END INTERFACE

   INTEGER :: I,J
   REAL (8), DIMENSION(0:NX,0:NY) :: TMP2D

  !$OMP PARALLEL DO PRIVATE(I, J)
  DO J = 0, NY; DO I = 0, NX
    TMP2D(I,J) = CPR(I,J) * CPR(I,J)
  END DO; END DO
  !$OMP END PARALLEL DO

  ZNORM = (INTEGRATE(TMP2D, DX, DY))
END SUBROUTINE NORMCHEM

SUBROUTINE RAD(P, RMS)
  ! Calculates the root mean square size RMS
  USE COMM_DATA, ONLY : NX, NY
  USE GPE_DATA, ONLY : DX, DY,  X2, Y2
  IMPLICIT NONE
  REAL (8), DIMENSION(0:NX, 0:NY), INTENT(IN) :: P
  REAL (8), DIMENSION(:), INTENT(OUT) :: RMS

  INTERFACE
    FUNCTION INTEGRATE(P2, DX, DY) RESULT(RES)
      USE COMM_DATA, ONLY : NX, NY
      IMPLICIT NONE
      REAL (8), DIMENSION(0:NX, 0:NY), INTENT(IN) :: P2
      REAL (8), INTENT(IN) :: DX, DY
      REAL (8) :: RES
    END FUNCTION INTEGRATE
  END INTERFACE

  INTEGER :: I,J
  REAL (8), DIMENSION(0:NX,0:NY) :: TMP2D

  FORALL(J=0:NY) TMP2D(:,J) = X2*P(:,J)
  RMS(1) = SQRT(INTEGRATE(TMP2D, DX, DY))

  FORALL(I=0:NX) TMP2D(I,:) = Y2*P(I,:)
  RMS(2) = SQRT(INTEGRATE(TMP2D, DX, DY))

  RMS(3)=  SQRT(RMS(1)**2+RMS(2)**2)
END SUBROUTINE RAD

SUBROUTINE CHEM(CP, MU, EN)
  ! Calculates the chemical potential MU and energy EN.  CP is the wave
  ! function, V is the potential and G is the nonlinearity.
  USE COMM_DATA, ONLY : NX, NY
  USE GPE_DATA, ONLY : DX, DY, V, G, OMEGA, X, Y,CI,XOP
  IMPLICIT NONE
  COMPLEX (8), DIMENSION(0:NX, 0:NY), INTENT(IN) :: CP
  REAL (8), INTENT(OUT) :: MU, EN

  INTERFACE
    PURE FUNCTION DIFF(P, DX) RESULT (DP)
      IMPLICIT NONE
      REAL (8), DIMENSION(0:), INTENT(IN) :: P
      REAL (8), INTENT(IN) :: DX
      REAL (8), DIMENSION(0:SIZE(P)-1) :: DP
    END FUNCTION DIFF

    FUNCTION INTEGRATE(P2, DX, DY) RESULT(RES)
      USE COMM_DATA, ONLY : NX, NY
      IMPLICIT NONE
      REAL (8), DIMENSION(0:NX, 0:NY), INTENT(IN) :: P2
      REAL (8), INTENT(IN) :: DX, DY
      REAL (8) :: RES
    END FUNCTION INTEGRATE

    SUBROUTINE NORMCHEM(CPR, ZNORM)
      USE COMM_DATA, ONLY : NX, NY
      IMPLICIT NONE
      REAL (8), DIMENSION(0:NX, 0:NY), INTENT(INOUT) :: CPR
      REAL (8), INTENT(OUT) :: ZNORM
    END  SUBROUTINE NORMCHEM
  END INTERFACE

  INTEGER :: I, J
  REAL (8), DIMENSION(0:NX, 0:NY) :: P2, GP2, DP2,P2R
  REAL (8), DIMENSION(0:NX,0:NY) ::  TMP2D, EMP2D, DPLZ, CPR,CPI, DPLY,DPLY2,DPLZ2
  REAL (8), DIMENSION(0:NX, 0:NY) :: DPXR, DPYR, DPXI, DPYI
  REAL (8) :: ZNORM

  MU = 0.0D0
  EN = 0.0D0

  CPR=CP
  CPI=AIMAG(CP)

  P2 = CP*CONJG(CP)

  DO I = 0, NX
    DPYR(I,0:NY) = DIFF(CPR(I,0:NY), DY)
    DPYI(I,0:NY) = DIFF(CPI(I,0:NY), DY)
  END DO

  DO J = 0, NY
    DPXR(0:NX,J) = DIFF(CPR(0:NX,J), DX)
    DPXI(0:NX,J) = DIFF(CPI(0:NX,J), DX)
  END DO

  DP2 = DPXR*DPXR + DPYR*DPYR

  DO J = 0, NY; DO I = 0, NX
    DPLZ(I,J) =  CPR(I,J)*( X(I)*DPYI(I,J) - Y(J)*DPXI(I,J) )
    DPLY(I,J) =  CPI(I,J)*( X(I)*DPYR(I,J) - Y(J)*DPXR(I,J) )
    DPLZ2(I,J) =  CPR(I,J)*( X(I)*DPYR(I,J) - Y(J)*DPXR(I,J) )
    DPLY2(I,J) =  CPI(I,J)*( X(I)*DPYI(I,J) - Y(J)*DPXI(I,J) )
  END DO; END DO

  P2R=CPR*CPR
  GP2 = G*P2
  TMP2D = (V + GP2)*P2R + DP2 - XOP*(OMEGA*DPLZ)
  EMP2D = (V + GP2/2.0D0)*P2R + DP2 - XOP*(OMEGA*DPLZ)

  CALL NORMCHEM(CPR,ZNORM)
  MU = INTEGRATE(TMP2D, DX, DY)/ZNORM
  EN = INTEGRATE(EMP2D, DX, DY)/ZNORM
END SUBROUTINE CHEM

FUNCTION INTEGRATE(U, DX, DY) RESULT(RES)
  USE COMM_DATA, ONLY : NX, NY
  IMPLICIT NONE
  REAL (8), DIMENSION(0:NX, 0:NY), INTENT(IN) :: U
  REAL (8), INTENT (IN) :: DX, DY
  REAL (8) :: RES

  INTERFACE
    PURE FUNCTION SIMP(F, DX)
      IMPLICIT NONE
      REAL (8), DIMENSION(0:), INTENT(IN) :: F
      REAL (8), INTENT(IN) :: DX
      REAL (8) :: SIMP
    END FUNCTION SIMP
  END INTERFACE

  REAL (8), DIMENSION(0:NX) :: TMP1D
  INTEGER :: I

  !$OMP PARALLEL DO PRIVATE(I)
  DO I=0,NX
    TMP1D(I) = SIMP(U(I,0:), DY)
  END DO
  !$OMP END PARALLEL DO

  RES = SIMP(TMP1D, DX)
END FUNCTION INTEGRATE

PURE FUNCTION SIMP(F, DX)
  ! Does the spatial integration with Simpson's rule.
  ! N refer to the number of integration points, DX space step, and
  ! F is the function to be integrated.
  IMPLICIT NONE
  REAL (8), DIMENSION(0:), INTENT(IN) :: F
  REAL (8), INTENT(IN) :: DX
  REAL (8) :: SIMP
  REAL (8) :: F1, F2
  INTEGER :: I, N

  N = SIZE(F) - 1
  F1 = F(1) + F(N-1) ! N EVEN
  F2 = F(2)
  DO I = 3, N-3, 2
    F1 = F1 + F(I)
    F2 = F2 + F(I+1)
  END DO
  SIMP = DX*(F(0) + 4.D0*F1 + 2.D0*F2 + F(N))/3.D0
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

  N = SIZE(P) - 1
  DP(0) = 0.D0
  DP(1) = (P(2) - P(0))/(2.D0*DX)
  FORALL(I=2:N-2)
    DP(I) = (P(I-2)-8.D0*P(I-1)+8.D0*P(I+1)-P(I+2))/(12.D0*DX)
  END FORALL
  DP(N-1) = (P(N) - P(N-2))/(2.D0*DX)
  DP(N) = 0.D0
END FUNCTION DIFF

SUBROUTINE WRITE_WAVE_FUNCTION(FUNIT, CP)
  ! Writes the value of wave function (poth real and imaginary parts)
  ! at every spatial point.
  USE COMM_DATA, ONLY : NX, NY
  USE GPE_DATA, ONLY : X, Y
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: FUNIT
  COMPLEX (8), DIMENSION(0:NX, 0:NY), INTENT(IN) :: CP

  INTEGER :: I, J

  DO J = 0, NY
    DO I = 0, NX
      WRITE (FUNIT, 998) REAL(CP(I,J)), AIMAG(CP(I,J))
    END DO
    WRITE(FUNIT,*)
  END DO
  998 FORMAT(ES16.9, 1X, ES16.9)
END SUBROUTINE WRITE_WAVE_FUNCTION

SUBROUTINE WRITE_DENSITY(FUNIT, U2)
  ! Writes the 2D density. The format is (Y,X,den(X,Y)).
  USE COMM_DATA, ONLY : NX, NY
  USE GPE_DATA, ONLY : X, Y
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: FUNIT
  REAL (8), DIMENSION(0:NX, 0:NY), INTENT(IN) :: U2

  INTEGER :: I, J

  DO J = 0, NY
    DO I = 0, NX
      WRITE(FUNIT, 999) Y(J), X(I), U2(I,J)
    END DO
    WRITE(FUNIT, *)
  END DO
  999 FORMAT(ES13.6, 1X, ES13.6, 1X, ES13.6)
END SUBROUTINE WRITE_DENSITY

SUBROUTINE DENSITY_1DX(FUNIT, U2, DE1DX)
  ! Writes the integrated 1D density over X dimension.
  USE COMM_DATA, ONLY : NX, NY
  USE GPE_DATA, ONLY : X, DY
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: FUNIT
  REAL (8), DIMENSION(0:NX, 0:NY), INTENT(IN) :: U2
  REAL (8), DIMENSION(0:NX), INTENT(OUT) :: DE1DX

  INTERFACE
    PURE FUNCTION SIMP(F, DX) RESULT (VALUE)
      IMPLICIT NONE
      REAL (8), DIMENSION(0:), INTENT(IN) :: F
      REAL (8), INTENT(IN) :: DX
      REAL (8) :: VALUE
    END FUNCTION SIMP
  END INTERFACE

  INTEGER :: I

  DO I = 0, NX
	  DE1DX(I) = SIMP(U2(I,:), DY)
    WRITE(FUNIT, 1001) X(I), DE1DX(I)
  END DO
  1001 FORMAT(ES13.6, 1X, ES13.6)
END SUBROUTINE DENSITY_1DX

SUBROUTINE DENSITY_1DY(FUNIT, U2, DE1DY)
  ! Writes the integrated 1D density over Y dimension.
  USE COMM_DATA, ONLY : NX, NY
  USE GPE_DATA, ONLY : Y, DX
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: FUNIT
  REAL (8), DIMENSION(0:NX, 0:NY), INTENT(IN) :: U2
  REAL (8), DIMENSION(0:NX), INTENT(OUT) :: DE1DY

  INTERFACE
    PURE FUNCTION SIMP(F, DX) RESULT (VALUE)
      IMPLICIT NONE
      REAL (8), DIMENSION(0:), INTENT(IN) :: F
      REAL (8), INTENT(IN) :: DX
      REAL (8) :: VALUE
    END FUNCTION SIMP
  END INTERFACE

  INTEGER :: J

  DO J = 0, NY
	  DE1DY(J) = SIMP(U2(:,J), DX)
    WRITE(FUNIT, 1001) Y(J), DE1DY(J)
  END DO
  1001 FORMAT(ES13.6, 1X, ES13.6)
END SUBROUTINE DENSITY_1DY
