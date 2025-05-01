! spin-SO-rot-imre2d-omp.f90
! Part of the BEC-GP-SPINOR-ROT-OMP software package
! Program for imaginary- and real-time propagation of rotating spin-orbit-coupled spin-1 Bose-Einstein condensates
! Numerical solution to the Gross-Pitaevskii equation in a two-dimensional harmonic trap
! (Fortran 90/95/2003 version)
! Developed by P. Muruganandam, A. Balaz, and S. K. Adhikari (21 January 2021)
!
!# To compile:
!# (1) Intel Fortran Compiler
! ifort -c -o util.o util.f90
! ifort -O3 -openmp -w -mcmodel medium -shared-intel spin-SO-rot-imre2d-omp.f90 util.o -o imre2d
!
!# (2) GNU Fortran (gfortran)
! gfortran -c -o util.o util.f90
! gfortran -O3 -fopenmp -w spin-SO-rot-imre1d-omp.f90 util.o -o imre2d
!
!# (3) PGI Fortran (pgfortran)
! pgfortran -c -o util.o util.f90
! pgfortran -O3 -fast -mp=allcores spin-SO-rot-imre2d-omp.f90 util.o -o imre2d
!
!# (4) Oracle/Sun Fortran (sunf90)
! sunf90 -c -o util.o util.f90
! sunf90 -fopenmp -fast spin-SO-rot-imre2d-omp.f90 util.o -o imre2d
!

! Begin selection of input parameters 
MODULE COMM_DATA
!*********************************************************************
! SELECT # OF SPACE POINTS NX, NY AND # OF TIME ITERATIONS NSTP, NPAS & NRUN
! USE NSTP = 0 TO READ FILE im-wave-fun-fin.txt
  INTEGER, PARAMETER :: NX = 160, NXX = NX-1, NX2 = NX/2
  INTEGER, PARAMETER :: NY = 160, NYY = NY-1, NY2 = NY/2
  INTEGER, PARAMETER :: NSTP = 10, NPAS = 50000, NRUN = 50000

!  Number of threads, the maximum possible number is the number of available cores
! NTHREADS = 0  takes all available threads 
  INTEGER, PARAMETER  :: NTHREADS = 0
  REAL (8), PARAMETER :: Pi = 3.14159265358979D0
END MODULE COMM_DATA 
!*********************************************************************
MODULE SPIN_PARS
  USE COMM_DATA, ONLY : PI
!---------------------------------------------------------------------
!******* SELECT  POLAR (C_2 > 0) OR FERROMAGNETIC (C_2 < 0) BEC ******
   REAL (8), PARAMETER :: C_0 = 482.D0, C_2 = 15.00D0
!   REAL (8), PARAMETER :: C_0 = 669.D0, C_2 = -3.100D0
 
END MODULE SPIN_PARS
!*********************************************************************
!
MODULE GPE_DATA
  USE COMM_DATA, ONLY : NX, NY, Pi
  USE SPIN_PARS, ONLY : C_0, C_2
!**********************************************************************
!*** SELECT OPTION FOR SO COUPLING AND STRENGTH GAMMA 
! INTEGER,PARAMETER :: OPT_SO = 0;REAL (8), PARAMETER :: GAM0 = .00D0  ! no SO coupling
  INTEGER,PARAMETER :: OPT_SO = 1;REAL (8), PARAMETER :: GAM0 = .500D0 !Sigma_x  p_y-Sigma_y p_x [Rashba] 
!  INTEGER,PARAMETER :: OPT_SO = 2;REAL (8), PARAMETER :: GAM0 = .500D0 !-Sigma_x p_y-Sigma_y p_x [Dresselhaus]
!********************************************************************
!*** SELECT RANDOM PHASE
   INTEGER, PARAMETER :: RANDOM=0 ! No random phase in the initial function
!  INTEGER, PARAMETER :: RANDOM=1 ! Random phase included in the initial function
!***************************************************************
!! SELECT SPACE STEP DX
  REAL (8), PARAMETER :: DX = 0.1D0, DY = DX
!! SELECT OPTION  FOR  PROPAGATION
   INTEGER, PARAMETER :: OPT_PROP = 1; REAL (8), PARAMETER :: DT = DX * DX * 0.1D0 ! IMAG TIME
!   INTEGER, PARAMETER :: OPT_PROP = 2; REAL (8), PARAMETER :: DT = DX * DX * 0.025D0 ! REAL TIME
!! SELECT PARAMETERS OF MODEL
  REAL (8), PARAMETER :: MAG_0 = 0.000000000D0, ACCUR = 1.0D-6
  REAL (8), PARAMETER :: OMEGA0 = .0000D0     !  RABI COUPLING
  REAL (8), PARAMETER :: OMEG  =  0.0000000d0   !   ROTATION
!----------------------------------------------------------------
!******  SELECT RING OR STRIPE STATES 
   REAL (8), PARAMETER :: OPT_ST =.73D0  
!****  RING-TYPE FOR GAM0 < OPT_ST, STRIPE FOR GAM0 > OPT_ST
!**********************************************************************
!! END OF PARAMETER SELECTION
!**********************************************************************
  REAL (8), PARAMETER :: LAMBDA = 1.0D0, LAMBDA2 = LAMBDA * LAMBDA

  REAL (8), PARAMETER :: OMEGA = OMEGA0/SQRT(2.0D0) - 1.0D-10, GAM = GAM0/SQRT(2.0D0)
  REAL (8), PARAMETER :: C0 = C_0, C2 = C_2
!
  COMPLEX (8) ::  CIJ
  COMPLEX (8), PARAMETER :: CI = (0.0D0, 1.0D0)
!
  REAL (8), DIMENSION(0:NX) :: X, X2
  REAL (8), DIMENSION(0:NY) :: Y, Y2
  REAL (8), DIMENSION(3,0:NX, 0:NY) :: V
  REAL (8), DIMENSION(0:NX, 0:NY) :: R2
  COMPLEX (8), DIMENSION(3,0:NX, 0:NY) :: CP
  REAL (8), DIMENSION(3) :: D, ZNORM
  REAL (8), DIMENSION(4) :: RAD 
  REAL (8) :: C0MC2, C0PC2, MAG, T
END MODULE GPE_DATA

MODULE CN_DATA
  USE COMM_DATA, ONLY : NX, NY, NXX, NYY 
!  COMPLEX (8) :: CAM, CBM
 COMPLEX (8), DIMENSION(0:NX) :: CBP, CBM
  COMPLEX (8), DIMENSION(0:NY) :: CAP, CAM
  COMPLEX (8), DIMENSION(0:NX,0:NY) :: CGAA, CGAB,  CALA,CALB
!  COMPLEX (8), DIMENSION(:,:):: CALA, CALB, CGAA, CGAB
  COMPLEX (8), DIMENSION(3,0:NX,0:NY) :: CBE
 COMPLEX (8) ::  CT0, CTMPX, CTMPY
  COMPLEX (8):: CC0
  COMPLEX (8), DIMENSION(3) :: CXX
  COMPLEX (8), DIMENSION(3) :: CYY

  REAL (8) :: CA0R, CB0R
END MODULE CN_DATA
!
!*******************************************************************************
PROGRAM SPIN_ONE_GPE2D_IMAG
!-------------------------------------------------------------------------------
  USE COMM_DATA
  USE GPE_DATA, ONLY : NX, NY, DT, X, Y, CP, C0, C2,MAG_0, &
       GAM0, OPT_PROP, RANDOM, DX, DY, RAD , CIJ, OPT_SO,OMEGA0,OMEG,OPT_ST
  USE OMP_LIB
  IMPLICIT NONE
!
  INTERFACE
    SUBROUTINE INITALIZE()
      IMPLICIT NONE
    END  SUBROUTINE INITALIZE
  END INTERFACE

  INTERFACE
    SUBROUTINE HERM(CP, DT)
      IMPLICIT NONE
      COMPLEX (8), DIMENSION(:,0:,0:), INTENT(INOUT) :: CP
      REAL (8), INTENT(IN) :: DT
    END SUBROUTINE HERM
  END INTERFACE

  INTERFACE 
    SUBROUTINE CALCNU(CP, DT)
      IMPLICIT NONE
      COMPLEX (8), DIMENSION(:,0:,0:), INTENT(INOUT) :: CP
      REAL (8), INTENT(IN) :: DT
    END SUBROUTINE CALCNU
  END INTERFACE
!
  INTERFACE 
    SUBROUTINE LUX(CP)
      IMPLICIT NONE
      COMPLEX (8), DIMENSION(:,0:,0:), INTENT(INOUT) :: CP
    END SUBROUTINE LUX
  END INTERFACE
!
  INTERFACE 
    SUBROUTINE LUY(CP)
      IMPLICIT NONE
      COMPLEX (8), DIMENSION(:,0:,0:), INTENT(INOUT) :: CP
    END SUBROUTINE LUY
  END INTERFACE
 
  INTERFACE
    SUBROUTINE RADIUS(CP, RAD )
      IMPLICIT NONE
      COMPLEX (8), DIMENSION(:,0:,0:), INTENT(INOUT) :: CP
   REAL (8), DIMENSION(4), INTENT(OUT) :: RAD 
    END  SUBROUTINE RADIUS
  END INTERFACE
!
  INTERFACE
    SUBROUTINE RENORM(CP,MAGNET)
      IMPLICIT NONE
      COMPLEX (8), DIMENSION(:,0:,0:), INTENT(INOUT) :: CP
      REAL (8) :: MAGNET 
    END  SUBROUTINE RENORM
  END INTERFACE
!
  INTERFACE
    SUBROUTINE ENERGY(CP,  EN)
      IMPLICIT NONE
      COMPLEX (8), DIMENSION(:,0:,0:), INTENT(IN) :: CP
      REAL (8), INTENT(OUT) :: EN
    END SUBROUTINE ENERGY
  END INTERFACE
!
  INTERFACE
    SUBROUTINE SO(CP,DT)
      IMPLICIT NONE
      COMPLEX (8), DIMENSION(:,0:,0:), INTENT(INOUT) :: CP
      REAL (8), INTENT(IN) :: DT
    END SUBROUTINE SO
  END INTERFACE
!-------------------------------------------------------------------------------
  REAL (8 ), DIMENSION(3,0:NX,0:NY) :: P, P2
  REAL (8) ::  EN , MAGNET
  INTEGER :: I, J, K, NO_OF_THREADS 
!
  REAL (8) :: T1, T2
  INTEGER :: CLCK_COUNTS_BEG, CLCK_COUNTS_END, CLCK_RATE
!
  CALL SYSTEM_CLOCK ( CLCK_COUNTS_BEG, CLCK_RATE )
  CALL CPU_TIME (T1)
!
 
  IF (NTHREADS /= 0) THEN
  CALL OMP_SET_NUM_THREADS(NTHREADS)
   END IF
  

  !$OMP PARALLEL
    !$OMP MASTER
      NO_OF_THREADS = OMP_GET_NUM_THREADS()
    !$OMP END MASTER
  !$OMP END PARALLEL

 
  903 FORMAT('     Nonlinearity C_0 = ', F0.6, ', C_2 = ',F0.6, ', MAG = ',F10.6) 
  905 FORMAT('     # of space steps: NX =   ', I0, ',  NY = ', I0)
  906 FORMAT('     SPACE AND TIME STEPS DX =',F9.5, ',  DY =',F9.5, ', DT =', E12.4)
  907 FORMAT('     # of Time steps:  NSTP = ', I0,', NPAS = ', I0,', NRUN = ',I0)
  1005 FORMAT('NSTP iter.:',x, 3F11.4, 2F11.5) 
  1006 FORMAT('NPAS iter.:',x, 3F11.4, 2F11.5) 
  1007 FORMAT('NRUN iter.:',x, 3F11.4, 2F11.5) 
  1008 FORMAT('NSTP iter.:',x, 3F11.3, 2F11.4) 
  1009 FORMAT('NPAS iter.:',x, 3F11.3, 2F11.4) 
  1010 FORMAT('NRUN iter.:',x, 3F11.3, 2F11.4) 
  1003 FORMAT (12X,'-----------------------------------------------------')
  1004 FORMAT (18X,'RAD(1)',5x, 'RAD(2)', 4X, 'RAD(3)', 4X, 'Energy', 6X, 'MAG')
!  1014 FORMAT (19X,'RAD (1)',7x, 'RAD (2)', 7X, 'rad(3)', 8X, '<LEN>', 6X) 
!  WRITE(7,*)
  SELECT CASE (OPT_PROP)
    CASE(1)
      CIJ = CMPLX(1.0D0, 0.0D0)
      OPEN(7, FILE = 'im-out.txt')
      WRITE(7,*)'    IMAGINARY-TIME PROPAGATION'
    CASE(2)
      CIJ = CMPLX(0.0D0, 1.0D0)
      OPEN(7, FILE = 're-out.txt')
      WRITE(7,*)'    REAL-TIME PROPAGATION'
  END SELECT  
  WRITE(7,*)
  WRITE(7,*)'    Number Threads   =', omp_get_max_threads()
  SELECT CASE (OPT_SO)
    CASE(1)
      WRITE(7,*) '    RASHBA SO Coupling,  ' , 'GAMMA = ', GAM0
    CASE(2)
      WRITE(7,*) '    DRESSELHAUS SO Coupling,  ', 'GAMMA = ', GAM0
    CASE DEFAULT
      WRITE(7,*) '    NO SO Coupling   ',   'GAMMA = ', GAM0
  END SELECT
  WRITE(7,*) '    RABI Coupling Omega = ', OMEGA0

    WRITE(7,*) '    Rotational vel = ', OMEG
  SELECT CASE (RANDOM)
    CASE(0)
       WRITE(7,*) '    NO RANDOM PHASE'
    CASE(1)
       WRITE(7,*) '    RANDOM PHASE INCLUDED'  
  END SELECT

  WRITE(7,903) C0, C2, MAG_0
  WRITE(7,*)
  WRITE(7,905) NX+1, NY+1
  WRITE(7,906) DX, DY, DT
  WRITE(7,907) NSTP, NPAS, NRUN
  WRITE(7,*)
 
  CALL INITIALIZE()
 
  CALL reNORM(CP,MAGNET)
 


 !$OMP PARALLEL DO PRIVATE(J, I)
  DO J = 0, NY; DO I = 0, NX
 !   DO L=1,3
    P(:,I,J) = ABS(CP(:,I,J))
    P2(:,I,J) = P(:,I,J) * P(:,I,J)
!    END DO 
  END DO; END DO
  !$OMP END PARALLEL DO

!  SELECT CASE (OPT_PROP)
!    CASE(1)
!      OPEN(11, FILE = 'im-den-ini.txt')
 !     OPEN(21, FILE = 'im-wave-fun-ini.txt')
!      OPEN(32, FILE = 'im-den-rad-ini.txt')
!    CASE(2)
!      OPEN(11, FILE = 're-den-ini.txt')
 !     OPEN(21, FILE = 're-wave-fun-ini.txt')
!      OPEN(32, FILE = 're-den-rad-ini.txt')
!  END SELECT  

!  DO I = 0, NX
!    WRITE(32, 1011) X(I),  P2(:,I,NX2), SUM(P2(:,I, NX2))
!    DO J = 0, NY
!      WRITE(11, 1000) X(I), Y(J), P2(:,I,J), SUM(P2(:,I, J))
 !     WRITE(21, 1001) X(I), Y(J), CP(:,I,J)
!    END DO
!    WRITE(11, *)
 !   WRITE(21, *)
!  END DO
!  CLOSE(32)
!  CLOSE(11)
 ! CLOSE(21)
 
  CALL COEF()
  CALL CALCULATE_TRAP()
  CALL RENORM(CP,MAGNET)
  IF (NSTP /= 0) THEN
    DO K = 1, NSTP
      CALL HERM(CP, DT)
      CALL CALCNU(CP, DT)
      CALL LUX(CP)
      CALL LUY(CP)
     IF(OPT_SO.NE.0) CALL SO(CP,  DT) 
      CALL RENORM(CP,MAGNET)
    END DO
  END IF
!
  CALL RADIUS(CP,RAD)    
  
 
   CALL ENERGY(CP,  EN)
  WRITE(7,1003)
  WRITE(7,1004)
  WRITE(7,1003)
  SELECT CASE (OPT_PROP)
    CASE(1)
      WRITE (7, 1005) RAD(1:3), EN, MAGNET
    CASE(2) 
      WRITE (7, 1008) RAD(1:3), EN, MAGNET 
  END SELECT   
!
  DO K = 1, NPAS
    CALL HERM(CP, DT)
    CALL CALCNU(CP, DT)
    CALL LUX(CP)
    CALL LUY(CP)
    IF(OPT_SO.NE.0) CALL SO(CP,  DT)
    CALL RENORM(CP,MAGNET)
  END DO
!   
  CALL RADIUS(CP, RAD)    
 
   CALL ENERGY(CP,  EN)
 
  SELECT CASE (OPT_PROP)
    CASE(1)
      WRITE (7,1006) RAD(1:3), EN , MAGNET
    CASE(2)
      WRITE (7,1009) RAD(1:3), EN , MAGNET
  END SELECT     
 
  DO K = 1, NRUN
    CALL HERM(CP, DT)
    CALL CALCNU(CP, DT)
    CALL LUX(CP)
    CALL LUY(CP)
    IF(OPT_SO.NE.0) CALL SO(CP,  DT)
    CALL RENORM(CP,MAGNET)
  END DO
!
  CALL RADIUS(CP,RAD)    
 
   CALL ENERGY(CP,  EN)
 
  SELECT CASE (OPT_PROP)
    CASE(1)
      WRITE (7, 1007) RAD(1:3), EN ,MAGNET
    CASE(2)
      WRITE (7, 1010) RAD(1:3), EN, MAGNET
  END SELECT   
  WRITE(7,1003)

  !$OMP PARALLEL DO PRIVATE(J, I)
  DO J = 0, NY; DO I = 0, NX
 
    P(:,I,J) = ABS(CP(:,I,J))
    P2(:,I,J) = P(:,I,J) * P(:,I,J)
 
  END DO; END DO
  !$OMP END PARALLEL DO

  SELECT CASE (OPT_PROP)
    CASE(1)
      OPEN(12, FILE = 'im-den-fin.txt')
      OPEN(32, FILE = 'im-den-rad-fin.txt')
      OPEN(22, FILE = 'im-wave-fun-fin.txt')
    CASE(2)
      OPEN(12, FILE = 're-den-fin.txt')
      OPEN(32, FILE = 're-den-rad-fin.txt')
      OPEN(22, FILE = 're-wave-fun-fin.txt')
  END SELECT  
 
  DO I = 0, NX
    WRITE(32, 1011)X(I),  P2(:,I,NX2), SUM(P2(:,I, NX2))
    DO J = 0, NY
      WRITE(12, 1000) X(I), Y(J), P2(:,I,J), SUM(P2(:,I, J))
      WRITE(22, 1001) X(I), Y(J), CP(:,I,J)
    END DO
    WRITE(12, *)
    WRITE(22, *)
  END DO
  CLOSE(32)
  CLOSE(12)
  CLOSE(22)


  998 FORMAT(6E20.4)
 SELECT CASE (OPT_PROP)
    CASE(1)
  OPEN(32, FILE = 'im-phase.txt')
  CASE(2)
 OPEN(32, FILE = 're-phase.txt')
  END SELECT  

  DO I = 1, NX-1  
  DO J = 1, NY-1      
 
        WRITE (32, 998)X(I),Y(J), ATAN2(AIMAG(CP(1,I,J)),real(CP(1,I,J)))  &
,ATAN2(AIMAG(CP(2,I,J)),real(CP(2,I,J))) ,ATAN2(AIMAG(CP(3,I,J)),real(CP(3,I,J))) 
 
 END DO  
 WRITE (32, 998)
 END DO
 CLOSE(32)



!
  CALL SYSTEM_CLOCK (CLCK_COUNTS_END, CLCK_RATE)
  CALL CPU_TIME(T2)
  WRITE (7, 1001)
  WRITE (7,*)
  WRITE (7,'(A,I7,A)') ' Clock Time: ', (CLCK_COUNTS_END - CLCK_COUNTS_BEG)/INT (CLCK_RATE,8), ' seconds'
  WRITE (7,'(A,I7,A)') '   CPU Time: ', INT(T2-T1), ' seconds'
  CLOSE (7)
!
  1011 FORMAT(F8.3, 4F16.8)
  1000 FORMAT(2F8.3, 4F16.8)
  1001 FORMAT(2F8.3, 6F16.8)
!-------------------------------------------------------------------------------
END PROGRAM SPIN_ONE_GPE2D_IMAG
!*******************************************************************************
  SUBROUTINE INITIALIZE()
  USE COMM_DATA, ONLY: NX, NX2, NY, NY2, NSTP
  USE GPE_DATA, ONLY: PI, DX, DY, X, Y, X2, Y2, R2, CP, C0, C2, C0MC2, C0PC2,   gam0, MAG_0,CI,RANDOM,OPT_PROP,OPT_ST, OPT_SO
  IMPLICIT NONE
  REAL (8), DIMENSION((1+NX)*(1+NY)) :: RANDNUM
  REAL (8) :: FAC_NORM, TX, TY, PSI_1R, PSI_1I, PSI_2R, PSI_2I, PSI_3R, PSI_3I,XJ
  COMPLEX (8), DIMENSION(1:3,0:200,0:200) :: TP
  INTEGER :: I, J, K
 INTEGER, DIMENSION (100) :: seed 
  C0MC2 = C0 - C2
  C0PC2 = C0 + C2 

     seed = 13
  IF(RANDOM.EQ.1)  CALL RANDOM_SEED(PUT=seed)
  
    IF(RANDOM.EQ.1)   call random_number(RANDNUM)
 
      IF(RANDOM.EQ.0) RANDNUM=0.d0 

  FORALL (J=0:NX) 
    X(J) = (J - NX2) * DX
    X2(J) = X(J) * X(J)
  END FORALL
  FORALL (J=0:NY) 
    Y(J) = (J - NY2) * DY
    Y2(J) = Y(J) * Y(J)
  END FORALL
  DO J = 0, NX
   R2(J,:) = X2(J) + Y2
  END DO
 
  IF (NSTP == 0) THEN
    IF (OPT_PROP.EQ.1) WRITE(*,'(a)') "Run the program using the input file to read. e.g.: ./imag2d < im2d-fin-wf.txt"
    IF (OPT_PROP.EQ.2) WRITE(*,'(a)') "Run the program using the input file to read. e.g.: ./real2d < im2d-fin-wf.txt"


  DO I = 0, NX; DO J = 0, NY
 
  K=(I+1)*(J+1)

      READ(*,*) TX, TY, PSI_1R, PSI_1I, PSI_2R, PSI_2I, PSI_3R, PSI_3I
      CP(1,I,J) = CMPLX(PSI_1R, PSI_1I, KIND(2))  * EXP (2.d0*PI*CI*RANDNUM(K))
      CP(2,I,J) = CMPLX(PSI_2R, PSI_2I, KIND(2))  * EXP (2.d0*PI*CI*RANDNUM(K))
      CP(3,I,J) = CMPLX(PSI_3R, PSI_3I, KIND(2))  * EXP (2.d0*PI*CI*RANDNUM(K))
!      tP(1,I,J) = CMPLX(PSI_1R, PSI_1I, KIND(2))  * EXP (2.d0*PI*CI*RANDNUM(K))
!      TP(2,I,J) = CMPLX(PSI_2R, PSI_2I, KIND(2))  * EXP (2.d0*PI*CI*RANDNUM(K))
!      tP(3,I,J) = CMPLX(PSI_3R, PSI_3I, KIND(2))  * EXP (2.d0*PI*CI*RANDNUM(K))
    END DO; END DO


!    DO I = 0, NX; DO J = 0, NY
!  CP(1,I,J) =tP(1,20+I,20+J)
!    CP(2,I,J) =  TP(2,20+I,20+J)
!  CP(3,I,J) = TP(3,20+I,20+J)
!    END DO; END DO


  ELSE

    FAC_NORM = SQRT(PI) ! (PI)^(1/4)


  IF(OPT_SO.EQ.0)THEN   ! NO SO coupling
  
 DO I = 0, NX; DO J = 0, NY 
   K=(I+1)*(J+1)
    CP(1,I,J) =(1.d0+MAG_0)* 0.5D0 * EXP(-R2(I,J)/8.0D0) / FAC_NORM  * EXP (2.d0*PI*CI*RANDNUM(K))
    CP(2,I,J) =SQRT((1.d0-MAG_0**2)* 0.5D0) * EXP(-R2(I,J)/8.0D0) / FAC_NORM   * EXP (2.d0*PI*CI*RANDNUM(K)) 
    CP(3,I,J) = (1.d0-MAG_0)* 0.5D0 *EXP(-R2(I,J)/8.0D0) / FAC_NORM   * EXP (2.d0*PI*CI*RANDNUM(K))
  END DO; END DO
  END IF
 
 

         IF(OPT_SO.EQ.1)THEN   ! Rashba coupling

   IF (GAM0.GT.0.0001d0.AND.C2.GT.0.000D0)  THEN   !! Ferromagnetic

 
 

 
!  The following choice is for (-1,0,+1) state  
  DO I = 0, NX; DO J = 0, NY
    K=(I+1)*(J+1)
    CP(1,I,J) =(X(I)-CI*Y(J))*(1.d0+MAG_0)* 0.5D0 * EXP(-R2(I,J)/8.0D0) / FAC_NORM   * EXP (2.d0*PI*CI*RANDNUM(K))
    CP(2,I,J) = SQRT((1.d0-MAG_0**2)* 0.5D0)* EXP(-R2(I,J)/8.0D0) / FAC_NORM   * EXP (2.d0*PI*CI*RANDNUM(K))
    CP(3,I,J) =(X(I)+CI*Y(J))* (1.d0-MAG_0)* 0.5D0 *EXP(-R2(I,J)/8.0D0) / FAC_NORM   * EXP (2.d0*PI*CI*RANDNUM(K))
  END DO; END DO
   END IF
   




  IF (GAM0.GT.0.0001d0.AND.C2.LT.0.0000D0)  THEN    !! Antiferromagnetic

 
   

!  The following choice is for (0,+1,+2) state    
 DO I = 0, NX; DO J = 0, NY
    CP(1,I,J) =(X(I)+CI*Y(J))**0*(1.d0+MAG_0)* 0.5D0 * EXP(-R2(I,J)/8.0D0) / FAC_NORM   * EXP (2.d0*PI*CI*RANDNUM(K))
    CP(2,I,J) =(X(I)+CI*Y(J))**1*SQRT((1.d0-MAG_0**2)* 0.5D0) * EXP(-R2(I,J)/8.0D0) / FAC_NORM   * EXP (2.d0*PI*CI*RANDNUM(K))
    CP(3,I,J) =(X(I)+CI*Y(J))**2* (1.d0-MAG_0)* 0.5D0 *EXP(-R2(I,J)/8.0D0) / FAC_NORM   * EXP (2.d0*PI*CI*RANDNUM(K))
  END DO; END DO
  END IF 

  

  END IF


   IF(OPT_SO.EQ.2)THEN   ! Dresselhaus coupling

   IF (GAM0.GT.0.0001d0.AND.C2.GT.0.000D0)  THEN   !! Ferromagnetic

 
 

 
!  The following choice is for (+1,0,-1) state for small gamma 
  DO I = 0, NX; DO J = 0, NY
    K=(I+1)*(J+1)
    CP(1,I,J) =(X(I)+CI*Y(J))*(1.d0+MAG_0)* 0.5D0 * EXP(-R2(I,J)/8.0D0) / FAC_NORM   * EXP (2.d0*PI*CI*RANDNUM(K))
    CP(2,I,J) = SQRT((1.d0-MAG_0**2)* 0.5D0)* EXP(-R2(I,J)/8.0D0) / FAC_NORM   * EXP (2.d0*PI*CI*RANDNUM(K))
    CP(3,I,J) =(X(I)-CI*Y(J))* (1.d0-MAG_0)* 0.5D0 *EXP(-R2(I,J)/8.0D0) / FAC_NORM   * EXP (2.d0*PI*CI*RANDNUM(K))
  END DO; END DO
   END IF
  



  IF (GAM0.GT.0.0001d0.AND.C2.LT.0.0000D0)  THEN    !! Antiferromagnetic

 


!  The following choice is for (0,-1,-2) state    
 DO I = 0, NX; DO J = 0, NY
    CP(1,I,J) =(X(I)-CI*Y(J))**0*(1.d0+MAG_0)* 0.5D0 * EXP(-R2(I,J)/8.0D0) / FAC_NORM   * EXP (2.d0*PI*CI*RANDNUM(K))
    CP(2,I,J) =(X(I)-CI*Y(J))**1*SQRT((1.d0-MAG_0**2)* 0.5D0) * EXP(-R2(I,J)/8.0D0) / FAC_NORM   * EXP (2.d0*PI*CI*RANDNUM(K))
    CP(3,I,J) =(X(I)-CI*Y(J))**2* (1.d0-MAG_0)* 0.5D0 *EXP(-R2(I,J)/8.0D0) / FAC_NORM   * EXP (2.d0*PI*CI*RANDNUM(K))
  END DO; END DO
  END IF 
 


  END IF



  END IF 
 

END SUBROUTINE INITIALIZE

SUBROUTINE CALCULATE_TRAP()
  USE GPE_DATA, ONLY : NX, LAMBDA2, X2, Y2, V
  IMPLICIT NONE
  INTEGER :: J, L 

    DO L = 1, 3
  !$OMP PARALLEL DO  
  DO J = 0, NX
 
      V(L,J,:) = LAMBDA2 * (X2(J) + Y2) / 2.0D0  
    END DO
 
 !$OMP END PARALLEL DO 
  END DO
END SUBROUTINE CALCULATE_TRAP

SUBROUTINE COEF()
  USE COMM_DATA, ONLY :  NX, NY,NXX,NYY
  USE GPE_DATA, ONLY : DT, DX, DY,  CIJ,OMEG,CI,X,Y
  USE CN_DATA
  IMPLICIT NONE
  INTEGER :: J,I
 
  COMPLEX (8) :: CA0,  CB0, DXX, DYY
  COMPLEX (8) :: CDT 
 
  DXX = 1.0D0 / (DX * DX * 2.0D0)
  DYY = 1.0D0 / (DY * DY * 2.0D0)
  CDT = CIJ * DT
  CA0 = 1.0D0 + CDT * DXX
  CA0R = 1.0D0 - CDT * DXX
  CB0 = 1.0D0 + CDT * DYY
  CB0R = 1.0D0 - CDT * DYY
!

  CC0 =CIJ* CI*DT*OMEG/2.0D0   !   added   *2  is for XOP



  CTMPX = CC0/(2.0D0*DX)
  CT0 = CDT*DXX/2.d0
  CAM = -(CT0 - CTMPX*Y)
  CAP = -(CT0 + CTMPX*Y)

 !$OMP PARALLEL DO  PRIVATE(I,J)
  DO J = 0, NY
    CALA(NXX, J) = 0.0D0
    CGAA(NXX, J) = -1.0D0/CA0
    DO I = NXX, 1, -1
      CALA(I-1,J) = CAM(J)*CGAA(I,J)
      CGAA(I-1,J) = -1.0D0/(CA0+CAP(J)*CALA(I-1,J))
    END DO
  END DO
 !$OMP END PARALLEL DO 




  CTMPY = CC0/(2.0D0*DY)
  CT0 = cDT*DYY/2.d0
  CBM = -(CT0 + CTMPY*X)
  CBP = -(CT0 - CTMPY*X)


  !$OMP PARALLEL DO  PRIVATE(I,J)
 DO I = 0, NX
    CALB(I, NYY) = 0.0D0
    CGAB(I, NYY) = -1.0D0/CB0
    DO J = NYY, 1, -1
      CALB(I,J-1) = CBM(I)*CGAB(I,J)
      CGAB(I,J-1) = -1.0D0/(CB0+CBP(I)*CALB(I,J-1))
    END DO
  END DO
 !$OMP END PARALLEL DO 
! 
END SUBROUTINE COEF

SUBROUTINE SO(CP,  DT)   
  USE COMM_DATA, ONLY : NX, NY
  USE GPE_DATA, ONLY : DX, DY, CI, CIJ, GAM,   C0PC2, C0MC2,  OPT_SO
  USE UTIL, ONLY :  DIFF
  IMPLICIT NONE
  COMPLEX (8), DIMENSION(:,0:,0:), INTENT(INOUT) :: CP
  REAL (8), INTENT(IN) :: DT
!-----------------------------------------------------------  
  REAL (8), DIMENSION(SIZE(CP,1),0:SIZE(CP,2)-1,0:SIZE(CP,3)-1) :: CPR, CPI, DCPR_X, DCPI_X, DCPR_Y, DCPI_Y 
  COMPLEX (8), DIMENSION(SIZE(CP,1),0:SIZE(CP,2)-1,0:SIZE(CP,3)-1) :: DPX, DPY

  INTEGER :: I, J, K
!
   !$OMP PARALLEL DO PRIVATE(I,J)
  DO I = 0, NX
  DO J = 0, NY
  CPR(:,I,J) = REAL(CP(:,I,J))
  CPI(:,I,J) = IMAG(CP(:,I,J))
 END DO
 END DO
  !$OMP END PARALLEL DO 
!
    DO K = 1, 3
   !$OMP PARALLEL DO  
 
    DO J = 0, NY    
      DCPR_X(K,:,J) = DIFF(CPR(K,:,J), DX)
      DCPI_X(K,:,J) = DIFF(CPI(K,:,J), DX)
    END DO   
  
   !$OMP END PARALLEL DO 
    END DO

    DO K = 1, 3
   !$OMP PARALLEL DO  
 
    DO I = 0, NX
      DCPR_Y(K,I,:) = DIFF(CPR(K,I,:), DY)
      DCPI_Y(K,I,:) = DIFF(CPI(K,I,:), DY)
    END DO
 
   !$OMP END PARALLEL DO 
  END DO

     !$OMP PARALLEL DO PRIVATE(I,J)
  DO I = 0, NX
  DO J = 0, NY
  DPX(:,I,J) = CMPLX(DCPR_X(:,I,J), DCPI_X(:,I,J))
  DPY(:,I,J) = CMPLX(DCPR_Y(:,I,J), DCPI_Y(:,I,J))
   END DO
 END DO
    !$OMP END PARALLEL DO 
  SELECT CASE (OPT_SO)
    CASE (1) ! Rashba

   !$OMP PARALLEL DO  PRIVATE(I,J)

      DO J = 0, NY; DO I = 0, NX
        CP(1,I,J) = CP(1,I,J) + CIJ * (CI * DPY(2,I,J) - DPX(2,I,J)) * DT * GAM  
        CP(3,I,J) = CP(3,I,J) + CIJ * (CI * DPY(2,I,J) + DPX(2,I,J)) * DT * GAM  
        CP(2,I,J) = CP(2,I,J) + CIJ * (CI * DPY(1,I,J) + DPX(1,I,J)  & 
                                   + CI * DPY(3,I,J) - DPX(3,I,J)) * DT * GAM  
      END DO; END DO 
 
   !$OMP END PARALLEL DO 
 
   CASE (2) ! Dresselhaus

   !$OMP PARALLEL DO  PRIVATE(I,J)

      DO J = 0, NY; DO I = 0, NX
        CP(1,I,J) = CP(1,I,J) + CIJ * (-CI * DPY(2,I,J) - DPX(2,I,J)) * DT * GAM  
        CP(3,I,J) = CP(3,I,J) + CIJ * (-CI * DPY(2,I,J) + DPX(2,I,J)) * DT * GAM  
        CP(2,I,J) = CP(2,I,J) + CIJ * (-CI * DPY(1,I,J) + DPX(1,I,J)  & 
                                   - CI * DPY(3,I,J) - DPX(3,I,J)) * DT * GAM  
      END DO; END DO 
   !$OMP END PARALLEL DO 

    CASE (3) ! Alternate

   !$OMP PARALLEL DO  PRIVATE(I,J)

      DO J = 0, NY; DO I = 0, NX
        CP(1,I,J) = CP(1,I,J) + CIJ * (CI * DPX(2,I,J) + DPY(2,I,J)) * DT * GAM  
        CP(3,I,J) = CP(3,I,J) + CIJ * (CI * DPX(2,I,J) - DPY(2,I,J)) * DT * GAM  
        CP(2,I,J) = CP(2,I,J) + CIJ * (CI * DPX(1,I,J) - DPY(1,I,J)  & 
                                   + CI * DPX(3,I,J) + DPY(3,I,J)) * DT * GAM  
      END DO; END DO 

   !$OMP END PARALLEL DO 
  END SELECT
 
END SUBROUTINE SO

SUBROUTINE LUX(CP)
  USE COMM_DATA
  USE CN_DATA, ONLY :  CA0R,   CALA, CGAA, CBE, CXX,CAM,CAP
  USE GPE_DATA, ONLY :OPT_PROP
  IMPLICIT NONE
  COMPLEX (8), DIMENSION(:,0:,0:), INTENT(INOUT) :: CP
  INTEGER :: I, J, L
!
  !$OMP PARALLEL DO  PRIVATE(I, J,CXX)
 
    DO J = 0, NY
!      CBE(:,NXX,J) = 0.0D0 !CP(L,NX,J)  
   IF (opt_PROP.eq.1) CBE(:,NXX,J) = 0.d0
    IF (opt_PROP.eq.2) CBE(:,NXX,J) = CP(:,NX,J)
      DO I = NXX, 1, -1
        CXX(:) = -CAP(J) * CP(:,I+1,J) + CA0R*CP(:,I,J) - CAM(J)*CP(:,I-1,J)        
        CBE(:,I-1,J) = CGAA(I,J)*(CAP(J)*CBE(:,I,J)-CXX(:)) 
      END DO 
 

      CP(:,0,J) = 0.0D0
      DO I = 0, NXX
        CP(:,I+1,J) = CALA(I,J)*CP(:,I,J)+CBE(:,I,J)
      END DO       
      CP(:,NX,J) = 0.0D0
    END DO
 
  !$OMP END PARALLEL DO
END SUBROUTINE LUX
! 
SUBROUTINE LUY(CP)
  USE OMP_LIB
  USE COMM_DATA  
  USE GPE_DATA, ONLY :OPT_PROP
  USE CN_DATA, ONLY :  CB0R,  CALB, CGAB, CBE, CYY,CBM,CBP
  IMPLICIT NONE
  COMPLEX (8), DIMENSION(:,0:,0:), INTENT(INOUT) :: CP
  INTEGER ::  I, J, L
  !$OMP PARALLEL DO   PRIVATE( I,J,CYY)
  
    DO I = 0, NX
  
  IF (opt_PROP.eq.1) CBE(:,I,NYY) = 0.d0
    IF (opt_PROP.eq.2) CBE(:,I,NYY) =   CP(:,I,NY)

      DO J = NYY, 1, -1
        CYY(:) = -CBP(I)*CP(:,I,J+1)+CB0R*CP(:,I,J)-CBM(I)*CP(:,I,J-1)      
        CBE(:,I,J-1) = CGAB(I,J)*(CBP(I)*CBE(:,I,J)-CYY(:))
      END DO
 
      CP(:,I,0) = 0.0D0
      DO J = 0, NYY
        CP(:,I,J+1) = CALB(I,J)*CP(:,I,J)+CBE(:,I,J)
      END DO     
      CP(:,I,NY) = 0.0D0
    END DO
 ! END DO
  !$OMP END PARALLEL DO
END SUBROUTINE LUY

SUBROUTINE HERM(CP, DT)
  USE GPE_DATA, ONLY : NX, NY, OMEGA,  C2, CIJ
  IMPLICIT NONE
  COMPLEX (8), DIMENSION(:,0:,0:), INTENT(INOUT) :: CP
  REAL (8), INTENT(IN) :: DT

  COMPLEX (8), DIMENSION(SIZE(CP,1),SIZE(CP,1)) :: V
  COMPLEX (8) :: A, B, CA, CB, CTMP(3)
  REAL (8) :: Cx, C2x
  INTEGER :: I, J
!
  !$OMP PARALLEL DO PRIVATE(J, I, V, A, B, CA, CB, CTMP, Cx, C2x)
  DO J = 0, NY; DO I = 0, NX
    A = OMEGA + C2 * CP(2,I,J) * CONJG(CP(3,I,J)) 
    B = OMEGA + C2 * CP(1,I,J) * CONJG(CP(2,I,J)) 
    C2x = ABS(A) ** 2 + ABS(B) ** 2
    Cx = SQRT(C2x)
    CA = CONJG(A)
    CB = CONJG(B)
    V(1,1) = A
    V(1,2) = B
    V(1,3) = A
    V(2,1) = Cx
    V(2,2) = 0.0D0
    V(2,3) = -Cx
    V(3,1) = CB
    V(3,2) = - CA
    V(3,3) = CB
    CTMP(1) = EXP(-Cx * CIJ * DT) * (CA * CP(1,I,J) + Cx * CP(2,I,J) &
              + B * CP(3,I,J)) / (2.0D0 * C2x)
    CTMP(2) = (CB * CP(1,I,J) - A * CP(3,I,J)) / C2x
    CTMP(3) = EXP(Cx * CIJ * DT) * (CA * CP(1,I,J) - Cx * CP(2,I,J)  &
              + B * CP(3,I,J)) / (2.0D0 * C2x)
    CP(:,I,J) = MATMUL(V, CTMP)
  END DO; END DO
  !$OMP END PARALLEL DO
END SUBROUTINE HERM

SUBROUTINE CALCNU(CP, DT)
  USE GPE_DATA, ONLY : NX, NY, V,   C0, C0MC2, C0PC2,CIJ
  IMPLICIT NONE
  COMPLEX (8), DIMENSION(:,0:,0:), INTENT(INOUT) :: CP
  REAL (8), INTENT(IN) :: DT

  REAL (8), DIMENSION(SIZE(CP,1),0:SIZE(CP,2)-1,0:SIZE(CP,3)-1) :: P, P2
  COMPLEX (8), DIMENSION(3) :: CTMP
  INTEGER :: I, J,L
!
 !$OMP PARALLEL DO PRIVATE(J, I)
  DO J = 0, NY; DO I = 0, NX
 
    P(:,I,J) = ABS(CP(:,I,J))
    P2(:,I,J) = P(:,I,J) * P(:,I,J)
 
  END DO; END DO
  !$OMP END PARALLEL DO

  !$OMP PARALLEL DO PRIVATE(J, I,CTMP)
  DO J = 0, NY; DO I = 0, NX
  
    CTMP(1) = V(1,I,J) + C0PC2 * (P2(1,I,J) + P2(2,I,J)) + C0MC2 * P2(3,I,J)
    CTMP(2) = V(2,I,J) + C0PC2 * (P2(3,I,J) + P2(1,I,J)) + C0 * P2(2,I,J)
    CTMP(3) = V(3,I,J) + C0PC2 * (P2(2,I,J) + P2(3,I,J)) + C0MC2 * P2(1,I,J)
    CP(:,I,J) = EXP(-CIJ * CTMP * DT) * CP(:,I,J)
  END DO; END DO
  !$OMP END PARALLEL DO
!
END SUBROUTINE CALCNU
 
FUNCTION INTEGRATE(U, DX, DY) RESULT(RES)
  USE COMM_DATA, ONLY : NX 
  USE UTIL, ONLY : SIMP
  IMPLICIT NONE
  REAL (8), DIMENSION(0:, 0:), INTENT(IN) :: U
  REAL (8), INTENT (IN) :: DX, DY
  REAL (8) :: RES
  REAL (8), DIMENSION(0:NX) :: TMP1D
  INTEGER :: I
!
  !$OMP PARALLEL DO PRIVATE(I)
  DO I = 0, NX
    TMP1D(I) = SIMP(U(I,0:), DY)
  END DO 
  !$OMP END PARALLEL DO
  RES = SIMP(TMP1D, DX)
END FUNCTION INTEGRATE

SUBROUTINE RADIUS(CP, RAD)
  USE COMM_DATA, ONLY : NX, NY
  USE GPE_DATA, ONLY : DX, DY, R2, ACCUR
  USE UTIL, ONLY : SIMP
  IMPLICIT NONE
  COMPLEX (8), DIMENSION(:,0:,0:), INTENT(INOUT) :: CP
  REAL (8), DIMENSION(4), INTENT(OUT) :: RAD
!----------------------------------------------------------------------
  INTERFACE
    FUNCTION INTEGRATE(P2, DX, DY) RESULT(RES)
      IMPLICIT NONE
      REAL (8), DIMENSION(0:,0:), INTENT(IN) :: P2
      REAL (8), INTENT(IN) :: DX, DY
      REAL (8) :: RES
    END FUNCTION INTEGRATE
  END INTERFACE
!----------------------------------------------------------------------
  REAL (8), DIMENSION(SIZE(CP,1),0:SIZE(CP,2)-1,0:SIZE(CP,3)-1) :: P, P2
  REAL (8), DIMENSION(0:SIZE(CP,2)-1,0:SIZE(CP,3)-1) :: TMP2D
  INTEGER :: I, J
 !$OMP PARALLEL DO PRIVATE(J, I)
  DO J = 0, NY; DO I = 0, NX
 
    P(:,I,J) = ABS(CP(:,I,J))
    P2(:,I,J) = P(:,I,J) * P(:,I,J)
 
  END DO; END DO
  !$OMP END PARALLEL DO
  IF (MAXVAL(P2(1,:,:)) < ACCUR) THEN
    rad(1) = 0.0D0
  ELSE
    rad(1) = INTEGRATE(P2(1,:,:) * R2, DX, DY) / INTEGRATE(P2(1,:,:), DX, DY)
  END IF
  IF(MAXVAL(P2(2,:,:)) < ACCUR) THEN
    RAD(2)=0.D0
  ELSE
    RAD(2) = INTEGRATE(P2(2,:,:) * R2, DX, DY) / INTEGRATE(P2(2,:,:), DX, DY)
  END IF
  IF (MAXVAL(P2(3,:,:)) < ACCUR) THEN
    RAD(3)=0
  ELSE
    RAD(3) = INTEGRATE(P2(3,:,:) * R2, DX, DY) / INTEGRATE(P2(3,:,:), DX, DY)
  END IF 
!$OMP PARALLEL DO PRIVATE(J, I)
  DO J = 0, NY; DO I = 0, NX
    TMP2D(I,J) = SUM(P2(:,I,J))
  END DO; END DO
  !$OMP END PARALLEL DO
  RAD(4) = INTEGRATE(TMP2D * R2, DX, DY) / INTEGRATE(TMP2D, DX, DY)
  RAD  = SQRT(RAD)
END SUBROUTINE RADIUS

!SUBROUTINE RENORM(CP)
!  USE COMM_DATA, ONLY : NX, NY
!  USE GPE_DATA, ONLY : DX, DY, D, ZNORM, MAG, MAG_0
!  IMPLICIT NONE
!  COMPLEX (8), DIMENSION(:,0:,0:), INTENT(INOUT) :: CP
!!----------------------------------------------------------------------
!  INTERFACE
!    FUNCTION INTEGRATE(P2, DX, DY) RESULT(RES)
!      IMPLICIT NONE
!      REAL (8), DIMENSION(0:,0:), INTENT(IN) :: P2
!      REAL (8), INTENT(IN) :: DX, DY
!      REAL (8) :: RES
!    END FUNCTION INTEGRATE
!  END INTERFACE
!!----------------------------------------------------------------------
!  REAL (8), DIMENSION(SIZE(CP,1),0:SIZE(CP,2)-1,0:SIZE(CP,3)-1) :: P, P2
 
!  REAL (8), DIMENSION(SIZE(CP,1)) :: D2 
!  REAL (8) :: MAG2, NUM, DEN
!  INTEGER :: I, J,L
!  !$OMP PARALLEL DO PRIVATE(J, I)
!  DO J = 0, NY; DO I = 0, NX
 
!    P(:,I,J) = ABS(CP(:,I,J))
!    P2(:,I,J) = P(:,I,J) * P(:,I,J)
 
!  END DO; END DO
!  !$OMP END PARALLEL DO
!  ZNORM(1) = INTEGRATE(P2(1,:,:), DX, DY)
!  ZNORM(2) = INTEGRATE(P2(2,:,:), DX, DY)
!  ZNORM(3) = INTEGRATE(P2(3,:,:), DX, DY)
 
!  MAG = MAG_0
!  MAG2 = MAG * MAG
!  NUM = 1.0D0 - MAG2
!  DEN = ZNORM(2) + SQRT(4.0D0 * NUM * ZNORM(1) * ZNORM(3) + MAG2 * ZNORM(2) * ZNORM(2))
!  D2(2) = NUM/DEN
!  D2(1) = (1.0D0 + MAG - D2(2) * ZNORM(2)) / (2.0D0 * ZNORM(1)) 
!  D2(3) = (1.0D0 - MAG - D2(2) * ZNORM(2)) / (2.0D0 * ZNORM(3)) 
!  IF (MINVAL(D2) < 0.0D0) THEN
     
!     STOP
!  END IF
!  D = SQRT(D2)
!   !$OMP PARALLEL DO PRIVATE(J, I) 
!   DO J = 0, NY; DO I = 0, NX
!  CP(1,i,j) = D(1) * CP(1,i,j)
!  CP(2,i,j) = D(2) * CP(2,I,j)
!  CP(3,i,j) = D(3) * CP(3,I,j)
!  END DO; END DO
!  !$OMP END PARALLEL DO
!END SUBROUTINE RENORM

SUBROUTINE RENORM(CP,MAGNET)
  USE COMM_DATA, ONLY : NX, NY
  USE GPE_DATA, ONLY : DX, DY, D, ZNORM, MAG, MAG_0,  OPT_SO
  IMPLICIT NONE
  COMPLEX (8), DIMENSION(:,0:,0:), INTENT(INOUT) :: CP
!----------------------------------------------------------------------
  INTERFACE
    FUNCTION INTEGRATE(P2, DX, DY) RESULT(RES)
      IMPLICIT NONE
      REAL (8), DIMENSION(0:,0:), INTENT(IN) :: P2
      REAL (8), INTENT(IN) :: DX, DY
      REAL (8) :: RES
    END FUNCTION INTEGRATE
  END INTERFACE
!----------------------------------------------------------------------
  REAL (8), DIMENSION(SIZE(CP,1),0:SIZE(CP,2)-1,0:SIZE(CP,3)-1) :: P, P2
 
  REAL (8), DIMENSION(SIZE(CP,1)) :: D2 
  REAL (8) :: MAG2, NUM, DEN,magnet
  INTEGER :: I, J,L
  !$OMP PARALLEL DO PRIVATE(J, I)
  DO J = 0, NY; DO I = 0, NX
 
    P(:,I,J) = ABS(CP(:,I,J))
    P2(:,I,J) = P(:,I,J) * P(:,I,J)
 
  END DO; END DO
  !$OMP END PARALLEL DO
  ZNORM(1) = INTEGRATE(P2(1,:,:), DX, DY)
  ZNORM(2) = INTEGRATE(P2(2,:,:), DX, DY)
  ZNORM(3) = INTEGRATE(P2(3,:,:), DX, DY)
 
  IF(OPT_SO.EQ.0)  THEN
  MAG = MAG_0
  MAG2 = MAG * MAG
  NUM = 1.0D0 - MAG2
  DEN = ZNORM(2) + SQRT(4.0D0 * NUM * ZNORM(1) * ZNORM(3) + MAG2 * ZNORM(2) * ZNORM(2))
  D2(2) = NUM/DEN
  D2(1) = (1.0D0 + MAG - D2(2) * ZNORM(2)) / (2.0D0 * ZNORM(1)) 
  D2(3) = (1.0D0 - MAG - D2(2) * ZNORM(2)) / (2.0D0 * ZNORM(3)) 
  magnet=mag_0
  IF (MINVAL(D2) < 0.0D0) THEN
     
     STOP
  END IF
  END  IF

  IF(OPT_SO.NE.0)  THEN
  DEN = 1.d0/(ZNORM(1)+ZNORM(2)+ZNORM(3))
  magnet= (ZNORM(1)-ZNORM(3))*DEN
  D2(1)=DEN
  D2(2)=DEN
  D2(3)=DEN
  END IF

  D = SQRT(D2)
   !$OMP PARALLEL DO PRIVATE(J, I) 
   DO J = 0, NY; DO I = 0, NX
  CP(1,i,j) = D(1) * CP(1,i,j)
  CP(2,i,j) = D(2) * CP(2,I,j)
  CP(3,i,j) = D(3) * CP(3,I,j)
  END DO; END DO
  !$OMP END PARALLEL DO
END SUBROUTINE RENORM
 
  SUBROUTINE ENERGY(CP,  EN) 
  USE COMM_DATA, ONLY : NX, NY
  USE GPE_DATA, ONLY : DX, DY, CI,  V, OMEGA, GAM, C0, C2, C0PC2, C0MC2,  OPT_SO,OMEG,X,Y
  USE UTIL, ONLY : DIFF
  IMPLICIT NONE
  COMPLEX (8), DIMENSION(:,0:,0:), INTENT(IN) :: CP
  REAL (8), INTENT(OUT) :: EN
!----------------------------------------------------------------------
  INTERFACE
    FUNCTION INTEGRATE(P2, DX, DY) RESULT(RES)
      IMPLICIT NONE
      REAL (8), DIMENSION(0:,0:), INTENT(IN) :: P2
      REAL (8), INTENT(IN) :: DX, DY
      REAL (8) :: RES
    END FUNCTION INTEGRATE
  END INTERFACE
!-----------------------------------------------------------  
  REAL (8), DIMENSION(SIZE(CP,1), 0:NX, 0:NY) :: CPR, CPI, CPR2, CPI2, P2,DPLZ
  REAL (8), DIMENSION(SIZE(CP,1), 0:NX, 0:NY) :: DCPR_X, DCPR_X2, DCPI_X, DCPI_X2
  REAL (8), DIMENSION(SIZE(CP,1), 0:NX, 0:NY) :: DCPR_Y, DCPR_Y2, DCPI_Y, DCPI_Y2
  COMPLEX (8), DIMENSION(SIZE(CP,1), 0:NX, 0:NY) :: CDPX
  COMPLEX (8), DIMENSION(SIZE(CP,1), 0:NX, 0:NY) :: CDPY
  COMPLEX (8), DIMENSION(SIZE(CP,1), 0:NX, 0:NY) :: DPX, DPY
  
  REAL (8), DIMENSION(0:NX, 0:NY) :: TMP2D   
 
 
 
  INTEGER :: I, J, K
!
  

  CPR = REAL(CP)
  CPI = AIMAG(CP)
  CPR2 = CPR * CPR
  CPI2 = CPI * CPI
  P2 = CPR2 + CPI2
 
  DO K = 1, 3
    DO J = 0, NY
      DCPR_X(K,:,J) = DIFF(CPR(K,:,J), DX)
      DCPI_X(K,:,J) = DIFF(CPI(K,:,J), DX)  !
      CDPX(K,:,J)  =   DCPR_X(K,:,J)+ CI *  DCPI_X(K,:,J)
      DCPR_X2(K,:,J) = DCPR_X(K,:,J) * DCPR_X(K,:,J) / 2.0D0
      DCPI_X2(K,:,J) = DCPI_X(K,:,J) * DCPI_X(K,:,J) / 2.0D0
    END DO
    DO I = 0, NX
      DCPR_Y(K,I,:) = DIFF(CPR(K,I,:), DY)
      DCPI_Y(K,I,:) = DIFF(CPI(K,I,:), DY)  !
      CDPY(K,I,:)  =   DCPR_Y(K,I,:)+ CI *  DCPI_Y(K,I,:)
      DCPR_Y2(K,I,:) = DCPR_Y(K,I,:) * DCPR_Y(K,I,:) / 2.0D0
      DCPI_Y2(K,I,:) = DCPI_Y(K,I,:) * DCPI_Y(K,I,:) / 2.0D0
    END DO
 
  END DO
  DPX = CMPLX(DCPR_X, DCPI_X)
  DPY = CMPLX(DCPR_Y, DCPI_Y)
 
  DO K = 1, 3
     DO J = 0, NY; DO I = 0, NX
     DPLZ(K,I,J) =  CPR(K,I,J)*( X(I)*DCPI_Y(K,I,J) - Y(J)*DCPI_X(K,I,J) )
   END DO; END DO
       END DO

      DO J = 0, NY; DO I = 0, NX
    TMP2D(I,J) =   -OMEG*X(I)*CI *(CONJG(CP(1,I,J)) * CDPY(1,I,J))  &
         -OMEG*X(I)*CI *(CONJG(CP(2,I,J)) * CDPY(2,I,J))  &
        -OMEG*X(I)*CI *(CONJG(CP(3,I,J)) * CDPY(3,I,J))   &
        + OMEG*Y(J)*CI*(CONJG(CP(1,I,J)) * CDPX(1,I,J))    &
         + OMEG*Y(J)*CI* (CONJG(CP(2,I,J)) * CDPX(2,I,J))    &
              + OMEG*Y(J)*CI* (CONJG(CP(3,I,J)) * CDPX(3,I,J))   
     END DO; END DO  
  TMP2D = V(1,:,:)*P2(1,:,:) +  V(2,:,:)*P2(2,:,:) + V(3,:,:)*P2(3,:,:)  &
         + C0*0.5D0*(P2(1,:,:)+P2(2,:,:)+P2(3,:,:))*(P2(1,:,:)+P2(2,:,:)+ P2(3,:,:)) &
         + C2*0.5D0 * (P2(1,:,:) + P2(2,:,:) - P2(3,:,:)) * P2(1,:,:) &
         + C2*0.5D0 * (P2(1,:,:) + P2(3,:,:)) * P2(2,:,:)  & 
         + C2 * 0.5D0 * (-P2(1,:,:) + P2(2,:,:)+P2(3,:,:))*P2(3,:,:) & 
         + C2*( CONJG(CP(3,:,:)) *CP(2,:,:) *CP(2,:,:) * CONJG(CP(1,:,:))  &
         + CP(3,:,:) * CONJG(CP(2,:,:))*CONJG(CP(2,:,:))*CP(1,:,:) ) &
         + DCPR_X2(1,:,:) + DCPR_X2(3,:,:) + DCPR_X2(2,:,:) + DCPI_X2(1,:,:)   & 
         + DCPI_X2(3,:,:) + DCPI_X2(2,:,:) + DCPR_Y2(1,:,:) + DCPR_Y2(3,:,:)    &
         + DCPR_Y2(2,:,:) + DCPI_Y2(1,:,:) + DCPI_Y2(3,:,:) + DCPI_Y2(2,:,:)  &
         + OMEGA *(CONJG(CP(1,:,:))* CP(2,:,:) + CONJG(CP(2,:,:))* CP(1,:,:) &
         + REAL(CONJG(CP(3,:,:))* CP(2,:,:) + CONJG(CP(2,:,:))* CP(3,:,:) ) ) &
         -TMP2D  
     EN = INTEGRATE(TMP2D, DX, DY)
!
  SELECT CASE (OPT_SO)
 
   CASE  (0)   !None
        TMP2D = 0.d0
    CASE (1)     !Rashba
      TMP2D = - GAM*CI  * (CONJG(CP(3,:,:)) + CONJG(CP(1,:,:)))*(DPY(2,:,:)) &  
              - GAM  * (-CONJG(CP(1,:,:)) + CONJG(CP(3,:,:)))*(DPX(2,:,:))  &
              - GAM*CI  * (CONJG(CP(2,:,:))) * (DPY(1,:,:) + DPY(3,:,:))     &
              + GAM  * (CONJG(CP(2,:,:))) * (-DPX(1,:,:) + DPX(3,:,:))
    CASE (2)     !Dresselhaus
      TMP2D = - GAM*CI *(CONJG(CP(3,:,:)) + CONJG(CP(1,:,:)))*(-DPY(2,:,:)) &  
              - GAM  * (-CONJG(CP(1,:,:)) + CONJG(CP(3,:,:)))*(DPX(2,:,:))  &
              + GAM*CI  * (CONJG(CP(2,:,:))) * (DPY(1,:,:) + DPY(3,:,:))     &
              + GAM  * (CONJG(CP(2,:,:))) * (-DPX(1,:,:) + DPX(3,:,:))
   CASE (3)     !Alternate 
      TMP2D = - GAM*CI  * (CONJG(CP(3,:,:)) + CONJG(CP(1,:,:)))*(DPX(2,:,:)) &  
              - GAM * (CONJG(CP(1,:,:)) - CONJG(CP(3,:,:)))*(DPY(2,:,:))  &
              - GAM*CI  * (CONJG(CP(2,:,:))) * (DPX(1,:,:) + DPX(3,:,:))     &
              + GAM  * (CONJG(CP(2,:,:))) * (DPY(1,:,:) - DPY(3,:,:))

  
  END SELECT
  EN = EN + INTEGRATE(TMP2D, DX, DY)
 
END SUBROUTINE ENERGY

FUNCTION SIMP(F, DX) RESULT (VALUE)
  IMPLICIT NONE
  REAL (8), DIMENSION(0:), INTENT(IN) :: F
  REAL (8), INTENT(IN) :: DX
  REAL (8) :: VALUE
  REAL (8) :: F1, F2
  INTEGER :: I, N
  N = SIZE(F) - 1
  F1 = F(1) + F(N-1) ! N EVEN 
  F2 = F(2)
 !$OMP PARALLEL DO PRIVATE(I,F1,F2)
  

  DO I = 3, N-3, 2
     F1 = F1 + F(I)
     F2 = F2 + F(I+1)
  END DO
  !$OMP END PARALLEL DO
  VALUE = DX * (F(0) + 4.0D0 * F1 + 2.0D0 * F2 + F(N))/3.0D0
END FUNCTION SIMP

  
