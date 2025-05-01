! util.f90
! Part of the BEC-GP-SPINOR-ROT-OMP software package 
! Program with utility subroutines
! (Fortran 90/95/2003 version)
! Developed by P. Muruganandam, A. Balaz, and S. K. Adhikari (21 January 2021)

MODULE UTIL

  INTERFACE SIMP
    MODULE PROCEDURE SIMP_R, SIMP_C
  END INTERFACE SIMP

  INTERFACE DIFF
    MODULE PROCEDURE DIFF_R, DIFF_C
  END INTERFACE DIFF

CONTAINS

PURE FUNCTION SIMP_R(F, DX) RESULT (VALUE)
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
END FUNCTION SIMP_R

PURE FUNCTION SIMP_C(F, DX) RESULT (VALUE)
  IMPLICIT NONE
  COMPLEX (8), DIMENSION(0:), INTENT(IN) :: F
  REAL (8), INTENT(IN) :: DX
  COMPLEX (8) :: VALUE
  COMPLEX (8) :: F1, F2
  INTEGER :: I, N

  N = SIZE(F) - 1
  F1 = F(1) + F(N-1) ! N EVEN 
  F2 = F(2)
  DO I = 3, N-3, 2
     F1 = F1 + F(I)
     F2 = F2 + F(I+1)
  END DO
  VALUE = DX*(F(0) + 4.0D0*F1 + 2.0D0*F2 + F(N))/3.0D0
END FUNCTION SIMP_C

PURE FUNCTION DIFF_R(P,DX) RESULT (DP)
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
END FUNCTION DIFF_R

PURE FUNCTION DIFF_C(P,DX) RESULT (DP)
! Computes the first derivative DP of P using
! Richardsonextrapolation formula. The derivative at the  
! boundaries are assumed to be zero
  IMPLICIT NONE
  COMPLEX (8), DIMENSION(0:), INTENT(IN) :: P
  REAL (8), INTENT(IN) :: DX
  COMPLEX (8), DIMENSION(0:SIZE(P)-1) :: DP
  INTEGER :: I, N

  N = SIZE(P) - 1
  DP(0) = 0.0D0
  DP(1) = (P(2) - P(0))/(2.0D0*DX)
  FORALL(I=2:N-2)
    DP(I) = (P(I-2)-8.0D0*P(I-1)+8.0D0*P(I+1)-P(I+2))/(12.0D0*DX)
  END FORALL
  DP(N-1) = (P(N) - P(N-2))/(2.0D0*DX)
  DP(N) = 0.0D0
END FUNCTION DIFF_C

END MODULE UTIL


