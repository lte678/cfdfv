MODULE MOD_Limiter
!===================================================================================================================================
! Module containing the limiter
!===================================================================================================================================
! MODULES
!-----------------------------------------------------------------------------------------------------------------------------------
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
INTERFACE Limiter
   MODULE PROCEDURE Limiter
END INTERFACE
INTERFACE Limiter_BarthJespersen
   MODULE PROCEDURE Limiter_BarthJespersen
END INTERFACE
INTERFACE Limiter_Venkatakrishnan
   MODULE PROCEDURE Limiter_Venkatakrishnan
END INTERFACE
!-----------------------------------------------------------------------------------------------------------------------------------
PUBLIC   :: Limiter
PRIVATE  :: Limiter_Venkatakrishnan, Limiter_BarthJespersen
!===================================================================================================================================

CONTAINS

SUBROUTINE Limiter(aElem)
!===================================================================================================================================
! Select Limiter
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Reconstruction_Vars,  ONLY:intLimiter
USE MOD_Mesh_Vars,            ONLY:tElem
!-----------------------------------------------------------------------------------------------------------------------------------
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
TYPE(tElem), POINTER    :: aElem
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
!===================================================================================================================================

! Determine uMax and uMin
SELECT CASE(intLimiter)
CASE (BARTHJESPERSEN)
  CALL Limiter_BarthJespersen(aElem)
CASE (VENKATAKRISHNAN)
  CALL Limiter_Venkatakrishnan(aElem)
CASE DEFAULT
  WRITE (*,*) ' ERROR in Limiter.f90: Limiter function unknown.'
  STOP
END SELECT
!-----------------------------------------------------------------------------------------------------------------------------------
END SUBROUTINE Limiter



SUBROUTINE Limiter_BarthJespersen(aElem)
!===================================================================================================================================
! Limiter after Barth & Jespersen
! 2D, unstructured Limiter
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Mesh_Vars,  ONLY:tElem,tSide
!-----------------------------------------------------------------------------------------------------------------------------------
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
TYPE(tElem), POINTER    :: aElem
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
REAL                    :: phi(NVAR), phiLoc(NVAR)
REAL                    :: uMax(NVAR), uMin(NVAR)
REAL                    :: MaxDiff(NVAR), MinDiff(NVAR), uDiff(NVAR)
REAL                    :: MaxDiff_sq(NVAR), MinDiff_sq(NVAR), uDiff_sq(NVAR)
INTEGER                 :: iVar
TYPE(tSide), POINTER    :: aSide
!===================================================================================================================================

  ! this routine computes the limited gradients aElem%u_x and aElem%u_y
  ! Get the min and max values of the neighbouring cells
  uMin = aElem%pvar
  uMax = aElem%pvar
  aSide => aElem%firstSide
  DO WHILE(ASSOCIATED(aSide))
    uMax(:) = MAX(uMax, aSide%connection%Elem%pvar)
    uMin(:) = MIN(uMin, aSide%connection%Elem%pvar)
    aSide => aSide%nextElemSide
  END DO

  MaxDiff(:) = uMax - aElem%pvar
  MinDiff(:) = uMin - aElem%pvar

  phi(:) = 1.0
  aSide => aElem%firstSide
  DO WHILE(ASSOCIATED(aSide))
    uDiff = aElem%u_x * aSide%GP(X_DIR) + aElem%u_y * aSide%GP(Y_DIR) 
    phiLoc(:) = 0
    DO iVar = 1, NVAR
      IF(uDiff(iVar) > 0.0) THEN
        phiLoc(iVar) = MIN(1.0, MaxDiff(iVar)/uDiff(iVar))
      ELSEIF(uDiff(iVar) < 0.0) THEN
        phiLoc(iVar) = MIN(1.0, MinDiff(iVar)/uDiff(iVar))
      END IF
    END DO

    ! Take the minimum scaler
    phi(:) = MIN(phi, phiLoc)
    aSide => aSide%nextElemSide
  END DO

  ! Apply limit
  aElem%u_x(:) = aElem%u_x * phi
  aElem%u_y(:) = aElem%u_y * phi

END SUBROUTINE Limiter_BarthJespersen


SUBROUTINE Limiter_Venkatakrishnan(aElem)
!===================================================================================================================================
! Venkatakrishnan limiter, additionally a limiting parameter k
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Mesh_Vars,ONLY:tElem,tSide
!-----------------------------------------------------------------------------------------------------------------------------------
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
TYPE(tElem), POINTER    :: aElem
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
REAL                    :: phi(NVAR), phiLoc(NVAR)
REAL                    :: uMax(NVAR), uMin(NVAR)
REAL                    :: MaxDiff(NVAR), MinDiff(NVAR), uDiff(NVAR)
REAL                    :: MaxDiff_sq(NVAR), MinDiff_sq(NVAR), uDiff_sq(NVAR)
INTEGER                 :: iVar
TYPE(tSide), POINTER    :: aSide
!===================================================================================================================================

! this routine computes the limited gradients aElem%u_x and aElem%u_y

  ! variables with reqired data:
  ! epsilon^2: aElem%venk_epsilon_sq

  ! Get the min and max values of the neighbouring cells
  uMin = aElem%pvar
  uMax = aElem%pvar
  aSide => aElem%firstSide
  DO WHILE(ASSOCIATED(aSide))
    uMax(:) = MAX(uMax, aSide%connection%Elem%pvar)
    uMin(:) = MIN(uMin, aSide%connection%Elem%pvar)
    aSide => aSide%nextElemSide
  END DO

  MaxDiff(:) = uMax - aElem%pvar
  MinDiff(:) = uMin - aElem%pvar
  MaxDiff_sq(:) = MaxDiff * MaxDiff
  MinDiff_sq(:) = MinDiff * MinDiff

  phi(:) = 1.0
  aSide => aElem%firstSide
  DO WHILE(ASSOCIATED(aSide))
    uDiff(:) = aElem%u_x * aSide%GP(X_DIR) + aElem%u_y * aSide%GP(Y_DIR)
    uDiff_sq(:) = uDiff * uDiff
    phiLoc(:) = 0
    DO iVar = 1, NVAR
      IF(uDiff(iVar) > 0.0) THEN
        phiLoc(iVar) = 1. / uDiff(iVar) * &
          (((MaxDiff_sq(iVar) + aElem%venk_epsilon_sq) * uDiff(iVar) + &
          2. * uDiff_sq(iVar) * MaxDiff(iVar))                       / &
          (MaxDiff_sq(iVar) + 2. * uDiff_sq(iVar) + uDiff(iVar)      * &
          MaxDiff(iVar) + aElem%venk_epsilon_sq))
      ELSEIF(uDiff(iVar) < 0.0) THEN
        phiLoc(iVar) = 1. / uDiff(iVar) * &
          (((minDiff_sq(iVar) + aElem%venk_epsilon_sq) * uDiff(iVar) + &
          2. * uDiff_sq(iVar) * minDiff(iVar))                       / &
          (minDiff_sq(iVar) + 2. * uDiff_sq(iVar) + uDiff(iVar)      * &
          minDiff(iVar) + aElem%venk_epsilon_sq))
      END IF
    END DO

    ! Take the minimum scaler
    phi(:) = MIN(phi, phiLoc)
    aSide => aSide%nextElemSide
  END DO

  ! Apply limit
  aElem%u_x(:) = aElem%u_x * phi
  aElem%u_y(:) = aElem%u_y * phi
!-----------------------------------------------------------------------------------------------------------------------------------
END SUBROUTINE Limiter_Venkatakrishnan

END MODULE MOD_Limiter
