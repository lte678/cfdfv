MODULE MOD_Boundary
!===================================================================================================================================
! Boundaries
!===================================================================================================================================
! MODULES
!-----------------------------------------------------------------------------------------------------------------------------------
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
INTERFACE InitBoundary
   MODULE PROCEDURE InitBoundary
END INTERFACE
INTERFACE SetBCatSides
   MODULE PROCEDURE SetBCatSides
END INTERFACE
INTERFACE SetBCatBarys
   MODULE PROCEDURE SetBCatBarys
END INTERFACE
INTERFACE Boundary
   MODULE PROCEDURE Boundary
END INTERFACE
!-----------------------------------------------------------------------------------------------------------------------------------
PUBLIC  :: InitBoundary, Boundary, SetBCatSides, SetBCatBarys
!===================================================================================================================================

CONTAINS


SUBROUTINE InitBoundary()
!===================================================================================================================================
! Init Boundaries
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Readintools  ,ONLY:GETINT,GETLOGICAL,GETREAL,GETREALARRAY,CNTSTR
USE MOD_Boundary_Vars,ONLY:tBoundary
USE MOD_Boundary_Vars,ONLY:nBC,FirstBC
USE MOD_Equation_Vars,ONLY:gamma
!-----------------------------------------------------------------------------------------------------------------------------------
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE 
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
INTEGER                  :: iBC, Int_in
REAL                     :: M, c, v
REAL                     :: alpha
TYPE(tBoundary), POINTER :: BC
!===================================================================================================================================

WRITE(*,*)
WRITE(*,*) '-[Boundary Conditions]--------------------------------------'
! number of BC
nBC = GETINT('nBC')
!WRITE(*,'(a, I2)') '  Number of Boundary Conditions: ', nBC

FirstBC => NULL()

DO iBC = 1, nBC
  !WRITE(*,'(a, I2, a)') '    Boundary Condition # ', iBC, ':'
! Create new Boundary Condition
  ALLOCATE(BC)
  BC%nextBC    => FirstBC
  FirstBC => BC
! Read boundary type and ID
  Int_in = GETINT('BCtype')
  BC%BCType = Int_in / 100
  BC%BCID   = MOD(Int_in, 100)
  !WRITE(*,'(a, I4)') '      BC Code: ', BC%BCID + BC%BCType*100
! Read additional information
  SELECT CASE(BC%BCType)
  CASE(SLIPWALL)
    WRITE(*,'(a)') '|      BC Type: Slipwall'
  ! Nothing to do
  CASE(INFLOW)
    WRITE(*,'(a, I2, a)') '|      BC Type: Inflow'
  ! Density
    BC%pvar(RHO) = GETREAL('rho')
    !WRITE(*,'(a, F10.4)') '        Density: ', BC%pvar(RHO)
  ! Mach number
    M = GETREAL('Mach')
    !WRITE(*,'(a, F10.4)') '        Mach Number: ', M
  ! angle of attack
    alpha = GETREAL('alpha')
    !WRITE(*,'(a, F10.4)') '        Angle of Attack: ', alpha
  ! pressure
    BC%pvar(P) = GETREAL('pressure')
    !WRITE(*,'(a, F10.4)') '        Pressure: ', BC%pvar(P)

    c = SQRT(gamma * BC%pvar(P) / BC%pvar(RHO))
    v = M * c
    BC%pvar(V1) = v * COS(alpha * ACOS(-1.) / 180.)
    BC%pvar(V2) = v * SIN(alpha * ACOS(-1.) / 180.)
  CASE(CHARACTERISTIC)
    WRITE(*,'(a, I2, a)') '|      BC Type: Characteristic'
  ! Density
    BC%pvar(RHO) = GETREAL('rho')
    !WRITE(*,'(a, F10.4)') '        Density: ', BC%pvar(RHO)
  ! Mach number
    M = GETREAL('Mach')
    !WRITE(*,'(a, F10.4)') '        Mach Number: ', M
  ! angle of attack
    alpha = GETREAL('alpha')
    !WRITE(*,'(a, F10.4)') '        Angle of Attack: ', alpha
  ! pressure
    BC%pvar(P) = GETREAL('pressure')
    !WRITE(*,'(a, F10.4)') '        Pressure: ', BC%pvar(P)
    c = SQRT(gamma * BC%pvar(P) / BC%pvar(RHO))
    v = M * c
    BC%pvar(V1) = v * COS(alpha * ACOS(-1.) / 180.)
    BC%pvar(V2) = v * SIN(alpha * ACOS(-1.) / 180.)
  CASE(OUTFLOW)
    WRITE(*,'(a, I2, a)') '|      BC Type: Outflow'
  ! Nothing to do
  CASE(PRESSURE_OUT)
    WRITE(*,'(a, I2, a)') '|      BC Type: Pressure Out'
  ! pressure
    BC%pvar(P) = GETREAL('pressure')
    !WRITE(*,'(a, F10.4)') '        Pressure: ', BC%pvar(P)
  CASE(EXACTSOL)
    WRITE(*,'(a, I2, a)') '|      BC Type: Exact Function'
  ! read exact function
    BC%ExactFunc = GETINT('BCExactFunc')
  CASE(PERIODIC)
    WRITE(*,'(a, I2, a)') '|      BC Type: Periodic'
  ! read exact function
    BC%connection(:) = GETREALARRAY('connection',2)
  END SELECT
ENDDO
END SUBROUTINE InitBoundary


SUBROUTINE SetBCatSides(time)
!===================================================================================================================================
! sets the ghost values at sides or elems
!===================================================================================================================================
! MODULES
USE MOD_Boundary_Vars,ONLY:
USE MOD_Mesh_Vars,    ONLY: nBCSides,BCSides
USE MOD_Mesh_Vars,    ONLY: tElem,tSide
!-----------------------------------------------------------------------------------------------------------------------------------
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)             :: time
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
TYPE(tElem), POINTER        :: aElem, gElem
TYPE(tSide), POINTER        :: aSide, gSide
INTEGER                     :: iSide
!===================================================================================================================================
!$omp parallel do private(gSide,gElem,aSide,aElem)

DO iSide = 1,nBCSides
  gSide => BCSides(iSide)%Side  ! Edge of the ghost cell
  aSide => gSide%connection     ! Edge of the solution cell
  gElem => gSide%Elem           ! Ghost cell
  CALL Boundary(gSide, time, aSide%pvar, gSide%pvar, gElem%Bary + gSide%GP)
ENDDO

!$omp end parallel do
!-----------------------------------------------------------------------------------------------------------------------------------
END SUBROUTINE


SUBROUTINE SetBCatBarys(time)
!===================================================================================================================================
! sets the ghost values at sides or elems
!===================================================================================================================================
! MODULES
USE MOD_Mesh_Vars,    ONLY: nBCSides,BCSides
USE MOD_Mesh_Vars,    ONLY: tElem,tSide
!-----------------------------------------------------------------------------------------------------------------------------------
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)             :: time
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
TYPE(tElem), POINTER        :: aElem, gElem
TYPE(tSide), POINTER        :: aSide, gSide
INTEGER                     :: iSide
!===================================================================================================================================
!$omp parallel do private(gSide,gElem,aSide,aElem)

DO iSide = 1,nBCSides
  gSide => BCSides(iSide)%Side  ! Edge of the ghost cell
  aSide => gSide%connection     ! Edge of the solution cell
  gElem => gSide%Elem           ! Ghost cell
  aElem => aSide%Elem           ! Solution cell
  CALL Boundary(gSide, time, aElem%pvar, gElem%pvar, gElem%Bary)
ENDDO

!$omp end parallel do
!---------------------------------------------------------------------------!
END SUBROUTINE


SUBROUTINE Boundary(aSide, time, int_pvar, ghost_pvar, x)
!===================================================================================================================================
! set boundary value at x
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_ExactFunc
USE MOD_EOS
USE MOD_Mesh_Vars    ,ONLY:tSide
USE MOD_Equation_Vars,ONLY:gamma
!-----------------------------------------------------------------------------------------------------------------------------------
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)                 :: time
REAL,INTENT(IN)                 :: x(2)
REAL,INTENT(IN)                 :: int_pvar(NVAR)
TYPE(tSide),POINTER,INTENT(IN)  :: aSide
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)                :: ghost_pvar(NVAR)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
REAL                            :: v_loc(2)
REAL                            :: c, v, pres, n(2)
REAL                            :: int_cvar(NVAR), ghost_cvar(NVAR)
REAL                            :: int_pvar_loc(NVAR), ghost_pvar_loc(NVAR)
REAL                            :: int_charvar(3), ghost_charvar(3) 
!===================================================================================================================================

!-----------------------------------------------------------------------------------------------------------------------------------
! Extract normal vector
n = aSide%connection%n
!-----------------------------------------------------------------------------------------------------------------------------------
! Determine type of boundary condition
SELECT CASE(aSide%BC%BCType)
!-----------------------------------------------------------------------------------------------------------------------------------
! Slipwall Boundary
CASE(SLIPWALL)
  ! Mirror the velocity normal to the cell boundary
  ghost_pvar(2:3) = int_pvar(2:3) - 2*n*DOT_PRODUCT(int_pvar(2:3), n)
  ! Set the scalar quantities to their original value
  ghost_pvar(1) = int_pvar(1)
  ghost_pvar(4) = int_pvar(4)
!-----------------------------------------------------------------------------------------------------------------------------------
! Supersonic Inflow Boundary
CASE(INFLOW)
  ! Use the provided state
  ghost_pvar(:) = aSide%BC%pvar(:)
!-----------------------------------------------------------------------------------------------------------------------------------
! Supersonic Outflow Boundary
CASE(OUTFLOW)
  ! Reusing the internal state, makes sure that the state is supersonic in all parts of this flow
  ! (Otherwise the intermediate states could potentially become subsonic!)
  ghost_pvar(:) = int_pvar(:)
!-----------------------------------------------------------------------------------------------------------------------------------
! Pressure Outflow 
CASE(PRESSURE_OUT)
  c  = SQRT(gamma * int_pvar(P) / int_pvar(RHO))
  v  = n(1) * int_pvar(V1) + n(2) * int_pvar(V2)
  IF (v/c < 1)  THEN
    pres = aSide%BC%pvar(P)
  ELSE
    pres = int_pvar(P)
  END IF 
  ghost_pvar      = int_pvar
  ghost_pvar(P)   = pres
  ghost_pvar(RHO) = int_pvar(RHO)*pres/int_pvar(P)
!-----------------------------------------------------------------------------------------------------------------------------------
! Characteristic Boundary Condition
CASE(CHARACTERISTIC)
! Compute Eigenvalues of the ghost cell:
  c = SQRT(gamma * aSide%BC%pvar(P) / aSide%BC%pvar(RHO))
  v = n(1) * aSide%BC%pvar(V1) + n(2) * aSide%BC%pvar(V2)
! Rotate primitve state into local coordinate system
  int_pvar_loc(:)  = int_pvar(:)
  int_pvar_loc(V1) =  n(1) * int_pvar(V1) + n(2) * int_pvar(V2)
  int_pvar_loc(V2) = -n(2) * int_pvar(V1) + n(1) * int_pvar(V2) 
  ghost_pvar_loc(:)  = aSide%BC%pvar
  ghost_pvar_loc(V1) =  n(1) * aSide%BC%pvar(V1) + n(2) * aSide%BC%pvar(V2)
  ghost_pvar_loc(V2) = -n(2) * aSide%BC%pvar(V1) + n(1) * aSide%BC%pvar(V2) 
! Compute conservative variables of both cells
  CALL PrimCons(int_pvar_loc, int_cvar)
  CALL PrimCons(ghost_pvar_loc, ghost_cvar)
! Compute characteristic variables of inner and ghost cell
  CALL ConsChar(int_charvar, int_cvar, int_pvar)
  CALL ConsChar(ghost_charvar, ghost_cvar, int_pvar)
! Determine characteristic state at boundary
  IF (v+c > 0.) THEN
    ghost_charvar(3) = int_charvar(3)
  END IF
  IF (v > 0.) THEN
    ghost_charvar(2) = int_charvar(2)
  END IF
  IF (v-c > 0.) THEN
    ghost_charvar(1) = int_charvar(1)
  END IF
! Determine the conservative state of the ghost cell
  CALL CharCons(ghost_charvar, ghost_cvar, int_pvar)
  IF (v > 0.) THEN
    ghost_cvar(M2) = int_cvar(M2) 
  END IF
! Determine the primitive state of the ghost cell
  CALL ConsPrim(ghost_pvar, ghost_cvar)
! Rotate the primitive state into the global coordinate system
  v_loc(1) = ghost_pvar(V1)
  v_loc(2) = ghost_pvar(V2)
  ghost_pvar(V1) = n(1) * v_loc(1) - n(2) * v_loc(2)
  ghost_pvar(V2) = n(2) * v_loc(1) + n(1) * v_loc(2)
!-----------------------------------------------------------------------------------------------------------------------------------
! Exact Solution
CASE(EXACTSOL)
  CALL ExactFunc(aSide%BC%ExactFunc, x(:), time, ghost_pvar(RHO:P))
!-----------------------------------------------------------------------------------------------------------------------------------
! Error Handler
CASE DEFAULT
  WRITE (*,*) 'ERROR in SetBoundaryConditions.f90:'
  WRITE (*,*) 'Boundary Condition Type ', aSide%BC%BCType, ' unknown.'
  WRITE (*,*) 'Aborting...'
  STOP
END SELECT
!-----------------------------------------------------------------------------------------------------------------------------------
END SUBROUTINE Boundary

END MODULE MOD_Boundary
