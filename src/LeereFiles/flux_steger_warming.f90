MODULE MOD_flux_steger_warming
!===================================================================================================================================
! Steger Warming
!===================================================================================================================================
! MODULES
IMPLICIT NONE
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
INTERFACE flux_steger_warming
   MODULE PROCEDURE flux_steger_warming
END INTERFACE
!-----------------------------------------------------------------------------------------------------------------------------------
PUBLIC  :: flux_steger_warming
!===================================================================================================================================

CONTAINS

SUBROUTINE flux_steger_warming( rho_l, rho_r, &
                                v1_l, v1_r,   &
                                v2_l, v2_r,   &
                                p_l, p_r,     &
                                flux_side     )
!===================================================================================================================================
! Steger Warming
!===================================================================================================================================
! MODULES
USE MOD_Equation_Vars, ONLY: gamma,gamma1,gamma1q
!-----------------------------------------------------------------------------------------------------------------------------------
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)             :: rho_l, rho_r
REAL,INTENT(IN)             :: v1_l , v1_r
REAL,INTENT(IN)             :: v2_l , v2_r
REAL,INTENT(IN)             :: p_l  , p_r
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)            :: flux_side(4)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
    
    !local variable declaration
    
!===================================================================================================================================

!-----------------------------------------------------------------------------------------------------------------------------------
   
   ! Insert your Code here
   
!-----------------------------------------------------------------------------------------------------------------------------------

END SUBROUTINE flux_steger_warming

END MODULE MOD_flux_steger_warming
