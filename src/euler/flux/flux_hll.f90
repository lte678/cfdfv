MODULE MOD_flux_hll
!===================================================================================================================================
! HLL flux
!===================================================================================================================================
! MODULES
!-----------------------------------------------------------------------------------------------------------------------------------
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
INTERFACE flux_hll
   MODULE PROCEDURE flux_hll
END INTERFACE
!-----------------------------------------------------------------------------------------------------------------------------------
PUBLIC  :: flux_hll
!===================================================================================================================================

CONTAINS

SUBROUTINE flux_hll( rho_l, rho_r, &
                      v1_l, v1_r,   &
                      v2_l, v2_r,   &
                      p_l, p_r,     &
                      flux_side     )
!===================================================================================================================================
! Computation of the HLL flux
!===================================================================================================================================
! MODULES
USE MOD_Equation_Vars, ONLY: gamma
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
REAL                        :: a_r_plus, a_l_minus
REAL                        :: c_l, c_r, e_tot_l, e_tot_r, h_l, h_r
REAL                        :: v1, v2, v_sq, h, c
REAL                        :: rho_l_root, rho_r_root
REAL                        :: k1
REAL, DIMENSION(4)          :: flux_l, flux_r
REAL, DIMENSION(4)          :: delta_u
!===================================================================================================================================

!-----------------------------------------------------------------------------------------------------------------------------------

c_l = SQRT(gamma*p_l/rho_l)
c_r = SQRT(gamma*p_r/rho_r)
e_tot_l = p_l/(gamma-1) + 0.5*rho_l*(v1_l*v1_l+v2_l*v2_l)
e_tot_r = p_r/(gamma-1) + 0.5*rho_r*(v1_r*v1_r+v2_r*v2_r)
h_l = (e_tot_l + p_l) / rho_l
h_r = (e_tot_r + p_r) / rho_r

! Calculate the necessary roe averages
rho_l_root = SQRT(rho_l)
rho_r_root = SQRT(rho_r)
k1 = rho_l_root + rho_r_root
v1 = (rho_r_root*v1_r + rho_l_root*v1_l) / k1
v2 = (rho_r_root*v2_r + rho_l_root*v2_l) / k1
v_sq = (v1*v1 + v2*v2)
h = (rho_r_root*h_r + rho_l_root*h_l) / k1
c = SQRT((gamma-1)*(h - 0.5*v_sq))

! Calculate the left and right eigen values
a_r_plus = MAX(0.0, v1_r+c_r, v1+c)
a_l_minus = MIN(0.0, v1_l-c_l, v1-c)

! Calculate left and right fluxes
flux_l(1) = rho_l*v1_l
flux_l(2) = rho_l*v1_l*v1_l + p_l
flux_l(3) = rho_l*v1_l*v2_l
flux_l(4) = v1_l*(e_tot_l + p_l)

flux_r(1) = rho_r*v1_r
flux_r(2) = rho_r*v1_r*v1_r + p_r
flux_r(3) = rho_r*v1_r*v2_r
flux_r(4) = v1_r*(e_tot_r + p_r)

! Calculate edge flux
delta_u(1) = rho_r - rho_l
delta_u(2) = rho_r*v1_r - rho_l*v1_l
delta_u(3) = rho_r*v2_r - rho_l*v2_l
delta_u(4) = e_tot_r - e_tot_l

flux_side(:) = (a_r_plus*flux_l - a_l_minus*flux_r)/(a_r_plus-a_l_minus) + ((a_r_plus*a_l_minus)/(a_r_plus-a_l_minus))*delta_u
!-----------------------------------------------------------------------------------------------------------------------------------

END SUBROUTINE flux_hll

END MODULE MOD_flux_hll
