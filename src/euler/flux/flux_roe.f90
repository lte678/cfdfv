MODULE MOD_flux_roe
!===================================================================================================================================
! Roe flux
!===================================================================================================================================
! MODULES
!-----------------------------------------------------------------------------------------------------------------------------------
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
INTERFACE flux_roe
   MODULE PROCEDURE flux_roe
END INTERFACE
!-----------------------------------------------------------------------------------------------------------------------------------
PUBLIC  :: flux_roe
!===================================================================================================================================

CONTAINS

SUBROUTINE flux_roe( rho_l, rho_r, &
                     v1_l, v1_r,   &
                     v2_l, v2_r,   &
                     p_l, p_r,     &
                     flux_side     )
!===================================================================================================================================
! Calculation of Roe flux
!===================================================================================================================================
! MODULES
USE MOD_Equation_Vars, ONLY: gamma,gamma1q
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
REAL                        :: v1, v2, v_sq, h, c      ! Not the values at the cell boundary. Rather: Roe averages
REAL                        :: e_tot_l, e_tot_r        ! Left and right total energy
REAL                        :: h_l, h_r                ! Left and right enthalpy
REAL                        :: rho_l_root, rho_r_root  ! Precomputed sqrts, since they are expensive
REAL                        :: a(4)                    ! Eigen values of linearized equations
REAL                        :: r(4, 4)                 ! Eigen vectors of linearized equations
real                        :: g(4)                    ! Eigenvector coefficients
REAL                        :: k1, k2, k3, k4, k5      ! Generic precomputed constants
!===================================================================================================================================

!-----------------------------------------------------------------------------------------------------------------------------------
rho_l_root = SQRT(rho_l)
rho_r_root = SQRT(rho_r)
k1 = rho_l_root + rho_r_root
k2 = gamma-1
k3 = rho_r - rho_l         ! Delta rho
k4 = rho_r*v1_r-rho_l*v1_l ! Delta m_1
k5 = rho_r*v2_r-rho_l*v2_l ! Delta m_2

e_tot_l = p_l/k2 + 0.5*rho_l*(v1_l*v1_l+v2_l*v2_l)
e_tot_r = p_r/k2 + 0.5*rho_r*(v1_r*v1_r+v2_r*v2_r)
h_l = (e_tot_l + p_l) / rho_l
h_r = (e_tot_r + p_r) / rho_r

! Calculate average quantities
v1 = (rho_r_root*v1_r + rho_l_root*v1_l) / k1
v2 = (rho_r_root*v2_r + rho_l_root*v2_l) / k1
v_sq = (v1*v1 + v2*v2)
h = (rho_r_root*h_r + rho_l_root*h_l) / k1
c = SQRT(k2*(h - 0.5*v_sq))

! Calculate the eigen values
a(1) = v1 - c
a(2) = v1
a(3) = v1
a(4) = v1 + c

! Calculate the eigen vectors
r(:,1) = (/ 1.0 , v1-c, v2  , h-v1*c   /)
r(:,2) = (/ 1.0 , v1  , v2  , 0.5*v_sq /)
r(:,3) = (/ 0.0 , 0.0 , 1.0 , v2       /)
r(:,4) = (/ 1.0 , v1+c, v2  , h+v1*c   /)

! Calculate the eigen vector coefficients
g(2) = -(k2/(c*c)) * (k3*(v1*v1-h) - v1*k4 + (e_tot_r - e_tot_l) - (k5-v2*k3)*v2)
g(1) = -(1/(2*c))*(k4-k3*(v1+c)) - 0.5*g(2)
g(3) = k5 - v2*k3
g(4) = k3 - g(1) - g(2)

! Calculate the final flux
flux_side(1) = 0.5*rho_l*v1_l + 0.5*rho_r*v1_r
flux_side(2) = 0.5*(rho_l*v1_l*v1_l + p_l) + 0.5*(rho_r*v1_r*v1_r + p_r)
flux_side(3) = 0.5*rho_l*v1_l*v2_l + 0.5*rho_r*v1_r*v2_r
flux_side(4) = 0.5*v1_l*(e_tot_l + p_l) + 0.5*v1_r*(e_tot_r + p_r)

flux_side(:) = flux_side(:) - 0.5*MATMUL(g*ABS(a),TRANSPOSE(r))
!-----------------------------------------------------------------------------------------------------------------------------------

END SUBROUTINE flux_roe

END MODULE MOD_flux_roe
