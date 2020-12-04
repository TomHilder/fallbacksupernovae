










c     ============================================================

      MODULE hydro_var

c     ============================================================

      IMPLICIT NONE

      SAVE

      real, dimension(-1:800+2,-1:1+2,-1:1+2,0:4) :: u,fl_u
      real, dimension(-2:800+2,-2:1+2,-2:1+2,1:3) :: b
      real, dimension(-1:800+2,-1:1+2,-1:1+2) :: c_s,p,t,eps
      real, dimension(1:800,1:1,1:1) :: s

      real :: dt,zeit=0.0d0
      real :: cfl_factor=0.4d0
      real :: rho_max,rho_min1,rho_min2,eps_rand

      real :: m_core

      real :: verlust_mtot=0.0d0,verlust_etot=0.0d0,
     &     mloss_inner=0.0d0

      END MODULE hydro_var

c     ============================================================


c     ============================================================

      MODULE param_var

c     ============================================================

      IMPLICIT NONE

      real, parameter :: k_b=1.38d-16
      real, parameter :: r_gas=8.31d7
      real, parameter :: mu=1.0d0
      real, parameter :: pi=3.141592653589793238462643383279d0
      real, parameter :: g_grav=6.6742d-8
      real, parameter :: zwei_pi_g=2.0d0*pi*g_grav
      real, parameter :: vier_pi_g=4.0d0*pi*g_grav
      real, parameter :: msun=1.98847d33

      integer, parameter :: l_max=10
      integer, parameter :: l_schritt=1


      real, parameter :: gamma_e  = 5.0d0/3.0d0
      real, parameter :: k_poly_e = 4.89d14

      END MODULE param_var

c     ============================================================


c     ============================================================

      MODULE gitter_var

c     ============================================================

      IMPLICIT NONE

      SAVE

      real, dimension(-3:800+3) :: r_if
      real, dimension(-2:800+3) :: r,r1
      real, dimension(-2:800+2) :: dr_if,dr_if1
      real, dimension(-2:800+3) :: dr,dr1

      real, dimension(-3:1+3) :: theta_if
      real, dimension(-2:1+3) :: theta
      real, dimension(-2:1+2) :: dtheta_if,dtheta_if1
      real, dimension(-2:1+3) :: dtheta,dtheta1
      real, dimension(-2:1+3) :: cot_theta

      real, dimension(-3:1+3) :: phi_if
      real, dimension(-2:1+3) :: phi
      real, dimension(-2:1+2) :: dphi_if,dphi_if1
      real, dimension(-2:1+3) :: dphi,dphi1

      real, dimension(-1:800+2,-1:1+2,-1:1+2) :: dv,dv1
      real, dimension(-1:800+2,-1:1+2,-1:1+2,3) :: da 
      real, dimension(-1:800+2,-1:1+2,-1:1+2,3) :: da1 

      real, dimension(-1:800+2,-1:1+2,-1:1+2) :: r_sin_theta,
     +     r_sin_theta1

c     Symmetriefaktor fuer Integration ueber das Gesamtvolumen
      real, parameter :: sym_fac=1.0d0

      END MODULE gitter_var

c     ============================================================



c     ============================================================

      MODULE hydro_sav_var

c     ============================================================

      IMPLICIT NONE

      SAVE

      real, dimension(-1:800+2,-1:1+2,-1:1+2,0:4) :: u_sav

      END MODULE hydro_sav_var

c     ============================================================
