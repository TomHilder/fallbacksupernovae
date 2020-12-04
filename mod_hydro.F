c     ============================================================

      MODULE hydro_var

c     ============================================================

      IMPLICIT NONE

      SAVE

      real, dimension(-1:NX+2,-1:NY+2,-1:NZ+2,0:4) :: u,fl_u
      real, dimension(-2:NX+2,-2:NY+2,-2:NZ+2,1:3) :: b
      real, dimension(-1:NX+2,-1:NY+2,-1:NZ+2) :: c_s,p,t,eps
      real, dimension(1:NX,1:NY,1:NZ) :: s

      real :: dt,zeit=0.0d0
      real :: cfl_factor=CFL_KRIT
      real :: rho_max,rho_min1,rho_min2,eps_rand

      real :: m_core

#ifdef CONS_TEST
      real :: verlust_mtot=0.0d0,verlust_etot=0.0d0,
     &     mloss_inner=0.0d0
#endif

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
#ifndef AEQUATOR
      integer, parameter :: l_schritt=1
#else
      integer, parameter :: l_schritt=2
#endif


      real, parameter :: gamma_e  = 5.0d0/3.0d0
      real, parameter :: k_poly_e = 4.89d14
#if EOS_TYP==1
      real, parameter :: gamma_b  = 2.5d0
      real, parameter :: k_poly_b = 1.0d-2
      real, parameter :: gamma_th = 1.5d0
      real, parameter :: gamma_e_m1_r  = 1.0d0/(gamma_e -1.0d0)
      real, parameter :: gamma_b_m1_r  = 1.0d0/(gamma_b -1.0d0)
      real, parameter :: gamma_th_m1   = gamma_th-1.0d0
#endif
#if EOS_TYP==2
      real, parameter :: gamma_b=2.5d0
      real, parameter :: gamma_th=1.34d0
      real, parameter :: gamma_e_m1_r  = 1.0d0/(gamma_e -1.0d0)
      real, parameter :: gamma_b_m1_r  = 1.0d0/(gamma_b -1.0d0)
      real, parameter :: gamma_th_m1   = gamma_th-1.0d0
      real, parameter :: rho_nuc=2.0d14
      real :: p_nuc!=(k_poly_e*rho_nuc**gamma_e)
      real :: k_poly_b!=p_nuc/rho_nuc**gamma_b
      real :: eps_nuc_rho!=p_nuc*(gamma_e_m1_r-gamma_b_m1_r)/
!     +     rho_nuc
#endif

      END MODULE param_var

c     ============================================================


c     ============================================================

      MODULE gitter_var

c     ============================================================

      IMPLICIT NONE

      SAVE

      real, dimension(-3:NX+3) :: r_if
      real, dimension(-2:NX+3) :: r,r1
      real, dimension(-2:NX+2) :: dr_if,dr_if1
      real, dimension(-2:NX+3) :: dr,dr1

      real, dimension(-3:NY+3) :: theta_if
      real, dimension(-2:NY+3) :: theta
      real, dimension(-2:NY+2) :: dtheta_if,dtheta_if1
      real, dimension(-2:NY+3) :: dtheta,dtheta1
      real, dimension(-2:NY+3) :: cot_theta

      real, dimension(-3:NZ+3) :: phi_if
      real, dimension(-2:NZ+3) :: phi
      real, dimension(-2:NZ+2) :: dphi_if,dphi_if1
      real, dimension(-2:NZ+3) :: dphi,dphi1

      real, dimension(-1:NX+2,-1:NY+2,-1:NZ+2) :: dv,dv1
      real, dimension(-1:NX+2,-1:NY+2,-1:NZ+2,3) :: da 
      real, dimension(-1:NX+2,-1:NY+2,-1:NZ+2,3) :: da1 

      real, dimension(-1:NX+2,-1:NY+2,-1:NZ+2) :: r_sin_theta,
     +     r_sin_theta1

c     Symmetriefaktor fuer Integration ueber das Gesamtvolumen
#ifdef OKTANT
      real, parameter :: sym_fac=8.0d0
#else
#ifdef AEQUATOR
      real, parameter :: sym_fac=2.0d0
#else
      real, parameter :: sym_fac=1.0d0
#endif
#endif

      END MODULE gitter_var

c     ============================================================



c     ============================================================

      MODULE hydro_sav_var

c     ============================================================

      IMPLICIT NONE

      SAVE

      real, dimension(-1:NX+2,-1:NY+2,-1:NZ+2,0:4) :: u_sav

      END MODULE hydro_sav_var

c     ============================================================
