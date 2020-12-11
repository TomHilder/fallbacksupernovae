      MODULE mod_eos

      IMPLICIT NONE

      CONTAINS

      SUBROUTINE eos(rho,eps,p,c_s,modus)

      USE param_var

      IMPLICIT NONE

      real, intent(in) :: rho
      real, intent(in) :: eps
      real, intent(out) :: p,c_s
      integer, intent(in) :: modus

#if EOS_TYP==1
      real :: p_e,p_b,p_th,eps_e,eps_b
#endif

#if EOS_TYP==2
      real :: p_e,p_th,eps_e,gamma_tmp
#endif

#if EOS_TYP==0
c      p=2.0d0/3.0d0*eps
c      c_s=sqrt(5.0d0/3.0d0*p/rho)
      p=(gamma_e-1.0d0)*eps
      c_s=sqrt(gamma_e*p/rho)
#endif
#if EOS_TYP==1
      p_e=k_poly_e*rho**gamma_e
      eps_e=p_e*gamma_e_m1_r
      p_b=k_poly_b*rho**gamma_b
      eps_b=p_b*gamma_b_m1_r
      p_th=max(eps-eps_e-eps_b,0.0d0)*gamma_th_m1
      p=p_e+p_b+p_th
      c_s=sqrt((gamma_e*p_e+gamma_b*p_b+gamma_th*p_th)/rho)
#endif
#if EOS_TYP==2
      if (rho.lt.rho_nuc) then
         p_e=k_poly_e*rho**gamma_e
         eps_e=p_e*gamma_e_m1_r
         gamma_tmp=gamma_e
      else
         p_e=k_poly_b*rho**gamma_b
         eps_e=p_e*gamma_b_m1_r+eps_nuc_rho*rho
         gamma_tmp=gamma_b
      end if
      p_th=max(eps-eps_e,0.0d0)*gamma_th_m1
      p=p_e+p_th
      c_s=sqrt((gamma_tmp*p_e+gamma_th*p_th)/rho)
#endif


      END SUBROUTINE eos

      END MODULE mod_eos
