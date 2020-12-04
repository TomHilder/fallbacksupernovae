










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



c      p=2.0d0/3.0d0*eps
c      c_s=sqrt(5.0d0/3.0d0*p/rho)
      p=(gamma_e-1.0d0)*eps
      c_s=sqrt(gamma_e*p/rho)


      END SUBROUTINE eos

      END MODULE mod_eos
