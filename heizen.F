c     ===============================================================

      MODULE mod_heizen

c     ===============================================================

      IMPLICIT NONE

      CONTAINS

#if EOS_TYP==2

c     ===============================================================

      SUBROUTINE heizen

c     ===============================================================
      
      USE hydro_var
      USE gitter_var
      USE param_var

      IMPLICIT NONE

      real :: rho,p_e,p_th,eps_e,gamma_tmp,lum,dlum,q

      integer :: i,j,k

      do k=1,NZ
#if NDIM==2
!$OMP PARALLEL DO
!$OMP& PRIVATE(i,j)
!$OMP& SHARED(k) 
#endif
         do j=1,NY
            lum=0.0d0
            do i=1,NX
               rho=u(i,j,k,0)
               if (rho.lt.rho_nuc) then
                  p_e=k_poly_e*rho**gamma_e
                  eps_e=p_e*gamma_e_m1_r
               else
                  p_e=k_poly_b*rho**gamma_b
                  eps_e=p_e*gamma_b_m1_r+eps_nuc_rho*rho
               end if
               if (rho.gt.1.0d7) then 
                  dlum=dt/5.0d-3*max(u(i,j,k,4)-eps_e,0.0d0)*dv(i,j,k)
               else
                  dlum=0.0d0
               endif
               q=(1.0d0-exp(-dr(i)*rho/1.0d17))*lum
               u(i,j,k,4)=u(i,j,k,4)+(q-dlum)*dv1(i,j,k)
               lum=lum-q+dlum
            end do
         end do
      end do

      return

      END SUBROUTINE heizen

c     ===============================================================
#else
      SUBROUTINE heizen
      END SUBROUTINE heizen
#endif


      END MODULE mod_heizen

c     ===============================================================

