










c     ===============================================================

      MODULE poisson

c     ===============================================================

      IMPLICIT NONE

      CONTAINS



c     =============================================================== 

      SUBROUTINE gravitation_1d(modus)

c     =============================================================== 

      USE hydro_var
      USE hydro_sav_var
      USE gitter_var
      USE param_var
      USE grav_var

      IMPLICIT NONE

      integer, intent(in) :: modus

      real, dimension(0:800) :: m

      real :: px

      integer :: i,j,k

      phi_grav_alt(:,:,:)=phi_grav(:,:,:)

c      m(0)=0.0d0
      m(0)=m_core
      do i=1,800
         m(i)=m(i-1)+sum(u(i,1:1,1:1,0)*dv(i,1:1,1:1))*sym_fac
      end do

      phi_grav(800,1,1)=-m(800)/r_if(800)
      do i=800-1,1,-1
         phi_grav(i,1,1)=phi_grav(i+1,1,1)-m(i)/r_if(i)**2*dr_if(i)
      end do
      phi_grav(0,1,1)=phi_grav(1,1,1)-(phi_grav(2,1,1)-phi_grav(1,1,1))
      phi_grav(-1,1,1)=phi_grav(0,1,1)-(phi_grav(2,1,1)-phi_grav(1,1,1))
      phi_grav(800+1,1,1)=phi_grav(800,1,1)
      phi_grav(800+2,1,1)=phi_grav(800,1,1)
      phi_grav(:,1,1)=g_grav*phi_grav(:,1,1)

c      a_grv(1,1,1,1)=(-0.5d0*m(1)/r_if(1)**2*g_grav)*r(1)/r_if(1)
!$OMP PARALLEL DO
      do i=1,800
c         a_grv(i,1,1,1)=-(m(i)/r_if(i)**2*
c     +        (r_if(i)-r(i))/(r_if(i)-r_if(i-1))+
c     +        m(i-1)/r_if(i-1)**2*
c     +        (r(i)-r_if(i-1))/(r_if(i)-r_if(i-1)))*
c     +        g_grav
         a_grv(i,1,1,1)=(phi_grav(i-1,1,1)-phi_grav(i+1,1,1))/
     +        (r(i+1)-r(i-1))
      end do


      if (modus.eq.0) then
c     Korrekturterm zur Energiebilanz -> 2.Ordnung fuer zeitliche Diskretisierung
         do i=1,800
            u(i,1,1,4)=u(i,1,1,4)+0.5d0*
     +           (phi_grav(i,1,1)*u_sav(i,1,1,0)+
     +           phi_grav_alt(i,1,1)*u(i,1,1,0)-
     +           phi_grav_alt(i,1,1)*u_sav(i,1,1,0)-
     +           phi_grav(i,1,1)*u(i,1,1,0))
         end do
      end if

      return

      END SUBROUTINE gravitation_1d

c     ===============================================================


      END MODULE poisson

c     ===============================================================

