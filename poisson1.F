c     ===============================================================

      MODULE poisson

c     ===============================================================

      IMPLICIT NONE

      CONTAINS

c     ===============================================================

      SUBROUTINE comp_legendre

c     ===============================================================
c     Berechnet die Legendre-Polynome P_l(cos(theta_if)) an den 
c     Zellraendern, den Zellzentren, sowie die Integrale von P_l
c     ueber eine Zelle in theta-Richtung.
c     ---------------------------------------------------------------

      USE size_cfc
      USE legendre
      USE phycon

      IMPLICIT NONE

      integer :: j,l

      p_leg(:,0)=1.0d0
      p_leg(:,0)=1.0d0
      p_leg(0:n,1)=cos(theta_if(0:n))
      p_leg_c(0:n,1)=cos(theta_if(0:n))
      if (theta_if(n).lt.0.5*pi+1.0d-6) p_leg(n,1)=0.0d0 !equatorial symmetry

      p1_leg(0:n,1)=-sin(theta_if(0:n))
      p1_leg(0:n,2)=-3*sin(theta_if(0:n))*cos(theta_if(0:n))

c     Rekursionsformel fuer die Legendre-Polynome:
c     (Jackson, Classical Electrodynamics, aeq. 3.29):
      do l=1,l_max
         p_leg  (:,l+1)=(p_leg  (:,l)*p_leg  (:,1)*(2*l+1)-
     +        l*p_leg  (:,l-1))/dble(l+1)
         p_leg_c(:,l+1)=(p_leg_c(:,l)*p_leg_c(:,1)*(2*l+1)-
     +        l*p_leg_c(:,l-1))/dble(l+1)
      end do

      do l=2,l_max
         p1_leg  (:,l+1)=(p1_leg  (:,l)*p_leg  (:,1)*(2*l+1)-
     +        (l+1)*p1_leg  (:,l-1))/dble(l)
         p1_leg_c(:,l+1)=(p1_leg_c(:,l)*p_leg_c(:,1)*(2*l+1)-
     +        (l+1)*p1_leg_c(:,l-1))/dble(l)
      end do


c     Exakte Werte an den Raendern setzen (nicht unbedingt noetig):
      do l=0,l_max+1
         if (mod(l,2).eq.0) then
            p_leg(0,l)=1.0d0
         else
            p_leg(0,l)=-1.0d0
         end if
#ifndef AEQUATOR
         p_leg(n,l)=1.0d0
#endif
      end do

c     Integrale ueber theta-Zelle mit Hilfe der Formel fuer die
c     Ableitungen der Legendre-Polynome (ibid., aeq.3.28)
      p_leg_int(1:n,0)=p_leg(1:n,1)-p_leg(0:n-1,1)
      do l=1,l_max
         p_leg_int(1:n,l)=
     +        ((p_leg(1:n  ,l+1)-p_leg(1:n  ,l-1))-
     +         (p_leg(0:n-1,l+1)-p_leg(0:n-1,l-1)))/
     +        dble(2*l+1)
#ifdef AEQUATOR
         p_leg_int(1:n,l)=2.0d0*p_leg_int(1:n,l)
#endif
         print *,l,sum(p_leg_int(1:n,l)*cos(theta(1:n))**1)
      end do

      p1_leg_int(1:n,1)=
     +     -(sin(2*theta_if(1:n))-sin(2*theta_if(0:n-1)))/4.0+
     +     (theta_if(1:n)-theta_if(0:n-1))/2.0
      p1_leg_int(1:n,2)=
     +     sin(theta_if(1:n))**3-sin(theta_if(0:n-1))**3

      do l=3,l_max
         p1_leg_int(1:n,l)=
     +        (l-2)*l/dble((l+1)*(l-1))*
     +        p1_leg_int(1:n,l-2)-
     +        (2*l-1)/dble((l+1)*(l-1))*
     +        (sin(theta_if(1:n  ))**2*p1_leg(1:n  ,l-1)-
     +         sin(theta_if(0:n-1))**2*p1_leg(0:n-1,l-1))
      end do

      return

      END SUBROUTINE comp_legendre

c     ===============================================================



c     ===============================================================

      SUBROUTINE hilfsgroessen_r

c     ===============================================================

      USE param_var
      USE gitter_var
      USE mod_hilfsgroessen_r

      IMPLICIT NONE

      real :: lp3_r,lm2_r
      integer :: i,l,lp3,lm2

      r_if_hoch_l(0,0)=1.0d0
      r_if_hoch_l(0,1:l_max+1)=0.0d0
      do l=0,l_max+2
         r_if_hoch_l(1:m,l)=r_if(1:m)**l
         r_if_hoch_l_r(1:m,l)=1.0d0/r_if_hoch_l(1:m,l)
         r_hoch_l(1:m,;)=r(1:m)**l
         r_hoch_lp1_r(1:m,;)=1.0d0/(r(1:m)**(l+1))
      end do

      do l=0,l_max
         lp3=l+3
         lp3_r=1.0d0/dble(lp3)

         do i=1,m
            gew_r_in(i,l)=(r_if(i)**lp3-r_if(i-1)**lp3)*lp3_r
         end do

         lm2=l-2

         if (l.ne.2) then
            lm2_r=1.0d0/dble(lm2)
            do i=2,m
               gew_r_ex(i,l)=lm2_r*
     +              ((1.0d0/r_if(i-1)**lm2)-(1.0d0/r_if(i)**lm2))
            end do
         else
            do i=2,m
               gew_r_ex(i,l)=log(r_if(i)/r_if(i-1))
            end do
         end if
         
          if (l.eq.0) then
             gew_r_ex(1,l)=r(1)*dr(1)
          else
             gew_r_ex(1,l)=0.0d0
          end if

      end do

      return

      END SUBROUTINE hilfsgroessen_r

c     ===============================================================




c     ===============================================================

      SUBROUTINE loese_poisson (rho,phi)

c     ===============================================================
c
c     Solves Poisson's equation div grad Phi = 4 pi rho.
c
c     ---------------------------------------------------------------

      USE size_cfc, ONLY: m,n,o
      USE legendre
      USE mod_hilfsgroessen_r

      IMPLICIT NONE

      real, intent (in) :: rho (1:m,1:n,1:o)
      real, intent (out) :: phi (1:m,1:n,1:o)

      real, dimension(1:m  ,0:l_max) :: mom_rho
      real, dimension(0:m  ,0:l_max) :: m_in,m_ex
      real, dimension(1:m+1,0:l_max) :: phi_l

      integer :: i,j,l
      integer :: lp1

      integer, parameter :: k=1 !Loesung nur in Axialsymmetrie

!$OMP PARALLEL DO
!$OMP&PRIVATE(i,j,l,lp1)
c!$OMP&SHARED(mom_rho,u,p_leg_int,phi_in,r_if_hoch,gew_r
!$OMP&DEFAULT(NONE)
      do l=0,l_max,l_schritt

         lp1=l+1
         lp1_r=1.0d0/dble(lp1)

         mom_rho(:,l)=0.0d0
         do j=1,n
            do i=1,m
               mom_rho(i,l)=mom_rho(i,l)+rho(i,j,k,0)*p_leg_int(j,l)
            end do
         end do

         m_in(0,l)=0.0d0
         do i=1,m
            m_in(i,l)=m_in(i-1,l)+mom_rho(i,l)*r_hoch_l(i,l)  *dvxtot(i)
         end do
         m_ex(m,l)=0.0
         do i=m-1,1,-1
            m_ex(i,l)=m_ex(i+1,l)+mom_rho(i,l)*r_hoch_lp1(i,l)*dvxtot(i)
         end do

         phi_l(m,l)=-m_in(m,l)/r(i)**(l+1) !oder doch anders???
         do i=m-1,1,-1
            phi_l(i,l)=phi_in(i+1,l)+
     +           (m_ex(i,l)*l  *r_if_hoch_l  (i,l+1)
     +           -m_in(i,l)*lp1/r_if_hoch_l_r(i,l+2))*dr_if(i)
         end do

      end do
!$OMP END PARALLEL DO

!$OMP PARALLEL DO
!$OMP&PRIVATE(i,j,l)
!$OMP&DEFAULT(NONE)
      do j=1,n
         phi(:,j,k)=0.0d0
         do l=0,l_max,l_schritt
            do i=1,m
               phi(i,j,k)=phi(i,j,k)+phi_l(i,l)*p_leg_c(j,l)
            end do
         end do
!$OMP END PARALLEL DO

      return
      
      END SUBROUTINE loese_poisson

c     ===============================================================



c     ===============================================================

      SUBROUTINE loese_vektor_poisson (s_vec,beta)

c     ===============================================================
c
c     Solves the vector Poisson equation div grad beta = 4 pi S.
c
c     ---------------------------------------------------------------

      USE param_var
      USE legendre
      USE gitter_var
      USE grav_var
      USE mod_hilfsgroessen_r

      IMPLICIT NONE

      real, intent (in) :: s_vec (1:m,1:n,1:o,3)
      real, dimension(1:m,1:n,1:o,3), intent (out) :: beta

      real, dimension(1:m,1:n,1:o) :: s_x,s_y,s_z

      real, dimension(1:m  ,0:l_max,3) :: mom_s
      real, dimension(0:m  ,0:l_max,3) :: m_in,m_ex
      real, dimension(1:m+1,0:l_max,3) :: phi_l

      real, dimension(1:m,1:n,1:o) :: beta_x,beta_y,beta_z

      real :: lp1_r,llp1_r

      integer :: i,j,l
      integer :: lp1

      integer, parameter :: k=1 !Loesung nur in Axialsymmetrie

!$OMP PARALLEL
!$OMP&PRIVATE(i,j,l,lp1,lp1_r,llp1_r)
!$OMP&DEFAULT(NONE)

!$OMP DO
      do j=1,n
         do i=1,m
            s_x(i,j,k)=
     +           sinus_theta(j)     *s_vec(i,j,k,1)+
     +           cosin_theta(j)*r(i)*s_vec(i,j,k,2)
            s_y(i,j,k)=
     +           sinus_theta(j)*r(i)*s_vec(i,j,k,3)
            s_z(i,j,k)=
     +           cosin_theta(j)     *s_vec(i,j,k,1)-
     +           sinus_theta(j)*r(i)*s_vec(i,j,k,2)
         end do
      end do
!$OMP END DO


!$OMP DO
      do l=1,l_max,l_schritt

         lp1=l+1
         lp1_r=1.0d0/dble(lp1)
         llp1_r=1.0d0/dble(lp1*l)

         mom_rho(:,l,1)=0.0d0
         mom_rho(:,l,2)=0.0d0
         do j=1,n
            do i=1,m
               mom_rho(i,l,1)=mom_rho(i,l,1)+s_x(i,j,k)*p1_leg_int(j,l)
               mom_rho(i,l,2)=mom_rho(i,l,2)+s_y(i,j,k)*p1_leg_int(j,l)
            end do
         end do

         m_in(0,l,1)=0.0d0
         m_in(0,l,2)=0.0d0
         do i=1,m
            m_in(i,l,1)=m_in(i-1,l,1)+mom_rho(i,l,1)*r_hoch_l(i,l)*
     +           dvxtot(i)
            m_in(i,l,2)=m_in(i-1,l,2)+mom_rho(i,l,2)*r_hoch_l(i,l)*
     +           dvxtot(i)
         end do
         m_ex(m,l,1)=0.0
         m_ex(m,l,2)=0.0
         do i=m-1,1,-1
            m_ex(i,l,1)=m_ex(i+1,l,1)+mom_rho(i,l,1)*r_hoch_lp1(i,l)*
     +           dvxtot(i)
            m_ex(i,l,2)=m_ex(i+1,l,2)+mom_rho(i,l,2)*r_hoch_lp1(i,l)*
     +           dvxtot(i)
         end do

         phi_l(m,l,1)=-m_in(m,l,1)/r(i)**(l+1) !oder doch anders???
         phi_l(m,l,2)=-m_in(m,l,2)/r(i)**(l+1) !oder doch anders???
         do i=m-1,1,-1
            phi_l(i,l,1)=phi_in(i+1,l,1)+
     +           (m_ex(i,l,1)*l  *r_if_hoch_l  (i,l+1)
     +           -m_in(i,l,1)*lp1/r_if_hoch_l_r(i,l+2))*dr_if(i)
            phi_l(i,l,2)=phi_in(i+1,l,2)+
     +           (m_ex(i,l,2)*l  *r_if_hoch_l  (i,l+1)
     +           -m_in(i,l,2)*lp1/r_if_hoch_l_r(i,l+2))*dr_if(i)
         end do

      end do
!$OMP END DO


!$OMP DO
      do l=0,l_max,l_schritt

         lp1=l+1
         lp1_r=1.0d0/dble(lp1)
         llp1_r=1.0d0/dble(lp1*l)

         mom_rho(:,l,3)=0.0d0
         do j=1,n
            do i=1,m
               mom_rho(i,l,3)=mom_rho(i,l,3)+s_z(i,j,k)*p_leg_int(j,l)
            end do
         end do


         m_in(0,l,3)=0.0d0
         do i=1,m
            m_in(i,l,3)=m_in(i-1,l,3)+mom_rho(i,l,3)*r_hoch_l(i,l)*
     +           dvxtot(i)
         end do
         m_ex(m,l,3)=0.0
         do i=m-1,1,-1
            m_ex(i,l,3)=m_ex(i+1,l,3)+mom_rho(i,l,3)*r_hoch_lp1(i,l)*
     +           dvxtot(i)
         end do

         phi_l(m,l,3)=-m_in(m,l,3)/r(i)**(l+1) !oder doch anders???
         do i=m-1,1,-1
            phi_l(i,l,3)=phi_in(i+1,l,3)+
     +           (m_ex(i,l,3)*l  *r_if_hoch_l  (i,l+1)
     +           -m_in(i,l,3)*lp1/r_if_hoch_l_r(i,l+2))*dr_if(i)
         end do

      end do
!$OMP END DO
!$OMP END PARALLEL


!$OMP PARALLEL DO
!$OMP&PRIVATE(i,j,l)
!$OMP&DEFAULT(NONE)
      do j=1,n
         beta_x(:,j,k)=0.0d0
         beta_y(:,j,k)=0.0d0
         beta_z(:,j,k)=0.0d0
         do l=1,l_max,l_schritt
            do i=1,m
               beta_x(i,j,k)=beta_x(i,j,k)+phi_l(i,l,1)*p1_leg_c(j,l)
               beta_y(i,j,k)=beta_y(i,j,k)+phi_l(i,l,2)*p1_leg_c(j,l)
            end do
         end do

         do l=1,l_max,l_schritt
            do i=0,m
               beta_z(i,j,k)=beta_z(i,j,k)+phi_l(i,l,3)* p_leg_c(j,l)
            end do
         end do

#ifndef AEQUATOR
         beta_x(:,j,k)=beta_x(:,j,k)*zwei_pi_g
         beta_y(:,j,k)=beta_y(:,j,k)*zwei_pi_g
         beta_z(:,j,k)=beta_z(:,j,k)*zwei_pi_g
#else
         beta_x(:,j,k)=beta_x(:,j,k)*vier_pi_g
         beta_y(:,j,k)=beta_y(:,j,k)*vier_pi_g
         beta_z(:,j,k)=beta_z(:,j,k)*vier_pi_g
#endif

      end do
!$OMP END PARALLEL DO



!$OMP DO
      do j=0,n
         do i=0,m
            beta(i,j,k,1)=
     +            sinus_theta(j)*beta_x(i,j,k)+
     +            cosin_theta(j)*beta_z(i,j,k)
            beta(i,j,k,2)=
     +           (cosin_theta(j)*beta_x(i,j,k)-
     +            sinus_theta(j)*beta_y(i,j,k))*r_inv(i)
            beta(i,j,k,3)=
     +           beta_y(i,j,k)*r_inv(i)*sinus_theta_inv(j)
         end do
      end do
!$OMP END DO
!$OMP END PARALLEL



      return
      
      END SUBROUTINE loese_vektor_poisson

c     ===============================================================


      END MODULE poisson

c     ===============================================================

