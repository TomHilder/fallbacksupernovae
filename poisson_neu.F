c     ===============================================================

      MODULE poisson

c     ===============================================================

      IMPLICIT NONE

      CONTAINS

c     ===============================================================

      SUBROUTINE comp_legendre

c     ===============================================================
c     Berechnet die Legendre-Polynome P_l(cos(theta_if)) an den 
c     Zellraendern und die Integrale von P_l ueber eine Zelle in
c     theta-Richtung.
c     ---------------------------------------------------------------

      USE legendre
      USE gitter_var

      IMPLICIT NONE

      integer :: j,l

      p_leg(:,0)=1.0d0
      p_leg(:,0)=1.0d0
      p_leg(0:NY,1)=cos(theta_if(0:NY))
#ifdef AEQUATOR
      p_leg(NY,1)=0.0d0
#endif

      p1_leg(0:NY,1)=-sin(theta_if(0:NY))
      p1_leg(0:NY,2)=-3*sin(theta_if(0:NY))*cos(theta_if(0:NY))

c     Rekursionsformel fuer die Legendre-Polynome:
c     (Jackson, Classical Electrodynamics, aeq. 3.29):
      do l=1,l_max
         p_leg(:,l+1)=(p_leg(:,l)*p_leg(:,1)*(2*l+1)-
     +        l*p_leg(:,l-1))/dble(l+1)
      end do

      do l=2,l_max
         p1_leg(:,l+1)=(p1_leg(:,l)*p_leg(:,1)*(2*l+1)-
     +        (l+1)*p1_leg(:,l-1))/dble(l)
      end do


c     Exakte Werte an den Raendern setzen (nicht unbedingt noetig):
      do l=0,l_max+1
         if (mod(l,2).eq.0) then
            p_leg(0,l)=1.0d0
         else
            p_leg(0,l)=-1.0d0
         end if
#ifndef AEQUATOR
         p_leg(NY,l)=1.0d0
#endif
      end do

c     Integrale ueber theta-Zelle mit Hilfe der Formel fuer die
c     Ableitungen der Legendre-Polynome (ibid., aeq.3.28)
      p_leg_int(1:NY,0)=p_leg(1:NY,1)-p_leg(0:NY-1,1)
      do l=1,l_max
         p_leg_int(1:NY,l)=
     +        ((p_leg(1:NY  ,l+1)-p_leg(1:NY  ,l-1))-
     +         (p_leg(0:NY-1,l+1)-p_leg(0:NY-1,l-1)))/
     +        dble(2*l+1)
#ifdef AEQUATOR
         p_leg_int(1:NY,l)=2.0d0*p_leg_int(1:NY,l)
#endif
         print *,l,sum(p_leg_int(1:NY,l)*cos(theta(1:NY))**1)
      end do

      p1_leg_int(1:NY,1)=
     +     -(sin(2*theta_if(1:NY))-sin(2*theta_if(0:NY-1)))/4.0+
     +     (theta_if(1:NY)-theta_if(0:NY-1))/2.0
      p1_leg_int(1:NY,2)=
     +     sin(theta_if(1:NY))**3-sin(theta_if(0:NY-1))**3

      do l=3,l_max
         p1_leg_int(1:NY,l)=
     +        (l-2)*l/dble((l+1)*(l-1))*
     +        p1_leg_int(1:NY,l-2)-
     +        (2*l-1)/dble((l+1)*(l-1))*
     +        (sin(theta_if(1:NY  ))**2*p1_leg(1:NY  ,l-1)-
     +         sin(theta_if(0:NY-1))**2*p1_leg(0:NY-1,l-1))
      end do

      do l=1,l_max
         print *,l,sum(p1_leg_int(1:NY,l)*sin(theta(1:NY))**1)
c         do j=1,NY
c            print *,j,0.5*(p1_leg(j,l)+p1_leg(j-1,l))
c            print *,j,p1_leg_int(j,l)/
c     +           (cos(theta_if(j))-cos(theta_if(j-1))),
c     +           sin(theta(j))
c         end do
         print *
c         if (l.eq.5) exit
      end do

      stop
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
      do l=0,l_max+1
         r_if_hoch_l(1:NX,l)=r_if(1:NX)**l
         r_if_hoch_l_r(1:NX,l)=1.0d0/r_if_hoch_l(1:NX,l)
      end do

      do l=0,l_max
         lp3=l+3
         lp3_r=1.0d0/dble(lp3)

         do i=1,NX
            gew_r_in(i,l)=(r_if(i)**lp3-r_if(i-1)**lp3)*lp3_r
         end do

         lm2=l-2

         if (l.ne.2) then
            lm2_r=1.0d0/dble(lm2)
            do i=2,NX
               gew_r_ex(i,l)=lm2_r*
     +              ((1.0d0/r_if(i-1)**lm2)-(1.0d0/r_if(i)**lm2))
            end do
         else
            do i=2,NX
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

      SUBROUTINE loese_poisson (rho,phi_if,phi)

c     ===============================================================
c
c     Solves Poisson's equation div grad Phi = 4 pi rho.
c
c     ---------------------------------------------------------------

      USE param_var
      USE legendre
      USE gitter_var
      USE grav_var
      USE mod_hilfsgroessen_r

      IMPLICIT NONE

      real, intent (in) :: rho (1:NX,1:NY,1:NZ)
      real, intent (out) :: phi_if (0:NX,0:NY,0:NZ)
      real, intent (out) :: phi (1:NX,1:NY,1:NZ)

      real, dimension(0:NX,0:l_max) :: mom_rho,phi_in_l,phi_ex_l

      real :: lp1_r,l2p1_r

      integer :: i,j,l
      integer :: lp1

      integer, parameter :: k=1 !Loesung nur in Axialsymmetrie

!$OMP PARALLEL DO
!$OMP&PRIVATE(i,j,l,lp1,lp1_r,l2p1_r)
c!$OMP&SHARED(mom_rho,u,p_leg_int,phi_in_l,r_if_hoch,gew_r
!$OMP&DEFAULT(NONE)
      do l=0,l_max,l_schritt

         lp1=l+1
         lp1_r=1.0d0/dble(lp1)
         l2p1_r=1.0d0/dble(l*2+1)

         mom_rho(:,l)=0.0d0
         do j=1,NY
            do i=1,NX
               mom_rho(i,l)=mom_rho(i,l)+rho(i,j,k,0)*p_leg_int(j,l)
            end do
         end do

         phi_in_l(0,l)=0.0d0
         do i=1,NX
            phi_in_l(i,l)=(phi_in_l(i-1,l)*r_if_hoch_l(i-1,lp1)+
     +           mom_rho(i,l)*gew_r_in(i,l))*r_if_hoch_l_r(i,lp1)
         end do
         phi_ex_l(NX,l)=0.0d0
         do i=NX,1,-1
            phi_ex_l(i-1,l)=(phi_ex_l(i,l)*r_if_hoch_l_r(i,l)+
     +           mom_rho(i,l)*gew_r_ex(i,l))*r_if_hoch_l(i-1,l)
         end do

      end do
!$OMP END PARALLEL DO

!$OMP PARALLEL DO
!$OMP&PRIVATE(i,j,l)
!$OMP&DEFAULT(NONE)
      do j=0,NY
         phi_if(:,j,k)=0.0d0
         do l=0,l_max,l_schritt
            do i=0,NX
               phi_if(i,j,k)=phi_if(i,j,k)-
     +              (phi_in_l(i,l)+phi_ex_l(i,l))*p_leg(j,l)
            end do
         end do

#ifndef AEQUATOR
         phi_if(:,j,k)=phi_if(:,j,k)*zwei_pi_g
#else
         phi_if(:,j,k)=phi_if(:,j,k)*vier_pi_g
#endif

      end do
!$OMP END PARALLEL DO


!$OMP PARALLEL DO
!$OMP&PRIVATE(i,j)
      do j=1,NY
         do i=1,NX
            phi(i,j,k)=0.25d0*
     +           (phi_if(i,  j,k)+phi_if(i-1,j  ,k)+
     +            phi_if(i,j-1,k)+phi_if(i-1,j-1,k))
         end do
      end do
!$OMP END PARALLEL DO

      return
      
      END SUBROUTINE loese_poisson

c     ===============================================================



c     ===============================================================

      SUBROUTINE loese_vektor_poisson (s_vec,beta,beta_if)

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

      real, intent (in) :: s_vec (1:NX,1:NY,1:NZ,3)
      real, dimension(0:NX,0:NY,0:NZ,3), intent (out) :: beta_if
      real, dimension(1:NX,1:NY,1:NZ,3), intent (out) :: beta

      real, dimension(1:NX,1:NY,1:NZ) :: s_x,s_y,s_z
      real, dimension(0:NX,0:NY,0:NZ) :: beta_x,beta_y,beta_z

      real, dimension(0:NX,0:l_max,3) :: mom_s,phi_in_l,phi_ex_l

      real :: lp1_r,llp1_r

      integer :: i,j,l
      integer :: lp1

      integer, parameter :: k=1 !Loesung nur in Axialsymmetrie

!$OMP PARALLEL
!$OMP&PRIVATE(i,j,l,lp1,lp1_r,llp1_r)
!$OMP&DEFAULT(NONE)

!$OMP DO
      do j=1,NY
         do i=1,NX
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
         do j=1,NY
            do i=1,NX
               mom_rho(i,l,1)=mom_rho(i,l,1)+s_x(i,j,k)*p1_leg_int(j,l)
               mom_rho(i,l,2)=mom_rho(i,l,2)+s_y(i,j,k)*p1_leg_int(j,l)
            end do
         end do

         phi_in_l(0,l,1)=0.0d0
         phi_in_l(0,l,2)=0.0d0
         do i=1,NX
            phi_in_l(i,l,1)=(phi_in_l(i-1,l,1)*r_if_hoch_l(i-1,lp1)+
     +           mom_rho(i,l,1)*gew_r_in(i,l))*r_if_hoch_l_r(i,lp1)
            phi_in_l(i,l,2)=(phi_in_l(i-1,l,2)*r_if_hoch_l(i-1,lp1)+
     +           mom_rho(i,l,2)*gew_r_in(i,l))*r_if_hoch_l_r(i,lp1)
         end do

         phi_ex_l(NX,l,1)=0.0d0
         phi_ex_l(NX,l,2)=0.0d0
         do i=NX,1,-1
            phi_ex_l(i-1,l,1)=(phi_ex_l(i,l,1)*r_if_hoch_l_r(i,l)+
     +           mom_rho(i,l,1)*gew_r_ex(i,l))*r_if_hoch_l(i-1,l)
            phi_ex_l(i-1,l,2)=(phi_ex_l(i,l,2)*r_if_hoch_l_r(i,l)+
     +           mom_rho(i,l,2)*gew_r_ex(i,l))*r_if_hoch_l(i-1,l)
         end do

         phi_in_l(:,l,1)=phi_in_l(:,l,1)*llp1_r
         phi_ex_l(:,l,2)=phi_ex_l(:,l,2)*llp1_r        

      end do
!$OMP END DO


!$OMP DO
      do l=0,l_max,l_schritt

         lp1=l+1
         lp1_r=1.0d0/dble(lp1)
         llp1_r=1.0d0/dble(lp1*l)

         mom_rho(:,l,3)=0.0d0
         do j=1,NY
            do i=1,NX
               mom_rho(i,l,3)=mom_rho(i,l,3)+s_z(i,j,k)*p_leg_int(j,l)
            end do
         end do

         phi_in_l(0,l,3)=0.0d0
         do i=1,NX
            phi_in_l(i,l,3)=(phi_in_l(i-1,l,3)*r_if_hoch_l(i-1,lp1)+
     +           mom_rho(i,l,3)*gew_r_in(i,l))*r_if_hoch_l_r(i,lp1)
         end do

         phi_ex_l(NX,l,3)=0.0d0
         do i=NX,1,-1
            phi_ex_l(i-1,l,3)=(phi_ex_l(i,l,3)*r_if_hoch_l_r(i,l)+
     +           mom_rho(i,l,3)*gew_r_ex(i,l))*r_if_hoch_l(i-1,l)
         end do

      end do
!$OMP END DO
!$OMP END PARALLEL


!$OMP PARALLEL DO
!$OMP&PRIVATE(i,j,l)
!$OMP&DEFAULT(NONE)
      do j=0,NY
         beta_x(:,j,k)=0.0d0
         beta_y(:,j,k)=0.0d0
         beta_z(:,j,k)=0.0d0
         do l=1,l_max,l_schritt
            do i=0,NX
               beta_x(i,j,k)=beta_x(i,j,k)-
     +              (phi_in_l(i,l,1)+phi_ex_l(i,l,1))*p1_leg(j,l)
               beta_y(i,j,k)=beta_y(i,j,k)-
     +              (phi_in_l(i,l,2)+phi_ex_l(i,l,2))*p1_leg(j,l)
            end do
         end do

         do l=1,l_max,l_schritt
            do i=0,NX
               beta_z(i,j,k)=beta_z(i,j,k)-
     +              (phi_in_l(i,l,3)+phi_ex_l(i,l,3))*p_leg(j,l)
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
      do j=0,NY
         do i=0,NX
            beta_if(i,j,k,1)=
     +            sinus_theta(j)*beta_x(i,j,k)+
     +            cosin_theta(j)*beta_z(i,j,k)
            beta_if(i,j,k,2)=
     +           (cosin_theta(j)*beta_x(i,j,k)-
     +            sinus_theta(j)*beta_y(i,j,k))*r_if_inv(i)
            beta_if(i,j,k,3)=
     +           beta_y(i,j,k)*r_if_inv(i)*sinus_theta_if_inv(j)
         end do
      end do
!$OMP END DO


!$OMP DO
      do j=1,NY
         do i=1,NX
            beta(i,j,k,1)=0.25d0*
     +           (beta_if(i-1,j-1,k,1)+beta_if(i-1,j,k,1)+
     +            beta_if(i  ,j-1,k,1)+beta_if(i  ,j,k,1))
            beta(i,j,k,2)=0.25d0*
     +           (beta_if(i-1,j-1,k,2)+beta_if(i-1,j,k,2)+
     +            beta_if(i  ,j-1,k,2)+beta_if(i  ,j,k,2))
            beta(i,j,k,3)=0.25d0*
     +           (beta_if(i-1,j-1,k,3)+beta_if(i-1,j,k,3)+
     +            beta_if(i  ,j-1,k,3)+beta_if(i  ,j,k,3))
         end do
      end do
!$OMP END DO


!$OMP END PARALLEL



      return
      
      END SUBROUTINE loese_vektor_poisson

c     ===============================================================


      END MODULE poisson

c     ===============================================================

