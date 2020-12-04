c     ===============================================================

      MODULE poisson

c     ===============================================================

      IMPLICIT NONE

      CONTAINS

#if NDIM>1
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


#if GEOM==1
c     ===============================================================

      SUBROUTINE loese_poisson

c     ===============================================================

      USE param_var
      USE legendre
      USE hydro_var
      USE gitter_var
      USE grav_var
      USE mod_hilfsgroessen_r

      IMPLICIT NONE

      real, dimension(0:NX,0:l_max) :: mom_rho,phi_in_l,phi_ex_l

      real :: lp1_r,l2p1_r

      integer :: i,j,l
      integer :: lp1
c#if NDIM<3
      integer, parameter :: k=1 !Loesung nur in Axialsymmetrie
c#endif

!$OMP PARALLEL DO
!$OMP&PRIVATE(i,j,l,lp1,lp1_r,l2p1_r)
c!$OMP&SHARED(mom_rho,u,p_leg_int,phi_in_l,r_if_hoch,gew_r
      do l=0,l_max,l_schritt

         lp1=l+1
         lp1_r=1.0d0/dble(lp1)
         l2p1_r=1.0d0/dble(l*2+1)

         mom_rho(:,l)=0.0d0
         do j=1,NY
            do i=1,NX
               mom_rho(i,l)=mom_rho(i,l)+u(i,j,k,0)*p_leg_int(j,l)
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
c        if (l.eq.0) then
c           phi_ex_l(0,l)=phi_ex_l(1,l)+mom_rho(1,l)*r(1)*dr(1)
c        else
c           phi_ex_l(0,l)=0.0d0
c        end if

      end do
!$OMP END PARALLEL DO

!$OMP PARALLEL DO
!$OMP&PRIVATE(i,j,l)
      do j=0,NY
         phi_grav_if(:,j,k)=0.0d0
         do i=0,NX
            phi_grav_if(i,j,k) = - m_core / r_if (i)
         end doxs
         do l=0,l_max,l_schritt
            do i=0,NX
               phi_grav_if(i,j,k)=phi_grav_if(i,j,k)-
     &              (phi_in_l(i,l)+phi_ex_l(i,l))*p_leg(j,l)
            end do
         end do

#ifndef AEQUATOR
         phi_grav_if(:,j,k)=phi_grav_if(:,j,k)*zwei_pi_g
#else
         phi_grav_if(:,j,k)=phi_grav_if(:,j,k)*vier_pi_g
#endif

      end do
!$OMP END PARALLEL DO


!$OMP PARALLEL DO
!$OMP&PRIVATE(i,j)
      do j=1,NY
         do i=1,NX
            a_grv(i,j,k,1)=
     +           ((phi_grav_if(i-1,  j,k)-phi_grav_if(i,  j,k))+
     +             phi_grav_if(i  ,j-1,k)-phi_grav_if(i,j-1,k))*
     +           dr1(i)
            a_grv(i,j,k,2)=
     +           (phi_grav_if(i  ,j-1,k)-phi_grav_if(i  ,j,k)+
     +            phi_grav_if(i-1,j-1,k)-phi_grav_if(i-1,j,k))*
     +           dtheta1(j)*r1(i)
            phi_grav(i,j,k)=0.25d0*
     +           (phi_grav_if(i,  j,k)+phi_grav_if(i-1,j  ,k)+
     +            phi_grav_if(i,j-1,k)+phi_grav_if(i-1,j-1,k))
         end do
      end do
!$OMP END PARALLEL DO

      return
      
      END SUBROUTINE loese_poisson

c     ===============================================================
#endif



c     =============================================================== 

      SUBROUTINE gravitation_2d(modus)

c     =============================================================== 

      USE grav_var
      USE hydro_var
      USE hydro_sav_var
      
      IMPLICIT NONE

      integer, intent(in) :: modus
      integer :: i,j

!$OMP PARALLEL DO
!$OMP& PRIVATE(j)      
      do j=1,NY
         phi_grav_alt(:,j,1)=phi_grav(:,j,1)
      end do
!$OMP END PARALLEL DO

#if GEOM==1
      call loese_poisson
#endif

#ifdef ENE_CONS_2
      if (modus.eq.0) then
c     Korrekturterm zur Energiebilanz -> 2.Ordnung fuer zeitliche Diskretisierung
!$OMP PARALLEL DO
!$OMP& PRIVATE(i,j)
         do j=1,NY
            do i=1,NX
               u(i,j,1,4)=u(i,j,1,4)+0.5d0*
     +              (phi_grav(i,j,1)*u_sav(i,j,1,0)+
     +              phi_grav_alt(i,j,1)*u(i,j,1,0)-
     +              phi_grav_alt(i,j,1)*u_sav(i,j,1,0)-
     +              phi_grav(i,j,1)*u(i,j,1,0))
c               u(i,j,1,4)=u(i,j,1,4)+
c     +              u_sav(i,j,1,0)*(phi_grav(i,j,1)-phi_alt(i,j,1))-
c     +              u    (i,j,,1,)*(phi_grav(i,j,1)-phi_alt(i,j,1))
            end do
         end do
!$OMP END PARALLEL DO
      end if
#endif

      return

      END SUBROUTINE gravitation_2d

c     ===============================================================



#endif


#if NDIM==1
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

      real, dimension(0:NX) :: m

      real :: px

      integer :: i,j,k

      phi_grav_alt(:,:,:)=phi_grav(:,:,:)

#if GEOM==1
c      m(0)=0.0d0
      m(0)=m_core
      do i=1,NX
         m(i)=m(i-1)+sum(u(i,1:NY,1:NZ,0)*dv(i,1:NY,1:NZ))*sym_fac
      end do

      phi_grav(NX,1,1)=-m(NX)/r_if(NX)
      do i=NX-1,1,-1
         phi_grav(i,1,1)=phi_grav(i+1,1,1)-m(i)/r_if(i)**2*dr_if(i)
      end do
      phi_grav(0,1,1)=phi_grav(1,1,1)-(phi_grav(2,1,1)-phi_grav(1,1,1))
      phi_grav(-1,1,1)=phi_grav(0,1,1)-(phi_grav(2,1,1)-phi_grav(1,1,1))
      phi_grav(NX+1,1,1)=phi_grav(NX,1,1)
      phi_grav(NX+2,1,1)=phi_grav(NX,1,1)
      phi_grav(:,1,1)=g_grav*phi_grav(:,1,1)

c      a_grv(1,1,1,1)=(-0.5d0*m(1)/r_if(1)**2*g_grav)*r(1)/r_if(1)
!$OMP PARALLEL DO
      do i=1,NX
c         a_grv(i,1,1,1)=-(m(i)/r_if(i)**2*
c     +        (r_if(i)-r(i))/(r_if(i)-r_if(i-1))+
c     +        m(i-1)/r_if(i-1)**2*
c     +        (r(i)-r_if(i-1))/(r_if(i)-r_if(i-1)))*
c     +        g_grav
         a_grv(i,1,1,1)=(phi_grav(i-1,1,1)-phi_grav(i+1,1,1))/
     +        (r(i+1)-r(i-1))
      end do
#endif


#ifdef ENE_CONS_2
      if (modus.eq.0) then
c     Korrekturterm zur Energiebilanz -> 2.Ordnung fuer zeitliche Diskretisierung
         do i=1,NX
            u(i,1,1,4)=u(i,1,1,4)+0.5d0*
     +           (phi_grav(i,1,1)*u_sav(i,1,1,0)+
     +           phi_grav_alt(i,1,1)*u(i,1,1,0)-
     +           phi_grav_alt(i,1,1)*u_sav(i,1,1,0)-
     +           phi_grav(i,1,1)*u(i,1,1,0))
         end do
      end if
#endif

      return

      END SUBROUTINE gravitation_1d

c     ===============================================================
#endif


      END MODULE poisson

c     ===============================================================

