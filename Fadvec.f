










c     ===============================================================

      MODULE advection

c     ===============================================================

      CONTAINS

c     ===============================================================

      SUBROUTINE advec_x(dt_modus)

c     ===============================================================

      USE hydro_var
      USE gitter_var
      USE grav_var
      USE poisson

      IMPLICIT NONE

      integer, intent(in) :: dt_modus

      real, dimension(-1:800+2,1:1,1:1) :: c
      real, dimension(-1:800+2,1:1,1:1,0:4) :: w,
     +     slope_r,slope_l,u_r,u_l
      real, dimension(-1:800+1,1:1,1:1,0:4) :: fl
      real, dimension(-1:800+2,1:1,1:1,0:4) :: u_alt
      real :: vx,vy,vz,rho

      integer :: i,j,k,i_rk,l
c     ---------------------------------------------------------------
c     Erster Runge-Kutta-Schritt
c     ---------------------------------------------------------------
      
c     Randbedingungen setzen:
      call randbed(1)
      call fluesse(1)

      do k=1,1
         do j=1,1
            do i=-1,800+2
               u_alt(i,j,k,0)=u(i,j,k,0)
               u_alt(i,j,k,1)=u(i,j,k,1)
               u_alt(i,j,k,2)=u(i,j,k,2)
               u_alt(i,j,k,3)=u(i,j,k,3)
               u_alt(i,j,k,4)=u(i,j,k,4)
c     Flux-Freezing-Geschwindigkeit:
               c(i,j,k)=c_s(i,j,k)+
     +              sqrt(u(i,j,k,1)**2+u(i,j,k,2)**2+u(i,j,k,3)**2)/
     +              u(i,j,k,0)
               w(i,j,k,:)=fl_u(i,j,k,:)/c(i,j,k)
            end do
c            c(:,j,k)=max(1.0d-2*maxval(c(:,j,k)),c(:,j,k))
c            c(:,j,k)=maxval(c(:,j,k))
         end do
      end do

      if (dt_modus.eq.0) then
         call zeitschritt(c(1:800,1:1,1:1))
      else
         dt=0.5d0*dt
      end if

      do k=1,1
         do j=1,1

            do i=-1,800+1
               fl(i,j,k,:)=
     +              0.25d0*(c(i+1,j,k)+c(i,j,k))*
     +              (u(i,j,k,:)+w(i,j,k,:))
            end do

            do i=-1,800+1
               fl(i,j,k,:)=fl(i,j,k,:)-
     +              0.25d0*(c(i+1,j,k)+c(i,j,k))*
     +              (u(i+1,j,k,:)-w(i+1,j,k,:))
            end do

c     Updates der erhaltenene Groessen:
            do i=0,800+1
               rho=u(i,j,k,0)
               vx=u(i,j,k,1)/rho
               vy=u(i,j,k,2)/rho
               vz=u(i,j,k,3)/rho

               u(i,j,k,0)=u_alt(i,j,k,0)+
     +              dt*(fl(i-1,j,k,0)*da(i-1,j,k,1)-
     +              fl(i,j,k,0)*da(i,j,k,1))*dv1(i,j,k)
               u(i,j,k,4)=u_alt(i,j,k,4)+
     +              dt*(fl(i-1,j,k,4)*da(i-1,j,k,1)-
     +              fl(i,j,k,4)*da(i,j,k,1))*dv1(i,j,k)-
c     Subtrahiere potentielle Energie:
     +              u(i,j,k,0)*phi_grav(i,j,k)
c               print *,i,u(i,j,k,:)
               u(i,j,k,1)=u_alt(i,j,k,1)+
     +              dt*((fl(i-1,j,k,1)-fl(i,j,k,1))*dr1(i)+
     +              rho*(
c     Scheinkraefte:
     +              r1(i)*(vy**2+vz**2-2.0d0*vx**2)+
c     Gravitationskraft:
     +              a_grv(i,j,k,1)))
               u(i,j,k,2)=u_alt(i,j,k,2)+
     +              dt*(fl(i-1,j,k,2)*da(i-1,j,k,1)-
     +              fl(i,j,k,2)*da(i,j,k,1))*dv1(i,j,k)
               u(i,j,k,3)=u_alt(i,j,k,3)+0.0d0*
     +              dt*(fl(i-1,j,k,3)*da(i-1,j,k,1)-
     +              fl(i,j,k,3)*da(i,j,k,1))*dv1(i,j,k)*
     +              r_sin_theta1(i,j,k)
            end do
         end do
      end do


c     ---------------------------------------------------------------
c     Zweiter Runge-Kutta-Schritt
c     ---------------------------------------------------------------

      dt=2.0d0*dt

c     Randbedingungen setzen:
      call randbed(1)
      call fluesse(1)

c     Flux-Freezing-Geschwindigkeit:
      do k=1,1
         do j=1,1
            do i=-1,800+2
               c(i,j,k)=c_s(i,j,k)+
     +              sqrt(u(i,j,k,1)**2+u(i,j,k,2)**2+u(i,j,k,3)**2)/
     +              u(i,j,k,0)
               w(i,j,k,:)=fl_u(i,j,k,:)/c(i,j,k)
               u_r(i,j,k,:)=0.5d0*(u(i,j,k,:)+w(i,j,k,:))
               u_l(i,j,k,:)=0.5d0*(u(i,j,k,:)-w(i,j,k,:))
            end do

            do l=0,4
               do i=0,800+1
                  slope_r(i,j,k,l)=van_leer(
     +                 (u_r(i,j,k,l)-u_r(i-1,j,k,l))*dr_if1(i-1),
     +                 (u_r(i+1,j,k,l)-u_r(i,j,k,l))*dr_if1(i))
                  slope_l(i,j,k,l)=van_leer(
     +                 (u_l(i,j,k,l)-u_l(i-1,j,k,l))*dr_if1(i-1),
     +                 (u_l(i+1,j,k,l)-u_l(i,j,k,l))*dr_if1(i))
               end do
            end do

            do i=0,800
               fl(i,j,k,:)=
     +              0.5d0*(c(i+1,j,k)+c(i,j,k))*
     +              (u_r(i,j,k,:)+0.5d0*slope_r(i,j,k,:)*dr(i))
            end do

            do i=0,800
               fl(i,j,k,:)=fl(i,j,k,:)-
     +              0.5d0*(c(i+1,j,k)+c(i,j,k))*
     +              (u_l(i+1,j,k,:)-0.5d0*slope_l(i+1,j,k,:)*dr(i+1))
            end do

c     Updates der erhaltenene Groessen:
            do i=1,800
               rho=u(i,j,k,0)
               vx=u(i,j,k,1)/rho
               vy=u(i,j,k,2)/rho
               vz=u(i,j,k,3)/rho

               u(i,j,k,0)=u_alt(i,j,k,0)+
     +              dt*(fl(i-1,j,k,0)*da(i-1,j,k,1)-
     +              fl(i,j,k,0)*da(i,j,k,1))*dv1(i,j,k)
               u(i,j,k,4)=u_alt(i,j,k,4)+
     +              dt*(fl(i-1,j,k,4)*da(i-1,j,k,1)-
     +              fl(i,j,k,4)*da(i,j,k,1))*dv1(i,j,k)-
c     Subtrahiere potentielle Energie:
     +              u(i,j,k,0)*phi_grav(i,j,k)
               u(i,j,k,1)=u_alt(i,j,k,1)+
     +              dt*((fl(i-1,j,k,1)-fl(i,j,k,1))*dr1(i)+
     +              rho*(
c     Scheinkraefte:
     +              r1(i)*(vy**2+vz**2-2.0d0*vx**2)+
c     Gravitationskraft:
     +              a_grv(i,j,k,1)))
               u(i,j,k,2)=u_alt(i,j,k,2)+
     +              dt*(fl(i-1,j,k,2)*da(i-1,j,k,1)-
     +              fl(i,j,k,2)*da(i,j,k,1))*dv1(i,j,k)
               u(i,j,k,3)=u_alt(i,j,k,3)+0.0d0*
     +              dt*(fl(i-1,j,k,3)*da(i-1,j,k,1)-
     +              fl(i,j,k,3)*da(i,j,k,1))*dv1(i,j,k)*
     +              r_sin_theta1(i,j,k)
            end do
         end do
      end do


      do k=1,1
         do j=1,1
            verlust_mtot=verlust_mtot+dt*fl(800,j,k,0)*da(800,j,k,1)
            verlust_etot=verlust_etot+dt*fl(800,j,k,4)*da(800,j,k,1)
            mloss_inner=mloss_inner-dt*fl(0,j,k,0)*da(0,j,k,1)
            m_core=m_core-dt*fl(0,j,k,0)*da(0,j,k,1)
         end do
      end do

      if (dt_modus.eq.0) zeit=zeit+dt

      return

      END SUBROUTINE advec_x

c     ===============================================================






c     ===============================================================

      REAL FUNCTION van_leer(a,b)

c     ===============================================================

      IMPLICIT NONE
      
      real :: a,b
      
      if (a*b.gt.0.0d0) then
         van_leer=2.0d0*a*b/(a+b) !van Leer
c         van_leer=sign(min(abs(a),abs(b)),a) !Minmod
      else
         van_leer=0.0d0
      end if

      return

      END FUNCTION van_leer

c     ===============================================================


c     ===============================================================

      SUBROUTINE randbed(richtung)

c     ===============================================================    

      USE grav_var
      USE hydro_var
      USE gitter_var

      IMPLICIT NONE

      integer, intent(in) :: richtung
      
      integer :: i,j,k
      
      select case(richtung)

c     ---------------------------------------------------------------
c     r-Richtung:
      case(1)
         do k=1,1
            do j=1,1
c     Hydro-Groessen:
c     Innenrand:               
c               u( 0,j,k,:)=u(1,j,k,:)
c               u(-1,j,k,:)=u(2,j,k,:)
c               u( 0,j,k,1)=-u(1,j,k,1)
c               u(-1,j,k,1)=-u(2,j,k,1)
               u( 0,j,k,0)=u(1,j,k,0)
               u( 0,j,k,1)=u(1,j,k,1)
               u( 0,j,k,2)=u(2,j,k,2)*(r(1)/r(0))
               u( 0,j,k,3)=u(2,j,k,3)*(r(1)/r(0))
               u( 0,j,k,4)=u(1,j,k,4)
               u(-1,j,k,0)=u(1,j,k,0)
               u(-1,j,k,1)=u(1,j,k,1)
               u(-1,j,k,2)=u(2,j,k,2)*(r(1)/r(-1))
               u(-1,j,k,3)=u(2,j,k,3)*(r(1)/r(-1))
               u(-1,j,k,4)=u(1,j,k,4)

c     Reflektierende Randbedingungen:
               u(800+1,j,k,0)=u(800-1,j,k,0)
               u(800+2,j,k,0)=u(800-2,j,k,0)
c         u(800+1,j,k,0)=rho_min1!u(800-1,j,k,0)
c         u(800+2,j,k,0)=rho_min1!u(800-2,j,k,0)

               u(800+1,j,k,1:3)=-u(800-1,j,k,1:3)*(r(800-1)/r(800+1))**2
               u(800+2,j,k,1:3)=-u(800-2,j,k,1:3)*(r(800-2)/r(800+2))**2
c         u(800+1,j,k,1:3)=0.0d0
c         u(800+2,j,k,1:3)=0.0d0

               u(800+1,j,k,4)=u(800-1,j,k,4)
               u(800+2,j,k,4)=u(800-2,j,k,4)
c         u(800+1,j,k,4)=eps_rand!u(800-1,j,k,4)
c         u(800+2,j,k,4)=eps_rand!u(800-2,j,k,4)


c     Gravitationsfeld:
               phi_grav(0,j,k)=phi_grav(1,j,k)
               phi_grav(-1,j,k)=phi_grav(2,j,k)
               phi_grav(800+1,j,k)=phi_grav(800,j,k)
               phi_grav(800+2,j,k)=phi_grav(800,j,k)
               
            end do
         end do




      end select

      return

      END SUBROUTINE randbed

c     ===============================================================




c     ===============================================================

      SUBROUTINE fluesse(richtung)

c     =============================================================== 

      USE hydro_var
      USE grav_var
      USE gitter_var
      USE param_var
      USE mod_eos

      IMPLICIT NONE

      integer, intent(in) :: richtung

      integer :: i_min,i_max,j_min,j_max,k_min,k_max

      integer :: i,j,k
      real :: v,rho


      select case(richtung)
      case(1)
         i_min=-1
         i_max=800+2
      case default
         stop 'FLUESSE: Richtung falsch gewaehlt!'
      end select


c      print *,u(30,1/2+1,1/2+1,3)**2,u(30,1/2+1,1/2+1,2)**2,
c     +     u(30,1/2+1,1/2+1,1)**2
      k=1
         j=1
            do i=i_min,i_max
               eps(i,j,k)=u(i,j,k,4)-
     +              0.5d0*(u(i,j,k,1)**2+u(i,j,k,2)**2+
     +              u(i,j,k,3)**2)/u(i,j,k,0)
               eps(i,j,k)=max(eps(i,j,k),0.0d0)
               call eos(u(i,j,k,0),eps(i,j,k),p(i,j,k),c_s(i,j,k),2)
c     Addiere potentielle Energie:
               u(i,j,k,4)=u(i,j,k,4)+u(i,j,k,0)*phi_grav(i,j,k)
               v=u(i,j,k,richtung)/u(i,j,k,0)
               fl_u(i,j,k,0)=u(i,j,k,richtung)
               fl_u(i,j,k,1)=v*u(i,j,k,1)
               fl_u(i,j,k,2)=v*u(i,j,k,2)
               fl_u(i,j,k,3)=v*u(i,j,k,3)
               fl_u(i,j,k,richtung)=fl_u(i,j,k,richtung)+p(i,j,k)
     +              +0.0d0*phi_grav(i,j,k)
               fl_u(i,j,k,4)=v*(u(i,j,k,4)+p(i,j,k))
            end do

c     Phi-Impulsfluss zu Drehimpulsfluss konvertieren:
               fl_u(:,j,k,3)=fl_u(:,j,k,3)*r_sin_theta(:,j,k)


      return

      END SUBROUTINE fluesse

c     =============================================================== 



 



c     =============================================================== 

      SUBROUTINE zeitschritt(c)

c     =============================================================== 

      USE hydro_var
      USE gitter_var

      IMPLICIT NONE

      real,dimension(1:800,1:1,1:1),intent(in) :: c

      integer i,j,k

      dt=1.0d9
      do i=1,800
         dt=min(dt,dr(i)/c(i,1,1))
      end do
      dt=cfl_factor*dt

      END SUBROUTINE zeitschritt 

c     =============================================================== 




c     =============================================================== 

      SUBROUTINE hydro_speichern

c     =============================================================== 

      USE hydro_var
      USE hydro_sav_var

      IMPLICIT NONE

      integer j,k

      do k=1,1
         do j=1,1
            u_sav(1:800,j,k,0)=u(1:800,j,k,0)
         end do
      end do

      END SUBROUTINE hydro_speichern

c     =============================================================== 




c     =============================================================== 

      SUBROUTINE courant_factor

c     =============================================================== 

      USE hydro_var
      USE hydro_sav_var
      USE mod_ausgabe, ONLY : ausgabe

      IMPLICIT NONE

      real :: drho_rho(1:800,1:1,1:1),de_e(1:800,1:1,1:1)

      real :: drho_rho_max,de_e_max
      integer :: i,j,k

      do k=1,1
         do j=1,1
            drho_rho(1:800,j,k)=
     +           abs(u(1:800,j,k,0)-u_sav(1:800,j,k,0))/
     +           u(1:800,j,k,0)
         end do
      end do

      drho_rho_max=0.0d0
      do k=1,1
         do j=1,1
            do i=1,800
               drho_rho_max=max(drho_rho_max,drho_rho(i,j,k))
            end do
         end do
      end do


      if (drho_rho_max.gt.6.0d-3) then
         cfl_factor=cfl_factor*4.0d-3/drho_rho_max
      else if(drho_rho_max.lt.3.0d-3) then
         cfl_factor=min(cfl_factor*(3.0d-3/drho_rho_max),0.4d0)
      end if

      cfl_factor=max(cfl_factor,1.0d-4)

      if (cfl_factor.lt.1.0d-3) then
         CALL ausgabe(-999)
c         stop 'CFL-Faktor zu klein!'
      end if

      END SUBROUTINE courant_factor

c     =============================================================== 



      END MODULE advection

c     ===============================================================
