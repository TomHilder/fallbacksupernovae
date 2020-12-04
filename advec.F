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

      real, dimension(-1:NX+2,1:NY,1:NZ) :: c
      real, dimension(-1:NX+2,1:NY,1:NZ,0:4) :: w,
     +     slope_r,slope_l,u_r,u_l
      real, dimension(-1:NX+1,1:NY,1:NZ,0:4) :: fl
      real, dimension(-1:NX+2,1:NY,1:NZ,0:4) :: u_alt
      real :: vx,vy,vz,rho

      integer :: i,j,k,i_rk,l
c     ---------------------------------------------------------------
c     Erster Runge-Kutta-Schritt
c     ---------------------------------------------------------------
      
c     Randbedingungen setzen:
      call randbed(1)
      call fluesse(1)

#if NDIM==3
!$OMP PARALLEL DO
!$OMP& PRIVATE(i,j,k)
#endif
      do k=1,NZ
#if NDIM==2
!$OMP PARALLEL DO
!$OMP& PRIVATE(i,j)
#endif
         do j=1,NY
            do i=-1,NX+2
#ifdef ATMOSPHAERE
               if(u(i,j,k,0).le.rho_min1) then
c                  u(i,j,k,0)=rho_min2
                  u(i,j,k,1)=0.0d0
c                  u(i,j,k,4)=eps_rand
               end if
#endif
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
         call zeitschritt(c(1:NX,1:NY,1:NZ))
      else
         dt=0.5d0*dt
      end if

#if 1
#if NDIM==3
!$OMP PARALLEL DO
!$OMP& PRIVATE(i,j,k,l,rho,vx,vy,vz)         
#endif
      do k=1,NZ
#if NDIM==2
!$OMP PARALLEL DO
!$OMP& PRIVATE(i,j,l,rho,vx,vy,vz)         
#endif
         do j=1,NY

            do i=-1,NX+1
               fl(i,j,k,:)=
     +              0.25d0*(c(i+1,j,k)+c(i,j,k))*
     +              (u(i,j,k,:)+w(i,j,k,:))
            end do

            do i=-1,NX+1
               fl(i,j,k,:)=fl(i,j,k,:)-
     +              0.25d0*(c(i+1,j,k)+c(i,j,k))*
     +              (u(i+1,j,k,:)-w(i+1,j,k,:))
            end do

c     Updates der erhaltenene Groessen:
            do i=0,NX+1
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
#if GEOM==1
c     Scheinkraefte:
     +              r1(i)*(vy**2+vz**2-2.0d0*vx**2)+
#endif
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

#endif

#if 1
c     ---------------------------------------------------------------
c     Zweiter Runge-Kutta-Schritt
c     ---------------------------------------------------------------

      dt=2.0d0*dt

c     Randbedingungen setzen:
      call randbed(1)
      call fluesse(1)

c     Flux-Freezing-Geschwindigkeit:
#if NDIM==3
!$OMP PARALLEL DO
!$OMP& PRIVATE(i,j,k,l,rho,vx,vy,vz)         
#endif
      do k=1,NZ
#if NDIM==2
!$OMP PARALLEL DO
!$OMP& PRIVATE(i,j,l,rho,vx,vy,vz)         
#endif
         do j=1,NY
            do i=-1,NX+2
               c(i,j,k)=c_s(i,j,k)+
     +              sqrt(u(i,j,k,1)**2+u(i,j,k,2)**2+u(i,j,k,3)**2)/
     +              u(i,j,k,0)
               w(i,j,k,:)=fl_u(i,j,k,:)/c(i,j,k)
               u_r(i,j,k,:)=0.5d0*(u(i,j,k,:)+w(i,j,k,:))
               u_l(i,j,k,:)=0.5d0*(u(i,j,k,:)-w(i,j,k,:))
            end do

            do l=0,4
               do i=0,NX+1
                  slope_r(i,j,k,l)=van_leer(
     +                 (u_r(i,j,k,l)-u_r(i-1,j,k,l))*dr_if1(i-1),
     +                 (u_r(i+1,j,k,l)-u_r(i,j,k,l))*dr_if1(i))
                  slope_l(i,j,k,l)=van_leer(
     +                 (u_l(i,j,k,l)-u_l(i-1,j,k,l))*dr_if1(i-1),
     +                 (u_l(i+1,j,k,l)-u_l(i,j,k,l))*dr_if1(i))
               end do
            end do

            do i=0,NX
               fl(i,j,k,:)=
     +              0.5d0*(c(i+1,j,k)+c(i,j,k))*
     +              (u_r(i,j,k,:)+0.5d0*slope_r(i,j,k,:)*dr(i))
            end do

            do i=0,NX
               fl(i,j,k,:)=fl(i,j,k,:)-
     +              0.5d0*(c(i+1,j,k)+c(i,j,k))*
     +              (u_l(i+1,j,k,:)-0.5d0*slope_l(i+1,j,k,:)*dr(i+1))
            end do

c     Updates der erhaltenene Groessen:
            do i=1,NX
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
#if GEOM==1
c     Scheinkraefte:
     +              r1(i)*(vy**2+vz**2-2.0d0*vx**2)+
#endif
c     Gravitationskraft:
     +              a_grv(i,j,k,1)))
               u(i,j,k,2)=u_alt(i,j,k,2)+
     +              dt*(fl(i-1,j,k,2)*da(i-1,j,k,1)-
     +              fl(i,j,k,2)*da(i,j,k,1))*dv1(i,j,k)
               u(i,j,k,3)=u_alt(i,j,k,3)+0.0d0*
     +              dt*(fl(i-1,j,k,3)*da(i-1,j,k,1)-
     +              fl(i,j,k,3)*da(i,j,k,1))*dv1(i,j,k)*
     +              r_sin_theta1(i,j,k)
#ifdef ATMOSPHARE
               if(u(i,j,k,0).le.rho_min1) then
c                  u(i,j,k,0)=rho_min2
                  u(i,j,k,1)=0.0d0
c                  u(i,j,k,4)=eps_rand
               end if
#endif
            end do
         end do
      end do
#endif

#ifdef N_R_SPHAER
      do i=1,N_R_SPHAER
         u(i,1:NY,1:NZ,0)=sum(u(i,1:NY,1:NZ,0)*dv(i,1:NY,1:NZ))/
     +        sum(dv(i,1:NY,1:NZ))
         u(i,1:NY,1:NZ,1)=sum(u(i,1:NY,1:NZ,1)*dv(i,1:NY,1:NZ))/
     +        sum(dv(i,1:NY,1:NZ))
         u(i,1:NY,1:NZ,4)=sum(u(i,1:NY,1:NZ,4)*dv(i,1:NY,1:NZ))/
     +        sum(dv(i,1:NY,1:NZ))
         u(i,1:NY,1:NZ,2)=0.0d0
         u(i,1:NY,1:NZ,3)=0.0d0
      end do
#endif

#ifdef CONS_TEST
      do k=1,NZ
#if NDIM==2
!$OMP PARALLEL DO
!$OMP& PRIVATE(j)
!$OMP& REDUCTION(+:verlust_mtot)
!$OMP& REDUCTION(+:verlust_etot)
#endif
         do j=1,NY
            verlust_mtot=verlust_mtot+dt*fl(NX,j,k,0)*da(NX,j,k,1)
            verlust_etot=verlust_etot+dt*fl(NX,j,k,4)*da(NX,j,k,1)
            mloss_inner=mloss_inner-dt*fl(0,j,k,0)*da(0,j,k,1)
            m_core=m_core-dt*fl(0,j,k,0)*da(0,j,k,1)
         end do
#if NDIM==2
!$OMP END PARALLEL DO
#endif
      end do
#endif

      if (dt_modus.eq.0) zeit=zeit+dt

      return

      END SUBROUTINE advec_x

c     ===============================================================


#if NDIM>=2
c     ===============================================================

      SUBROUTINE advec_y(dt_modus)

c     =============================================================

      USE grav_var
      USE hydro_var
      USE gitter_var

      IMPLICIT NONE

      integer, intent(in) :: dt_modus

      real, dimension(1:NX,-1:NY+2,1:NZ) :: c
      real, dimension(1:NX,-1:NY+2,1:NZ,0:4) :: w,
     +     slope_r,slope_l,u_r,u_l
      real, dimension(1:NX,-1:NY+1,1:NZ,0:4) :: fl
      real, dimension(1:NX,-1:NY+2,1:NZ,0:4) :: u_alt
      real :: vx,vy,vz,rho

      integer :: i,j,k,i_rk,l

c     ---------------------------------------------------------------
c     Erster Runge-Kutta-Schritt
c     ---------------------------------------------------------------
      
c     Randbedingungen setzen:
      call randbed(2)
      call fluesse(2)
#if NDIM==3
!$OMP PARALLEL DO
!$OMP& PRIVATE(i,j,k,l,rho,vx,vy,vz)         
#endif
      do k=1,NZ
#if NDIM==2
!$OMP PARALLEL DO
!$OMP& PRIVATE(i,j,l,rho,vx,vy,vz)         
#endif
         do j=-1,NY+2
            do i=1,NX
#ifdef ATMOSPHAERE
               if(u(i,j,k,0).le.rho_min1) then
c                  u(i,j,k,0)=rho_min2
                  u(i,j,k,1)=0.0d0
c                  u(i,j,k,4)=eps_rand
               end if
#endif
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
         end do
      end do

      if (dt_modus.eq.0) then
         call zeitschritt(c(1:NX,1:NY,1:NZ))
      else
         dt=0.5d0*dt
      end if

#if NDIM==3
!$OMP PARALLEL DO
!$OMP& PRIVATE(i,j,k,l,rho,vx,vy,vz)         
#endif
      do k=1,NZ
#if NDIM==2
!$OMP PARALLEL DO
!$OMP& PRIVATE(i,j,l,rho,vx,vy,vz)         
#endif
         do j=-1,NY+1
            do i=1,NX
               fl(i,j,k,:)=
     +              0.25d0*(c(i,j+1,k)+c(i,j,k))*
     +              (u(i,j,k,:)+w(i,j,k,:))
               fl(i,j,k,:)=(fl(i,j,k,:)-
     +              0.25d0*(c(i,j+1,k)+c(i,j,k))*
     +              (u(i,j+1,k,:)-w(i,j+1,k,:)))
            end do
         end do

c     Updates der erhaltenene Groessen:
#if NDIM==2
!$OMP PARALLEL DO
!$OMP& PRIVATE(i,j,l,rho,vx,vy,vz)         
#endif
         do j=-1,NY+1
            do i=1,NX
               rho=u(i,j,k,0)
               vx=u(i,j,k,1)/rho
               vy=u(i,j,k,2)/rho
               vz=u(i,j,k,3)/rho

               u(i,j,k,0)=u_alt(i,j,k,0)+
     +              dt*(fl(i,j-1,k,0)*da(i,j-1,k,2)-
     +              fl(i,j,k,0)*da(i,j,k,2))*dv1(i,j,k)
               u(i,j,k,4)=u_alt(i,j,k,4)+
     +              dt*(fl(i,j-1,k,4)*da(i,j-1,k,2)-
     +              fl(i,j,k,4)*da(i,j,k,2))*dv1(i,j,k)-
c     Subtrahiere potentielle Energie:
     +              u(i,j,k,0)*phi_grav(i,j,k)
               u(i,j,k,1)=u_alt(i,j,k,1)+
     +              dt*(fl(i,j-1,k,1)*da(i,j-1,k,2)-
     +              fl(i,j,k,1)*da(i,j,k,2))*dv1(i,j,k)
               u(i,j,k,2)=u_alt(i,j,k,2)+dt*r1(i)*
     +              ((fl(i,j-1,k,2)-fl(i,j,k,2))*dtheta1(j)+
     +              rho*(
#if GEOM==1
c     Scheinkraefte:
     +              (vz**2-vy**2)*cot_theta(j)-vx*vy)+
#endif
c     Gravitationskraft:
     +              a_grv(i,j,k,2))
               u(i,j,k,3)=u_alt(i,j,k,3)+
     +              dt*(fl(i,j-1,k,3)*da(i,j-1,k,2)-
     +              fl(i,j,k,3)*da(i,j,k,2))*dv1(i,j,k)*
     +              r_sin_theta1(i,j,k)

#ifdef ATMOSPHARE
               if(u(i,j,k,0).le.rho_min1) then
c                  u(i,j,k,0)=rho_min2
                  u(i,j,k,1)=0.0d0
c                  u(i,j,k,4)=eps_rand
               end if
#endif
            end do
         end do
      end do

#if 1
c     ---------------------------------------------------------------
c     Zweiter Runge-Kutta-Schritt
c     ---------------------------------------------------------------

      dt=2.0d0*dt

c     Randbedingungen setzen:
      call randbed(2)
      call fluesse(2)

c     Flux-Freezing-Geschwindigkeit:
#if NDIM==3
!$OMP PARALLEL DO
!$OMP& PRIVATE(i,j,k,rho,vx,vy,vz)         
#endif
      do k=1,NZ
#if NDIM==2
!$OMP PARALLEL DO
!$OMP& PRIVATE(i,j,l,rho,vx,vy,vz)         
#endif
         do j=-1,NY+2
            do i=1,NX
               c(i,j,k)=c_s(i,j,k)+
     +              sqrt(u(i,j,k,1)**2+u(i,j,k,2)**2+u(i,j,k,3)**2)/
     +              u(i,j,k,0)
               w(i,j,k,:)=fl_u(i,j,k,:)/c(i,j,k)
               u_r(i,j,k,:)=0.5d0*(u(i,j,k,:)+w(i,j,k,:))
               u_l(i,j,k,:)=0.5d0*(u(i,j,k,:)-w(i,j,k,:))
            end do
         end do

#if NDIM==2
!$OMP PARALLEL DO
!$OMP& PRIVATE(i,j,l)         
#endif
         do j=0,NY+1
            do l=0,4
               do i=1,NX
                  slope_r(i,j,k,l)=van_leer(
     +                 (u_r(i,j,k,l)-u_r(i,j-1,k,l))*dtheta_if1(j-1),
     +                 (u_r(i,j+1,k,l)-u_r(i,j,k,l))*dtheta_if1(j))
                  slope_l(i,j,k,l)=van_leer(
     +                 (u_l(i,j,k,l)-u_l(i,j-1,k,l))*dtheta_if1(j-1),
     +                 (u_l(i,j+1,k,l)-u_l(i,j,k,l))*dtheta_if1(j))
               end do
            end do
         end do

#if NDIM==2
!$OMP PARALLEL DO
!$OMP& PRIVATE(i,j,l,rho,vx,vy,vz)         
#endif
         do j=0,NY
            do i=1,NX
               fl(i,j,k,:)=
     +              0.5d0*(c(i,j+1,k)+c(i,j,k))*
     +              (u_r(i,j,k,:)+0.5d0*slope_r(i,j,k,:)*dtheta(j))
               fl(i,j,k,:)=fl(i,j,k,:)-
     +              0.5d0*(c(i,j+1,k)+c(i,j,k))*
     +              (u_l(i,j+1,k,:)-0.5d0*slope_l(i,j+1,k,:)*
     +              dtheta(j+1))
            end do
         end do

c     Updates der erhaltenene Groessen:
#if NDIM==2
!$OMP PARALLEL DO
!$OMP& PRIVATE(i,j,l,rho,vx,vy,vz)         
#endif
         do j=1,NY
            do i=1,NX
               rho=u(i,j,k,0)
               vx=u(i,j,k,1)/rho
               vy=u(i,j,k,2)/rho
               vz=u(i,j,k,3)/rho

               u(i,j,k,0)=u_alt(i,j,k,0)+
     +              dt*(fl(i,j-1,k,0)*da(i,j-1,k,2)-
     +              fl(i,j,k,0)*da(i,j,k,2))*dv1(i,j,k)
               u(i,j,k,4)=u_alt(i,j,k,4)+
     +              dt*(fl(i,j-1,k,4)*da(i,j-1,k,2)-
     +              fl(i,j,k,4)*da(i,j,k,2))*dv1(i,j,k)-
c     Subtrahiere potentielle Energie:
     +              u(i,j,k,0)*phi_grav(i,j,k)
               u(i,j,k,1)=u_alt(i,j,k,1)+
     +              dt*(fl(i,j-1,k,1)*da(i,j-1,k,2)-
     +              fl(i,j,k,1)*da(i,j,k,2))*dv1(i,j,k)
               u(i,j,k,2)=u_alt(i,j,k,2)+dt*r1(i)*
     +              ((fl(i,j-1,k,2)-fl(i,j,k,2))*dtheta1(j)+
     +              rho*(
#if GEOM==1
c     Scheinkraefte:
     +              (vz**2-vy**2)*cot_theta(j)-vx*vy)+
#endif
c     Gravitationskraft:
     +              a_grv(i,j,k,2))
               u(i,j,k,3)=u_alt(i,j,k,3)+
     +              dt*(fl(i,j-1,k,3)*da(i,j-1,k,2)-
     +              fl(i,j,k,3)*da(i,j,k,2))*dv1(i,j,k)*
     +              r_sin_theta1(i,j,k)

            end do
         end do
      end do
#endif

#ifdef N_R_SPHAER
      do i=1,N_R_SPHAER
         u(i,1:NY,1:NZ,0)=sum(u(i,1:NY,1:NZ,0)*dv(i,1:NY,1:NZ))/
     +        sum(dv(i,1:NY,1:NZ))
         u(i,1:NY,1:NZ,1)=sum(u(i,1:NY,1:NZ,1)*dv(i,1:NY,1:NZ))/
     +        sum(dv(i,1:NY,1:NZ))
         u(i,1:NY,1:NZ,4)=sum(u(i,1:NY,1:NZ,4)*dv(i,1:NY,1:NZ))/
     +        sum(dv(i,1:NY,1:NZ))
         u(i,1:NY,1:NZ,2)=0.0d0
         u(i,1:NY,1:NZ,3)=0.0d0
      end do
#endif

      if (dt_modus.eq.0) zeit=zeit+dt

      return

      END SUBROUTINE advec_y

c     ===============================================================
#endif


#if NDIM==3
c     ===============================================================

      SUBROUTINE advec_z(dt_modus)

c     ===============================================================

      USE grav_var
      USE hydro_var
      USE gitter_var

      IMPLICIT NONE

      integer, intent(in) :: dt_modus

      real, dimension(1:NX,1:NY,-1:NZ+2) :: c
      real, dimension(1:NX,1:NY,-1:NZ+2,0:4) :: w,
     +     slope_r,slope_l,u_r,u_l
      real, dimension(1:NX,1:NY,-1:NZ+1,0:4) :: fl
      real, dimension(-1:NX+2,-1:NY+2,-1:NZ+2,0:4) :: u_alt
      real :: vx,vy,vz,rho

      integer :: i,j,k,i_rk,l

c     ---------------------------------------------------------------
c     Erster Runge-Kutta-Schritt
c     ---------------------------------------------------------------
      
c     Randbedingungen setzen:
      call randbed(3)
      call fluesse(3)


!$OMP PARALLEL DO
!$OMP& PRIVATE(i,j,k)
      do k=-1,NZ+1
         do j=1,NY
            do i=1,NX
#ifdef ATMOSPHAERE
               if(u(i,j,k,0).le.rho_min1) then
c                  u(i,j,k,0)=rho_min2
                  u(i,j,k,1)=0.0d0
c                  u(i,j,k,4)=eps_rand
               end if
#endif
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
         end do
      end do

      if (dt_modus.eq.0) then
         call zeitschritt(c(1:NX,1:NY,1:NZ))
      else
         dt=0.5d0*dt
      end if

!$OMP PARALLEL DO
!$OMP& PRIVATE(i,j,k)
      do k=-1,NZ+1
         do j=1,NY
            do i=1,NX
               fl(i,j,k,:)=
     +              0.25d0*(c(i,j,k+1)+c(i,j,k))*
     +              (u(i,j,k,:)+w(i,j,k,:))
c            end do
c         end do
c      end do

c      do k=-1,NZ+1
c         do j=1,NY
c            do i=1,NX
               fl(i,j,k,:)=fl(i,j,k,:)-
     +              0.25d0*(c(i,j,k+1)+c(i,j,k))*
     +              (u(i,j,k+1,:)-w(i,j,k+1,:))
            end do
         end do
      end do
         

c     Updates der erhaltenene Groessen:
!$OMP PARALLEL DO
!$OMP& PRIVATE(i,j,k)
      do k=-1,NZ+1
         do j=1,NY
            do i=1,NX

               u(i,j,k,0)=u_alt(i,j,k,0)+
     +              dt*(fl(i,j,k-1,0)*da(i,j,k-1,3)-
     +              fl(i,j,k,0)*da(i,j,k,3))*dv1(i,j,k)
               u(i,j,k,4)=u_alt(i,j,k,4)+
     +              dt*(fl(i,j,k-1,4)*da(i,j,k-1,3)-
     +              fl(i,j,k,4)*da(i,j,k,3))*dv1(i,j,k)-
c     Subtrahiere potentielle Energie:
     +              u(i,j,k,0)*phi_grav(i,j,k)
               u(i,j,k,1)=u_alt(i,j,k,1)+
     +              dt*(fl(i,j,k-1,1)*da(i,j,k-1,3)-
     +              fl(i,j,k,1)*da(i,j,k,3))*dv1(i,j,k)
               u(i,j,k,2)=u_alt(i,j,k,2)+
     +              dt*(fl(i,j,k-1,2)*da(i,j,k-1,3)-
     +              fl(i,j,k,2)*da(i,j,k,3))*dv1(i,j,k)
               u(i,j,k,3)=u_alt(i,j,k,3)+
     +              dt*((fl(i,j,k-1,3)-fl(i,j,k,3))*dphi1(k)*
     +              r_sin_theta1(i,j,k)+
c     Gravitationskraft:
     +              a_grv(i,j,k,3)*u(i,j,k,0))

            end do
         end do
      end do


#if 1
c     ---------------------------------------------------------------
c     Zweiter Runge-Kutta-Schritt
c     ---------------------------------------------------------------

      dt=2.0d0*dt

c     Randbedingungen setzen:
      call randbed(3)
      call fluesse(3)

c     Flux-Freezing-Geschwindigkeit:
!$OMP PARALLEL DO
!$OMP& PRIVATE(i,j,k)
      do k=-1,NZ+2
         do j=1,NY
            do i=1,NX
               c(i,j,k)=c_s(i,j,k)+
     +              sqrt(u(i,j,k,1)**2+u(i,j,k,2)**2+u(i,j,k,3)**2)/
     +              u(i,j,k,0)
               w(i,j,k,:)=fl_u(i,j,k,:)/c(i,j,k)
               u_r(i,j,k,:)=0.5d0*(u(i,j,k,:)+w(i,j,k,:))
               u_l(i,j,k,:)=0.5d0*(u(i,j,k,:)-w(i,j,k,:))
            end do
         end do
      end do

!$OMP PARALLEL DO      
!$OMP& PRIVATE(i,j,k,l)
      do k=0,NZ+1
         do j=1,NY
            do l=0,4
               do i=1,NX
                  slope_r(i,j,k,l)=van_leer(
     +                 (u_r(i,j,k,l)-u_r(i,j,k-1,l))/dphi_if(k-1),
     +                 (u_r(i,j,k+1,l)-u_r(i,j,k,l))/dphi_if(k))
                  slope_l(i,j,k,l)=van_leer(
     +                 (u_l(i,j,k,l)-u_l(i,j,k-1,l))/dphi_if(k-1),
     +                 (u_l(i,j,k+1,l)-u_l(i,j,k,l))/dphi_if(k))
               end do
            end do
         end do
      end do

!$OMP PARALLEL DO
!$OMP& PRIVATE(i,j,k)
      do k=0,NZ
         do j=1,NY
            do i=1,NX
               fl(i,j,k,:)=
     +              0.5d0*(c(i,j,k+1)+c(i,j,k))*
     +              (u_r(i,j,k,:)+0.5d0*slope_r(i,j,k,:)*dphi(k))
c            end do
c         end do
c      end do

c      do k=0,NZ
c         do j=1,NY
c            do i=1,NX
               fl(i,j,k,:)=fl(i,j,k,:)-
     +              0.5d0*(c(i,j,k+1)+c(i,j,k))*
     +              (u_l(i,j,k+1,:)-0.5d0*slope_l(i,j,k+1,:)*dphi(k+1))
            end do
         end do
      end do

c     Updates der erhaltenene Groessen:
!$OMP PARALLEL DO
!$OMP& PRIVATE(i,j,k,rho,vx,vy,vz)
      do k=1,NZ
         do j=1,NY
            do i=1,NX
               rho=u(i,j,k,0)
               vx=u(i,j,k,1)/rho
               vy=u(i,j,k,2)/rho
               vz=u(i,j,k,3)/rho

               u(i,j,k,0)=u_alt(i,j,k,0)+
     +              dt*(fl(i,j,k-1,0)*da(i,j,k-1,3)-
     +              fl(i,j,k,0)*da(i,j,k,3))
               u(i,j,k,4)=u_alt(i,j,k,4)+
     +              dt*(fl(i,j,k-1,4)*da(i,j,k-1,3)-
     +              fl(i,j,k,4)*da(i,j,k,3))-
c     Subtrahiere potentielle Energie:
     +              u(i,j,k,0)*phi_grav(i,j,k)
               u(i,j,k,1)=u_alt(i,j,k,1)+
     +              dt*(fl(i,j,k-1,1)*da(i,j,k-1,3)-
     +              fl(i,j,k,1)*da(i,j,k,3))
               u(i,j,k,2)=u_alt(i,j,k,2)+
     +              dt*(fl(i,j,k-1,2)*da(i,j,k-1,3)-
     +              fl(i,j,k,2)*da(i,j,k,3))*dv1(i,j,k)
               u(i,j,k,3)=u_alt(i,j,k,3)+
     +              dt*((fl(i,j,k-1,3)-fl(i,j,k,3))*dphi1(k)*
     +              r_sin_theta1(i,j,k)+
c     Gravitationskraft:
     +              a_grv(i,j,k,3)*u(i,j,k,0))

            end do
         end do
      end do
#endif

      if (dt_modus.eq.0) zeit=zeit+dt

      return

      END SUBROUTINE advec_z

c     ===============================================================
#endif


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
         do k=1,NZ
#if NDIM==2
!$OMP PARALLEL DO
!$OMP& PRIVATE(j)
!$OMP& SHARED(k,u,r)
#endif
            do j=1,NY
c     Hydro-Groessen:
#if RAND_MODUS<2
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
#endif

#if RAND_MODUS==0
c     Inflow-Randbedingungen:
               u(NX+1,j,k,:)=((r(NX+1)-r(NX-1))*u(NX-1,j,k,:)+
     +              (r(NX)-r(NX+1))*u(NX,j,k,:))/
     +              (r(NX)-r(NX-1))
               u(NX+2,j,k,:)=((r(NX+2)-r(NX-1))*u(NX-1,j,k,:)+
     +              (r(NX)-r(NX+2))*u(NX,j,k,:))/
     +              (r(NX)-r(NX-1))
               u(NX+1,j,k,:)=u(NX,j,k,:)
               u(NX+2,j,k,:)=u(NX,j,k,:)
               u(NX+1,j,k,1)=u(NX,j,k,1)*r(NX)**2/r(NX+1)**2
               u(NX+2,j,k,1)=u(NX,j,k,1)*r(NX)**2/r(NX+2)**2
#endif
#if RAND_MODUS==1
c     Reflektierende Randbedingungen:
               u(NX+1,j,k,0)=u(NX-1,j,k,0)
               u(NX+2,j,k,0)=u(NX-2,j,k,0)
c         u(NX+1,j,k,0)=rho_min1!u(NX-1,j,k,0)
c         u(NX+2,j,k,0)=rho_min1!u(NX-2,j,k,0)

               u(NX+1,j,k,1:3)=-u(NX-1,j,k,1:3)*(r(NX-1)/r(NX+1))**2
               u(NX+2,j,k,1:3)=-u(NX-2,j,k,1:3)*(r(NX-2)/r(NX+2))**2
c         u(NX+1,j,k,1:3)=0.0d0
c         u(NX+2,j,k,1:3)=0.0d0

               u(NX+1,j,k,4)=u(NX-1,j,k,4)
               u(NX+2,j,k,4)=u(NX-2,j,k,4)
c         u(NX+1,j,k,4)=eps_rand!u(NX-1,j,k,4)
c         u(NX+2,j,k,4)=eps_rand!u(NX-2,j,k,4)
#endif
#if RAND_MODUS==2
c     Periodische Randbedingungen
               u( 0  ,j,k,:)=u(NX  ,j,k,:)
               u(-1  ,j,k,:)=u(NX-1,j,k,:)
               u(NX+1,j,k,:)=u(1   ,j,k,:)
               u(NX+2,j,k,:)=u(2   ,j,k,:)

#endif


c     Gravitationsfeld:
               phi_grav(0,j,k)=phi_grav(1,j,k)
               phi_grav(-1,j,k)=phi_grav(2,j,k)
               phi_grav(NX+1,j,k)=phi_grav(NX,j,k)
               phi_grav(NX+2,j,k)=phi_grav(NX,j,k)
               
            end do
         end do

#if NDIM>=2
c     ---------------------------------------------------------------
c     theta-Richtung:
      case(2)
         do k=1,NZ
#if NDIM==2
!$OMP PARALLEL DO
!$OMP& PRIVATE(i)
!$OMP& SHARED(j,k,u)
#endif
            do i=1,NX
c     Hydro-Groessen:
               u(i, 0,k,0)= u(i,1,k,0)
               u(i,-1,k,0)= u(i,2,k,0)
               u(i, 0,k,1)= u(i,1,k,1)
               u(i,-1,k,1)= u(i,2,k,1)
               u(i, 0,k,2)=-u(i,1,k,2)
               u(i,-1,k,2)=-u(i,2,k,2)
               u(i, 0,k,3)= u(i,1,k,3)
               u(i,-1,k,3)= u(i,2,k,3)
               u(i, 0,k,4)= u(i,1,k,4)
               u(i,-1,k,4)= u(i,2,k,4)

               u(i,NY+1,k,0)= u(i,NY  ,k,0)
               u(i,NY+2,k,0)= u(i,NY-1,k,0)
               u(i,NY+1,k,1)= u(i,NY  ,k,1)
               u(i,NY+2,k,1)= u(i,NY-1,k,1)
#ifdef AEQUATOR
               u(i,NY+1,k,2)= u(i,NY  ,k,2)
               u(i,NY+2,k,2)= u(i,NY-1,k,2)
#else 
               u(i,NY+1,k,2)=-u(i,NY  ,k,2)
               u(i,NY+2,k,2)=-u(i,NY-1,k,2)
#endif
               u(i,NY+1,k,3)= u(i,NY  ,k,3)
               u(i,NY+2,k,3)= u(i,NY-1,k,3)
               u(i,NY+1,k,4)= u(i,NY  ,k,4)
               u(i,NY+2,k,4)= u(i,NY-1,k,4)

c     Gravitationsfeld:
               phi_grav(i,   0,k)=phi_grav(i, 1,k)
               phi_grav(i,  -1,k)=phi_grav(i, 2,k)
               phi_grav(i,NY+1,k)=phi_grav(i,NY,k)
               phi_grav(i,NY+2,k)=phi_grav(i,NY,k)

            end do
         end do
#endif


#if NDIM==3
c     ---------------------------------------------------------------
c     phi-Richtung:
      case(3)
!$OMP PARALLEL DO
!$OMP& PRIVATE(i,j)
         do j=1,NY
            do i=1,NX
c     Hydro-Groessen:
               u(i,j, 0,0)= u(i,j,NZ  ,0)
               u(i,j,-1,0)= u(i,j,NZ-1,0)
               u(i,j, 0,1)= u(i,j,NZ  ,1)
               u(i,j,-1,1)= u(i,j,NZ-1,1)
               u(i,j, 0,1)=-u(i,j,NZ  ,1)
               u(i,j,-1,1)=-u(i,j,NZ-1,1)
               u(i,j, 0,3)= u(i,j,NZ  ,3)
               u(i,j,-1,3)= u(i,j,NZ-1,3)
               u(i,j, 0,4)= u(i,j,NZ  ,4)
               u(i,j,-1,4)= u(i,j,NZ-1,4)
               
               u(i,j,NZ+1,0)= u(i,j,1,0)
               u(i,j,NZ+2,0)= u(i,j,2,0)
               u(i,j,NZ+1,1)= u(i,j,1,1)
               u(i,j,NZ+2,1)= u(i,j,2,1)
               u(i,j,NZ+1,2)= u(i,j,1,2)
               u(i,j,NZ+2,2)= u(i,j,2,2)
               u(i,j,NZ+1,3)= u(i,j,1,3)
               u(i,j,NZ+2,3)= u(i,j,2,3)
               u(i,j,NZ+1,4)= u(i,j,1,4)
               u(i,j,NZ+2,4)= u(i,j,2,4)

c     Graviationsfeld:
c     nocht irrelevant

            end do
         end do
#endif

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
         i_max=NX+2
#if NDIM>=2
         j_min=1
         j_max=NY
#endif
#if NDIM==3
         k_min=1
         k_max=NZ
#endif
#if NDIM>=2
      case(2)

         i_min=1
         i_max=NX
         j_min=-1
         j_max=NY+2
#if NDIM==3
         k_min=1
         k_max=NZ
#endif
#endif
#if NDIM==3
      case(3)
         i_min=1
         i_max=NX
         j_min=1
         j_max=NY
         k_min=-1
         k_max=NZ+2
#endif
      case default
         stop 'FLUESSE: Richtung falsch gewaehlt!'
      end select


c      print *,u(30,NY/2+1,NZ/2+1,3)**2,u(30,NY/2+1,NZ/2+1,2)**2,
c     +     u(30,NY/2+1,NZ/2+1,1)**2
#if NDIM>=3
      do k=k_min,k_max
#else
      k=1
#endif
#if NDIM>=2
!$OMP PARALLEL DO
!$OMP&PRIVATE(i,j,v)
!$OMP&SHARED(k)
         do j=j_min,j_max
#else
         j=1
#endif
            do i=i_min,i_max
#ifdef NEC
               u(i,j,k,0)=max(u(i,j,k),rho_min2)
#endif
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

#if NDIM==3
            if (richtung.lt.3) then
#endif
c     Phi-Impulsfluss zu Drehimpulsfluss konvertieren:
               fl_u(:,j,k,3)=fl_u(:,j,k,3)*r_sin_theta(:,j,k)
#if NDIM==3
            end if
#endif

#if NDIM>=2
         end do
#endif
#if NDIM==3
      end do
#endif

      return

      END SUBROUTINE fluesse

c     =============================================================== 



 



c     =============================================================== 

      SUBROUTINE zeitschritt(c)

c     =============================================================== 

      USE hydro_var
      USE gitter_var

      IMPLICIT NONE

      real,dimension(1:NX,1:NY,1:NZ),intent(in) :: c

      integer i,j,k

      dt=1.0d9
#if NDIM==3
      do k=1,NZ
!$OMP PARALLEL DO
!$OMP& PRIVATE(i,j)
!$OMP& REDUCTION(min:dt)
         do j=1,NY
            do i=1,NX
               dt=min(dt,min(dr(i),r(i)*dtheta(j),
     +              r_sin_theta(i,j,k)*dphi(k))/c(i,j,k))
            end do
         end do
      end do
#endif
#if NDIM==2
!$OMP PARALLEL DO
!$OMP& PRIVATE(i,j)
!$OMP& REDUCTION(min:dt)
      do j=1,NY
         do i=1,NX
            dt=min(dt,min(dr(i),r(i)*dtheta(j))/c(i,j,1))
         end do
      end do
!$OMP END PARALLEL DO
#endif
#if NDIM==1
      do i=1,NX
         dt=min(dt,dr(i)/c(i,1,1))
      end do
#endif
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

      do k=1,NZ
#if NDIM==2
!$OMP PARALLEL DO
!$OMP& PRIVATE(j)
!$OMP& SHARED(k)         
#endif
         do j=1,NY
            u_sav(1:NX,j,k,0)=u(1:NX,j,k,0)
#ifdef ALLES_SPEICHERN
            u_sav(1:NX,j,k,1)=u(1:NX,j,k,1)
            u_sav(1:NX,j,k,2)=u(1:NX,j,k,2)
            u_sav(1:NX,j,k,3)=u(1:NX,j,k,3)
            u_sav(1:NX,j,k,4)=u(1:NX,j,k,4)
#endif
         end do
      end do

      END SUBROUTINE hydro_speichern

c     =============================================================== 


#ifdef ATMOSPHAERE
c     ===============================================================

      SUBROUTINE bestimme_rho_min

c     ===============================================================

      USE hydro_var

      IMPLICIT NONE

      integer :: i,j,k

      rho_max=0.0d0
      do k=1,NZ
!$OMP PARALLEL DO
!$OMP& PRIVATE(i,j)
!$OMP& REDUCTION(max:rho_max)
         do j=1,NY
            do i=1,NX
               rho_max=max(rho_max,u(i,j,k,0))
            end do
         end do
!$OMP END PARALLEL DO
      end do
      rho_min1=RHO_MIN_FAC*rho_max
      rho_min2=0.5d0*rho_min1
c      rho_min2=1.0d-2*rho_min1

      END SUBROUTINE bestimme_rho_min

c     ===============================================================
#endif


c     =============================================================== 

      SUBROUTINE courant_factor

c     =============================================================== 

      USE hydro_var
      USE hydro_sav_var
      USE mod_ausgabe, ONLY : ausgabe

      IMPLICIT NONE

      real :: drho_rho(1:NX,1:NY,1:NZ),de_e(1:NX,1:NY,1:NZ)

      real :: drho_rho_max,de_e_max
      integer :: i,j,k

      do k=1,NZ
#if NDIM==2
!$OMP PARALLEL DO
!$OMP& PRIVATE(j)
!$OMP& SHARED(k)         
#endif
         do j=1,NY
            drho_rho(1:NX,j,k)=
     +           abs(u(1:NX,j,k,0)-u_sav(1:NX,j,k,0))/
     +           u(1:NX,j,k,0)
#ifdef ALLES_SPEICHERN
            de_e(1:NX,j,k)=
     +           abs(u(1:NX,j,k,4)-u_sav(1:NX,j,k,1))/
     +           u(1:NX,j,k,4)
#endif
         end do
      end do

      drho_rho_max=0.0d0
      do k=1,NZ
#if NDIM==2
!$OMP PARALLEL DO
!$OMP& PRIVATE(j)
!$OMP& SHARED(k) 
!$OMP& REDUCTION(max:drho_rho_max)        
#endif
         do j=1,NY
            do i=1,NX
               drho_rho_max=max(drho_rho_max,drho_rho(i,j,k))
#ifdef ALLES_SPEICHERN
               de_e_max=max(de_e_max,de_e(i,j,k))
#endif
            end do
         end do
      end do

#ifdef ALLES_SPEICHERN
      drho_rho_max=max(drho_rho_max,de_e_max)
#endif

      if (drho_rho_max.gt.6.0d-3) then
         cfl_factor=cfl_factor*4.0d-3/drho_rho_max
      else if(drho_rho_max.lt.3.0d-3) then
         cfl_factor=min(cfl_factor*(3.0d-3/drho_rho_max),CFL_KRIT)
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
