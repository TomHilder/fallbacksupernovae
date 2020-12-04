#define RHO1 1.0d0
#define RHO4 0.125d0
#define DRU1 1.0d0
#define DRU4 0.1d0
#define ZENTRALDICHTE 1.0d7
#define CORE_MASS 1.5d0

c     ============================================================

      SUBROUTINE hydro_initialisieren

c     ============================================================

      USE hydro_var
      USE grav_var
      USE gitter_var
      USE param_var
      USE advection, ONLY : randbed

      IMPLICIT NONE

      character(256) strrd

      real, dimension(0:NX) :: m
      real :: dist2
      real :: gamma_tmp,k_poly_tmp
      real :: lambda !pressure scale height
     
      integer :: i,j,k,i_aus

      u=0.0d0
      eps_rand=0.0d0

#if EOS_TYP==0
      gamma_tmp=5.0d0/3.0d0
      k_poly_tmp=k_poly_e
#endif
#if EOS_TYP>=1
      gamma_tmp=4.00/3.0d0
      k_poly_tmp=k_poly_e
#endif

#if EOS_TYP==2
      p_nuc=(k_poly_e*rho_nuc**gamma_e)
      k_poly_b=p_nuc/rho_nuc**gamma_b
      eps_nuc_rho=p_nuc*(gamma_e_m1_r-gamma_b_m1_r)/
     +     rho_nuc
#endif

#if GEOM==0
      u(:,:,:,0)=RHO1
      u(NX/4:3*NX/4,:,:,0)=RHO4

      u(:,:,:,4)=1.5d0*DRU1
      u(NX/4:3*NX/4,:,:,4)=1.5d0*DRU4

      a_grv(:,:,:,:)=0.0d0
      phi_grav(:,:,:)=0.0d0
#else
c     Spherical grid

      rho_max=ZENTRALDICHTE
      rho_min1=RHO_MIN_FAC*rho_max
      rho_min2=0.5d0*rho_min1
      
      u(:,:,:,0)=rho_min2
      u(:,:,:,4)=1.5d0*u(:,:,:,0)

      u(1,:,:,0)=ZENTRALDICHTE
      m(0)=CORE_MASS*msun
      i_aus=NX+1
      lambda = r_if(1) * 0.25d0

      p(1,:,:,0) = u(1,:,:,0)*g_grav*m(0)/r_if(0)**2 ! P = rho*g * Lambda

c     solve for hydrostatic structure
      CYCLUS: do i=1,NX
         if (u(i,1,1,0).lt.rho_min1) then
            u(i,:,:,0)=u(i-1,1,1,0)!rho_min2
            eps_rand=u(i-1,1,1,4)
         end if
         
         p (i,:,:) = p (i-1,:,:)

         u  (i,:,:,4)=p(i,1,1)/(gamma_tmp-1.0d0)

         m(i)=m(i-1)+sum(u(i,1:NY,1:NZ,0)*dv(i,1:NY,1:NZ))*sym_fac
         u(i+1,:,:,0)=u(i,1,1,0)-
     +        dr_if(i)*g_grav*m(i)/r_if(i)**2/gamma_tmp/k_poly_tmp*
     +        u(i,1,1,0)**(2.0d0-gamma_tmp)

      end do CYCLUS


      if (i_aus.le.NX) then
         u(i_aus:NX,:,:,0)=u(i_aus-1,1,1,0)
         u(i_aus:NX,:,:,4)=u(i_aus-1,1,1,4)
#ifndef ATMOSPHARE
         stop 'hydro_initialisieren: Dichte zu gering!'
#endif
      else
         eps_rand=u(NX,1,1,4)
      end if
#ifndef ATMOSPHARE
      rho_min1=u(NX,1,1,0)
      rho_min2=rho_min1
#endif


#if EOS_TYP==1
      p(:,:,:)=k_poly_e*u(:,:,:,0)**gamma_e
      u(:,:,:,4)=k_poly_e*u(:,:,:,0)**gamma_e*gamma_e_m1_r+
     +     k_poly_b*u(:,:,:,0)**gamma_b*gamma_b_m1_r
      u(:,:,:,4)=1.00001d0*u(:,:,:,4)
#endif
#if EOS_TYP==2
      p(:,:,:)=k_poly_e*u(:,:,:,0)**gamma_e
      u(:,:,:,4)=p(:,:,:)*gamma_e_m1_r
c      u(:,:,:,4)=1.00001d0*u(:,:,:,4)
#endif


c     Rotationsgeschwindigkeit:
c      do k=1,NZ
c         do j=1,NY
c            do i=1,i_aus
c               u(i,j,k,3)=r(i)*sin(theta(j))*1.0d-2*u(i,j,k,0)
c               u(i,j,k,4)=u(i,j,k,4)+0.5d0*
c     +              u(i,j,k,3)**2/u(i,j,k,0)
c            end do
c         end do
c      end do


      call randbed(1)
#if NDIM>=2
      call randbed(2)
#endif
#if NDIM==3
      call randbed(3)
#endif
      do k=-1,NZ+2
         do j=-1,NY+2
            do i=-1,NX+2
               eps(i,j,k)=u(i,j,k,4)-
     +              0.5d0*(u(i,j,k,1)**2+u(i,j,k,2)**2+u(i,j,k,3)**2)/
     +              u(i,j,k,0)
               eps(i,j,k)=max(eps(i,j,k),0.0d0)
            end do
         end do
      end do

      print *,'Gesamtmasse: ', m(i_aus-1)

      open(11,file='modell.dat',form='unformatted',
     +     status='replace')
      write(11) NX
      write(11) r(1:NX)
      write(11) u(1:NX,1,1,:),eps(1:NX,1,1)
      close(11)

c      u(1:NX,1,1,1)=u(1:NX,1,1,0)*r(1:NX)

c      do j=1,NY
c         u(:,j,:,0)=(0.5d0+0.5d0*sin(theta(j)))**8*u(:,j,:,0)
c      end do

c      stop

#endif

      print *,'Anfangsmodell bestimmt'

      return

      END SUBROUTINE hydro_initialisieren

c     ============================================================


c     ============================================================

      SUBROUTINE modell_lesen

c     ============================================================

      USE hydro_var
      USE gitter_var

      IMPLICIT NONE

      real,allocatable :: r_tmp(:),u_tmp(:,:),eps_tmp(:)
      real :: xi,xi1
      integer :: nx_tmp
      
      integer :: i,j,k


      open(11,file='modell.dat',form='unformatted',
     +     status='old')
      read(11) nx_tmp
      allocate(r_tmp(0:nx_tmp+1),u_tmp(0:nx_tmp+1,0:4),
     +     eps_tmp(0:nx_tmp+1))
      read(11) r_tmp(1:nx_tmp)
      r_tmp(0)       =0.0d0
      r_tmp(nx_tmp+1)=10.0d0*r_tmp(nx_tmp)
      read(11) u_tmp(1:nx_tmp,0:4),eps_tmp(:)
      u_tmp(0       ,0:4)=u_tmp(1     ,0:4)
      u_tmp(nx_tmp+1,0:4)=u_tmp(nx_tmp,0:4)
      eps_tmp(0)       =eps_tmp(1)
      eps_tmp(nx_tmp+1)=eps_tmp(nx_tmp)
      close(11)

      j=0
      do i=1,NX
 300     continue
         if (r(i).lt.r_tmp(j+1).and.r(i).gt.r_tmp(j)) then
            xi=(r(i)-r_tmp(j))/(r_tmp(j+1)-r_tmp(j))
            xi1=1.0d0-xi
            u(i,1,1,0:4)=xi*u_tmp(j+1,0:4)+xi1*u_tmp(j,0:4)
c            print *,r(i),u(i,1,1,0)
         else
            j=j+1
            goto 300
         end if
      end do

      do k=1,NZ
         do j=1,NY
            u(:,j,k,:)=u(:,1,1,:)
         end do
      end do

      return

      END SUBROUTINE modell_lesen

c     ============================================================



c     ============================================================

      SUBROUTINE gitter_initialisieren

c     ============================================================

      USE param_var
      USE gitter_var

      IMPLICIT NONE

      real :: alpha,d,r_test,dr_test_dalpha,dalpha,delta_r
      integer :: i,j,k,i_it

#if GEOM==0
c     Kartesische Koordinaten
      do i=-3,NX+3
         r_if(i)=dble(i)/dble(NX)*RADIUS
      end do

      do i=-2,NX+3
         r(i)=0.5d0*(r_if(i)+r_if(i-1))
      end do
      r(-2:0)=-r(-2:0)
      do i=-2,NX+2
         dr_if(i)=abs(r(i+1)-r(i))
      end do
      do i=-2,NX+3
         dr(i)=r_if(i)-r_if(i-1)
      end do

      r1(:)=1.0d0
      dr_if1(:)=1.0d0/dr_if(:)
      dr1(:)=1.0d0/dr(:)

      do i=-3,NY+3
         theta_if(i)=dble(i)/dble(NY)*RADIUS
      end do
      do i=-2,NY+3
         theta(i)=0.5d0*(theta_if(i)+theta_if(i-1))
         cot_theta(i)=0.0d0
      end do
      do i=-2,NY+2
         dtheta_if(i)=theta(i+1)-theta(i)
      end do
      do i=-2,NY+3
         dtheta(i)=theta_if(i)-theta_if(i-1)
      end do

      dtheta_if1(:)=1.0d0/dtheta_if(:)
      dtheta1(:)=1.0d0/dtheta(:)

      do i=-3,NZ+3
         phi_if(i)=dble(i)/dble(NZ)*RADIUS
      end do
      do i=-2,NZ+3
         phi(i)=0.5d0*(phi_if(i)+phi_if(i-1))
      end do
      do i=-2,NZ+2
         dphi_if(i)=phi(i+1)-phi(i)
         dphi_if1(i)=1.0d0/dphi_if(i)
      end do
      do i=-2,NZ+3
         dphi(i)=phi_if(i)-phi_if(i-1)
         dphi1(i)=1.0d0/dphi(i)
      end do

      do k=-1,NZ+2
         do j=-1,NY+2
            do i=-1,NX+2
               dv(i,j,k)=dr(i)*dtheta(j)*dphi(k)
               dv1(i,j,k)=1.0d0/dv(i,j,k)
               r_sin_theta(i,j,k)=1.0d0
               r_sin_theta1(i,j,k)=1.0d0
            end do
         end do
      end do

      do k=-1,NZ+2
         do j=-1,NY+2
            do i=-1,NX+2
               da(i,j,k,1)=dtheta(j)*dphi(k)
               da1(i,j,k,1)=1.0d0/da(i,j,k,1)
               da(i,j,k,2)=dr(i)*dphi(k)
               da1(i,j,k,2)=1.0d0/da(i,j,k,2)
               da(i,j,k,3)=dr(i)*dtheta(j)
               da1(i,j,k,3)=1.0d0/da(i,j,k,2)
            end do
         end do
      end do
#endif /*GEOM==0*/

#if GEOM==1
c     Sphaerische Polarkoordinaten
#ifdef LOG_GITTER
      r_if(0)=0.0d0
      d=R_MIN
      r_if(1)=d
      alpha=(RADIUS/R_MIN)**(1.0d0/dble(NX))

      newton: do i_it=1,50
         r_test=d*(alpha**NX-1.0d0)/(alpha-1.0d0)
         dr_test_dalpha=
     +        (d*NX*alpha**(NX-1)+r_test)/(alpha-1.0d0)
         dalpha=-(r_test-RADIUS)/dr_test_dalpha
         alpha=alpha+dalpha
         if (abs(dalpha/alpha).lt.1.0d-15) exit newton
      end do newton

      if (i_it.ge.50) then
         stop 'gitter_initialisieren: Newton-Iteration!'
      end if

      delta_r=d
      do i=2,NX+3
         delta_r=alpha*delta_r
         r_if(i)=r_if(i-1)+delta_r
      end do
      r_if(-1)=-r_if(1)
      r_if(-2)=-r_if(2)
      r_if(-3)=-r_if(3)
#else
      do i=-3,NX+3
         r_if(i)=dble(i)/dble(NX)*RADIUS
      end do
#endif

      do i=-2,NX+3
         r(i)=sqrt((r_if(i)**2+abs(r_if(i)*r_if(i-1))+
     +        r_if(i-1)**2)/3.0d0)
c         r(i)=0.5d0*(r_if(i)+r_if(i-1))
      end do
      r(-2:0)=-r(-2:0)
      do i=-2,NX+2
         dr_if(i)=abs(r(i+1)-r(i))
      end do
c      dr_if(0)=r(0)+r(1)
      do i=-2,NX+3
         dr(i)=r_if(i)-r_if(i-1)
      end do

      r1(:)=1.0d0/r(:)
      dr_if1(:)=1.0d0/dr_if(:)
      dr1(:)=1.0d0/dr(:)

      do i=-3,NY+3
#ifdef AEQUATOR
         theta_if(i)=dble(i)/dble(NY)*pi*0.5d0
#else
         theta_if(i)=dble(i)/dble(NY)*pi
#endif
      end do
      do i=-2,NY+3
         theta(i)=0.5d0*(theta_if(i)+theta_if(i-1))
         cot_theta(i)=cos(theta(i))/sin(theta(i))
      end do
      do i=-2,NY+2
         dtheta_if(i)=theta(i+1)-theta(i)
      end do
      do i=-2,NY+3
         dtheta(i)=theta_if(i)-theta_if(i-1)
      end do

      dtheta_if1(:)=1.0d0/dtheta_if(:)
      dtheta1(:)=1.0d0/dtheta(:)

      do i=-3,NZ+3
#ifndef OKTANT
         phi_if(i)=dble(i)/dble(NZ)*2.0d0*pi
#else
         phi_if(i)=dble(i)/dble(NZ)*0.5d0*pi
#endif
      end do
      do i=-2,NZ+3
         phi(i)=0.5d0*(phi_if(i)+phi_if(i-1))
      end do
      do i=-2,NZ+2
         dphi_if(i)=phi(i+1)-phi(i)
         dphi_if1(i)=1.0d0/dphi_if(i)
      end do
      do i=-2,NZ+3
         dphi(i)=phi_if(i)-phi_if(i-1)
         dphi1(i)=1.0d0/dphi(i)
      end do

      do k=-1,NZ+2
         do j=-1,NY+2
            do i=-1,NX+2
               dv(i,j,k)=(r_if(i)**3-r_if(i-1)**3)/3.0d0*dphi(k)*
     +              abs(cos(theta_if(j-1))-cos(theta_if(j)))
               dv1(i,j,k)=1.0d0/dv(i,j,k)
               r_sin_theta(i,j,k)=r(i)*sin(theta(j))
               r_sin_theta1(i,j,k)=1.0d0/r_sin_theta(i,j,k)
c               print *,i,j,k,dv(i,j,k)
c               print *,i,j,k,r_sin_theta(i,j,k)
            end do
         end do
      end do
      r_sin_theta1(:,:,:)=1.0d0/r_sin_theta(:,:,:)


      do k=-1,NZ+2
         do j=-1,NY+2
            do i=-1,NX+2
               da(i,j,k,1)=r_if(i)**2*dphi(k)*
     +              abs(cos(theta_if(j-1))-cos(theta_if(j)))
               da1(i,j,k,1)=1.0d0/da(i,j,k,1)
               da(i,j,k,2)=0.5d0*sin(theta_if(j))*
     +              (r_if(i)**2-r_if(i-1)**2)*dphi(k)
               da1(i,j,k,2)=1.0d0/da(i,j,k,2)
               da(i,j,k,3)=0.5d0*(r_if(i)**2-r_if(i-1)**2)*
     +              abs(cos(theta_if(j-1))-cos(theta_if(j)))
               da1(i,j,k,3)=1.0d0/da(i,j,k,2)
            end do
         end do
      end do
#endif /*GEOM==1*/


      print *,'Gitter initialisiert.'
      print *

      return

      END SUBROUTINE gitter_initialisieren

c     ============================================================





