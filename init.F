#define RHO1 1.0d0
#define RHO4 0.125d0
#define DRU1 1.0d0
#define DRU4 0.1d0
#define ZENTRALDICHTE 1.0d10

c     ============================================================

      SUBROUTINE hydro_initialisieren

c     ============================================================

      USE hydro_var
      USE grav_var
      USE gitter_var
      USE param_var
      USE advection, ONLY : randbed

      IMPLICIT NONE

      real, dimension(0:NX) :: m
      real :: dist2
      real :: gamma_tmp,k_poly_tmp

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

c      do j=1,NY
c         do i=1,NX
c            u(i,j,:,3)=10.0d0*r(i)
c         end do
c      end do

c      u(:,:,:,0)=RHO1
c      u(1:NX/2,:,:,0)=RHO4

c      u(:,:,:,4)=1.5d0*DRU1
c      u(1:NX/2,:,:,4)=1.5d0*DRU4
      a_grv(:,:,:,:)=0.0d0
      phi_grav(:,:,:)=0.0d0
#else

      rho_max=ZENTRALDICHTE
      rho_min1=RHO_MIN_FAC*rho_max
      rho_min2=0.5d0*rho_min1

      u(:,:,:,0)=rho_min2
      u(:,:,:,4)=1.5d0*u(:,:,:,0)

      u(1,:,:,0)=ZENTRALDICHTE
      m(0)=0.0d0
      i_aus=NX+1
      CYCLUS: do i=1,NX
         if (u(i,1,1,0).lt.rho_min1) then
            u(i,:,:,0)=u(i-1,1,1,0)!rho_min2
            eps_rand=u(i-1,1,1,4)
c            u(i,1,1,4)=u(i-1,1,1,4)!eps_rand
c            i_aus=i+1
c            exit
         end if
         p  (i,:,:)  =k_poly_tmp*u(i,1,1,0)**gamma_tmp
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
      USE param_var

      IMPLICIT NONE


      real, parameter :: r_dump = 2.0d0 * R_MIN
      real, parameter :: ma_sh = 1.5
      real, parameter :: infall_v_frac = 1.0

      ! set to 0 for not initially collapsing, set to 1 for initially collapsing
      real, parameter :: infall = 1

      ! set to 0 to use infall v as fraction of esc_v, and 1 for fraction of sound speed
      real, parameter :: use_sound = 1

      real,allocatable :: r_tmp(:),m_tmp(:),rho_tmp(:),p_tmp(:)
      real :: xi,xi1,esc_v,sound_v,use_v
      integer :: nx_tmp,dump_i

      integer :: i,j,k


      open(11,file='s32',form='formatted',
     +     status='old')
      read(11,*) nx_tmp
      allocate(r_tmp(1:nx_tmp),m_tmp(1:nx_tmp),
     &     rho_tmp(1:nx_tmp),p_tmp(1:nx_tmp))
      do i = 1, nx_tmp
         read (11,*) r_tmp(i), m_tmp(i), rho_tmp(i), p_tmp(i)
c         print *,i,r_tmp(i),m_tmp(i),rho_tmp(i),p_tmp(i)
      end do
      close(11)

      j=1
      do i=1,NX
 300     continue
         if (r(i).lt.r_tmp(j+1).and.r(i).gt.r_tmp(j)) then
            xi=(r(i)-r_tmp(j))/(r_tmp(j+1)-r_tmp(j))
            xi1=1.0d0-xi
            u(i,1,1,0)=xi*rho_tmp(j+1)+xi1*rho_tmp(j)
            u(i,1,1,4)=xi*p_tmp(j+1)+xi1*p_tmp(j)
c     pressure to energy
            u(i,1,1,4) = u(i,1,1,4) / (gamma_e - 1.0d0)
            if (i .eq. 1) m_core = xi*m_tmp(j+1)+xi1*m_tmp(j)
c            print *,i,j,r(i),r_tmp(j:j+1),m_tmp(j),rho_tmp(j),p_tmp(j)
         else if (r(i).gt.r_tmp(nx_tmp)) then
            u (i,1,1,0) = (r (i-1) / r(i)) ** 3
            u (i,1,1,4) = (r (i-1) / r(i)) ** 3
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

c     find index where radius => r_dump
      do i=1,NX
        if (r(i) .ge. r_dump) then
          dump_i = i
          print*,'dump_i = ',i
          EXIT
        end if
      end do

c     set in-falling material velocity profile:
      if (infall .eq. 1) then
        esc_v = sqrt(6.674E-8 * 1.988E33 * 32 / r_dump)
        sound_v = sqrt(gamma_e*
     &       (gamma_e-1.0)*u(dump_i,1,1,4)/u(dump_i,1,1,0))
        print*,'esc_v = ', esc_v
        print*,'sound_v = ', sound_v
        print*,'sound_v/esc_v= ', sound_v/esc_v
        use_v = esc_v
        if (use_sound .eq. 1) then
          use_v = sound_v
        end if
        !print*,'use_v = ',use_v
        do k=1,NZ
          do j=1,NY
            do i=1,NX
              if (r(i) .gt. r_dump) then
                u(i,j,k,1)=-1 * u(i,j,k,0) * infall_v_frac * use_v *
     &             (r(i) / r_dump)**(-2.) ! move this line to before if statement for superposition profile
                u(i,j,k,4)=u(i,j,k,4) + 0.5 * u(i,j,k,1)**2 / u(i,j,k,0)
              end if
            end do
          end do
        end do
      end if

c     dump extra energy to explode model:
      do k=1,NZ
         do j=1,NY
            do i=1,NX
               if (r(i) .le. r_dump) then
                  u(i,j,k,1)=u(i,j,k,1) +
     &                 u(i,j,k,0) * (r(i)/r_dump) *
     &                 (ma_sh * sound_v +
     &                 u(dump_i,1,1,1)/u(dump_i,1,1,0))
c     &                 sqrt(gamma_e*(gamma_e-1.0)*u(i,j,k,4)/u(i,j,k,0))
                  u(i,j,k,4)=u(i,j,k,4) * (1.0+ma_sh) +
     &                 0.5d0 * u(i,j,k,1)**2 / u(i,j,k,0)
               end if
            end do
         end do
      end do

c     set in-falling material velocity profile
!      if (infall .eq. 1) then
!        esc_v = sqrt(6.674E-8 * 1.988E33 * 32 / r_dump)
!        sound_v = sqrt(5./3. * p_tmp(dump_i) / u(dump_i,1,1,0))
!        print*,'esc_v = ', esc_v
!        print*,'sound_v = ', sound_v
!        print*,'esc_v/sound_v= ', esc_v/sound_v
!        use_v = esc_v
!        if (use_sound .eq. 1) then
!          use_v = sound_v
!        end if
!        !print*,'use_v = ',use_v
!        do k=1,NZ
!          do j=1,NY
!            do i=1,NX
!              if (r(i) .gt. r_dump) then
!              u(i,j,k,1)=-1 * u(i,j,k,0) * infall_v_frac * use_v *
!     &             (r(i) / r_dump)**(-2.)
!              u(i,j,k,4)=u(i,j,k,4) + 0.5 * u(i,j,k,1)**2 / u(i,j,k,0)
!              end if
!            end do
!          end do
!        end do
!      end if

      print*,'Core mass (M_sun):', m_core / msun
c      stop

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
      alpha=(RADIUS/R_MIN)**(1.0d0/dble(NX+1))
      r_if(0)=R_MIN
      r_if(1)=R_MIN*alpha

      do i=2,NX+3
         r_if(i)=r_if(i-1)*alpha
      end do
      r_if(-1)=r_if(0)/alpha
      r_if(-2)=r_if(0)/alpha**2
      r_if(-3)=r_if(0)/alpha**3
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
c      r(-2:0)=-r(-2:0)
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
