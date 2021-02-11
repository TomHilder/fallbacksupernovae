











c     ============================================================

      SUBROUTINE hydro_initialisieren

c     ============================================================

      USE hydro_var
      USE grav_var
      USE gitter_var
      USE param_var
      USE advection, ONLY : randbed

      IMPLICIT NONE

      real, dimension(0:800) :: m
      real :: dist2
      real :: gamma_tmp,k_poly_tmp

      integer :: i,j,k,i_aus

      u=0.0d0
      eps_rand=0.0d0

      gamma_tmp=5.0d0/3.0d0
      k_poly_tmp=k_poly_e



      rho_max=1.0d10
      rho_min1=1.0d-14*rho_max
      rho_min2=0.5d0*rho_min1

      u(:,:,:,0)=rho_min2
      u(:,:,:,4)=1.5d0*u(:,:,:,0)

      u(1,:,:,0)=1.0d10
      m(0)=0.0d0
      i_aus=800+1
      CYCLUS: do i=1,800
         if (u(i,1,1,0).lt.rho_min1) then
            u(i,:,:,0)=u(i-1,1,1,0)!rho_min2
            eps_rand=u(i-1,1,1,4)
c            u(i,1,1,4)=u(i-1,1,1,4)!eps_rand
c            i_aus=i+1
c            exit
         end if
         p  (i,:,:)  =k_poly_tmp*u(i,1,1,0)**gamma_tmp
         u  (i,:,:,4)=p(i,1,1)/(gamma_tmp-1.0d0)

         m(i)=m(i-1)+sum(u(i,1:1,1:1,0)*dv(i,1:1,1:1))*sym_fac
         u(i+1,:,:,0)=u(i,1,1,0)-
     +        dr_if(i)*g_grav*m(i)/r_if(i)**2/gamma_tmp/k_poly_tmp*
     +        u(i,1,1,0)**(2.0d0-gamma_tmp)

      end do CYCLUS


      if (i_aus.le.800) then
         u(i_aus:800,:,:,0)=u(i_aus-1,1,1,0)
         u(i_aus:800,:,:,4)=u(i_aus-1,1,1,4)
         stop 'hydro_initialisieren: Dichte zu gering!'
      else
         eps_rand=u(800,1,1,4)
      end if
      rho_min1=u(800,1,1,0)
      rho_min2=rho_min1




c     Rotationsgeschwindigkeit:
c      do k=1,1
c         do j=1,1
c            do i=1,i_aus
c               u(i,j,k,3)=r(i)*sin(theta(j))*1.0d-2*u(i,j,k,0)
c               u(i,j,k,4)=u(i,j,k,4)+0.5d0*
c     +              u(i,j,k,3)**2/u(i,j,k,0)
c            end do
c         end do
c      end do


      call randbed(1)
      do k=-1,1+2
         do j=-1,1+2
            do i=-1,800+2
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
      write(11) 800
      write(11) r(1:800)
      write(11) u(1:800,1,1,:),eps(1:800,1,1)
      close(11)

c      u(1:800,1,1,1)=u(1:800,1,1,0)*r(1:800)

c      do j=1,1
c         u(:,j,:,0)=(0.5d0+0.5d0*sin(theta(j)))**8*u(:,j,:,0)
c      end do

c      stop


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


      real, parameter :: r_dump = 2.0d0 * 1.5d8
      real, parameter :: ma_sh = 1.50
      real, parameter :: infall_v_frac = 1.0

      ! set to 0 for not initially collapsing, set to 1 for initially collapsing, set to 2 for split infall profile
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
      do i=1,800
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

      do k=1,1
         do j=1,1
            u(:,j,k,:)=u(:,1,1,:)
         end do
      end do

c     find index where radius => r_dump
      do i=1,800
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
        do k=1,1
          do j=1,1
            do i=1,800
              u(i,j,k,1)=-1 * u(i,j,k,0) * infall_v_frac * use_v *
     &             (r(i) / r_dump)**(-2.)
              if (r(i) .gt. r_dump) then
                u(i,j,k,4)=u(i,j,k,4) + 0.5 * u(i,j,k,1)**2 / u(i,j,k,0)
              end if
            end do
          end do
        end do
      end if

c     set in-falling material split velocity profile
      if (infall .eq. 2) then
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
        do k=1,1
          do j=1,1
            do i=1,800
              u(i,j,k,1)=-1 * u(i,j,k,0) * infall_v_frac * use_v *
     &             (r(i) / r_dump)**(-2.)
              if (r(i) .gt. r_dump) then
                u(i,j,k,4)=u(i,j,k,4) + 0.5 * u(i,j,k,1)**2 / u(i,j,k,0)
              end if
            end do
          end do
        end do
      end if

c     dump extra energy to explode model:
      do k=1,1
         do j=1,1
            do i=1,800
               if (r(i) .le. r_dump) then
                  u(i,j,k,1)=u(i,j,k,1) +
     &                 u(i,j,k,0) * ma_sh * (r(i)/r_dump) *
     &                 sqrt(gamma_e*(gamma_e-1.0)*u(i,j,k,4)/u(i,j,k,0))
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
!        do k=1,1
!          do j=1,1
!            do i=1,800
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


c     Sphaerische Polarkoordinaten
      alpha=(1d13/1.5d8)**(1.0d0/dble(800+1))
      r_if(0)=1.5d8
      r_if(1)=1.5d8*alpha

      do i=2,800+3
         r_if(i)=r_if(i-1)*alpha
      end do
      r_if(-1)=r_if(0)/alpha
      r_if(-2)=r_if(0)/alpha**2
      r_if(-3)=r_if(0)/alpha**3

      do i=-2,800+3
         r(i)=sqrt((r_if(i)**2+abs(r_if(i)*r_if(i-1))+
     +        r_if(i-1)**2)/3.0d0)
c         r(i)=0.5d0*(r_if(i)+r_if(i-1))
      end do
c      r(-2:0)=-r(-2:0)
      do i=-2,800+2
         dr_if(i)=abs(r(i+1)-r(i))
      end do
c      dr_if(0)=r(0)+r(1)
      do i=-2,800+3
         dr(i)=r_if(i)-r_if(i-1)
      end do

      r1(:)=1.0d0/r(:)
      dr_if1(:)=1.0d0/dr_if(:)
      dr1(:)=1.0d0/dr(:)

      do i=-3,1+3
         theta_if(i)=dble(i)/dble(1)*pi
      end do
      do i=-2,1+3
         theta(i)=0.5d0*(theta_if(i)+theta_if(i-1))
         cot_theta(i)=cos(theta(i))/sin(theta(i))
      end do
      do i=-2,1+2
         dtheta_if(i)=theta(i+1)-theta(i)
      end do
      do i=-2,1+3
         dtheta(i)=theta_if(i)-theta_if(i-1)
      end do

      dtheta_if1(:)=1.0d0/dtheta_if(:)
      dtheta1(:)=1.0d0/dtheta(:)

      do i=-3,1+3
         phi_if(i)=dble(i)/dble(1)*2.0d0*pi
      end do
      do i=-2,1+3
         phi(i)=0.5d0*(phi_if(i)+phi_if(i-1))
      end do
      do i=-2,1+2
         dphi_if(i)=phi(i+1)-phi(i)
         dphi_if1(i)=1.0d0/dphi_if(i)
      end do
      do i=-2,1+3
         dphi(i)=phi_if(i)-phi_if(i-1)
         dphi1(i)=1.0d0/dphi(i)
      end do

      do k=-1,1+2
         do j=-1,1+2
            do i=-1,800+2
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


      do k=-1,1+2
         do j=-1,1+2
            do i=-1,800+2
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


      print *,'Gitter initialisiert.'
      print *

      return

      END SUBROUTINE gitter_initialisieren

c     ============================================================
