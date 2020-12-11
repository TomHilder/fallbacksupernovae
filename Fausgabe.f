










c     ============================================================

      MODULE mod_ausgabe

c     ============================================================

      CONTAINS

c     ============================================================

      SUBROUTINE ausgabe (istep)

c     ============================================================

      USE hydro_var
      USE gitter_var
      USE grav_var
      USE param_var
      USE mod_eos

      IMPLICIT NONE

      integer, intent(in) :: istep
      integer :: i
      real :: eps_tmp
      character(len=8) :: num

      WRITE(num,'(I8.8)') istep
      OPEN (1, FILE = 'output'//num//'.dat', FORM = 'FORMATTED')

      do i=1,800

      ! Recompute EOS so it may be outputted:

          eps(i,0,0)=u(i,0,0,4)-
     +        0.5d0*(u(i,0,0,1)**2+u(i,0,0,2)**2+
     +        u(i,0,0,3)**2)/u(i,0,0,0)
          eps(i,0,0)=max(eps(i,0,0),0.0d0)
      call eos(u(i,0,0,0),eps(i,0,0),p(i,0,0),c_s(i,0,0),2)

      ! output columns: radius, density, velocity (momentum density/density),
      !   ???, energy (energy density/density), grav potential, grav acc, time,
      !   pressure, sound speed,

         write(1,'(10(E16.7))') r(i),u(i,1/2+1,1,0),
     +        u(i,1/2+1,1,1)/u(i,1/2+1,1,0),
     +        u(i,1/2+1,1,3)/u(i,1/2+1,1,0),
     +        u(i,1,1,4),phi_grav(i,1,1),a_grv(i,1,1,1),zeit,
     +        p(i,0,0),c_s(i,0,0)
      end do
c      do i=1,1
c         write(1,'(7(D16.7))') theta(i),u(4,i,1,2),u(5,i,1,0),
c     +        u(5,i,1,1),u(5,i,1,4),u(5,i,1,3)
c      end do
      write(1,*)

      CLOSE(1)

      return

      END SUBROUTINE ausgabe

c     ============================================================




c     ============================================================

      SUBROUTINE protokoll(iter)

c     ============================================================

      USE hydro_var
      USE grav_var
      USE gitter_var
      USE param_var

      IMPLICIT NONE

      integer, intent(in) :: iter

      real :: m_ges,e_ges,p_z,rho_c,e_grav,e_pos

      integer :: i,j,k

      m_ges=sum(u(1:800,1:1,1:1,0)*dv(1:800,1:1,1:1))*sym_fac
      e_ges=(sum(u(1:800,1:1,1:1,4)*dv(1:800,1:1,1:1)))*sym_fac
      e_grav=0.5d0*sum(u(1:800,1:1,1:1,0)*dv(1:800,1:1,1:1)*
     +     phi_grav(1:800,1:1,1:1))*sym_fac
      rho_c=sum(u(1,1:1,1:1,0)*dv(1,1:1,1:1))/
     +     sum(dv(1,1:1,1:1))


      e_pos = 0.0d0

      do k=1,1
         do j=1,1
            do i=1,800
               IF (u(i,j,k,4)+u(i,j,k,0)*phi_grav(i,j,k)
     &              .GT. 0.0) THEN
                  e_pos = e_pos +
     &                 (u(i,j,k,4)+u(i,j,k,0)*phi_grav(i,j,k)) *
     &                 dv (i,j,k) * sym_fac
               END IF
            end do
         end do
      end do

c      print *,e_ges,e_grav

      write(*,'(I16,7(D16.7))') iter,zeit,2*cfl_factor,m_ges,0*e_ges+
     +     e_grav+0*verlust_etot*sym_fac,
     +     e_pos,m_core

      write(97,*) zeit,rho_c
c
      write(98,*) zeit,e_ges+e_grav+verlust_etot*sym_fac
      write(99,*) zeit,m_ges+verlust_mtot*sym_fac

c      if (zeit.gt.0.85d1) stop

      return

      END SUBROUTINE protokoll

c     ============================================================



      END MODULE mod_ausgabe

c     ============================================================