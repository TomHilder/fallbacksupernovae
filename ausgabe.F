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
      integer :: i, j, k
      real :: eps_tmp, e_ges, e_grav
      real, dimension(0:NX) :: m_enclosed
      character(len=8) :: num

      WRITE(num,'(I8.8)') istep
      OPEN (1, FILE = 'output'//num//'.dat', FORM = 'FORMATTED')

      e_ges=(sum(u(1:NX,1:NY,1:NZ,4)*dv(1:NX,1:NY,1:NZ)))*sym_fac
      e_grav=0.5d0*sum(u(1:NX,1:NY,1:NZ,0)*dv(1:NX,1:NY,1:NZ)*
     +     phi_grav(1:NX,1:NY,1:NZ))*sym_fac

      m_enclosed(0) = m_core ! mass coordinate calcs:
      do i=1,NX
        m_enclosed(i) = m_enclosed(i-1)
        do j = 1, NY
          do k = 1, NZ
            m_enclosed(i) = m_enclosed(i) + u(i,j,k,0)*dv(i,j,k)
     +        *sym_fac
          enddo
        enddo

#if EOS_TYP==1
         eps_tmp=max(1.0d-50,eps(i,NY/2+1,1)-
     +        (k_poly_e*u(i,NY/2,NZ/2,0)**gamma_e*gamma_e_m1_r+
     +        k_poly_b*u(i,NY/2,NZ/2,0)*gamma_b*gamma_b_m1_r))/
     +        (u(i,NY/2+1,1,0)*29979245800.0d0**2)
#endif

      ! Recompute EOS so it may be outputted:

          eps(i,NY/2+1,1)=u(i,NY/2+1,1,4)-
     +        0.5d0*(u(i,NY/2+1,1,1)**2+u(i,NY/2+1,1,2)**2+
     +        u(i,NY/2+1,1,3)**2)/u(i,NY/2+1,1,0)
          eps(i,NY/2+1,1)=max(eps(i,NY/2+1,1),0.0d0)
          call eos(u(i,NY/2+1,1,0),eps(i,NY/2+1,1),p(i,NY/2,1),
     +        c_s(i,NY/2,1),0)

      ! output columns: radius, density, velocity (momentum density/density),
      !   ???, energy (energy density/density), grav potential, grav acc, time,
      !   pressure, sound speed, mass coordinate, E_tot

         write(1,'(12(E16.7))') r(i),u(i,NY/2+1,1,0),
     +        u(i,NY/2+1,1,1)/u(i,NY/2+1,1,0),
     +        u(i,NY/2+1,1,3)/u(i,NY/2+1,1,0),
     +        u(i,1,1,4),phi_grav(i,1,1),a_grv(i,1,1,1),zeit,
     +        p(i,NY/2+1,1),c_s(i,NY/2+1,1),m_enclosed(i),
     +        0*e_ges+e_grav+0*verlust_etot*sym_fac
      end do
c      do i=1,NY
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

      m_ges=sum(u(1:NX,1:NY,1:NZ,0)*dv(1:NX,1:NY,1:NZ))*sym_fac
      e_ges=(sum(u(1:NX,1:NY,1:NZ,4)*dv(1:NX,1:NY,1:NZ)))*sym_fac
      e_grav=0.5d0*sum(u(1:NX,1:NY,1:NZ,0)*dv(1:NX,1:NY,1:NZ)*
     +     phi_grav(1:NX,1:NY,1:NZ))*sym_fac
      rho_c=sum(u(1,1:NY,1:NZ,0)*dv(1,1:NY,1:NZ))/
     +     sum(dv(1,1:NY,1:NZ))


      e_pos = 0.0d0

      do k=1,NZ
         do j=1,NY
            do i=1,NX
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

#if (NDIM>=2) && !(defined(AEQUATOR))
      p_z=0.0d0
      do k=1,NZ
         do j=1,NY
            do i=1,NX
               p_z=p_z+u(i,j,1,1)*cos(theta(j))+
     +              u(i,j,1,2)*sin(theta(j))*dv(i,j,k)
   end dozeit
         end do
      end do
      p_z=p_z*sym_fac

      write(*,'(I16,7(D16.7))') iter,zeit,2*cfl_factor,m_ges,e_ges,
     +     rho_c,p_z/m_ges
#else
      write(*,'(I16,7(D16.7))') iter,zeit,2*cfl_factor,m_ges,0*e_ges+
     +     e_grav+0*verlust_etot*sym_fac,
     +     e_pos,m_core
#endif

      write(97,*) zeit,rho_c
#ifdef CONS_TEST
c
      write(98,*) zeit,e_ges+e_grav+verlust_etot*sym_fac
      write(99,*) zeit,m_ges+verlust_mtot*sym_fac
#endif

c      if (zeit.gt.0.85d1) stop

      return

      END SUBROUTINE protokoll

c     ============================================================



      END MODULE mod_ausgabe

c     ============================================================
