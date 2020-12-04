










c     ============================================================

      PROGRAM TVD

c     ============================================================

      USE advection
      USE poisson
      USE mod_heizen
      USE mod_ausgabe

      IMPLICIT NONE

      integer :: iter,iter2

c      open(1,file='prot.dat',status='replace')
      
      CALL gitter_initialisieren
c      CALL hydro_initialisieren
      CALL modell_lesen
      CALL gravitation_1d(1)

 
      CALL ausgabe(0)
      print *, 'Ausgabe des Anfangsmodells.'
      print *

      write(*,'(7(A16))') 'Iteration','Zeit/s','CFL-Faktor',
     +     'M_tot/g','E_tot/erg','E_+/erg'
      write(*,'(96("-"))')
      CALL protokoll(0)

      do iter2=0,1200
         do iter=1,1000
c            print *,'Iteration',iter+iter2*1000

            CALL hydro_speichern

            CALL advec_x(0)
c            CALL heizen
            CALL gravitation_1d(0)
            CALL hydro_speichern
c            CALL heizen
            CALL advec_x(1)
            CALL gravitation_1d(0)

            CALL courant_factor
   
            if (mod(iter,100).eq.0) CALL protokoll(iter+1000*iter2)
         
c           CALL ausgabe

         end do
         CALL ausgabe(iter-1+1000*iter2)

      end do

c      close(1)


      END PROGRAM TVD

c     ============================================================
