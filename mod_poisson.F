c     ===============================================================

      MODULE legendre

c     ===============================================================

      USE param_var

      IMPLICIT NONE

      SAVE

      real, dimension(0:NY,0:l_max+1) :: p_leg
      real, dimension(1:NY,0:l_max) :: p_leg_int

      real, dimension(1:NY,1:l_max) :: p_leg_c

      real, dimension(0:NY,1:l_max+1) :: p1_leg
      real, dimension(1:NY,1:l_max) :: p1_leg_int
      real, dimension(1:NY,1:l_max) :: p1_leg_c

      END MODULE legendre

c     ===============================================================


c     ===============================================================

      MODULE mod_hilfsgroessen_r

c     ===============================================================

      USE param_var

      IMPLICIT NONE

      SAVE

      real, dimension(1:NX,0:l_max) :: r_hoch_l
      real, dimension(1:NX,0:l_max) :: r_hoch_lp1_r
      real, dimension(0:NX,0:l_max+1) :: r_if_hoch_l
      real, dimension(1:NX,0:l_max+1) :: r_if_hoch_l_r
      real, dimension(0:NX,0:l_max) :: gew_r_in,gew_r_ex

      END MODULE mod_hilfsgroessen_r

c     ===============================================================



