










c     ===============================================================

      MODULE grav_var

c     ===============================================================

      IMPLICIT NONE

      SAVE

      real, dimension(0:800,0:1,0:1) :: phi_grav_if
      real, dimension(-1:800+2,-1:1+2,-1:1+2) :: phi_grav,phi_grav_alt
      real, dimension(0:800+1,0:1+1,0:1+1,3) :: a_grv
      real, dimension(0:800+1,0:1+1,0:1+1,3) :: a_grv_alt=0.0d0

      END MODULE grav_var

c     ===============================================================

c     ===============================================================

      MODULE legendre

c     ===============================================================

      USE param_var

      IMPLICIT NONE

      SAVE

      real, dimension(0:1,0:l_max+1) :: p_leg
      real, dimension(1:1,0:l_max) :: p_leg_int


      real, dimension(0:1,1:l_max+1) :: p1_leg
      real, dimension(1:1,1:l_max) :: p1_leg_int

      END MODULE legendre

c     ===============================================================


c     ===============================================================

      MODULE mod_hilfsgroessen_r

c     ===============================================================

      USE param_var

      IMPLICIT NONE

      SAVE

      real, dimension(0:800,0:l_max+1) :: r_if_hoch_l
      real, dimension(1:800,0:l_max+1) :: r_if_hoch_l_r
      real, dimension(0:800,0:l_max) :: gew_r_in,gew_r_ex

      END MODULE mod_hilfsgroessen_r

c     ===============================================================



