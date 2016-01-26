! This script will compute the reflection coefficient taking into
! account the parametrization done by Simmons and Backus (1996)
!23456789     
       SUBROUTINE f2m(vp,vs,rho)

       REAL vp,vs,rho,vp_tmp

       !Velocities will be inserted in m/s mode.
       vs = (vp - 1360)/1.16

       !convert Vp to ft/s in order to use Gardner's relationship
       vp_tmp = vp/0.3048
       rho=0.2127*(vp_tmp**0.25)

       END