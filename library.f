c This script will compute the reflection coefficient taking into
c account the parametrization done by Simmons and Backus (1996)
!23456789
      MODULE library
      CONTAINS

       SUBROUTINE f2m(vp,vs,rho)
       REAL vp,vs,rho,vp_tmp
       !Velocities will be inserted in m/s mode.
       vs = (vp - 1360)/1.16

       !convert Vp to ft/s in order to use Gardner's relationship
       vp_tmp = vp/0.3048
       rho=0.2127*(vp_tmp**0.25)

       END SUBROUTINE

!! NEW SUBROUTINE RC - SIMMONS AND BACKUS FINAL       
       SUBROUTINE rc(vp1,vp2,vs1,vs2,rho1,rho2,angr,angr_inc,i_rock,
     1    RCF,dRrho,dRsh) !output

       REAL vp1,vp2,vs1,vs2,rho1,rho2,angr,angr_inc,tmp1,tmp2,RCG,Y
       REAL vp_avg,vs_avg,rho_avg,d_vp,d_rho,Ra,Rp,tmp3,tmp4,Rb,Rrho
       REAL dRb,dRrho,dRsh,RCF,k,Ro
       INTEGER i_rock

  
c Aux Variables to compute the eq.18
       vp_avg  = (vp1 + vp2)/2.0
       vs_avg  = (vs1 + vs2)/2.0
       Y       = 4*(vs_avg/vp_avg)**2.0
       k       = 0.86*(vp_avg/vs_avg)

c Variables to compute the normal-incidence compressional-wave
coefficient (Ro) based on the equations (6)
       rho_avg = (rho1 + rho2)/2.0
       d_vp    = (vp2 - vp1)
       d_rho   = (rho2 - rho1)
       Ra      = (1/2.0)*(d_vp/vp_avg)
       Rp      = (1/2.0)*(d_rho/rho_avg)
       Ro      = Ra + Rp

c RC equation without any assumptions (Gardner nor Castagna) eq.18
       tmp1    = (0.8 - 0.2*Y - 1.6*k*Y)*(sin(angr)**2)
       tmp2    = 0.8*((sin(angr)**2))*(tan(angr)**2)
       RCG     = (1 + tmp1 + tmp2)*Ro

c Computing dRsh - prediction error for shear wave
       Rb      = k*Ra
       Rrho    = 0.25*Ra 
       dRb     = Rb   - 0.8*k*Ro
       dRrho   = Rrho - 0.2*Ro 
       dRsh    = dRb  + dRrho


c RC final equation after Gardner and Castagna relationship and the
c prediction errors for hydrocarbon zones. eq. 23
       tmp3    = 2*Y*(sin(angr)**2)*dRsh
       tmp4    = (Y -1 - (tan(angr)**2))*(sin(angr)**2)*dRrho
       RCF      = RCG  - tmp3 + tmp4
       END SUBROUTINE


! NEW MODULE RCB       
       SUBROUTINE rcb(vp1,vp2,vs1,vs2,rho1,rho2,angr,angr_inc,i_rock,
     1    RC) !output

       REAL vp1,vp2,vs1,vs2,rho1,rho2,angr,angr_inc,tmp1,tmp2,RCG,Y
       REAL vp_avg,vs_avg,rho_avg,d_vp,d_rho,Ra,Rp,tmp3,tmp4,Rb,Rrho
       REAL dRb,dRrho,dRsh,RCF,RSH,Rho,k,Ro,d_vs,RC
       INTEGER i_rock

  
c Aux Variables to compute the eq.18
       vp_avg  = (vp1 + vp2)/2.0
       vs_avg  = (vs1 + vs2)/2.0
       Y       = 4*(vs_avg/vp_avg)**2.0
       k       = (vp_avg/vs_avg)*0.86
c Variables to compute the normal-incidence compressional-wave
coefficient (Ro) based on the equations (6)
       rho_avg = (rho1 + rho2)/2.0
       d_vp    = (vp2 - vp1)
       d_vs    = (vs2 - vs1)
       d_rho   = (rho2 - rho1)
       Ra      = (1/2.0)*(d_vp/vp_avg)
       Rb      = (1/2.0)*(d_vs/vs_avg)
       Rp      = (1/2.0)*(d_rho/rho_avg)
       Ro      = Ra + Rp
        
       tmp1 = (1+(sin(angr)**2)+(sin(angr)**2)*(tan(angr)**2))
       tmp2 = 2*Y*(sin(angr)**2)
       tmp3 = (Y-1-(tan(angr)**2))
       
       RCF =(Ra + Rp)
       RSH = (Rb + Rp)

       RC=tmp1*RCF - tmp2*RSH +tmp3*(sin(angr)**2)*Rp

       END SUBROUTINE


!! NEW SUBROUTINE AKI AND RICHARDS       
       SUBROUTINE aki(vp1,vp2,vs1,vs2,rho1,rho2,angr,angr_inc,i_rock,
     1    RC) !output

       REAL vp1,vp2,vs1,vs2,rho1,rho2,angr,angr_inc,tmp1,tmp2,Y
       REAL vp_avg,vs_avg,Ro,rho_avg,d_vp,d_rho,Ra,Rp,tmp3,Rb
       REAL k,RC,d_vs
       INTEGER i_rock

  
c Aux Variables to compute the eq.18
       vp_avg  = (vp1 + vp2)/2.0
       vs_avg  = (vs1 + vs2)/2.0
       Y       = 4*(vs_avg/vp_avg)**2.0
       k       = (vp_avg/vs_avg)*0.86
c Variables to compute the normal-incidence compressional-wave
coefficient (Ro) based on the equations (6)
       rho_avg = (rho1 + rho2)/2.0
       d_vp    = (vp2 - vp1)
       d_vs    = (vs2 - vs1)
       d_rho   = (rho2 - rho1)
       Ra      = (1/2.0)*(d_vp/vp_avg)
       Rb      = (1/2.0)*(d_vs/vs_avg)
       Rp      = (1/2.0)*(d_rho/rho_avg)
       Ro      = Ra + Rp
        
       tmp1 = Ra + Rp
       tmp2 = Ra - Y*Rp -2*Y*Rb
       tmp3 = Ra

       RC=tmp1+tmp2*(sin(angr)**2)+tmp3*(sin(angr)**2)*(tan(angr)**2)

       END SUBROUTINE


       END MODULE library