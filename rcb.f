c This script will compute the reflection coefficient taking into
c account the parametrization done by Simmons and Backus (1996)
!23456789     
       SUBROUTINE rcb(vp1,vp2,vs1,vs2,rho1,rho2,angr,angr_inc,i_rock,
     1    RC) !output

       REAL vp1,vp2,vs1,vs2,rho1,rho2,angr,angr_inc,tmp1,tmp2,RCG,Y
       REAL vp_avg,vs_avg,rho_avg,d_vp,d_rho,Ra,Rp,tmp3,tmp4,Rb,Rrho
       REAL dRb,dRrho,dRsh,RCF,RSH,Rho,k
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

       END