!2345678901234567890
       SUBROUTINE shuey(vp1,vp2,vs1,vs2
     1            ,rho1,rho2, angr, angr_inc,  i_rock,   			!input
     1            rc) 								!output
     
       REAL vp1,vp2,vs1,vs2,rho1,rho2,angr
       REAL rc,A,B,C,VPP,VSS,RHO, angr_inc, ang_m
       INTEGER i_rock
       
       VPP = (vp1 + vp2)/2
       VSS = (vs1 + vs2)/2
       RHO = (rho1 + rho2)/2
       
       VPPP = (vp2 - vp1)/VPP
       VSSS = (vs2 - vs1)/VSS
       RHOO = (rho2 - rho1)/RHO
       
       
       
       A = (0.5)*(VPPP + RHOO)
       
       B = (0.5)*VPPP -(2*(VSS/VPP)**2)*(2*VSSS + RHOO)
       C = (0.5)*VPPP
       
       ang_m = (angr + angr_inc)/2
       rc = A + B*(sin(angr))**2 + C*((sin(angr))**2)*((tan(angr))**2) 

       RETURN
       END      
       
       
       