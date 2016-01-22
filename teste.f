      PROGRAM plot_test

      REAL v1(90,2),angr_inc,pi,vp1,vp2,vs1,vs2,rho1,rho2,angr,RCF
      INTEGER i,j,f,i_rock
  
c Temp values
      angr_inc=0.0
      pi = 3.14567
      i_rock=0
c  Rock physics parameters in SI units (m/s, g/cm3)
      vp1 = 4000
      vp2 = 5000
      vs1 = 2300
      vs2 = 2600
      rho1= 1.6
      rho2= 1.9
!2345678      
      angr=0
      DO 20 i=1,50
       CALL rc(vp1,vp2,vs1,vs2,rho1,rho2,angr,angr_inc,i_rock,
     1    RCF)
       v1(i,2)=RCF
       angr=i*pi/180
20    ENDDO 

      DO j=1,50
        v1(j,1)=j
      ENDDO
      
      OPEN(unit=21,file='plot_test.txt')
      DO 11 f=1,50
        WRITE(21,97)(v1(f,i),i=1,2)
11    ENDDO        
97    FORMAT(2f15.5)    
      END
