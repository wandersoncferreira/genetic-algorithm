      PROGRAM plot_test
        USE library

      REAL v1(90,2),pi,vp1,vp2,vs1,vs2,rho1,rho2,angr,RCF
      REAL v2(90,2),v3(90,2),dRrho,dRsh,RCF_GARD,RCG,RCF_AKI
      INTEGER i,j,f
  
c Temp values
      pi = 3.141592653
      
      vp1 = (9000)*0.3048
      !vs2 = (6128.2)*0.3048
      !rho2= 2.4783
      
      vp2 = (9080)*0.3048
      vs2 = (5523.4)*0.3048
      rho2= 2.456
      !vp2 = (8889.7)*0.3048
      !vs2 = (5514.6)*0.3048
      !rho2= 2.463

      CALL  f2m(vp1,vs1,rho1)
      !CALL  f2m(vp2,vs2,rho2)
      write(*,*)vp1/0.3048,vs1/0.3048,rho1
      write(*,*)vp2/0.3048,vs2/0.3048,rho2
      
      !vp2 = 3094.6
      !vs2 = 1453.1
      !rho2= 2.1350
      
      !vp2=2979.9
      !vs2=1536.5  !gas with 30% porosity
      !rho2=1.9096
      !vp2=5336.8
      !vs2=2583.3
      !rho2=0.6755

c     TESTE COM MODELO PROPOSTO POR SIMMONS E BACKUS
!2345678 
      angr=0.0
      DO 20 i=1,80
       CALL sims(vp1,vp2,vs1,vs2,rho1,rho2,angr,
     1    RCF,dRrho,dRsh,RCG)
       v1(i,2)=RCF
       v2(i,2)=RCG
       angr=i*pi/180.0
20    ENDDO 
       write(*,*)"Drho: ", dRrho, "dRsh: ", dRsh

      DO j=1,80
        v1(j,1)=j
      ENDDO

      OPEN(unit=21,file='c1.txt')
      DO 11 f=1,80
        WRITE(21,97)(v1(f,i),i=1,2)
11    ENDDO        
97    FORMAT(2f15.5)    

c     TESTE COM MODELO IMPLEMENTADO NA TESE

      angr=0.0
      DO 30 i=1,80
       CALL aki(vp1,vp2,vs1,vs2,rho1,rho2,angr,
     1    RCF_AKI)
       angr=i*pi/180.0
30    ENDDO 

      DO j=1,80
        v2(j,1)=j
      ENDDO
      
      OPEN(unit=71,file='c2.txt')
      DO 41 f=1,80
        WRITE(71,98)(v2(f,i),i=1,2)
41    ENDDO        
98    FORMAT(2f15.5)    
 
c     TESTE COM MODELO PROPOSTO POR SIMMONS SEM BACKUS

      angr=0.0
      DO 40 i=1,80
       CALL gard(vp1,vp2,vs1,vs2,rho1,rho2,angr,
     1    RCF_GARD)
       v3(i,2)=RCF_GARD
       angr=i*pi/180.0
40    ENDDO 

      DO j=1,80
        v3(j,1)=j
      ENDDO
      
      OPEN(unit=91,file='c3.txt')
      DO 111 f=1,80
        WRITE(91,99)(v3(f,i),i=1,2)
111    ENDDO        
99    FORMAT(2f15.5)    
      
      END
