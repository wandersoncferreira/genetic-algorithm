      PROGRAM plot_test
        USE library

      REAL v1(90,2),angr_inc,pi,vp1,vp2,vs1,vs2,rho1,rho2,angr,RCF
      REAL v2(90,2),v3(90,2),dRrho,dRsh
      INTEGER i,j,f,i_rock
  
c Temp values
      angr_inc=0.0
      pi = 3.14567
      i_rock=0
      
      vp1 = (8200)*0.3048
      !vs1 = (3833)*0.3048
      !rho1= 2.30
      
      vp2 = (9800)*0.3048
      !vs2 = (5071)*0.3048
      !rho2= 2.27

      CALL  f2m(vp1,vs1,rho1)
      CALL  f2m(vp2,vs2,rho2)
      !write(*,*)vp1,vs1,rho1
      !write(*,*)vp2,vs2,rho2
      
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
       CALL rc(vp1,vp2,vs1,vs2,rho1,rho2,angr,angr_inc,i_rock,
     1    RCF,dRrho,dRsh)
       v1(i,2)=RCF
       angr=i*pi/180.0
       write(*,*)dRrho,dRsh
20    ENDDO 

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
       CALL aki(vp1,vp2,vs1,vs2,rho1,rho2,angr,angr_inc,i_rock,
     1    RCF)
       v2(i,2)=RCF
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
       CALL rcb(vp1,vp2,vs1,vs2,rho1,rho2,angr,angr_inc,i_rock,
     1    RCF)
       v3(i,2)=RCF
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
