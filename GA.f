       PROGRAM GA
       
       USE ga_methods
       
       
       PARAMETER (maxray=6500,no_traces=12,Nsel=10,npar=4,nbits=5) 
       INTEGER, PARAMETER :: popsize=5, nl=4  
       INTEGER  maxit,iga, Nt, keep, nparM,zz,nPPD(3*npar,2**(nbits)+1)
       INTEGER tt,kk,nindex,mo_p(Nsel/2),fa_p(Nsel/2),winner(1),k
       INTEGER seed,vv,nopt1,nopt2,indx(popsize),size_t,count
       INTEGER individuos(Nsel),u,choices(Nsel),N,indx_t(4),b,f
       
       REAL*4 teste_u(4,npar*nbits),thickness(nl),noise(6500,1)
       REAL*4 hi(3*npar),lo(3*npar),hi_s,lo_s,hi_p,lo_p,hi_rho,lo_rho       
       REAL*4 umax,umin,uavg,fitness(popsize),rr,fit,Sti(Nsel/2)
       REAL*4 PPD(3*npar,2**(nbits)+1),grc,grcc(popsize),Stp(Nsel/2)
       REAL*4 of1_fitness(Nsel/2),of2_fitness(Nsel/2),sumPPD,t1,t2
       REAL*4 model_final(maxray,no_traces),VP(1,nl),VS(1,nl),RHO(1,nl)
       REAL*4 synthetic(maxray,no_traces),teste_pu(4),reflex(popsize)

       REAL*4 mutrate, mincost,selection,dt,teste(16),M,Pc,Pm,Pu
       REAL*4, DIMENSION(:,:), ALLOCATABLE  :: pop_p,pop_rho
       REAL*4, DIMENSION(:,:), ALLOCATABLE  :: offsp1,offsp2
       REAL*4, DIMENSION(:,:), ALLOCATABLE  :: quant,pop,par
       REAL*4, DIMENSION(:,:), ALLOCATABLE  :: par_p
       REAL*4 par_s(popsize,npar),par_rho(popsize,npar)
       REAL*4 mo(Nsel/2,npar*nbits),fa(Nsel/2,npar*nbits)
       
       REAL*4, DIMENSION(:,:), ALLOCATABLE :: o1parp,o2parp,o1pars
       REAL*4, DIMENSION(:,:), ALLOCATABLE :: o2pars,o1parho,o2parho  
       REAL*4, DIMENSION(:,:), ALLOCATABLE :: offsp1_p,offsp1_s
       REAL*4, DIMENSION(:,:), ALLOCATABLE :: offsp1_rho,offsp2_p
       REAL*4, DIMENSION(:,:), ALLOCATABLE :: offsp2_s,offsp2_rho
       
       REAL*4, DIMENSION(:), ALLOCATABLE  :: vp_it
       REAL*4, DIMENSION(:), ALLOCATABLE  :: vs_it
       REAL*4, DIMENSION(:), ALLOCATABLE  :: rho_it
       REAL*4, DIMENSION(:), ALLOCATABLE  :: depth_it
       REAL*4, DIMENSION(:), ALLOCATABLE  :: time_series

 
c      The GA parameters is defined in here     
       CALL CPU_TIME(t1)
       maxit = 200
       Pm = 0.01         ! mutation rate
       Pc = 0.9           ! crossover rate 
       Pu = 0.4       ! update rate
       selection = 0.5   ! selection rate
       Nt = nbits*npar
       keep = FLOOR(selection*popsize)
c       iga = 0     !iteration number
       
       dt = 0.004  ! time sample rate
       N = 1
       size_t = 4
       count = 0
       sumPPD = 0
       

c     ----- CALLING ZEROS FOR SEVERAL ARRAYS -----
       CALL zero(synthetic,maxray)
       CALL zero(model_final,maxray)
       CALL zero(fitness,popsize)
       CALL zero(of1_fitness,Nsel/2)
       CALL zero(of2_fitness,Nsel/2)
       CALL zero(choices,Nsel)
       CALL zero(individuos,Nsel)
       CALL zero(indx,popsize)

c      -----------------------------------------------------------------
c      Creating the initial population P-WAVE
        hi_p = 1500.0     ! bounds of P-wave velocity
        lo_p = 4000.0        
        ALLOCATE(pop_p(popsize,Nt))
        ALLOCATE(offsp1_p(Nsel/2,Nt))
        ALLOCATE(offsp2_p(Nsel/2,Nt))
        CALL init_random_seed()
        CALL RANDOM_NUMBER(pop_p)   
        pop_p = ANINT(pop_p) 
        ALLOCATE(par_p(popsize,npar))
        ALLOCATE(o1parp(Nsel/2,npar))
        ALLOCATE(o2parp(Nsel/2,npar))
        par_p = gadecode(pop_p,hi_p,lo_p,nbits,popsize,npar)
       
       

c      Creating the initial population S-WAVE 
    
        ALLOCATE(offsp1_s(Nsel/2,Nt))
        ALLOCATE(offsp2_s(Nsel/2,Nt))
        ALLOCATE(o1pars(Nsel/2,npar))
        ALLOCATE(o2pars(Nsel/2,npar))



                    
             
c       Creating the initial population DENSITY 
     

        ALLOCATE(offsp1_rho(Nsel/2,Nt))
        ALLOCATE(offsp2_rho(Nsel/2,Nt))
        ALLOCATE(o1parho(Nsel/2,npar))
        ALLOCATE(o2parho(Nsel/2,npar))



c      Creating the initial population of thicknesses
        thickness = (/1000,50,50,500/) 
        
        
        ALLOCATE(pop(popsize,Nt))
        ALLOCATE(par(popsize,3*npar))
        ALLOCATE(offsp1(Nsel/2,Nt))
        ALLOCATE(offsp2(Nsel/2,Nt))
       
        pop=pop_p
       
        DO i=1,popsize
            par_rho(i,:) = 0.23*((par_p(i,:))**0.25)
            par_s(i,:) =(par_p(i,:) - 1360)/1.16
        ENDDO
       


        DO f=1,3*npar
            DO l=1,2**nbits+1
                nPPD(f,l)=0
                PPD(f,l)=0
            ENDDO
        ENDDO

        hi_s   = (hi_p - 1360)/1.16
        hi_rho = 0.23*((hi_p**0.25))
        lo_s   = (lo_p - 1360)/1.16
        lo_rho = 0.23*((lo_p**0.25))
		

		
       hi = (/ hi_p,hi_p,hi_p,hi_p,hi_s,hi_s,hi_s,hi_s,hi_rho,hi_rho,
     1       hi_rho, hi_rho /)
       lo = (/ lo_p,lo_p,lo_p,lo_p,lo_s,lo_s,lo_s,lo_s,lo_rho,lo_rho,
     1       lo_rho,lo_rho /)
c        lo = (/ lo_p,lo_p,lo_s,lo_s,lo_rho,lo_rho /)
c        hi = (/ hi_p,hi_p,hi_s,hi_s,hi_rho,hi_rho /)



Cc      ---------------------------------------------------------------
Cc      ---------------------------------------------------------------
Cc      ---------------------------------------------------------------
Cc      ---------------  TRUE MODEL I WANT TO FIND OUT   --------------
Cc      ---------------------------------------------------------------
Cc      ---------------------------------------------------------------
Cc      ---------------------------------------------------------------
      
        VP(1,:)  = (/ 2002.357, 2699.043, 2394.03, 2658.021 /)
        VS(1,:)  = (/ 778.59, 1327.48, 1106.25, 1358.98  /)
        RHO(1,:) = (/ 1.1408, 1.2932, 1.1466, 1.44164  /)
      
        CALL d2t(VP(1,:),VS(1,:),RHO(1,:),thickness,dt,nl,
     1       vp_it, vs_it, rho_it, depth_it,time_series,nopt1)
     

        CALL synt(vp_it,vs_it,rho_it,depth_it,
     1       nopt1,model_final,grc)


c        CALL init_random_seed()
c        CALL RANDOM_NUMBER(noise)
c        noise = noise/5 
c        model_final = model_final + noise
          

        OPEN (unit = 75, file = 'cmp_invertido4.txt') 
        DO vv=1,nopt1
        WRITE(75,995)(model_final(vv,j), j=1,no_traces)
        ENDDO
995     FORMAT(12F14.7)

        write(*,*) " fitness_max ","     fitness_min ","     fitness_avg
     1  ", "         iteration", "   time_elapsed"

        DEALLOCATE(vp_it)
        DEALLOCATE(vs_it)
        DEALLOCATE(rho_it)
        DEALLOCATE(depth_it)
        DEALLOCATE(time_series)
C       
Cc      ---------------------------------------------------------------
Cc      ---------------------------------------------------------------
Cc      ---------------------------------------------------------------
Cc      -----------------    END OF TRUE MODELING    ------------------
Cc      ---------------------------------------------------------------
Cc      ---------------------------------------------------------------
Cc      ---------------------------------------------------------------
C
C!23456
        DO 20 i=1,popsize !compute one synthetic seismogram for each elem.
Cc  
            CALL d2t(par_p(i,:),par_s(i,:),par_rho(i,:),thickness,dt,nl,
     1          vp_it, vs_it, rho_it, depth_it,time_series,nopt2)
       
            CALL synt(vp_it,vs_it,rho_it,depth_it,
     1          nopt2,synthetic,grc)
      
            count = count + 1

            DEALLOCATE(vp_it)
            DEALLOCATE(vs_it)
            DEALLOCATE(rho_it)
            DEALLOCATE(depth_it)
            DEALLOCATE(time_series)
       
     
            IF(nopt1.GT.nopt2)THEN
            CALL xcross2d(synthetic(1:nopt1,:),model_final(1:nopt1,:),
     1	    nopt1,no_traces,rr)
            ELSE
            
            CALL xcross2d(synthetic(1:nopt2,:),model_final(1:nopt2,:),
     1	    nopt2,no_traces,rr)
            ENDIF
             
            IF(rr.lt.0) THEN
                fitness(i) = 0
            ELSE
                fitness(i) = rr
            ENDIF
          

20      ENDDO


        umax = MAXVAL(fitness)
        umin = MINVAL(fitness)
        uavg = SUM(fitness)/popsize
      
    
        CALL scalepop(fitness,popsize,umax,uavg,umin) !scale the fitness landscape       
        CALL sort(popsize,fitness,indx)               !sorting the fitness values
        CALL reverse(indx,popsize)                   !reversing sorting higher to lower
        fitness=fitness(indx)
        pop=pop(indx,:)
       
       

Cc ---------------------------------------------------------------------------------------
Cc ---------------------------------------------------------------------------------------
Cc -------------------- Now we start the loop to find the best population ----------------
Cc ---------------------------------------------------------------------------------------
Cc ---------------------------------------------------------------------------------------

C
      DO 9999 u =1,maxit
C      
        IF(u.GT.maxit/2)then
            Pu = 0.5
        ENDIF
C
        CALL sswr(fitness(1:keep),Nsel,choices)
        CALL shuffle(choices,Nsel)	   
        individuos = choices(1:Nsel)
        individuos(1) = 1
c        write(*,*) individuos(1:5)
		
        DO 30 i=1,Nsel/2
            Stp(i) = 1 + (i-1)*2
30      ENDDO
        
        DO 40 i=1,Nsel/2
            Sti(i) = 2 + (i-1)*2
40	ENDDO

        mo_p = individuos(int(Stp))
        fa_p = individuos(int(Sti))
        
        mo   = pop(int(mo_p),:)
        fa   = pop(int(fa_p),:) 

CCc ---------------------------------------------------------------------------------------
CCc ---------------------------------------------------------------------------------------
CCc --------------------------- 		Crossover time 		----------------------------------
CCc ---------------------------------------------------------------------------------------
CCc ------------------------- Verifying update probability --------------------------------
CCc ---------------------------------------------------------------------------------------
CCc ---------------------------------------------------------------------------------------
        DO 50 i=1,(Nsel/2)
            CALL crossover(mo(i,:),fa(i,:),offsp1(i,:),
	1        offsp2(i,:),Nt,Pc,Pm)

            offsp1_p(i,:)   = offsp1(i,1:Nt)				
            offsp2_p(i,:)   = offsp2(i,1:Nt)

            o1parp=gadecode(offsp1_p(i,:),hi_p,lo_p,nbits,1,npar)
            o1pars  = (o1parp - 1360)/1.16
            o1parho = 0.23*((o1parp**0.25))

            o2parp=gadecode(offsp2_p(i,:),hi_p,lo_p,nbits,1,npar)
            o2pars  = (o2parp - 1360)/1.16
            o2parho = 0.23*((o2parp**0.25))   

CCc ------------------------ Calculo feito para o offspring 1 -----------------------------

        CALL d2t(o1parp,o1pars,o1parho,thickness,dt,nl,
     1       vp_it, vs_it, rho_it, depth_it,time_series,nopt2)
     
   		   
       
        CALL synt(vp_it,vs_it,rho_it,depth_it,
     1       nopt2,synthetic,grc)
     


        DEALLOCATE(vp_it)
        DEALLOCATE(vs_it)
        DEALLOCATE(rho_it)
        DEALLOCATE(depth_it)
        DEALLOCATE(time_series)
		   
        IF(nopt1.GT.nopt2)THEN
            CALL xcross2d(synthetic(1:nopt1,:),model_final(1:nopt1,:),
     1	    nopt1,no_traces,rr)
        ELSE
            CALL xcross2d(synthetic(1:nopt2,:),model_final(1:nopt2,:),
     1	   nopt2,no_traces,rr)
        ENDIF				
        
        
        IF(rr.lt.0) THEN
            of1_fitness(i) = 0
        ELSE
            of1_fitness(i) = rr
        ENDIF
       	   
        count = count + 1
       	    
CCc ----------------------- Calculo feito para o offspring 2 ------------------------------     

        CALL d2t(o2parp,o2pars,o2parho,thickness,dt,nl,
     1       vp_it, vs_it, rho_it, depth_it,time_series,nopt2)
       
        CALL synt(vp_it,vs_it,rho_it,depth_it,
     1       nopt2,synthetic,grc)


        DEALLOCATE(vp_it)
        DEALLOCATE(vs_it)
        DEALLOCATE(rho_it)
        DEALLOCATE(depth_it)
        DEALLOCATE(time_series)
		   
        
        IF(nopt1.GT.nopt2)THEN
            CALL xcross2d(synthetic(1:nopt1,:),model_final(1:nopt1,:),
     1	         nopt1,no_traces,rr)
        ELSE
            CALL xcross2d(synthetic(1:nopt2,:),model_final(1:nopt2,:),
     1	         nopt2,no_traces,rr)
        ENDIF				
        
        
        IF(rr.lt.0) THEN
            of2_fitness(i) = 0
        ELSE
            of2_fitness(i) = rr
        ENDIF
       	   
        count = count + 1
       	   
        teste_pu(1) = of1_fitness(i)
        teste_pu(2) = of2_fitness(i)       	    
        teste_pu(3) = fitness(int(fa_p(i)))
        teste_pu(4) = fitness(int(mo_p(i)))
      	          	   
        CALL sort(size_t,teste_pu,indx_t)        !sorting the fitness values
        CALL reverse(indx_t,size_t)
            
        teste_u(1,:) = offsp1(i,:)
        teste_u(2,:) = offsp2(i,:)
        teste_u(3,:) = fa(i,:)
        teste_u(4,:) = mo(i,:)
           
        teste_pu = indx_t
        teste_u = teste_u(indx_t,:)
   
   
        flag = 0
        IF(teste_pu(1).EQ.3 .OR. teste_pu(1).EQ.4)then
            CALL flip(N,Pu,winner)
                IF(winner(1).EQ.1)then
                    pop(keep+i,:) = teste_u(1,:)
                ELSE  
                    DO k=1,4
                        IF(teste_pu(k).EQ.1 .OR. teste_pu(k).EQ.2)then
                            pop(keep+i,:) = teste_u(k,:)
                            flag = k
                            EXIT
                        ENDIF
                    ENDDO    
                ENDIF 
        ENDIF
   		   
   		   
        IF(teste_pu(2).EQ.3 .OR. teste_pu(2).EQ.4)then
            CALL flip(N,Pu,winner)
                IF(winner(1).EQ.1)then
                    pop(keep+i+(Nsel/2),:) = teste_u(2,:)
                ELSE  
                    DO k=1,4
                        IF(teste_pu(k).EQ.1 .OR. teste_pu(k).EQ.2)then
                            IF(flag.NE.k)then
                                pop(keep+i,:) = teste_u(k,:)
                            ENDIF
                            EXIT
                        ENDIF
                    ENDDO    
                ENDIF 
        ENDIF
   		   
   		   
        IF(teste_pu(1).EQ.1.AND.teste_pu(2).EQ.2)then	   
            pop(keep+i,:) = teste_u(1,:)		
            pop(keep+i+(Nsel/2),:) = teste_u(2,:)
        ENDIF
        
        IF(teste_pu(1).EQ.2.AND.teste_pu(2).EQ.1)then
            pop(keep+i,:) = teste_u(1,:)		
            pop(keep+i+(Nsel/2),:) = teste_u(2,:)
        ENDIF
             
           
           
           
50	ENDDO

         
CCc ---------------------------------------------------------------------------------------
CCc ---------------------------------------------------------------------------------------
CCc ------------------------- End of verifying update probability -------------------------          
CCc ---------------------------------------------------------------------------------------
CCc --------------------------------------------------------------------------------


        par_p=gadecode(pop,hi_p,lo_p,nbits,popsize,npar)
    
        DO i=1,popsize
            par_rho(i,:) = 0.23*((par_p(i,:))**0.25)
            par_s(i,:) =(par_p(i,:) - 1360)/1.16
        ENDDO
    

        DO 70 i=1,popsize !compute one synthetic seismogram for each elem.

            CALL d2t(par_p(i,:),par_s(i,:),par_rho(i,:),thickness,dt,nl,
     1           vp_it, vs_it, rho_it, depth_it,time_series,nopt2)
       
            CALL synt(vp_it,vs_it,rho_it,depth_it,
     1           nopt2,synthetic,grc)
     
            reflex(i) = synthetic(147,1)

            DEALLOCATE(vp_it)
            DEALLOCATE(vs_it)
            DEALLOCATE(rho_it)
            DEALLOCATE(depth_it)
            DEALLOCATE(time_series)
		   
            IF(nopt1.GT.nopt2)THEN
            CALL xcross2d(synthetic(1:nopt1,:),model_final(1:nopt1,:),
     1	         nopt1,no_traces,rr)
            ELSE
            CALL xcross2d(synthetic(1:nopt2,:),model_final(1:nopt2,:),
     1	         nopt2,no_traces,rr)
            ENDIF				
        
        
            IF(rr.lt.0) THEN
                fitness(i) = 0
            ELSE
                fitness(i) = rr
            ENDIF
       	   
c    --- P wave velocity 
            grcc(i) = grc      	   
            count = count + 1

70      ENDDO

C		 
        umax = MAXVAL(fitness)
        umin = MINVAL(fitness)
        uavg = SUM(fitness)/popsize
        
c		   CALL scalepop(fitness,popsize,umax,uavg,umin) !scale the fitness landscape
        CALL sort(popsize,fitness,indx)               !sorting the fitness values
        CALL reverse(indx,popsize)                    !reversing sorting higher to lower
        
        fitness=fitness(indx)
        pop=pop(indx,:)

        CALL CPU_TIME(t2)
        write(*,*)umax,umin,uavg, u, t2-t1

        DO 10 i=1,popsize
            par(i,:) = (/ par_p(i,:),par_s(i,:),par_rho(i,:) /)
10      CONTINUE
 
c		   write(*,*)par(1,:)
        OPEN (unit = 333, file = 'checking.txt')   
        DO 354 b =1,popsize 
            WRITE(333,919)(par(b,f), f=1,3*npar)
354     ENDDO
919     FORMAT(12F14.7)

        DO 778 f=1,3*npar
            DO 777 l=1,popsize
                nindex = 1 + ((par(l,f)-hi(f))/(lo(f)-hi(f)))*(2**nbits) 
                nPPD(f,nindex) = nPPD(f,nindex) + 1
                PPD(f,nindex) = PPD(f,nindex) + fitness(l)		   
777         ENDDO 
778     ENDDO
	
        sumPPD = sumPPD + SUM(fitness)  		
	   
9999    ENDDO  

        
        DO f=1,3*npar
            DO i=1,2**nbits+1
                IF(PPD(f,i).gt.0.AND.nPPD(f,i).gt.0)THEN
                    PPD(f,i) = PPD(f,i)/nPPD(f,i)
                ENDIF	
            ENDDO
        ENDDO 
			
C
        
        
        OPEN (unit = 22, file = 'reflex.txt')   
        DO 13 b =1,2**nbits+1 
            WRITE(22,99)(PPD(f,b), f=1,3*npar)
13      ENDDO
99      FORMAT(12F14.7)
			
			
			
			DO f=1,3*npar
			DO i=1,2**nbits+1
			IF(PPD(f,i).gt.0.AND.nPPD(f,i).gt.0)THEN
				PPD(f,i) = (nPPD(f,i) * PPD(f,i))/sumPPD
			ENDIF	
			ENDDO
			ENDDO 
			

C
		OPEN (unit = 23, file = 'reflex_sumPPD.txt')   
		DO 12 b =1,2**nbits+1 
        	WRITE(23,98)(PPD(f,b), f=1,3*npar)
12      ENDDO
98     FORMAT(12F14.7)	
C
C
C			DO f=1,3*npar
C			DO i=1,2**nbits+1
C			IF(PPD(f,i).gt.0.AND.nPPD(f,i).gt.0)THEN
C				PPD(f,i) = nPPD(f,i)
C			ENDIF	
C			ENDDO
C			ENDDO 
C
C		OPEN (unit = 21, file = 'reflex_nPPD.txt')   
C		DO 11 b =1,2**nbits+1 
C        	WRITE(21,97)(PPD(f,b), f=1,3*npar)
C11      ENDDO
C97     FORMAT(6F14.7)
C		
C		
        write(*,*) "Numero de modelagens calculadas: ", count
        write(*,*) "Tempo gasto (min): ", (t2-t1)/60
C      
C 
C
Cc        OPEN (unit = 22, file = 'synthetic.txt') 
Cc        DO vv=1,popsize
Cc        WRITE(22,999)(par(vv,j),  j=1,3*npar)
Cc        ENDDO
Cc999     FORMAT(6F14.7)
C       		  
       END
       
     
c      - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  
c      sub function to compute the seed for the random number generation
!23456   
      subroutine init_random_seed()
      use iso_fortran_env, only: int64
        implicit none
        integer, allocatable :: seed(:)
        integer :: i, n, un, istat, dt(8), pid
        integer(int64) :: t
          
        call random_seed(size = n)
        allocate(seed(n))
       
!23456             ! First try if the OS provides a random number generator
        open(newunit=un, file="/dev/urandom", access="stream", 
     1   form="unformatted", action="read", status="old", iostat=istat)
            if (istat == 0) then
               read(un) seed
               close(un)
            else
               ! Fallback to XOR:ing the current time and pid. The PID is
               ! useful in case one launches multiple instances of the same
               ! program in parallel.
               call system_clock(t)
               if (t == 0) then
                  call date_and_time(values=dt)
                t = (dt(1) - 1970) * 365_int64 * 24 * 60 * 60 * 1000 
     1              + dt(2) * 31_int64 * 24 * 60 * 60 * 1000 
     1             + dt(3) * 24_int64 * 60 * 60 * 1000 
     1              + dt(5) * 60 * 60 * 1000 
     1             + dt(6) * 60 * 1000 + dt(7) * 1000 
     1              + dt(8)
               end if
               pid = getpid()
               t = ieor(t, int(pid, kind(t)))
               do i = 1, n
                  seed(i) = lcg(t)
               end do
            end if
            call random_seed(put=seed)
          contains
            ! This simple PRNG might not be good enough for real work, but is
            ! sufficient for seeding a better PRNG.
            function lcg(s)
              integer :: lcg
              integer(int64) :: s
              if (s == 0) then
                 s = 104729
              else
                 s = mod(s, 4294967296_int64)
              end if
              s = mod(s * 279470273_int64, 4294967291_int64)
              lcg = int(mod(s, int(huge(0), int64)), kind(0))
            end function lcg
          end subroutine init_random_seed
      

 

