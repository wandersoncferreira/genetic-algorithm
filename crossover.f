!2345678
      SUBROUTINE crossover(mother,father,offspring1,
     1            offspring2,Nt,Pc,Pm)
      PARAMETER (m = 6)
      REAL*4 mother(Nt),father(Nt),r,offspring1(Nt)
      REAL*4 offspring2(Nt),Pc,Pm
      INTEGER inxd(m),flag,j,bernoulli(1),N,bermut(Nt)
      INTEGER nbits,k(m)

      N = 1
      nbits = Nt/6
c     ---- I need to set one random crossover point for each model parameter
c     ---- Let's call m = number of model parameters.
c     ---- I need also to know the nbits size of each model parameters, through nbits
      
c     ---- 
      CALL flip(N,Pc,bernoulli)
	  IF(bernoulli(1).GT.0)THEN
	  		  DO 20 i =1,m
	  			 CALL init_random_seed()
	  			 CALL RANDOM_NUMBER(r)
	  			 k(i)=((i-1)*nbits+1)+FLOOR((i*nbits+1-((i-1)*nbits+1))*r)
20	  		  ENDDO	     
		
	  
      
c    ---- Computing the crossonver
	  		offspring1 = mother
	  		offspring2 = father
	  
	  		 DO 30 j=1,m
	  		 	offspring1(k(j):j*nbits) = father(k(j):j*nbits)
	  		 	offspring2(k(j):j*nbits) = mother(k(j):j*nbits)
30	  	     ENDDO
	   
c     ----     
      ELSE   !no crossover flip == 0
             offspring1 = mother;
             offspring2 = father;
      ENDIF       
      
c     Mutation part of the code. Low rates are usually applied.  
c     Often considered only as a background process in GA schemes    
      
c     Pm = 0.9     !probability of mutation
      
      CALL flip(Nt,Pm,bermut)
      DO 40 i=1,Nt
      	IF(bermut(i).gt.0) then
      		IF(offspring1(i).gt.0) then
      			offspring1(i) = 0
      		ELSE
      			offspring1(i) = 1
      		ENDIF		
      	ENDIF		
40    ENDDO 

      CALL flip(Nt,Pm,bermut)
      DO 50 i=1,Nt
      	IF(bermut(i).gt.0) then
      		IF(offspring2(i).gt.0) then
      			offspring2(i) = 0
      		ELSE
      			offspring2(i) = 1
      		ENDIF		
      	ENDIF		
50    ENDDO         
            
      
      RETURN
      END
      
      SUBROUTINE flip(N,p,bernoulli)
      INTEGER bernoulli(N)
      REAL*4 p,ran(N)
      
      CALL init_random_seed()
      CALL RANDOM_NUMBER(ran)
      
      DO 10 i=1,N
      	IF(ran(i).le.p) then
        	 bernoulli(i) = 1
      	ELSE
        	 bernoulli(i) = 0
      	ENDIF
10    ENDDO      
      RETURN
      END       
      
      