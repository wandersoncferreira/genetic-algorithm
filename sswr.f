        SUBROUTINE SSWR(fitness,popsize,choices)
        IMPLICIT NONE
        
        INTEGER popsize,winner(1),N,it
        REAL*4 fitness(popsize),expected,fractional(popsize)
        REAL*4 randx,uavg
        INTEGER k,jassign,j,choices(popsize)
      
        N = 1      
        j = 0
        k = 0
      
        CALL zero(choices,popsize)
        CALL zero(fractional,popsize)
      
        uavg = SUM(fitness)/popsize
c       write(*,*)uavg

        DO 1 j=1,popsize
c         j = j + 1
            expected    = fitness(j)/uavg
            jassign     = AINT(expected)
            fractional(j) = expected - jassign
c            if(expected > 1) then
c                write(*,*) "Expected:", expected
c            endif    
            
            DO WHILE(jassign.GT.0) 
                k = k + 1
                jassign = jassign - 1
                choices(k) = j
            ENDDO		 
1       ENDDO


        j = 0

        DO WHILE(k.LT.popsize)
            j = j + 1
         
            IF(j.GT.popsize) THEN
                j = 1
            ENDIF
        
            IF(fractional(j).GT.0)THEN
                CALL flip(N,fractional(j),winner(1))         
                IF(winner(1).EQ.1)THEN
                    k = k + 1
                    choices(k) = j
                    fractional(j) = fractional(j) - 1
                ENDIF
            ENDIF 
        ENDDO

        RETURN
        END 
