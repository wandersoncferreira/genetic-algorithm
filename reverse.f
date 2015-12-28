!2345678
       SUBROUTINE reverse(array,popsize)
       
       IMPLICIT NONE
       REAL array(*),temp
       INTEGER popsize,K
       
       DO 41 K=1,popsize/2
          temp=array(K)
          array(K)=array(popsize+1-K)
          array(popsize+1-K)=temp
41     ENDDO
       END     
          
       