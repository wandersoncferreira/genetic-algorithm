      subroutine convol3(trace,it1,it2,filter,if1,if2,nzero,result,
     1       maxn_raw, max_fil)
c-----------------------------------------------------------------
c  THIS ROUTINE CONVOLVES PORTIONS OF THE REAL ARRAYS TRACE AND
c  FILTER AND RETURNS THE ARRAY BACK IN TRACE.
c
c  ARGUMENT DEFINITIONS . . .
c
c     TRACE  = REAL ARRAY TO BE CONVOLVED WITH FILTER.
c     IT1      = THE FIRST WORD OF THE TRACE ARRAY TO BE USED FOR THE
c        CONVOLUTION.
c     IT2      = THE LAST WORD OF THE TRACE ARRAY TO BE USED FOR THE
c        CONVOLUTION.
c     FILTER   = REAL ARRAY TO BE CONVOLVED WITH THE TRACE ARRAY.
c     IF1      = THE FIRST WORD OF THE FILTER ARRAY TO BE USED FOR THE
c     IF2      = THE LAST WORD OF THE FILTER ARRAY TO BE USED FOR THE
c        CONVOLUTION.
c     NZERO    = SAMPLE WITHIN THE FILTER ARRAY WITH RESPECT TO THE
c                FIRST SAMPLE WHERE TIME ZERO IS LOCATED.
c                FOR SYMMETRICAL FILTERS, THIS IS THE CENTER SAMPLE
c
c---------------------------------------------------------------------
!234567    
		IMPLICIT NONE
		integer maxn_raw, max_fil, k
		real*4 trace(maxn_raw), result(maxn_raw),filter(max_fil)
	
		integer if1,if2,it1,it2,nfil, j, icur, nzero


c	  PAD TO AVOID EDGE EFFECTS . . .
			nfil = if2-if1+1
			DO 10 j=1,nfil
			result(j) = 0.
					                     
                                                 !
c	                              ! was TRACE(IT1)
10	        result(j+nfil+it2-it1+1) = 0.                                   !
c	          	                  !was TRACE(IT2)
			DO 15 j=it1,it2
15	        result(j+nfil-it1+1) = trace(j)

			icur = nfil+1 +nzero-1
			do 4 j=it1,it2
			trace(j) = 0.
			
			do 6 k=if1,if2
6           trace(j) = trace(j) +filter(k)*result(icur+if1-k)
4           icur = icur +1

				

		  RETURN
		  END
