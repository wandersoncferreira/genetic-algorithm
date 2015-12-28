       SUBROUTINE ricker(Fc,timeT,dtn,npoint,ricker_wavelet)
       
       INTEGER Fc, npoint
       REAL sumtime, pi, tt(npoint), p1(npoint), ricker_wavelet(npoint)
       
       pi = 3.14159265358979
       
       sumtime = 0
       
c      creating the time array
 
       DO 10 i = 1,npoint
c       sumtime = sumtime + dtn
       tt(i) = -timeT + sumtime
       sumtime = sumtime + dtn 

            
10     CONTINUE
    
       p1 = Fc*Fc*tt*tt
       ricker_wavelet = (1 - 2*pi*pi*p1)*EXP(-pi*pi*p1)
       
       RETURN
       END
       
