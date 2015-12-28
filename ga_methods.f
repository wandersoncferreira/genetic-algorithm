      MODULE ga_methods
      
          
      
      CONTAINS 
      
       FUNCTION gadecode (pop,hi,lo,nbits,popsize,npar)
    
       IMPLICIT NONE
       REAL*4 pop(popsize,nbits*npar),popT(nbits*npar,popsize)
       REAL*4 hi,lo
       INTEGER  nbits, nparM, popsize, npar,i
       REAL*4, DIMENSION(:,:), ALLOCATABLE  :: quant
       REAL*4, DIMENSION(:,:), ALLOCATABLE  :: ct
       REAL*4, DIMENSION(:,:), ALLOCATABLE  :: ctT 
       REAL*4, DIMENSION(:,:), ALLOCATABLE  :: par
       REAL*4, DIMENSION(:,:), ALLOCATABLE  :: f
       REAL*4, DIMENSION(:,:), ALLOCATABLE  :: fTT
       REAL*4, DIMENSION(:,:), ALLOCATABLE  :: gadecode
       
      
c      creating the decode function based on the number of bits and variables    
      
         
       ALLOCATE(quant(nbits,1))
       DO 20 i = 1,nbits
             quant(i,1) = 0.5**i
20     CONTINUE
       
       DO 30 i = 1,nbits
             quant(i,1) = quant(i,1)/SUM(quant(:,1))
30     CONTINUE
       
       popT = TRANSPOSE(pop)
       nparM = npar*popsize
       
c      allocating the size of the variable ct

       ALLOCATE(ct(nbits,nparM))
       
       ct = RESHAPE(popT, (/ nbits, nparM /) )
       
c      in case of non-squared matrices:
       ALLOCATE(ctT(nparM,nbits))
       ctT = TRANSPOSE(ct)
       
c      matrix multiplication

       ALLOCATE(par(nparM,1))
       par = ((MATMUL(ctT,quant))*(hi - lo) + lo)
       
       ALLOCATE(f(npar,popsize))
       f = RESHAPE(par,(/ npar, popsize /) )
       
       ALLOCATE(fTT(popsize,npar))
       fTT = TRANSPOSE(f)
       
       
       ALLOCATE(gadecode(popsize,npar))
       gadecode = fTT
       pop = TRANSPOSE(popT)
       
       DEALLOCATE(fTT)
       DEALLOCATE(f)
       DEALLOCATE(par)
       DEALLOCATE(ctT)
       DEALLOCATE(ct)
       DEALLOCATE(quant)
       
       RETURN
       END       
       
       
        SUBROUTINE scalepop(objective,popsize,umax,uavg,umin)
        IMPLICIT NONE
            INTEGER popsize,i
            REAL*4 sumfitness,a,b,umax,umin,uavg
            REAL*4 objective(popsize),fitness(popsize)

            CALL zero(fitness,popsize)
            CALL prescale(umax, uavg, umin, a, b)
            sumfitness = 0
       
            DO 10 i = 1, popsize
                CALL scale(objective(i),a,b)
                IF (objective(i) .LT. 0) then
                    objective(i) = 0
                ENDIF
                sumfitness = sumfitness + fitness(i)
10          CONTINUE

       RETURN
       END 
    

        SUBROUTINE prescale(umax, uavg, umin, a, b)
        IMPLICIT NONE
        REAL*4 umax, uavg, umin, flag, delta, fmultiple, a,b
       
        fmultiple = 2
        flag = (fmultiple*uavg - umax)/(fmultiple - 1.0 )
        
        IF (umin .GT. flag) then
            delta = umax - uavg
            a = (fmultiple - 1.0) * uavg/delta
            b = uavg * (1-a)
        ELSE
            delta = uavg - umin
            a = uavg/delta
            b = -umin*uavg/delta
        ENDIF
        
    
        IF (delta.LT.0.000001.AND.delta.GT.-0.000001)then
            a = 1
            b = 0
        ENDIF
       
        RETURN
        END
       
       
        SUBROUTINE scale(objective,a,b)
        IMPLICIT NONE 
        REAL objective, a, b      
            objective = a * objective + b
        RETURN
        END 
       

!23456
       SUBROUTINE d2t(Vp, Vs, rho, thickness, dt, nlayer,
     1            vp_it, vs_it, rho_it, depth_it,time_series,nopt)
       
       PARAMETER (maxx=6000)
       INTEGER nlayer,i,tpp(nlayer),flag,j,kk,nopt
       REAL*4 Vp(nlayer), Vs(nlayer), rho(nlayer), thickness(nlayer)
       REAL*4 dt, start_log_time
       REAL*4 tf,soma,tp(nlayer),sumt
       REAL*4, DIMENSION(:), ALLOCATABLE :: vp_it, vs_it, rho_it
       REAL*4, DIMENSION(:), ALLOCATABLE :: depth_it, time_series
       
c      defining parameters for starting log time and additional time (tf)       
       start_log_time = 0.0
       tf = 0.0
       
c      computing the traveltimes for each layer in order to discretize
c      the model
       tp = 2*thickness/Vp
       tp(1) = tp(1) + start_log_time
       tp(nlayer) = tp(nlayer) + tf
       
       soma= 0
       DO 10 i=1,nlayer
       soma = soma + tp(i)
       tp(i) = soma
10     CONTINUE  

c      Allocating space for the output variables
	   tpp = INT(FLOOR(tp/dt))
	   nopt = tpp(nlayer)
       ALLOCATE(vp_it(tpp(nlayer)))
       ALLOCATE(vs_it(tpp(nlayer)))
       ALLOCATE(rho_it(tpp(nlayer)))
       ALLOCATE(depth_it(tpp(nlayer)))
       
       flag = 1
    
       DO 20 j = 1,tpp(nlayer)
       
          IF (j .LE. tpp(flag)) then
              vp_it(j) = Vp(flag)
              vs_it(j) = Vs(flag)
              rho_it(j)= rho(flag)
              depth_it(j) = (Vp(flag) * dt)/2
              
          ELSE
              flag = flag + 1
              vp_it(j)    = Vp(flag)
              vs_it(j)    = Vs(flag)
              rho_it(j)   = rho(flag)
              depth_it(j) = (Vp(flag) * dt)/2
          END IF
          
20    CONTINUE

c      computing the traveltime series
       ALLOCATE(time_series(tpp(nlayer)))
       DO 30 j = 1,tpp(nlayer)
       time_series(j) = j
30     CONTINUE            
       time_series = time_series*dt
       
       sumt = 0
       DO 40 j = 1,tpp(nlayer)
       depth_it(j) = depth_it(j)+sumt
       sumt = depth_it(j)
40     CONTINUE 

       
c       DEALLOCATE(vp_it)
c       DEALLOCATE(vs_it)
c       DEALLOCATE(rho_it)
c       DEALLOCATE(depth_it)
c       DEALLOCATE(time_series)
       
                     
       RETURN
       END    
          
          
       SUBROUTINE synt(vp_it,vs_it,rho_it,depth_it,
     1            nopt,synthetic,grc)
       
       INTEGER nopt
       PARAMETER (maxn_raw=24000,max_ray=6500,npoint=251,no_traces=12)
       REAL*4 :: n,m,tmax_raw, ray(36,3,max_ray),samp_rate
       REAL*4 vp_it(max_ray),vs_it(max_ray),rho_it(max_ray)
       REAL*4 depth_it(max_ray)
       REAL*4 temp1(maxn_raw),temp2(maxn_raw),xoff_beg,xinc,pi
       REAL*4 pray_sin(36), zzz,px(max_ray),tx(max_ray),ang_cut
       REAL*4 rc(max_ray),ricker_wavelet(npoint),dtn,timeT, ang_inc
       
       INTEGER nang, ips,ttt
       INTEGER i_beg_tx, ierr,i_rock,i_mute_top,i_mute_bot,ifirst
       INTEGER Fc,uu,index_ang(max_ray),vv
       INTEGER nout, max_fil, xspread, t_mute
       INTEGER i_nmo, nosamp_out, toptime, vrms, in_mini
       REAL*4 tr_out(max_ray),xcur,raw_d(max_ray)
       REAL*4 synthetic(max_ray,no_traces),grc
       REAL*4 t1, t2, model_final(max_ray,no_traces)
       

c      defining variables
       pi = 3.141592654
       samp_rate = 0.004
       nang = 36        
       ang_inc = 2
       ips = 0
       i_rock = 0
       i_mute_top = 1
       i_mute_bot = nopt
       in_mini = 0
       tmax_raw = float(max_ray)*samp_rate
       

       xoff_beg = 0
       xinc = 100
       
       
c      wavelet parameters
       Fc = 30
       timeT = 0.5
       dtn = 0.004
       
       
c      nmo parameters
       ang_cut = 89.
       nout    = nopt + .1
       max_fil = 301
       xspread = 0.
       t_mute = .020
       i_nmo  = 4
       nosamp_out = 10
       toptime = 0.004
       vrms = 0
       
    
        
!2345678901234567890
       CALL ray_tab(vp_it,depth_it,nopt,tmax_raw,nang,ang_inc,  
     1              max_ray,vs_it,ips, 
     +              ray)
    
        
c      preparing the angle index for the function ray_extract     
       DO 30 i = 1,nopt
       index_ang(i) = 1    
30     CONTINUE
     
c      table of sin squared for ray_extract            
       DO 40  i=1,nang
	   zzz = (float(i-1)*ang_inc)*(pi/180.)
	   pray_sin(i) =(sin(zzz))**2.
40     CONTINUE

c      start creating the synthetic trace

       DO 60 i = 1,no_traces
       xcur = xoff_beg + float(i-1) * xinc  ! current offset distance

       
       CALL ray_extract(ray, nopt, xcur, index_ang,	!input
     1        max_ray, nang, samp_rate , ang_inc,pray_sin,      !input
     1       px,tx, i_beg_tx, ierr  )
     

     
       CALL ray_rc(px,nopt,vp_it,vs_it,rho_it,ang_inc, 	!input
     1	     max_ray, ips, i_rock, i_mute_top, i_mute_bot,	!input  
     1     rc, ifirst,grc)						!output

      
c    need to adjust the rc samples by one     
       DO 50  j=1,nopt-1
	     rc(j) = rc(j+1)
                 
50     CONTINUE
	     rc(nopt) = 0.
  
            


c      creating the ricker wavelet	     
	   CALL ricker(Fc,timeT,dtn,npoint,ricker_wavelet)



       CALL ray_nmo(rc,tx,px,vp_it,nopt,ang_inc,ang_cut, !input
     1      samp_rate, ifirst,nout,tmax_raw, max_ray,maxn_raw, 	!input
     1      max_fil, ricker_wavelet, npoint, xspread, t_mute,						!input
     1	    temp1, temp2, toptime, i_nmo, xcur,	vrms,nosamp_out,!temp
     1      i_beg_tx,tr_out,in_mini)						!output	
   

       
       DO 110 vv = 1, nopt
       synthetic(vv,i) = tr_out(vv)
110    CONTINUE       


60    CONTINUE  !end of the number of traces



	  model_final = synthetic
c	  write(*,*) model_final(1:196,2)
	  



c       OPEN (unit = 22, file = 'synthetic.txt')       
c       DO kk = 1,nopt
c       WRITE(22,999)(model_final(kk,j), j =1,no_traces)
c       ENDDO
c999    FORMAT(15F14.7)

c       OPEN (unit = 55, file = 'raw_data.txt')
c       OPEN (unit = 66, file = 'synthetic_seismogram.txt')


       RETURN
       END subroutine
       
       
                 
       
       
       END MODULE ga_methods
