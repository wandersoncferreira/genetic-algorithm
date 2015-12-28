	subroutine ray_nmo(rc,tx,px,pwavet,notd,ang_inc,ang_cut, !input
     1      samp_rate, ifirst,nout,tmax_raw, max_ray,maxn_raw, 	!input
     1      max_fil, wave, nwave, xspread, t_mute,						!input
     1	    temp1, temp2, toptime, i_nmo, xcur,	vrms,nosamp_out,!temp
     1      i_beg_tx,tr_out,in_mini  )						!output		
c----------------------------------------------------------------
c  This routine applies a geophone array and filter in the raw time
c  state and then nnmo corrects the trace.
c
c  INPUT
c 	rc()	= reflection coefficients in tzero time
c  	tx()    = tx times at tzero time.
c       px()    = first layer angle for computing layer angle 
c       pwavet()= pwave velocity of model
c       notd    = number samples in tzero arrays
c	wave()  = seismic wavelet
c       nwave   = number samples in seismic wavelet
c       ang_inc = angle increment in incident angle table
c       ang_cut = incident angle to mute at each refletion
c	samp_rate = samp_rate in seconds
c	ifirst  = first live sample in tzero array to map
c       nout    = number of samples in the nmo-corrected trace.
c                 This is notd plus a defined padding time.
c	tmax_raw= maximum arrival time in seconds for raw tr_out array
c	xspread  = length of geophone array in feet
c       t_mute   = time to taper mute front of synthetic in sec. 
c       i_nmo 	= flag for type of nmo
c		= 0 ray trace
c		= 1 VRMS 
c		= 2 static shift
c               = 3 no NMO
c               = 4 flat events - no stretch
c	xcur    = offset distance
c       vrms    = vrms velocity
c       nosamp_out = number of output samples to write if i_nmo =3
c       in_mini = flag=1 for minimum phase wavelet; =0 for zero phase
c
c  OUPUT
c      tr_out() = nmo corrected trace with nout samples
c      
c  TEMPORARY
c      temp1() = time trace holding the raw field trace 
c      temp2() = second temp array 
c
c  ARRAY DIMENSIONS
c	max_ray	= number samples in time arrays
c	maxn_raw= number samples in temporary raw trace 
c	max_fil = filter array samples
c	   
c
c-----------------------------------------------------------------
c Hilterman  -  September 2004
c-----------------------------------------------------------------
	implicit none 
	integer max_ray,max_fil,maxn_raw
	real*4 rc(max_ray), tx(max_ray), px(max_ray), pwavet(max_ray)
	real*4 wave(max_fil), tr_out(max_ray), temp1(maxn_raw)
	real*4 temp2(maxn_raw), tspread, pray, t_array(15), tdel_ph
	real*4 t_inc, tfac, xphone, tsamp, znil, t_mute, angr_cut
	real*4 ang_inc, ang_cut,samp_rate, tmax_raw,xspread,ang_cur
	real*4 toptime, tshift, xcur, xx, txx, vrms(max_ray)
	
	integer notd, nwave, ifirst, nout,nphone, noph_half,j,i,mt
	integer  ibeg,  iend_raw, nzero,  mtzero,istar
	integer i_nmo, mt_max, nosamp_out,  i_beg_tx
	integer in_mini,kk
c	integer m_mute,istart_cut,istart

c------------------
c  preset parameters ....
c-----------------
    
        
	nphone = 15	! number of geophones in array
	noph_half = 8
	angr_cut = sin(ang_cut* (0.017453293)) !sine of specified cutoff angle 
	xphone = nphone
	znil = 1.e-10
c--------------------
c  zero temp arrays ...
c--------------------

	call zero(temp1,maxn_raw)
	call zero(temp2,maxn_raw)
	

c--------------------
c  apply geophone array from end to first live sample ...
c--------------------
	ibeg = ifirst + 1
	if(ibeg.lt.2) ibeg =2
	
	do 100 j=notd, ibeg,-1
	if(abs(rc(j)).lt.znil) go to 100  ! no contribution because rc=0
	if(tx(j-1).gt.tmax_raw) go to 100  ! too great arrival time
	

c  compare incident angle to cutoff angle ...
	pray = sin((px(j-1)-1.)*ang_inc *(0.017453293))/ pwavet(1) ! ray parameter
	ang_cur = pwavet(j-1) * pray ! sine of current incident angle
	if(ang_cur.gt.angr_cut) then
	  go to 120  !zeroing data above cutoff angle
	endif
	

	
c  designing arrival times for each geophone
 	tsamp = tx(j-1) / samp_rate  !time in samples  
 	if(i_nmo.eq.4) tsamp = float(j)
	tspread = (xspread * pray)/ samp_rate  ! total traveltime across geophone array in samples 	
	t_array(noph_half) = tsamp  !arrival at mid phone that is considered at offset
	tdel_ph = tspread/float(nphone-1)  ! arrival time incement between geophones
	  do 20 i=1,noph_half-1
	  t_inc = float(i) * tdel_ph    ! time from center of geophone array
	  t_array(noph_half-i) = tsamp - t_inc
20	  t_array(noph_half+i) = tsamp + t_inc
	
		

c  placing geophone response along with reflection coefficient at arrtival
c  time in temp array ...

	  do 40 i=1,nphone
	  mt = t_array(i)
	  tfac = t_array(i) - mt  
	  temp1(mt) = temp1(mt) + (1.-tfac)*rc(j)/xphone	  
 	  temp1(mt+1) = temp1(mt+1) + (tfac)*rc(j)/xphone
 	  
40	  continue




100	continue    ! end of loop over tzero times
120	continue   
        

c----------------------------	  
c  taper mute front end of trace ...
c----------------------------------
c  Remove mute May 25, 2009
c 	istart_cut = tx(j-2)/samp_rate    !March 2008
c	istart = tx(ifirst) / samp_rate
c 	if(istart.lt.istart_cut) istart = istart_cut
c	m_mute = t_mute / samp_rate  ! taper mute in samples
c	do 50 i=1,m_mute
c	tfac = float(i)/float(m_mute+1)	  
c50	temp1(istart+i-1) = tfac * temp1(istart+i-1) ! un-nmo trace

c-------------------------------
c convolve with wavelet ...
c-------------------------------
     
	iend_raw = tx(notd)/samp_rate
c	write(*,*)samp_rate,tx(notd),iend_raw
c	write(*,*)(maxn_raw-3*max_fil)
	if(iend_raw.lt.0) then  ! no live samples  
	  call zero(tr_out,notd+1)
	  return
	endif	
	
	if(iend_raw.gt.(maxn_raw-3*max_fil)) 
     1              iend_raw=maxn_raw-3*max_fil
	nzero = nwave/2
	nzero = nzero + 1 ! time zero postion for wavelet
	
	       
	if(in_mini.eq.1) nzero=1
c    	call convol3(temp1,ifirst,iend_raw,wave,1,nwave,nzero,temp2,
      	call convol3(temp1,0,iend_raw,wave,0,nwave,nzero,temp2,
     1       maxn_raw, max_fil)
     
		

         

c------------------------------	
c  do nmo ... ray-trace method
c------------------------------
	call zero(tr_out,notd+1)
    
 	
c check type of NMO ...

c ray-trace NMO
	if(i_nmo.eq.0) then   ! conventional NMO
c   	  do 300 i=notd,ibeg, -1
c   	  do 300 i=notd,2, -1
    	  do 300 i=notd,i_beg_tx, -1   	   	  
 	    mt = tx(i)/samp_rate 
 	    tfac = (tx(i)/samp_rate) - float(mt)
300 	  tr_out(i) = tr_out(i) + temp1(mt)*(1.-tfac)+temp1(mt+1)*tfac
 	  return
 	endif

    

c Vrms NMO
	if(i_nmo.eq.1) then   ! VRMS NMO
	  mt_max = tx(notd) / samp_rate  ! max live sample in raw trace
 
c  	  do 350 i=notd,ibeg, -1
c   	  do 350 i=notd,2, -1 
    	  do 350 i=notd,i_beg_tx, -1
 	    xx = (xcur/samp_rate)**2  ! 
 	    txx = sqrt(float(i*i) + (xx/(vrms(i)**2)))
 	    mt = txx   ! sample for tx
 	    if(mt.gt.mt_max) then
 	      tr_out(i) = 0.
 	      go to 350
 	    endif
 	    tfac = (txx) - float(mt)
	  tr_out(i) = tr_out(i) + temp1(mt)*(1.-tfac)+temp1(mt+1)*tfac
350	  continue
 	  return
 	endif
 	
c  static shift NMO ...
	if(i_nmo.eq.2) then  !static shot NMO   Stage A
	  mtzero = toptime/samp_rate    ! tzero sample to flatten to
	  if(mtzero.gt.notd) then	  
	    write(3,*) 'Flattening depth greater than well depth' 
	    return
	  endif	
	  if(mtzero.le.ibeg) then	  
 	    write(3,*) 'Flattening depth past criticaal angle cutoff' 
	    return
	  endif		
		
	  tshift  = tx(mtzero) - (float(mtzero) * samp_rate)  ! trace shift (sec)
 	  mt =  tshift/samp_rate	      ! trace shift (samples)
 	  tfac =  tshift/samp_rate - float(mt) 	!fraction of shift
	  istar = notd
	  if((istar+mt+1).gt.iend_raw) then
	     istar = iend_raw - mt-1
	     if(istar.le.ibeg) return
	  endif 
	
 	  do 400 i=istar,ibeg, -1		 	
400 	  tr_out(i)=tr_out(i)+temp1(i+mt)*(1.-tfac)+temp1(i+mt+1)*tfac 		
	endif	! stage A

c  no nmo ...
	if(i_nmo.eq.3) then
	  do 500 i=1,nosamp_out
500	  tr_out(i) = temp1(i)
	endif

c no mo ior nmo corrections ...
	if(i_nmo.eq.4) then
	  do 600 i=1,notd
600	  tr_out(i) = temp1(i)
	endif

c	    do 812, j=1,nosamp_out
c	    write(21,799) j, tr_out(j)
c812	continue
c799	format (i8,5x, e10.4)

c 	  write(21) (tr_out(i), i=1,nosamp_out)
 

c	   OPEN (unit = 44, file = 'generated_data.txt')
c       
c       DO kk = 1, maxn_raw
c       WRITE(44,99)(tr_out(kk))
c       ENDDO
c99     FORMAT(1F20.5)

	
 	return
 	end
