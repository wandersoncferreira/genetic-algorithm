	subroutine ray_extract(ray, nopt, offset, index_ang,	!input
     1        max_ray, nang, samp_rate , ang_inc,pray_sin,      !input
     1       px,tx, i_beg_tx, ierr  )					!output
c-------------------------------------------------------------------
c  This routine finds the traveltime and ray parameter for each time
c  sample associated with the a particular offset
c
c  INPUT:
c    ray(i,j,k)	= i=35 different angles,from 0-70 degrees in 1st layer
c 		  j=1 traveltime
c		  j=2 offset distance
c                 k = number of samples
c     nopt = number of samples in active tzero time array
c     max_ray = parameter for indexing the array ray
c     offset   = offset of current AVO trace
c     nang     = number of angles to search
c     index_ang()  = index specifying where to start in angle loop ... 
c     samp_rate  = sample rate in sec.
c     ang_inc    = increment in incident angle
c     pray_sin() = array of sin square incident angles searched
c
c  OUTPUT
c     tx()	= traveltimes for each trace
c     px()	= location of offset distance in angle loop is isotropic
c		= angle in radians if anisotropic
c     i_beg_tx  = sample for first possible NMO traveltime
c     ierr      = 0 success
c               = 1 problems 
c     
c-------------------------------------------------------------------
c  Hilterman September 2004
c-------------------------------------------------------------------
	implicit none
	integer max_ray, nang
	real*4 ray(36,3,max_ray), offset, tx(max_ray), px(max_ray)
	real*4 offsq , xnil, pray_sin(nang)
	real*4 w1, w2, samp_rate , fac,  ang_inc
	real*4 zz
c	real*4 x_bug(36), t_bug(36), p_bug(36)   !remove tdebug

	integer nopt,index_ang(max_ray),i,j, ierr
	integer i_beg_tx
c	integer kk

c  Presets
	 i_beg_tx = 1
	ierr = 0
c	fac = (3.141592654 / 180.)
	 fac = 1.    ! angles in ray(1,3,j) are radians
	 xnil = 1.
	 offsq = offset**2	
c-----------------------
c  test for zero offset ...
c-----------------------

    
    
	if(offset.lt.xnil)then
	  tx(1) = samp_rate
	  px(1) = 1.    ! refers to angle 0 degrees
	  do 50 i=2,nopt
	  tx(i) = tx(i-1) + samp_rate
	  px(i) = 1.
50	continue
          go to 9999
        endif
		


c---------------------------------
c  loop layers ..
c---------------------------------	
	do 900 i=1,nopt   	
c------------------------------	
c   Loop over angles ...
c--------------------------------
    
 	do 800 j= index_ang(i), nang 
	if(ray(j,2,i).gt.offset) go to 700   !distance is gt than desired offset
	
800	continue 	! loop over all incident angles  

       
c  No solution found for layer i ..
	index_ang(i) = nang 
	
	tx(i) = -999.
	px(i) = -999.
	i_beg_tx = i
	go to 900
	
c  Solution offset found ...
700	continue	! solution for offset found
			! interpolate with x**2, t**2	
	index_ang(i) = j 
 	w1 =(ray(j,2,i)**2-offsq)/(ray(j,2,i)**2-ray(j-1,2,i)**2)
 	w2 = 1. - w1
	tx(i) = sqrt(w2*(ray(j,1,i))**2 + w1*(ray(j-1,1,i))**2)
	
	
c  isotropic case ...

	  zz = sqrt(w2*pray_sin(j) + w1*pray_sin(j-1))
c	  write(*,*) zz
	  if(zz.gt.1.) then  
	       ierr = 1
	       go to 9999
	  endif
	  zz = (180./3.141592654) * asin(zz) 
  	  px(i) = zz/ang_inc + 1. 
       

900	continue  	!loop over number of time samples 

c-------------------------
c  end of program
c-------------------------
	 
9999	continue
	i_beg_tx = i_beg_tx + 1
	return
	end

