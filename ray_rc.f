	subroutine ray_rc(px, nopt,vel,swave,den,ang_inc, max_ray,ips, !input
     1                    i_rock,i_mute_top,i_mute_bot,			!input  
     1     rc, ifirst,grc)					  !output
c-------------------------------------------------------------------
c  This routine computes the approximate Zoeppritz rc
c
c  INPUT:
c     px()     = array with location points in angle array
c     nopt     = number of samples in active tzero time array
c     max_ray  = parameter for indexing the array ray
c     nang     = number of angles 
c     vel()    = pwave velocity time array
c     swave()  = shear wave velocity array
c     den()    = density array
c     ang_inc  = angle increment for incident angle
c     ips      = 0 PP path
c    	       = 1 PS path
c     i_rock   = indicates total or partial AVO response
c                 = 0  Plot total Shuey AVO response
c		  = 1  Plot Shuey AI contribution   if ips=0
c		       Plot Bort Shear contribution if ips=1 
c                 = 2  Plot Shuey PR contribution   if ips=0
c		       Plot Bort density contribu   if ips=1 
c                 = 3  Plot Shuey delV/V contribu   if ips=0
c    i_mute_top = mute trace to this sample to avoid reflections from
c                 top of upper layer, "topmed".
c    i_mute_bot = mute trace to this sample to avoid reflections from
c                 bottom of lower layer "botmed".
c
c  OUTPUT
c     rc()	= refelction coefficients ... -999. is a null value
c     ifirst    = ifirst non-zero rc in array ....
c     
c-------------------------------------------------------------------
c  Hilterman    September 2004
c-------------------------------------------------------------------
	implicit none
	integer max_ray
	real*4 px(max_ray),rc(max_ray),vel(max_ray),den(max_ray),
     1          swave(max_ray)
	real*4 ang_inc, test, zclose, angr, ang, ang_IN, angr_inc
	real*4 grc
 
	
	integer nopt,j, ifirst, kk , ips, i_rock
	integer i_mute_top, i_mute_bot
	
c----------------
c  preset 
c----------------
	rc(1) = 0.      ! note first sample in array is for time zero	
	zclose = .999

c---------------------------------
c  loop of time samples ....
c---------------------------------	
     
	do 900 j=1,nopt-1
	if(px(j).eq.-999.) then
	   rc(j+1) = -999.
	   go to 900
	endif
	
c	ang_IN = px(j) * ang_inc
	ang = (px(j)-1.) * ang_inc
	angr = ang * .017453293  ! radians
c	angr_inc = ang_IN * .017453293 

c critical angle test	
c	test = (vel(j+1)/vel(1)) * sin(angr_t) ! FEB 14, 2011
  	test = (vel(j)/vel(1)) * sin(angr)  ! original code	
	
		
c  new test code ...
c  	test2 = (vel(j+1)/vel(j)) * sin(angr) ! feb 14, 2011	
	if(test.gt.zclose) then  ! incident at boundary feb 14, 2011
c	if(test.gt.zclose.or.test2.gt.zclose) then  ! incident at boundary	
	    rc(j+1) = -999.
	    go to 900
	endif 
	
	  
	angr = asin(test)  ! incident at boundary
c	  WRITE(*,*) angr
	
c  calculate reflection coefficient ...
	if(ips.eq.0) then
c	  call aki_single(vel(j),vel(j+1),swave(j),swave(j+1) !input
c     1   , den(j),den(j+1), angr,     			!input
c     1    rc(j+1))					!output
     
      
	  call shuey(vel(j),vel(j+1),swave(j),swave(j+1) !input
     1   , den(j),den(j+1), angr, angr_inc,  i_rock,   			!input
     1    rc(j+1))					!output   
     
     	IF(rc(j+1).GT.0)then
     		grc = rc(j+1)
     	ENDIF	

c       WRITE(*,*) j 
     
c	else
c	  call aki_sing_ps(vel(j),vel(j+1),swave(j),swave(j+1) !input
c     1   , den(j),den(j+1), angr,  i_rock,   			!input
c     1    rc(j+1))					!output	
	
	
	endif
900	continue			

c---------------------------------
c  find deepest -999. and zero rc upward ...
c---------------------------------
	ifirst = 1
	do 950 j=1,nopt
	kk = nopt- j +1
	if(rc(kk).eq.-999.) then
	   ifirst = kk
	   go to 960
	endif
950	continue  

960	continue
     

c  zero rc array from 1 to ifirst
	if(ifirst.lt.i_mute_top) ifirst = i_mute_top
	do 980 j=1,ifirst
980	rc(j) = 0.

c  zero end of rc array ...
	do 990 j= i_mute_bot,nopt
990	rc(j) = 0.  

c-------------------------
c  End of program
c-------------------------	 
	return
	end

