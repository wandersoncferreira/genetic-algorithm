	subroutine ray_tab(vel,dt_samp, nopt, tmax_raw,nang,ang_inc,	!input
     1     max_ray, shear, ips,					! input
     1       ray  )					!output
c-------------------------------------------------------------------
c  This routine computes a ray table for building AVO traces.
c
c  INPUT:
c     vel  = pwave velocity time array
c     dt_samp = depth sample value for each t0 time
c     nopt = number of samples in time array
c     max_ray = parameter for indexing ray
c     tmax_raw = maximum time allowed for tx.
c     nang = number of incident angles to search
c     ang_inc = angle in degrees between incident angles
c     shear() = shear wave veloicty time array,
c     ips    = 0 compute pp ray path
c            = 1 compute ps ray path
c
c  OUTPUT
c    ray(i,j,k)	= i=36 different angles , from 0 to 70 in first layer
c 		  j=1 traveltime
c		  j=2 offset distance
c                 j=3 incident angle in degrees
c                 k = number of samples needed in 
c     
c-------------------------------------------------------------------
c  Hilterman August 2004
c-------------------------------------------------------------------
	implicit none
	integer max_ray 
	real*4 ray(36,3,max_ray), vel(max_ray), dt_samp(max_ray)
	real*4 ang, ang_inc, pi, zclose, p, delz, rad2, tmax_raw 
	real*4 rad_to_ang, shear(max_ray), ang_ps, rad3
	
	integer nopt, nang, i, j, k, ips
		
	pi = 3.141592654
    	zclose = .999	! sin(86.4) test value for reaching critical angle
  	rad_to_ang = 180./pi
c  	ang_inc = 2.
c  	nang = 36	!number of angles in table
  	
c---------------------------------
c  loop over incident angles from 0 to 70 degrees ...
c---------------------------------	
	do 200 i=1,nang
	ang =  ang_inc * float(i-1)
	ray(i,3,1) = ang  !1st layer angle in degrees    	 		
	ang = ang * (pi /180.)	!convert to radians
	p = sin(ang)/ vel(1)  ! ray parameter ...
	
c  build first layer traveltime and offset ...
	ray(i,2,1) = 2. * tan(ang) * dt_samp(1) !1st layer offset
	ray(i,1,1) =(2.*dt_samp(1))/( cos(ang)*vel(1))!1st layer tx
   	if(ips.eq.1) then
  	  ang_ps = asin(p*shear(1)) 
  	  ray(i,1,1) = (ray(i,1,1)/2.) +   !p time
     1	       (dt_samp(1))/( cos(ang_ps)*shear(1))!s time
  	  ray(i,2,1) = (tan(ang)+tan(ang_ps)) * dt_samp(1) !1st layer offset
  	endif
c----------------------
c  loop over time layers ...
c---------------------
	delz = dt_samp(1)
	do 300 j=2,nopt	
     	delz = dt_samp(j) - dt_samp(j-1)  ! thickness of layer 
     	rad2 =  p*p * vel(j) * vel(j)  ! pwave sin**2 in layer
     	
c test for critical angle ...
	if(rad2.gt.zclose) then
50 	   continue   
	   do 100 k=j, nopt
	     ray(i,1,k) = -999.
	     ray(i,2,k) = -999. 
	     ray(i,3,k) = -999.	         	
100	   continue  
	   go to 200
	endif   
	
c  test is PP or PS travel path		   			
	  
	if(ips.eq.0) then  
	  rad2 = sqrt(1-rad2)  ! cos in layer 
	  ray(i,1,j) = ray(i,1,j-1) + 2.* (delz/(vel(j) * rad2))  ! PP t	 	  	  		 	
	  ray(i,2,j) = ray(i,2,j-1) + 2.*( (delz*p * vel(j))/(rad2))  ! offset
	else
	  rad2 = sqrt(1-rad2)  ! cos in layer 	
     	  rad3 = p*p * shear(j) * shear(j)	
	  rad3 = sqrt(1-rad3)  ! shear cos in layer 	   
	  ray(i,1,j) = ray(i,1,j-1) + (delz/(vel(j) * rad2))  +
     1      		(delz/(shear(j) * rad3))                   ! PS t	 	
	  ray(i,2,j) = ray(i,2,j-1) + ((delz*p * vel(j))/(rad2))  
     1           		+ ((delz*p * shear(j))/(rad3))! PS offset	
	endif	  


        ray(i,3,j) = rad_to_ang * asin(vel(j) * p) ! pwave angle in degrees	
       
300	continue	! end of number of time layers loop

200	continue	! end of number of angles loop	

c----------------------
c  end of subroutine
c---------------------
9999	continue
	return
	end

