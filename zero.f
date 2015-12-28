	subroutine zero(wave, nopts)	
c-------------------------------------------------------------------
c   Routine zeros array wave that has nopts samples 
c  
c-------------------------------------------------------------------
c  Hilterman September 2004
c-------------------------------------------------------------------
	integer  nopts, i
	real*4 wave(nopts) 
	
	do 100 i=1,nopts
100	wave(i) = 0.

	return
	end
	
  
	
