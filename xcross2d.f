!23456   
       SUBROUTINE xcross2D(array1,array2,
     1            n,m,r)
       INTEGER n,m
       REAL*4 array1(n,m),array2(n,m),tempy(m,n),numerator(m)
       REAL*4 tempx(m,n),num,dem1,dem2,pax(m),pay(m),temp1(n,1)
       REAL*4 u(1),uu(1),uuu(1),temp2(n,1),r
       
       tempy = TRANSPOSE(array2)
       tempx = TRANSPOSE(array1)
       

       DO 1 i=1,m
            temp1(:,1) = array1(:,i)
            temp2(:,1) = array2(:,i)
      		u   =  MATMUL(tempy(i,:),temp1)
      		uu  =  MATMUL(tempy(i,:),temp2)
      		uuu =  MATMUL(tempx(i,:),temp1)
      		numerator(i) = u(1)
       		pax(i)       = uu(1)
       		pay(i)       = uuu(1)
1      ENDDO
       

       num  = SUM(numerator)
       dem1 = SUM(pax)
       dem2 = SUM(pay)
       
       r = (2*num)/(dem1 + dem2) 			
       RETURN
       END
       