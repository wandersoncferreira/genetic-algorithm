       SUBROUTINE sort(n,arr,indx)
       INTEGER n,indx(n),M,NSTACK
       REAL arr(n)
       PARAMETER (M=7,NSTACK=50)
       INTEGER i,indxt,ir,itemp,j,jstack,k,l,istack(NSTACK)
       REAL a
       
       DO j=1,n
          indx(j)=j
       ENDDO
       jstack=0
       l=1
       ir=n
       
!2345678       
1      IF(ir-l.lt.M)then
          DO j=l+1,ir
              indxt=indx(j)
              a=arr(indxt)
              DO i=j-1,1,-1
                  IF(arr(indx(i)).le.a)GOTO 2
                  indx(i+1)=indx(i)
              ENDDO 
              i=l-1
2             indx(i+1)=indxt
          ENDDO
          IF(jstack.eq.0)RETURN
          ir=istack(jstack)
          l=istack(jstack-1)
          jstack=jstack-2
!2345678          
       ELSE
          k=(l+ir)/2
          itemp=indx(k)
          indx(k)=indx(l+1)
          indx(l+1)=itemp
          IF(arr(indx(l)).gt.arr(indx(ir)))then
              itemp=indx(l)
              indx(l)=indx(ir)
              indx(ir)=itemp
          ENDIF
          IF(arr(indx(l+1)).gt.arr(indx(ir)))then
              itemp=indx(l+1)
              indx(l+1)=indx(ir)
              indx(ir)=itemp
          ENDIF
          IF(arr(indx(l)).gt.arr(indx(l+1)))then
              itemp=indx(l)
              indx(l)=indx(l+1)
              indx(l+1)=itemp
          ENDIF
          i=l+1
          j=ir
          indxt=indx(l+1)
          a=arr(indxt)
3         CONTINUE
              i=i+1
          IF(arr(indx(i)).lt.a)GOTO 3
4         CONTINUE
              j=j-1
          IF(arr(indx(j)).gt.a)GOTO 4
          IF(j.lt.i)GOTO 5
          itemp=indx(i)
          indx(i)=indx(j)
          indx(j)=itemp
          GOTO 3
5         indx(l+1)=indx(j)
          indx(j)=indxt
          jstack=jstack+2
          IF(jstack.gt.NSTACK) then
               stop
          ENDIF     
          IF(ir-i+1.ge.j-1)then
              istack(jstack)=ir
              istack(jstack-1)=i
              ir=j-1
          ELSE
              istack(jstack)=j-1
              istack(jstack-1)=l
              l=i
          ENDIF
!2345678                                           
       ENDIF
       GOTO 1
       END                 