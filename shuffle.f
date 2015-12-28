      subroutine Shuffle(a,Nsel)
      integer, intent(inout) :: a(Nsel)
      integer :: i, randpos, temp
      real :: r
      
      do i = Nsel, 2, -1
         call random_number(r)
         randpos = int(r * i) + 1
         temp = a(randpos)
         a(randpos) = a(i)
         a(i) = temp
      end do
 
      end subroutine Shuffle