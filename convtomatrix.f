

         Subroutine conv2mat(a,b,m)
         
         
         integer m, i, j, k
         double precision a(2*m), b(m,m) 

         do i=1,m
          k=0
           do j=1,m
           k=((m*(j-1))+i)
           b(i,j)=a(k)
           end do
          end do

          end 
