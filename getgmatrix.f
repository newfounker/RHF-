
c          Get G-matrix or the two-electron matrix
c          define
c          G(i,j)=Sumn(k,l)[pmatrix(k,l)*{
c           g2matrix(i,j,l,k)-0.5*g2matrix(i,k,l,j)}]

       Subroutine getgmatrix(pmatrix,glinmatrix,m,gmatrix) 

       implicit double precision (a-h,o-z)
       integer m, i, j, k
       double precision pmatrix(m,m), glinmatrix(m,m,m,m),
     + interm(m,m), gmatrix(m,m),sum

       

       do i=1,m
        do j=1,m
         sum=0
          do k=1,m
          do l=1,m
          interm(k,l)=(pmatrix(k,l))*(glinmatrix(i,j,l,k)
     +  -((0.5)*glinmatrix(i,k,l,j)))
          sum=sum + interm(k,l)
          end do
          end do
          gmatrix(i,j)= sum
        end do
       end do          










c       do i=1,m
c       do j=1,m
c        gmatrix(i,j)=0
c         do k=1,m
c         do l=1,m
c         interm(k,l)=(pmatrix(k,l)*(glinmatrix(i,j,l,k)
c     +   -((0.5)*glinmatrix(i,k,l,j))))
c         gmatrix(i,j)=gmatrix(i,j)+interm(k,l)
c         end do
c         end do
        
        
 
c      end do
c       end do


       end 
