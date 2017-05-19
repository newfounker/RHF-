          
       Subroutine maddn(m, a, b, c)

       integer m, i, j, k
       double precision a(m,m), b(m,m), c(m,m)

c      Matrix addition
       do i=1,m
       do j=1,m
       c(i,j)=a(i,j) + b(i,j)
       end do
       end do
       

       end   


