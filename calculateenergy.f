            

c          calculating the energy
c          E= 0.5 *[sum(i)sum(j)[p(j,i)*{hmatrix(i,j)+fmatrix(i,j)}]]




       Subroutine calculateenergy (hmatrix,fmatrix,pmatrix,m,energy)


       implicit double precision (a-h,o-z)
       integer m, i, j, k
       double precision hmatrix(m,m), fmatrix(m,m), 
     + pmatrix(m,m)

       do i=1,m
        sum=0.0d0
        do j=1,m
        sum=sum +(0.5*(pmatrix(j,i)*(hmatrix(i,j)+fmatrix(i,j))))
        end do
       energy=energy+sum
       end do
        
       end 

