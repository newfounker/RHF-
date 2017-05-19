   

       Subroutine formdensitymat (cmatrix,occ,m,dmatrix)

       implicit double precision (a-h,o-z)
       integer m, i, j, k
       double precision cmatrix(m,m), occ(m,m), 
     + int(m,m), dmatrix(m,m), kmatrix(m,m)

       parameter (alpha1=1.0d0)
       parameter (alpha=1.0d0)
       parameter (beta=0.0d0)
       
       do i=1,m
       do j=1,m
         kmatrix(i,j)=0.0d0
       end do
       end do

       do i=1,0.5*m
       do j=1,0.5*m
         kmatrix(i,j)=cmatrix(i,j)
       end do
       end do 



       call dgemm('N', 'N', m, m, m, alpha1,cmatrix, m,
     &    occ,m,beta,kmatrix,m)

       call dgemm('N', 'T', m, m, m, alpha,kmatrix, m,
     &   cmatrix,m,beta,dmatrix,m)
       
       end


