    
             program test
             implicit double precision (a-h,o-z)
             
            double precision, allocatable, dimension(:,:)  :: 
     +      hmatrix, smatrix, 
     +      cmatrix, ctranspose, pmatrix
     
c            double precision, allocatable, dimension(:,:)  ::
c     +       dmatrix

c            double precision dmatrix(1000,1000), 
c    +     feigenmatrix(1000,1000), seigenmatrix(1000,1000),        
c    +     xmatrix(2,2)   

           double precision, allocatable, dimension(:,:)  ::
     +     dmatrix, feigenmatrix, seigenmatrix,
     +     xmatrix                      

           double precision, allocatable, dimension(:,:)  :: 
     +     g3matrix,gmatrix, fmatrix, intmatrix, occ
            

           double precision, allocatable, dimension(:,:)  ::
     +     ovrlpmatrix, trovrlpmatrix 

c          double precision g2matrix(100,100,100,100), energy(50)
           double precision energy(100)  

          double precision, allocatable, dimension(:,:,:,:) ::
     +    g2matrix  

        double precision, allocatable, dimension(:) :: feigenvs,
     +  seigenvs   
      
           integer La(100), Ma(100) 

c        double precision interm(1000,1000), intematrix(100,100) 
         double precision, allocatable, dimension(:,:)  ::
     +   interm, intematrix      
  

         double precision alpha, beta, alpha1, nuclearrep, finalenergy
 
            parameter (alpha=1.0d0)
            parameter (beta=0.0d0)

            parameter (alpha1=2.0d0)
            

c            real, allocatable, dimension(:) :: oneh,ovrlp
          integer n, m, o,i, j, ldim, ldim1, ldim2,k,l, nop
             dimension oneh(2000),ovrlp(2000) 
             dimension buf(600),ibuf(600)
             logical FileExist
             call aces_init_rte
             call aces_ja_init

c            call getrec(4,"JOBARC", "NAOBASFN", 1, m)
            
c             ldim= m*(m+1)          
c             allocate (oneh(ldim))
c             allocate (ovrlp(ldim))   
c            dimension oneh(m),ovrlp(m), buf(600),ibuf(600) 
c             integer n, m, o, i, ldim
c            real k
c             parameter (ldim=m(m+1)/2) 
c            call aces_init_rte
c            call aces_ja_init
 
c             ldim1= (m*(m+1))/2
              
             energyprevious=0.0d0  
            
             call getrec(4,"JOBARC", "NATOMS", 1, n)
             print *," Number of atoms in this calculation", n
          
             call getrec(4,"JOBARC", "COMPNIRR", 1, o)
             print *," Number of Irreps in this calculation", o
            
             call getrec(4,"JOBARC", "NAOBASFN", 1, m)
             print *," Number of Atomic orbitals in basis function", m

             call getrec(4,"JOBARC", "NBASTOT", 1, k)
             print *," Number of Atomic orbitals NBASTOT", k

             call getrec(4,"JOBARC", "NMPROTON", 1, nop)
             print *," Number of protons", nop

             call getrec(4,"JOBARC", "NUCREP", 1, nuclearrep)
             print *," Nuclear repulsion energy", nuclearrep
                 

             ldim1= (m*(m+1))/2            

             l=nop   

             allocate (hmatrix(m,m))
             allocate (smatrix(m,m))

             allocate (cmatrix(m,m)) 
 
             allocate (ctranspose(m,m))   
             allocate (pmatrix(m,m))            

c             allocate (dmatrix(m,m))             


c             allocate (gmatrix(m,m))  
              allocate (g3matrix(ldim1,ldim1))
              allocate (g2matrix(m,m,m,m))


              allocate (dmatrix(m,m))
              allocate (feigenmatrix(m,m))
              allocate (seigenmatrix(m,m))
              allocate (xmatrix(m,m)) 

              allocate (ovrlpmatrix(m,m))
              allocate (trovrlpmatrix(m,m)) 

              allocate (interm(m,m))
              allocate (intematrix(m,m))


c             call getrec(4,"JOBARC", "LDIM", 1, m)
c             print *, m 
c             call get1eint(oneh,ovrlp,buf,ibuf, ldim)
c            call Get2EInt(G2,G3,buf,ibuf,naobasfn,ldim2) 

c             allocate (occ(m,m)) 
              
c             ldim1= (m*(m+1))/2
 
              print *, m    

c             allocate (oneh(ldim))
c             allocate (ovrlp(ldim))
             
            

              call get1eint(oneh,ovrlp,buf,ibuf,ldim1)
              
              
 
              do j=1, ldim1
              print *, oneh(j)
              end do

         
c              do j=1, ldim1
c              if (ovrlp(j).lt.0.0d0) then
c              ovrlp(j)=0.0d0
c              end if
c              end do

      
              
              print *, "The overlap matrix is made from here "
              do j=1,ldim1
              print *, ovrlp(j)
              end do 
              
          
             
              
              print *, "h-matrix" 
              
c               do i=1,m
c               k=0
c               do j=1,i
c               k= ((m*(j-1))-((j-1)*(j-2)/2))+(i-j+1)
c              hmatrix(i,j)= oneh(k)
c               smatrix(i,j)=ovrlp(k)
c               ovrlpmatrix(i,j)=ovrlp(k)
c               enddo
c               enddo


                do i=1,m
                 do j=1,i
                  k=(((i*(i-1))/2)+j)
                  hmatrix(i,j)= oneh(k)
                  smatrix(i,j)=ovrlp(k)
                  ovrlpmatrix(i,j)=ovrlp(k)  
                 end do
                end do    

               do i=1,m
                do j=(i+1),m 
                hmatrix(i,j)=hmatrix(j,i)
                smatrix(i,j)=smatrix(j,i)
                ovrlpmatrix(i,j)=ovrlpmatrix(j,i)
                enddo 
               enddo


c            The 1-e h-matrix is printed below    
              
             call output(hmatrix,1,m,1,m,m,m,2)               
c             do i= 1,m
c             do j=1,m
c             print *, hmatrix(i,j)
c             enddo
c             print *," "
c             enddo
             
 

c            The overlap matrix is printed below 

            print *, "The overlap-matrix elements are "
            call output(smatrix,1,m,1,m,m,m,2)  
            
c           do i= 1,m
c           do j=1,m
c           print *, smatrix(i,j)
c           enddo
c           print *," "
c           enddo
            
c            allocate (dmatrix(m,m))
            allocate (occ(m,m))
c            allocate (dmatrix(m,m)) 

c            Initial guess for Coefficient matrix
c        
             do i=1,m
             do j=1,i-1
             cmatrix(i,j)=0
             enddo
             enddo
 
             do i=1,nop
             cmatrix(i,i)=0
             enddo

             do i=(nop+1),m
             cmatrix(i,i)=0
             enddo

             do i=1,m
             do j=(i+1),m
             cmatrix(i,j)=cmatrix(j,i)
             enddo
             enddo

                          
             print *, " The initial coefficient matrix"  
             call output(cmatrix,1,m,1,m,m,m,2)  

c             do i=1,m 
c             do j=1,m
c            print *, cmatrix(i,j)
c             enddo
c             print *," "
c             enddo 



             
c            Getting the transpose of the coefficient matrix

             call transp(cmatrix,ctranspose,m,m)
             
c            do i=1,m
c            do j=1,m
c            print *, ctranspose(i,j)
c            enddo
c            print *," "
c            enddo

     
c            Calculating the density matrix P
c            P= 2*C(transpose)*C
c          Do we form the density matrix by just using 
c          the coeeficients of occupied orbitals.
c          How do you ensure that ?
c          Do you use an occupation matrix?
c          So my density matrix should be
c             P= 2*[C(transpose)*occ*C]          

c          So lets define an occupation matrix

             do i=1,m
              do j=1,m
              occ(i,j)=0.0d0
              enddo
             enddo

             do i=1,(0.5*nop)
              occ(i,i)=2.0d0
             enddo

           print *, " The occupation matrix"
           call output(occ,1,m,1,m,m,m,2)


c          Now forming the density matrix

           call formdensitymat(cmatrix,occ,m,pmatrix)            

c     &       call xgemm(n,n,m,m,m,1,cmatrix,
c             m,ctranspose,m,0,pmatrix,m)           

c         call dgemm('N', 'N', m, m, m, alpha1,ctranspose, m,
c    &    occ,m,beta,dmatrix,m)

c         call dgemm('N', 'N', m, m, m, alpha,dmatrix, m,
c    &    cmatrix,m,beta,pmatrix,m)
           print *, "The density matrix"  
           call output(pmatrix,1,m,1,m,m,m,2)

c             do i=1,m
c             do j=1,m
c             print *, pmatrix(i,j)
c             enddo
c             print *," "
c             enddo

c           Now we need to get the 2-electron integrals
c           call get2eint(g2matrix,g3matrix,buf,ibuf,m,ldim2)

c           allocating the dimensions for the matrix 
c            allocate (g2matrix(m,m,m,m)) 

             allocate (gmatrix(m,m))
c             allocate (g3matrix(m,m))    
             allocate (fmatrix(m,m)) 

             allocate (intmatrix(m,m))

             allocate (seigenvs(m))
             allocate (feigenvs(m))

             call get2eint1(g2matrix,g3matrix,buf,ibuf,m,ldim1)

c            print *, "The g2-matrix"  
c            call output(g2matrix,1,m,1,m,m,m,2)


             print *, "The g2-matrix"             
             do i=1,m
             do j=1,m
             do k=1,m
              do l=1,m
              print *, i,j,k,l, g2matrix(i,j,k,l)
              end do
             print *," "
             end do
             end do
             end do  

c             print *, "The wanted element" , g2matrix(2,4,1,1)

c            Now trying to multiply the density matrix 
c            with g2matrix to get G-matrix or the two-electron matrix 
c          define 
c          G(i,j)=Sumn(k,l)[pmatrix(k,l)*{
c           g2matrix(i,j,l,k)-0.5*g2matrix(i,k,l,j)}]  

          call getgmatrix(pmatrix,g2matrix,m,gmatrix)
          
c        do i=1,m
c        do j=1,m
c         gmatrix(i,j)=0
c          do k=1,m
c          do l=1,m
c          interm(k,l)=(pmatrix(k,l)*(g2matrix(i,j,l,k)
c    +     -((0.5)*g2matrix(i,k,l,j))))
c          gmatrix(i,j)=gmatrix(i,j)+interm(k,l)
c          end do
c          end do

c         end do
c        end do            

         print *, " The G-matrix"
         call output(gmatrix,1,m,1,m,m,m,2) 
 
c        print *, " The G-matrix"  
c         do i=1,m
c         do j=1,m
c         print *, gmatrix(i,j)
c         enddo
c         print *," "
c         enddo 

              


c           Making the complete fock matrix 
c           fmatrix=hmatrix + gmatrix

            call maddn(m,hmatrix,gmatrix,fmatrix)

c           do i=1,m
c           do j=1,m
c           fmatrix(i,j)= hmatrix(i,j) + gmatrix(i,j)
c           end do
c           end do

c           Printing the fock matrix
            print *, "The fock matrix"
            call output(fmatrix,1,m,1,m,m,m,2)


c              do i=1,m
c              do j=1,m
c              print *, fmatrix(i,j)  
c              end do
c              print *, " "
c              end do
 
c          Getting a transformation matrix by diagonalising 
c          the overlap matrix
           call eig(smatrix, seigenmatrix, 2, m, 0)

c          Just for now
c           seigenmatrix(2,2)=seigenmatrix(2,1)
c           seigenmatrix(2,1)=seigenmatrix(1,1)

c          printing the diagonalised s-matrix elements 
            print *, "The diagonalised s-matrix elements are "
            call output(smatrix,1,m,1,m,m,m,2)

c            do i= 1,m
c            do j=1,m
c            print *, smatrix(i,j)
c            enddo
c            print *," "
c            enddo                 


c          Taking the square root of the diagonalised overlap matrix
           do i=1,m
           smatrix(i,i)=sqrt(smatrix(i,i))
           end do

c          Also arranging the eigen vectors in the matrix form
          

c          do i=1,m
c          k=0
c           do j=1,m
c           k=((m*(i-1))+(j-i+1))
c           seigenmatrix(i,j)=seigenvs(k)
c          print *, seigenmatrix(i,j)
c           end do
c          end do

c          call conv2mat(seigenvs,seigenmatrix,m)

            print *, " The s-eigenvectors "   
            call output(seigenmatrix,1,m,1,m,m,m,2)

c           do i= 1,m
c           do j=1,m
c           print *, seigenmatrix(i,j)
c            enddo
c            print *," "
c            enddo
 



            
          print *, " The square-rooted diagonalised s-matrix "
          call output(smatrix,1,m,1,m,m,m,2) 

c           do i= 1,m
c           do j=1,m
c           print *, smatrix(i,j)
c           enddo
c           print *," "
c           enddo
  
c          Taking the inverse of the square-rooted diagonalised matrix
           call minv(smatrix,m,m,La, det, 1.0 d-8, 0, 2)

         
           print *, " Theinverseofsquareroot of diagonalized matrix "
           call output(smatrix,1,m,1,m,m,m,2) 

c            do i= 1,m
c            do j=1,m
c            print *, smatrix(i,j)
c           enddo
c            print *," "
c            enddo

 
c         Now using the inverted matrix to get the transformation matrix
c         X=seigenmatix*smatrix ; 
c         smatrix-diagonalized-squarerooted-inverse of overlap matrix
          call dgemm('N', 'N', m, m, m, alpha, seigenmatrix, m,
     &    smatrix,m,beta,xmatrix,m)

          print *, "The new transformation matrix"
          call output(xmatrix,1,m,1,m,m,m,2)
c         Just for now
c         xmatrix(1,1)=xmatrix(1,2) 
c         xmatrix(1,2)=xmatrix(2,1)
c         xmatrix(2,1)=xmatrix(1,1)
c         xmatrix(2,2)=-xmatrix(1,2) 

          call output(xmatrix,1,m,1,m,m,m,2) 

c         Checking if the transformation matrix really orthonormalise
c          the s-matrix 

          call dgemm('T', 'N', m, m, m, alpha, xmatrix, m,
     &    ovrlpmatrix,m,beta,intematrix,m)

          call dgemm('N', 'N', m, m, m, alpha, intematrix, m,
     &    xmatrix,m,beta,trovrlpmatrix,m)

          print *, "The transformed overlap matrix"
          call output(trovrlpmatrix,1,m,1,m,m,m,2)
          
            
c            do i= 1,m
c              do j=1,m
c               print *, xmatrix(i,1)
c              end do
c            print *," "
c            end do

           

c          transform the fock matrix - 
c          fmatrix= xmatrix(tr)*fmatrix*xmatrix   

              
          call dgemm('T', 'N', m, m, m, alpha, xmatrix, m,
     &    fmatrix,m,beta,intmatrix,m)  

          
            do i= 1,m
            do j=1,m
            print *, intmatrix(i,j)
            enddo
            print *," "
            enddo
           
           call dgemm('N', 'N', m, m, m, alpha, intmatrix, m,
     &    xmatrix,m,beta,fmatrix,m)

           print *, " The transformed fock-matrix "
           call output(fmatrix,1,m,1,m,m,m,2)

c            do i= 1,m
c            do j=1,m
c            print *, fmatrix(i,j)
c            enddo
c            print *," "
c            enddo

c          Now diagonalising the transformed fock matrix
          call eig(fmatrix,feigenmatrix,2, m,0)
              
           print *, " The diagonalised fock matrix " 
           call output(fmatrix,1,m,1,m,m,m,2)

c          do i= 1,m
c          do j=1,m
c          print *, fmatrix(i,j)
c          enddo
c           print *," "
c           enddo



c          Printing the eigenvalues and the eigenvectors
           
c           do j=1,2*m
c           print *, feigenvs(j)
c           enddo  
           
           do k=1,50 

c          first we need to arrange the eigenvectors which
c          are in the form of column to a matrix form


c          call conv2mat(feigenvs,feigenmatrix,m)

c           do i=1,m
c           k=0 
c           do j=1,m
c           k=((m*(i-1))+(j-i+1))
c            feigenmatrix(i,j)=feigenvs(k)
c            print *, feigenmatrix(i,j)
c           end do
c           end do

            print *, "feigenmatrix"
            call output(feigenmatrix,1,m,1,m,m,m,2)
                      
c            do i= 1,m
c            do j=1,m
c            print *, feigenmatrix(2,1)
c            end do
c            print *," "
c            end do



                  
c          Getting the new coefficient matrix
c          C=X(transformation matrix) * feigenmatrix(neweigenvalues)
          call dgemm('N', 'N', m, m, m, alpha,xmatrix, m,
     &    feigenmatrix,m,beta,cmatrix,m) 

           print *, "The new coefficient matrix"
           call output(cmatrix,1,m,1,m,m,m,2)
      
c          do i= 1,m
c          do j=1,m
c         print *, cmatrix(i,j)
c           enddo
c           print *," "
c           enddo

c          Calculate the new density matrix 
           call formdensitymat(cmatrix,occ,m,pmatrix)
           print *, "The new density matrix"
           call output(pmatrix,1,m,1,m,m,m,2)

c          form new gmatrix
           call getgmatrix(pmatrix,g2matrix,m,gmatrix)
           print *, "The new G-matrix"
           call output(gmatrix,1,m,1,m,m,m,2) 

c          form new fock matrix
           call maddn(m,hmatrix,gmatrix,fmatrix) 
           print *, "The new fock matrix"
           call output(fmatrix,1,m,1,m,m,m,2) 

c          calculating the energy
c          E= 0.5 *[sum(i)sum(j)[p(j,i)*{hmatrix(i,j)+fmatrix(i,j)}]]
c          call calculateenergy (hmatrix,fmatrix,pmatrix,m,energy)

           call calculateenergy (hmatrix,fmatrix,pmatrix,m,energy(k))   
           print *, "The energy is", energy(k)
          
            
            
           if (k.ge.2) then 
           diff=energy(k)-energy(k-1)
           end if

c          if (diff.ge.1.0d-5) then do
           
c           else do 
c           break
c           end if
           
c          Want to get this in iteration form
c          Again transforming the fock matrix 

           call dgemm('T', 'N', m, m, m, alpha, xmatrix, m,
     &    fmatrix,m,beta,intmatrix,m)
           
          call dgemm('N', 'N', m, m, m, alpha, intmatrix, m,
     &    xmatrix,m,beta,fmatrix,m)

c          Diagonalising the transformed fock matrix
 
           call eig(fmatrix,feigenmatrix,2, m,1)

c          Printing the eigenvalues
           print *, "The new eigenvalues"
           call output(fmatrix,1,m,1,m,m,m,2)


           end do

c          converting eigenvectors to matrix form
c          call conv2mat(feigenvs,feigenmatrix,m)  

c          Getting the new coefficient matrix
c          C=X(transformation matrix) * feigenmatrix(neweigenvalues)
          
c          call dgemm('N', 'N', m, m, m, alpha, xmatrix, m,
c    &     feigenmatrix,m,beta,cmatrix,m)
   




           do k=1,100
           print *, "The electronic energy for iteration", k,
     +     "is", energy (k)
           end do

           finalenergy=energy(50)+nuclearrep
           print *, "The total energy including nuclear repulsion",
     +     finalenergy, "Hartrees"

           print *, "The final orbital energies are "
           
           call output(fmatrix,1,m,1,m,m,m,2)
c          do i= 1,m
c          do j=1,m
c          print *, fmatrix(i,j)
c          enddo
c          print *," "
c          enddo

c          printing the final eigenvectors
           print *, "The final eigen vectors are "

           call output(cmatrix,1,m,1,m,m,m,2)           


              print *, "thats it"

              call aces_ja_fin
              call aces_fin
             

              end  



