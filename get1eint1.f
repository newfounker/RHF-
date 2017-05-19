C*************************************************************
      subroutine Get1EInt1(oneh,ovrlp,buf,ibuf,ldim)
C     
C     Get one electron integrals
C         
C*************************************************************
      implicit double precision (a-h,o-z)
c
      logical FileExist
c
      dimension oneh(ldim),ovrlp(ldim),
     &          buf(600),ibuf(600)
c
      ilnbuf=600
      FileExist=.false.
      inquire(file='IIII',exist=FileExist)
      if (FileExist) then
c         write(*,*) 'IIII file exists'
         open(unit=10,file='IIII',form='UNFORMATTED',
     &        access='SEQUENTIAL')
         rewind 10
         call locate(10,'ONEHAMIL')
         call zero(oneh,ldim)
         nut = ilnbuf
         do while (nut.eq.ilnbuf)
            read(10) buf, ibuf, nut
            do int = 1, nut
               oneh(ibuf(int)) = buf(int)
            end do
         end do
c
         call locate(10,'OVERLAP ')
         call zero(ovrlp,ldim)
         nut = ilnbuf
         do while (nut.eq.ilnbuf)
            read(10) buf, ibuf, nut
            do int = 1, ilnbuf
               ovrlp(ibuf(int)) = buf(int)
            end do
         end do
         close(10)
      else
         write(*,*) 'IIII file does not exist'
         stop
      end if
      return
      end
