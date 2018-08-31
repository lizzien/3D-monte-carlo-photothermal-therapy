! Subroutine to read in 3D array of optical properties of tissue

      subroutine struc(absCoeff)
c
      implicit none
      
      include 'grid.txt'

      integer i,j,k,iostat
      real*8 absCoeff(nxg, nyg, nzg)
        
      do i=1,nxg
         do j=1,nyg
            do k=1,nzg
               absCoeff(i,j,k)=0.d0
            end do
          end do
      end do

      ! Open file containing absorption coefficient structure
      open(30, FILE='data/abs_coeff_structure.dat', STATUS='unknown')
      ! Code adapted from:
      ! https://stackoverflow.com/questions/22584283/getting-fortran
      ! -runtime-error-end-of-file?rq=1
      k = nzg
      do while(k >= 1)
          read(30,*, IOSTAT=iostat) (absCoeff(i,1,k),i=1,nxg)
          if ( iostat < 0 ) then
           print*, 'Warning: File containts less than nrow entries'
           exit
          else if( iostat > 0 )then
           write(6,'(A)') 'Error: error reading file'
           stop
          end if
          k = k - 1
      end do
      close(30)

      do j=1, nyg
        absCoeff(:,j,:) = DBLE(absCoeff(:,1,:))
      end do
 
      return
      end