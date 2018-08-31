! set up arrays of zeros
!
! adapted from code from Kenneth Wood:
! http://www-star.st-and.ac.uk/~kw25/research/montecarlo/points/points.html
! accessed 04.07.2018

      subroutine iarray(xface,yface,zface,rhokap,
     +                  fimage,qimage,uimage,
     +                  jmean,den,absorb,flux)

      include 'grid.txt'
      include 'images.txt'

      integer i,j,k

      do i=1,nxg+1
        xface(i)=0.d0
      end do
      do i=1,nyg+1
        yface(i)=0.d0
      end do
      do i=1,nzg+1
        zface(i)=0.d0
      end do

      do i=1,nxg
          flux(i)=0.d0
          do j=1,nyg
            do k=1,nzg
               rhokap(i,j,k)=0.d0
               jmean(i,j,k)=0.d0
               den(i,j,k)=0.d0
               absorb(i,j,k)=0.d0
            end do
          end do
      end do

      do i=1,nxim
         do j=1,nyim
            fimage(i,j)=0.d0
            qimage(i,j)=0.d0
            uimage(i,j)=0.d0
         end do
      end do



      return
      end
