! sets source location
!
! adapted from code from Kenneth Wood:
! http://www-star.st-and.ac.uk/~kw25/research/montecarlo/points/points.html
! accessed 04.07.2018

      subroutine sources(xsource,ysource,zsource,lsource,lumtot,iseed,
     +                  lum)

      implicit none

      include 'sources.txt'

      integer i,iseed
      real*8 lum

c**** Set photon locations and luminosities
      ! Gaussian beam centre at centre of top surface
      xsource(1) = 0.d0
      ysource(1) = 0.d0
      zsource(1) = 0.999d0
      lsource(1) = lum

c**** Calculate total luminosity of all sources
      lumtot=0.d0
      do i=1,nsource
        lumtot=lumtot+lsource(i)
      end do

      return

      end

