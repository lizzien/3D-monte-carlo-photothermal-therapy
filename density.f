! get density of grid cell
!
! adapted from code from Kenneth Wood:
! http://www-star.st-and.ac.uk/~kw25/research/montecarlo/points/points.html
! accessed 04.07.2018

      subroutine density(x,y,z,rho,zmax,absCoeff,i,j,k)

      implicit none
      
      include 'grid.txt'

      real*8 x,y,z,rho,zmax

      real*8 w,w2,r,r2,h0,h
      
      real*8 absCoeff(nxg, nyg, nzg)
      
      integer i,j,k

c***** three lines of comments below from astro code
c***** calculate some distances for use in setting up density 
c***** structure. Note that distances are in units of xmax, ymax, and zmax 
c***** as called from the loop over cells in gridset.f

      w2=x*x+y*y
      w=dsqrt(w2)
      r2=w2+z*z
      r=dsqrt(r2)

c************** Tissue Geometry ************************************
      ! Get absorption coefficient from grid
      if ((z.gt.-zmax).and.(z.lt.zmax)) then
        rho=absCoeff(i,j,k)
      else
        rho = 0.d0
      endif
      
c*****************************************************************

      return
      end

