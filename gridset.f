! sets up grid faces and density grid
! calculates optical depths to edge of grid
!
! adapted from code from Kenneth Wood:
! http://www-star.st-and.ac.uk/~kw25/research/montecarlo/points/points.html
! accessed 04.07.2018

      subroutine gridset(xface,yface,zface,rhokap,xmax,ymax,zmax,
     +                  kappa,absCoeff)

      implicit none

      include 'grid.txt'

      real*8 xmax,ymax,zmax,kappa
      real*8 absCoeff(nxg, nyg, nzg)

      integer i,j,k
      real*8 x,y,z,rho,taueq,taupole
      !real rhosum,dV

      print *, 'Setting up density grid....'

c*********************************************************************
c**********  Linear Cartesian grid. Set up grid faces ****************
c*********************************************************************
      do i=1,nxg+1
         xface(i)=DBLE(i-1)*2.d0*xmax/DBLE(nxg)
      end do
      do i=1,nyg+1
         yface(i)=DBLE(i-1)*2.d0*ymax/DBLE(nyg)
      end do
      do i=1,nzg+1
         zface(i)=DBLE(i-1)*2.d0*zmax/DBLE(nzg)
      end do
c**********************************************************************

      !rhosum=0.
c**************  Loop through x, y, and z to set up grid density.  ****
      do i=1,nxg
       do j=1,nyg
        do k=1,nzg
           x=xface(i)-xmax+xmax/DBLE(nxg)
           y=yface(j)-ymax+ymax/DBLE(nyg)
           z=zface(k)-zmax+zmax/DBLE(nzg)

c**********************Call density setup subroutine *****************
           call density(x,y,z,rho,zmax,absCoeff,i,j,k)
           rhokap(i,j,k)=rho*kappa!*1.5e13 ! rho*kappa*R, R=1AU=1.5e13cm    !me!

           !rhosum=rhosum+rhokap(i,j,k)/(kappa*1.5e13)
c*********************************************************************

        end do
       end do
      end do
c*********************************************************************

      !dV=2.*xmax/nxg
      !dV=dV*2.*ymax/nyg
      !dV=dV*2.*zmax/nzg
      !dV=dV*1.5e13
      !dV=dV/2.e33  ! dV/M_sun
      !dV=dV*1.5e13
      !dV=dV*1.5e13

      !print *,'Mass (M_sun) = ',rhosum*dV

c****************** Calculate equatorial and polar optical depths ****
      taueq=0.d0
      taupole=0.d0
      do i=1,nxg
         taueq=taueq+rhokap(i,nyg/2,nzg/2)
      enddo
      do i=1,nzg
         taupole=taupole+rhokap(nxg/2,nyg/2,i)
      enddo
      taueq=taueq*2.d0*xmax/DBLE(nxg)
      taupole=taupole*2.d0*zmax/DBLE(nzg)
      print *,'taueq = ',taueq,'  taupole = ',taupole

c************** Write out cut through density grid *******************
c      open(unit=10,file='rhoxz.dat',status='unknown')
c        do i=1,nyg
c           write(10,*) (rhokap(k,nyg/2+1,i),k=1,nxg)
c        end do
c      close(10)

      return
      end
