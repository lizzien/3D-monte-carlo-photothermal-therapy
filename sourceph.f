! subroutine to emit photon
! sets direction and source location of photon
!
! adapted from code from Kenneth Wood:
! http://www-star.st-and.ac.uk/~kw25/research/montecarlo/points/points.html
! accessed 04.07.2018

      subroutine sourceph(xs,ys,zs,xp,yp,zp,nxp,nyp,nzp,
     +                    sint,cost,sinp,cosp,phi,fi,fq,fu,fv,
     +                    xmax,ymax,zmax,twopi,
     +                    xcell,ycell,zcell,nxg,nyg,nzg,iseed,
     +                    b,delta)

      implicit none

      include 'photon.txt'

      integer xcell,ycell,zcell,nxg,nyg,nzg,iseed
      real*8 xs,ys,zs,xmax,ymax,zmax,twopi
      real*8 ran2
      real*8 b,rphi,rr,delta

c***** emit photon isotropically from point source location ************

      ! Gaussian light source
      rphi = twopi * ran2(iseed)
      rr = b * dsqrt(-dlog(ran2(iseed)))
      xp = xs + (rr * dcos(rphi))
      yp = ys + (rr * dsin(rphi))
      zp = zs

      ! commented when two lines below were added
      !cost=2.d0*ran2(iseed)-1.d0
      !sint=(1.d0-cost*cost)
      !if(sint.le.0.d0)then
      !  sint=0.d0
      !else
      !  sint=dsqrt(sint)
      !endif
      
      ! Added to make beam parallel
      cost=dcos(twopi/2.d0)
      sint=dsin(twopi/2.d0)

      phi=twopi*ran2(iseed)
      cosp=dcos(phi)
      sinp=dsin(phi)

c***** Set photon direction cosines for direction of travel *********
      nxp=sint*cosp  
      nyp=sint*sinp
      nzp=cost

c***** Set Stokes fluxes ********************************************
      fi=1.d0
      fq=0.d0
      fu=0.d0
      fv=0.d0

c*************** Linear Grid *************************
      xcell=int(DBLE(nxg)*(xp+xmax)/(2.d0*xmax))+1
      ycell=int(DBLE(nyg)*(yp+ymax)/(2.d0*ymax))+1
      zcell=int(DBLE(nzg)*(zp+zmax)/(2.d0*zmax))+1
c*****************************************************

      return
      end

