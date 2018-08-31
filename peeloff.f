! Implements "peeling off" procedure
! A fraction of the photon's energy is binned in the direction of the
! observer every time the photon is scattered.
!
! code from Kenneth Wood:
! http://www-star.st-and.ac.uk/~kw25/research/montecarlo/points/points.html
! accessed 04.07.2018

      subroutine peeloff(xp,yp,zp,nxp,nyp,nzp,sint,cost,sinp,cosp,phi,
     +                   fi,fq,fu,fv,
     +                   obsx,obsy,obsz,sinto,costo,sinpo,cospo,phio,
     +                   xmax,ymax,zmax,xface,yface,zface,rhokap,
     +                   rimage,wgt1,albedo,hgg,g2,pl,pc,sc,
     +                   pi,twopi,fourpi,
     +                   fimage,qimage,uimage,
     +                   xcell,ycell,zcell,nscatt,delta)

      implicit none

      include 'grid.txt'
      include 'images.txt'

      integer xcell,ycell,zcell,nscatt
      real xp,yp,zp,nxp,nyp,nzp,fi,fq,fu,fv,obsx,obsy,obsz
      real sinto,costo,sinpo,cospo,phio
      real xmax,ymax,zmax,rimage,wgt1,albedo,pl,pc,sc
      real hgg,g2,fio,fqo,fuo,fvo,sint,cost,sinp,cosp,phi
      real pi,twopi,fourpi,delta

      integer xl,yl
      real phot,photq,photu,hgfac,ximage,yimage,tau2

      call scatt1(fi,fq,fu,fv,fio,fqo,fuo,fvo,nxp,nyp,nzp,
     +            sint,cost,sinp,cosp,phi,
     +            obsx,obsy,obsz,sinto,costo,sinpo,cospo,phio,
     +            hgg,g2,pl,pc,sc,hgfac,
     +            pi,twopi,fourpi)

      call taufind1(xp,yp,zp,obsx,obsy,obsz,xmax,ymax,zmax,
     +                    xface,yface,zface,rhokap,
     +                    tau2,xcell,ycell,zcell,delta)

      if(tau2.eq.0.) goto 666 ! photon is on edge of grid

c      phot=wgt1*hgfac*albedo**nscatt*exp(-tau2)*fio
c      photq=wgt1*hgfac*albedo**nscatt*exp(-tau2)*fqo
c      photu=wgt1*hgfac*albedo**nscatt*exp(-tau2)*fuo

      phot=wgt1*hgfac*albedo*exp(-tau2)*fio
      photq=wgt1*hgfac*albedo*exp(-tau2)*fqo
      photu=wgt1*hgfac*albedo*exp(-tau2)*fuo

c*********** Bin the photon into the image according to its position and 
c*********** direction of travel. 
      yimage=rimage+zp*sinto-yp*costo*sinpo-xp*costo*cospo
      ximage=rimage+yp*cospo-xp*sinpo
      xl=int(nxim*ximage/(2.*rimage))+1
      yl=int(nyim*yimage/(2.*rimage))+1

      if((xl.gt.nxim) .or. (xl.lt.1)) then
         print *,'peeloff: xl out of bounds',xl
         print *,xp,yp,zp,tau2
         goto 666
      endif
      if((yl.gt.nyim) .or. (yl.lt.1)) then
         print *,'peeloff: yl out of bounds',yl
         print *,xp,yp,zp,tau2
         goto 666
      endif

c*********** place weighted photon into image location ********************
c            if(nscatt.eq.1) fimage(xl,yl)=fimage(xl,yl)+phot
c            if(nscatt.eq.1) qimage(xl,yl)=qimage(xl,yl)+photq
c            if(nscatt.eq.1) uimage(xl,yl)=uimage(xl,yl)+photu
            fimage(xl,yl)=fimage(xl,yl)+phot
            qimage(xl,yl)=qimage(xl,yl)+photq
            uimage(xl,yl)=uimage(xl,yl)+photu

666   continue

      end



