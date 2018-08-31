! generates random optical depth 
! propagates photon along optical depth
! until photon reaches optical depth or leaves grid.
!
! adapted from code from Kenneth Wood:
! http://www-star.st-and.ac.uk/~kw25/research/montecarlo/points/points.html
! accessed 04.07.2018

      subroutine tauint2(xp,yp,zp,nxp,nyp,nzp,xmax,ymax,zmax,
     +                  xface,yface,zface,rhokap,
     +                  xcell,ycell,zcell,tflag,iseed,delta,
     +                  jmean,absorb,flux)

      implicit none

      include 'grid.txt'

      integer tflag,iseed,xcell,ycell,zcell
      real*8 xp,yp,zp,nxp,nyp,nzp,xmax,ymax,zmax
      real*8 ran2

      integer celli,cellj,cellk
      integer cellr,zflag
      real*8 tau,taurun,taucell,d,d1,dcell,xcur,ycur,zcur,dsx,dsy,dsz
      real*8 dx,dy,dz,smax,delta

c***** tflag=0 means photon is in envelope
      tflag=0
      zflag=0

c**** generate random optical depth tau
      tau=-dlog(ran2(iseed))

c***** set the cumulative distance and optical depth (d and taurun) 
c***** along the photon path to zero.  set the current photon coordinates.
c***** note that the origin of the (xcur,ycur,zcur) system is at the 
c***** bottom corner of the grid.
      taurun=0.d0
      d=0.d0
      xcur=xp+xmax
      ycur=yp+ymax
      zcur=zp+zmax

      celli=xcell
      cellj=ycell
      cellk=zcell

c***** calculate smax -- maximum distance photon can travel
      if(nxp.gt.0.d0) then
         dsx=(2.d0*xmax-xcur)/nxp
      elseif(nxp.lt.0.d0) then
         dsx=-xcur/nxp
      elseif(nxp.eq.0.d0) then
         dsx=1.d2*xmax
      endif

      if(nyp.gt.0.d0) then
         dsy=(2.d0*ymax-ycur)/nyp
      elseif(nyp.lt.0.d0) then
         dsy=-ycur/nyp
      elseif(nyp.eq.0.d0) then
         dsy=1.d2*ymax
      endif

      if(nzp.gt.0.d0) then
         dsz=(2.d0*zmax-zcur)/nzp
      elseif(nzp.lt.0.d0) then
         dsz=-zcur/nzp
      elseif(nzp.eq.0.d0) then
         dsz=1.d2*zmax
      endif

      smax=dmin1(dsx,dsy,dsz)
      if(smax.lt.delta) then
         tflag=1
         return
      endif
       
c***** integrate through grid
      dowhile((taurun.lt.tau).and.(d.lt.(.999d0*smax)))

c***** find distance to next x, y, and z cell walls.
c***** note that dx is not the x-distance, but the actual distance along 
c*****the direction of travel to the next x-face, and likewise for dy and dz.
         if(nxp.gt.0.d0) then
            dx=(xface(celli+1)-xcur)/nxp
            if(dx.lt.delta) then
               xcur=xface(celli+1)
               celli=celli+1
               dx=(xface(celli+1)-xcur)/nxp
            endif
         elseif(nxp.lt.0.d0) then
            dx=(xface(celli)-xcur)/nxp
            if(dx.lt.delta) then
               xcur=xface(celli)
               dx=(xface(celli-1)-xcur)/nxp
               celli=celli-1
            endif
         elseif(nxp.eq.0.d0) then
            dx=1.d2*xmax
         endif

         if(nyp.gt.0.d0) then
            dy=(yface(cellj+1)-ycur)/nyp
            if(dy.lt.delta) then
               ycur=yface(cellj+1)
               cellj=cellj+1
               dy=(yface(cellj+1)-ycur)/nyp
            endif
         elseif(nyp.lt.0.d0) then
            dy=(yface(cellj)-ycur)/nyp
            if(dy.lt.delta) then
               ycur=yface(cellj)
               dy=(yface(cellj-1)-ycur)/nyp
               cellj=cellj-1
            endif
         elseif(nyp.eq.0.d0) then
            dy=1.d2*ymax
         endif

         if(nzp.gt.0.d0) then
            dz=(zface(cellk+1)-zcur)/nzp
            if(dz.lt.delta) then
               zcur=zface(cellk+1)
               cellk=cellk+1
               dz=(zface(cellk+1)-zcur)/nzp
            endif
         elseif(nzp.lt.0.d0) then
            dz=(zface(cellk)-zcur)/nzp
            if(dz.lt.delta) then
               zcur=zface(cellk)
               dz=(zface(cellk-1)-zcur)/nzp
               cellk=cellk-1
            endif
         elseif(nzp.eq.0.d0) then
            dz=1.d2*zmax
         endif

c***** distances are only zero if photon is on cell wall.  if it is 
c***** on cell wall then set to arbitrary large distance, since we will
c***** in fact hit another wall
         if( (dx.eq.0.d0) .or. ((dabs(dx)).lt.(delta)) ) dx=1.d2*xmax
         if( (dy.eq.0.d0) .or. ((dabs(dy)).lt.(delta)) ) dy=1.d2*ymax
         if( (dz.eq.0.d0) .or. ((dabs(dz)).lt.(delta)) ) dz=1.d2*zmax

c***** find distance to next cell wall -- minimum of dx, dy, and dz
         dcell=dmin1(dx,dy,dz)
         if(dcell.le.0.d0) then
            print *,'tauint2: dcell < 0'
         endif
         if(dx.lt.0.d0) dcell=dmin1(dy,dz)
         if(dy.lt.0.d0) dcell=dmin1(dx,dz)
         if(dz.lt.0.d0) dcell=dmin1(dx,dy)

c***** optical depth to next cell wall is 
c***** taucell= (distance to cell)*(opacity of current cell)
         taucell=dcell*rhokap(celli,cellj,cellk)

c***** if taurun+taucell>tau then scatter at distance d+d1.  
c***** update photon position and cell.  
c***** if taurun+taucell<tau then photon moves distance dcell 
c***** (i.e. ends up on next cell wall) and update photon position
c***** and cell.
         if((taurun+taucell).ge.tau) then
            d1=(tau-taurun)/rhokap(celli,cellj,cellk)
            d=d+d1
            taurun=taurun+taucell
            xcur=xcur+d1*nxp
            ycur=ycur+d1*nyp
            zcur=zcur+d1*nzp
            
            jmean(celli,cellj,cellk) = jmean(celli,cellj,cellk) + d1

c*************** Linear Grid ************************
            celli=int(DBLE(nxg)*xcur/(2.d0*xmax))+1
            cellj=int(DBLE(nyg)*ycur/(2.d0*ymax))+1
            cellk=int(DBLE(nzg)*zcur/(2.d0*zmax))+1
c****************************************************

         else
            d=d+dcell
            taurun=taurun+taucell
            xcur=xcur+dcell*nxp
            ycur=ycur+dcell*nyp
            zcur=zcur+dcell*nzp
            
            jmean(celli,cellj,cellk) = jmean(celli,cellj,cellk) + dcell

c*************** Linear Grid ************************
            celli=int(DBLE(nxg)*xcur/(2.d0*xmax))+1
            cellj=int(DBLE(nyg)*ycur/(2.d0*ymax))+1
            cellk=int(DBLE(nzg)*zcur/(2.d0*zmax))+1
c****************************************************

          endif

      end do

c***** calculate photon final position.  if it escapes envelope then
c***** set tflag=1.  if photon doesn't escape leave tflag=0 and update 
c***** photon position.
      if((d.ge.(.999d0*smax))) then
          tflag=1
          !if (cellk == 202) then
          if (zcur > 2.d0*zmax-delta) then
              xcur = xcur - xmax
              ycur = ycur - ymax
              if (dsqrt(xcur**2 + ycur**2) < xmax) then
                  cellr=int(DBLE(nxg)*dsqrt(xcur**2+ycur**2)/xmax)+1
                  flux(cellr) = flux(cellr) + 1.d0
              endif
          endif
      else
         xp=xp+d*nxp
         yp=yp+d*nyp
         zp=zp+d*nzp
         xcell=int(DBLE(nxg)*(xp+xmax)/(2.d0*xmax))+1
         ycell=int(DBLE(nyg)*(yp+ymax)/(2.d0*ymax))+1
         zcell=int(DBLE(nzg)*(zp+zmax)/(2.d0*zmax))+1
         !absorb(xcell,ycell,zcell) = absorb(xcell,ycell,zcell) + 1.d0
      endif

      return
      end
