! adapted from code from Kenneth Wood:
! http://www-star.st-and.ac.uk/~kw25/research/montecarlo/points/points.html
! accessed 04.07.2018

      subroutine tauint(xp,yp,zp,nxp,nyp,nzp,xmax,ymax,zmax,
     +                  xface,yface,zface,rhokap,
     +                  tau,xcell,ycell,zcell,delta)

      implicit none

      include 'grid.txt'

      integer xcell,ycell,zcell
      real xp,yp,zp,nxp,nyp,nzp,xmax,ymax,zmax,tau

      real xcur,ycur,zcur,taurun,taucell,d,d1,dsx,dsy,dsz,dx,dy,dz
      real smax,dcell,delta


c***** set the cumulative distance and optical depth (d and taurun) 
c***** along the photon path to zero.  set the current photon coordinates.
c***** note that the origin of the (xcur,ycur,zcur) system is at the 
c***** bottom corner of the grid.
      taurun=0.
      taucell=0.
      d=0.
      xcur=xp+xmax
      ycur=yp+ymax
      zcur=zp+zmax

c***** calculate smax -- maximum distance photon can travel
      if(nxp.gt.0.) then
         dsx=(2.*xmax-xcur)/nxp 
      elseif(nxp.lt.0.) then
         dsx=-xcur/nxp
      elseif(nxp.eq.0.) then
         dsx=1.e2*xmax
      endif

      if(nyp.gt.0.) then
         dsy=(2.*ymax-ycur)/nyp
      elseif(nyp.lt.0.) then
         dsy=-ycur/nyp
      elseif(nyp.eq.0.) then
         dsy=1.e2*ymax
      endif

      if(nzp.gt.0.) then
         dsz=(2.*zmax-zcur)/nzp
      elseif(nzp.lt.0.) then
         dsz=-zcur/nzp
      elseif(nzp.eq.0.) then
         dsz=1.e2*zmax
      endif

      smax=amin1(dsx,dsy,dsz)
       
c***** integrate through grid
      dowhile((taurun.lt.tau).and.(d.lt.(.999*smax)))

c***** find distance to next x, y, and z cell walls.  
c***** note that dx is not the x-distance, but the actual distance along 
c*****the direction of travel to the next x-face, and likewise for dy and dz.
         if(nxp.gt.0.) then
            dx=(xface(xcell+1)-xcur)/nxp
            if(dx.lt.delta) then
               xcur=xface(xcell+1)
               xcell=xcell+1
               dx=(xface(xcell+1)-xcur)/nxp
            endif
         elseif(nxp.lt.0.) then
            dx=(xface(xcell)-xcur)/nxp
            if(dx.lt.delta) then
               xcur=xface(xcell)
               dx=(xface(xcell-1)-xcur)/nxp
               xcell=xcell-1
            endif
         elseif(nxp.eq.0.) then
            dx=1.e2*xmax
         endif

         if(nyp.gt.0.) then
            dy=(yface(ycell+1)-ycur)/nyp
            if(dy.lt.delta) then
               ycur=yface(ycell+1)
               ycell=ycell+1
               dy=(yface(ycell+1)-ycur)/nyp
            endif
         elseif(nyp.lt.0.) then
            dy=(yface(ycell)-ycur)/nyp
            if(dy.lt.delta) then
               ycur=yface(ycell)
               dy=(yface(ycell-1)-ycur)/nyp
               ycell=ycell-1
            endif
         elseif(nyp.eq.0.) then
            dy=1.e2*ymax
         endif

         if(nzp.gt.0.) then
            dz=(zface(zcell+1)-zcur)/nzp
            if(dz.lt.delta) then
               zcur=zface(zcell+1)
               zcell=zcell+1
               dz=(zface(zcell+1)-zcur)/nzp
            endif
         elseif(nzp.lt.0.) then
            dz=(zface(zcell)-zcur)/nzp
            if(dz.lt.delta) then
               zcur=zface(zcell)
               dz=(zface(zcell-1)-zcur)/nzp
               zcell=zcell-1
            endif
         elseif(nzp.eq.0.) then
            dz=1.e2*zmax
         endif

c***** distances are only zero if photon is on cell wall.  if it is 
c***** on cell wall then set to arbitrary large distance, since we will
c***** in fact hit another wall
         if( (dx.eq.0.) .or. ((abs(dx)).lt.(delta)) ) dx=1.e2*xmax
         if( (dy.eq.0.) .or. ((abs(dy)).lt.(delta)) ) dy=1.e2*ymax
         if( (dz.eq.0.) .or. ((abs(dz)).lt.(delta)) ) dz=1.e2*zmax

c***** find distance to next cell wall -- minimum of dx, dy, and dz
         dcell=amin1(dx,dy,dz)
         if(dcell.le.0.) then
            print *,'tauint: dcell < 0'
c            stop
         endif
         if(dx.lt.0.) dcell=amin1(dy,dz)
         if(dy.lt.0.) dcell=amin1(dx,dz)
         if(dz.lt.0.) dcell=amin1(dx,dy)

c***** optical depth to next cell wall is 
c***** taucell= (distance to cell)*(opacity of current cell)
         taucell=dcell*rhokap(xcell,ycell,zcell)

c***** if taurun+taucell>tau then scatter at distance d+d1.  
c***** update photon position and cell.  
c***** if taurun+taucell<tau then photon moves distance dcell 
c***** (i.e. ends up on next cell wall) and update photon position
c***** and cell.
         if((taurun+taucell).ge.tau) then
            d1=(tau-taurun)/rhokap(xcell,ycell,zcell)
            d=d+d1
            taurun=taurun+taucell
            xcur=xcur+d1*nxp
            ycur=ycur+d1*nyp
            zcur=zcur+d1*nzp

c           goto 100

c*************** Linear Grid ************************
            xcell=int(nxg*xcur/(2.*xmax))+1
            ycell=int(nyg*ycur/(2.*ymax))+1
            zcell=int(nzg*zcur/(2.*zmax))+1
c****************************************************

         else
            d=d+dcell
            taurun=taurun+taucell
            xcur=xcur+dcell*nxp
            ycur=ycur+dcell*nyp
            zcur=zcur+dcell*nzp

c*************** Linear Grid ************************
            xcell=int(nxg*xcur/(2.*xmax))+1
            ycell=int(nyg*ycur/(2.*ymax))+1
            zcell=int(nzg*zcur/(2.*zmax))+1
c****************************************************

         endif

      end do

      xp=xcur-xmax
      yp=ycur-ymax
      zp=zcur-zmax
      xcell=int(nxg*(xp+xmax)/(2.*xmax))+1
      ycell=int(nyg*(yp+ymax)/(2.*ymax))+1
      zcell=int(nzg*(zp+zmax)/(2.*zmax))+1

      return
      end
