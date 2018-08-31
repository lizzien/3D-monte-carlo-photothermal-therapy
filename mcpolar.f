! main program for implementing Monte Carlo 3D grid code
! simulates light transport through tissue
!
! adapted from code from Kenneth Wood:
! http://www-star.st-and.ac.uk/~kw25/research/montecarlo/points/points.html
! accessed 04.07.2018

      program mcpolar

      implicit none

      include 'grid.txt'
      include 'sources.txt'
      include 'photon.txt'
      include 'images.txt'

c***** Parameter declarations ****************************************
      integer nphotons,nph,iseed,is,j,jcount,nscatt,totscatt,xl,yl,i
      integer xcell,ycell,zcell,tflag
      integer k
      real*8 lum,b,dVi,jsum,totscattR,nphotonsR
      real*8 kappa,albedo,hgg,pl,pc,sc,xmax,ymax,zmax,rimage
      real*8 viewthet,viewphi,obsx,obsy,obsz
      real*8 sinto,costo,sinpo,cospo,phio
      real*8 pi,twopi,fourpi,g2,ximage,yimage,ph,wgt,tau1,imtot
      real*8 delta

      real*8 ran2

      real*8 absCoeff(nxg, nyg, nzg)
      real*8 r,dr
      real*8 w,wa,threshold,chance
      integer cellr,zflag
      
      character(48) filename
      character(41) filePath

c***** Read in parameters from params.par using namelist command *****
      namelist /params/ nphotons,iseed,
     +                  kappa,albedo,hgg,pl,pc,sc,
     +                  xmax,ymax,zmax,rimage,
     +                  viewthet,viewphi,
     +                  lum,b
      open(10,file='params.par',status='unknown')
          read(10,params)
      close(10)
c*********************************************************************

c***** Set up constants, pi and 2*pi  ********************************
      pi=4.d0*datan(1.d0)
      twopi=2.d0*pi
      fourpi=4.d0*pi

      iseed=-abs(iseed)  ! Random number seed must be negative for ran2

      g2=hgg*hgg  ! HG parameter, hgg^2
      
      threshold = 10.d-5
      chance = 0.1d0

c***** Observation direction
      sinto=dsin(viewthet*pi/180.d0)
      costo=dcos(viewthet*pi/180.d0)
      sinpo=dsin(viewphi*pi/180.d0)
      cospo=dcos(viewphi*pi/180.d0)
      phio=datan2(sinpo,cospo)
      obsx=sinto*cospo
      obsy=sinto*sinpo
      obsz=costo

c**** Initialize arrays to zero *************************************
      call iarray(xface,yface,zface,rhokap,
     +                  fimage,qimage,uimage,
     +                  jmean,den,absorb,flux)

      ! set up absorption coefficient structure
      call struc(absCoeff)

c***** Set up density grid *******************************************
      call gridset(xface,yface,zface,rhokap,xmax,ymax,zmax,
     +            kappa,absCoeff)

c***** Set small distance for use in optical depth integration routines 
c***** for roundoff effects when crossing cell walls
      delta=1.d-5*(2.d0*xmax/DBLE(nxg))

c***** Set up point source locations, luminosities, total luminosity
      call sources(xsource,ysource,zsource,lsource,lumtot,iseed,lum)

c      goto 666

c**** Scattered photon loop ******************************************
      totscatt=0
      jcount=0

c**** Loop over sources. nph=number of photons to release from each source
      do is=1,nsource

        nph=int(nphotons*lsource(is)/lumtot)

c**** Loop over nph photons from each source *************************
        do j=1,nph

          jcount=jcount+1
          if(mod(jcount,100000).eq.0)then
             print *, jcount,' scattered photons completed'
          end if

c***** Release photon from point source *******************************
          call sourceph(xsource(is),ysource(is),zsource(is),
     +                  xp,yp,zp,nxp,nyp,nzp,sint,cost,sinp,cosp,phi,
     +                  fi,fq,fu,fv,xmax,ymax,zmax,twopi,
     +                  xcell,ycell,zcell,nxg,nyg,nzg,iseed,b,delta)

          nscatt=0

! commented out to remove force first scattering
!c******* Find optical depth, tau1, to edge of grid
!          call taufind1(xp,yp,zp,nxp,nyp,nzp,xmax,ymax,zmax,
!     +                  xface,yface,zface,rhokap,
!     +                  tau1,xcell,ycell,zcell,delta)
!          if(tau1 .lt. 1.e-3) goto 100
!          wgt1=(1.-exp(-tau1))
!c******* Force photon to scatter at optical depth tau before edge of grid
!          tau=-alog(1.-ran2(iseed)*wgt1)
!c******* Find scattering location of tau
!          call tauint(xp,yp,zp,nxp,nyp,nzp,xmax,ymax,zmax,
!     +                  xface,yface,zface,rhokap,
!     +                  tau,xcell,ycell,zcell)

! Added code
c****** Find scattering location
          call tauint2(xp,yp,zp,nxp,nyp,nzp,xmax,ymax,zmax,
     +                  xface,yface,zface,rhokap,
     +                  xcell,ycell,zcell,tflag,iseed,delta,
     +                  jmean,absorb,flux)

c******** Photon scatters in grid until it exits (tflag=1) or number 
c******** of scatterings exceeds a set value (nscatt)
          nscatt=1
          tflag=0
          zflag=0
          ! set weight to 1
          w = 1.d0

c          dowhile((tflag.eq.0) .and. (nscatt.lt.5))
          dowhile(tflag.eq.0)
c************ Do peeling off and project weighted photons into image
           ! commented out to remove peeling off
!              if( (ran2(iseed).lt.albedo) ) then! .or. (nscatt.eq.1)) then
!              call peeloff(xp,yp,zp,nxp,nyp,nzp,sint,cost,sinp,cosp,phi,
!     +                   fi,fq,fu,fv,
!     +                   obsx,obsy,obsz,sinto,costo,sinpo,cospo,phio,
!     +                   xmax,ymax,zmax,xface,yface,zface,rhokap,
!     +                   rimage,wgt1,albedo,hgg,g2,pl,pc,sc,
!     +                   pi,twopi,fourpi,
!     +                   fimage,qimage,uimage,
!     +                   xcell,ycell,zcell,nscatt,delta)

c************ ! add photon weight to local bin
              wa = w*(1.d0-albedo)
              absorb(xcell,ycell,zcell) = absorb(xcell,ycell,zcell)+wa
              w = w*albedo
c************ Scatter photon into new direction and update Stokes parameters
              call stokes(nxp,nyp,nzp,sint,cost,sinp,cosp,phi,
     +                  fi,fq,fu,fv,pl,pc,sc,hgg,g2,pi,twopi,iseed)
              ! if weight goes below threshold roulette method
              ! to terminate or increase weight of photon
              if (w < threshold) then
                  if (ran2(iseed) <= chance) then
                      w = w/chance
                  else
                      goto 100
                  endif
              endif

c************ Find next scattering location
              call tauint2(xp,yp,zp,nxp,nyp,nzp,xmax,ymax,zmax,
     +                  xface,yface,zface,rhokap,
     +                  xcell,ycell,zcell,tflag,iseed,delta,
     +                  jmean,absorb,flux)
              nscatt=nscatt+1
          end do
          nscatt=nscatt-1
          totscatt=totscatt+nscatt

100      continue

        end do      ! end loop over nph photons

      end do        ! end loop over nsource sources

      print *,"jcount,nphotons,totscatt: ",jcount,nphotons,totscatt

        totscattR = DBLE(totscatt)
        nphotonsR = DBLE(nphotons)
        print*,totscattR,nphotonsR
        print*,'Avereage number of scatterings = ',(totscattR/nphotonsR)
        
        ! Photon absorption count
        open(10,file='data/absorption_count_xz.dat',
     +       status='unknown')
        do k=1,nzg
            write(10,*) (absorb(i,101,k),i=1,nxg)
        enddo
        close(10)
        open(10,file='data/absorption_count_yz.dat',
     +       status='unknown')
        do k=1,nzg
            write(10,*) (absorb(101,j,k),j=1,nyg)
        enddo
        close(10)
        open(10,file='data/absorption_count_xy.dat',
     +       status='unknown')
        do j=1,nyg
            write(10,*) (absorb(i,j,101),i=1,nxg)
        enddo
        close(10)
        
        !filePath = "data\"
        !write(filename,'(I3.3,".dat")')1
        !filename = filepath//filename
        
        ! all slices of cube
        !do i = 1, nxg
        !    !filePath = "data\"
        !    write(filename,'(I3.3,".dat")')i-1
        !    filename = filepath//filename	

        !   open(10,file=filename,status='unknown')
        !   do k=1,nzg
        !        write(10,*) (absorb(i,j,k),j=1,nyg)
        !    enddo
        !    close(10)
        !enddo

!666   continue     ! commented because not used

! commented out all below because don't need image plane
!c****** Direct photon loop.  Loop over sources and weight photons by 
!c****** W=ph*exp(-tau1)/4pi 
!      do is=1,nsource
!
!            ph=int(nphotons*lsource(is)/lumtot)

!  c****** Set photon location, grid cell, and direction of observation
!          xp=xsource(is)
!          yp=ysource(is)
!          zp=zsource(is)
!          xcell=int(nxg*(xp+xmax)/(2.*xmax))+1
!          ycell=int(nyg*(yp+ymax)/(2.*ymax))+1
!          zcell=int(nzg*(zp+zmax)/(2.*zmax))+1
!          nxp=obsx
!          nyp=obsy        
!          nzp=obsz
!
!  c****** Find optical depth, tau1, to edge of grid along viewing direction
!          call taufind1(xp,yp,zp,nxp,nyp,nzp,xmax,ymax,zmax,
!       +                xface,yface,zface,rhokap,
!       +                tau1,xcell,ycell,zcell,delta)
!
!  c****** direct photon weight is exp(-tau1)/4pi
!          wgt=ph*exp(-tau1)/fourpi
!
!  c****** bin the photon into the image according to its position and 
!  c****** direction of travel.  
!          yimage=rimage+zp*sinto-yp*costo*sinpo-xp*costo*cospo
!          ximage=rimage+yp*cospo-xp*sinpo
!          xl=int(nxim*ximage/(2.*rimage))+1
!          yl=int(nyim*yimage/(2.*rimage))+1
!
!          if((xl.gt.nxim) .or. (xl.lt.1)) then
!             print *,'mcpolar: xl out of bounds'
!             goto 111
!          endif
!
!          if((yl.gt.nyim) .or. (yl.lt.1)) then
!             print *,'mcpolar: yl out of bounds'
!             goto 111
!          endif
!
!  c*********** place weighted photon into image location ********************
!          fimage(xl,yl)=fimage(xl,yl)+wgt
!
!   111    continue
!
!        end do
!
!  c***** Normalize images by nphotons
!        do i=1,nxim
!          do j=1,nyim
!             fimage(i,j)=fimage(i,j)/nphotons
!             qimage(i,j)=qimage(i,j)/nphotons
!             uimage(i,j)=uimage(i,j)/nphotons
!          end do
!        end do
!
!        do i=1,nxim
!           do j=1,nyim
!              imtot=imtot+fimage(i,j)
!           end do
!        end do
!        print *,'4*pi*imtot = ',imtot*fourpi
!
!  c***** put results into output files
!        print *,size(fimage)
!
!  c      open(unit=10,file='fimage.dat',status='unknown',
!  c     +     form='unformatted')
!  c           write(10) fimage
!  c      close(10)
!  c      open(unit=10,file='qimage.dat',status='unknown',
!  c     +     form='unformatted')
!  c           write(10) qimage
!  c      close(10)
!  c      open(unit=10,file='uimage.dat',status='unknown',
!  c     +     form='unformatted')
!  c           write(10) uimage
!  c      close(10)
!
!        open(unit=10,file='data/fimage.dat',status='unknown')
!          do i=1,nyim
!             write(10,*) (fimage(j,i),j=1,nxim)
!          end do
!        close(10)
!
!        open(unit=10,file='data/qimage.dat',status='unknown')
!          do i=1,nyim
!             write(10,*) (qimage(j,i),j=1,nxim)
!          end do
!        close(10)
!
!        open(unit=10,file='data/uimage.dat',status='unknown')
!          do i=1,nyim
!             write(10,*) (uimage(j,i),j=1,nxim)
!          end do
!        close(10)
          
      

        stop
        end
