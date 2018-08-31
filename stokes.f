! adapted from code from Kenneth Wood:
! http://www-star.st-and.ac.uk/~kw25/research/montecarlo/points/points.html
! accessed 04.07.2018

      subroutine stokes(nxp,nyp,nzp,sint,cost,sinp,cosp,phi,
     +                  fi,fq,fu,fv,pl,pc,sc,hgg,g2,pi,twopi,iseed)

      implicit none

      integer iseed
      real*8 nxp,nyp,nzp,sint,cost,sinp,cosp,phi,fi,fq,fu,fv
      real*8 pl,pc,sc,hgg,g2,pi,twopi,wght


      real*8 p1,p2,p3,p4,a11,a12,a13,a21,a22,a23,a24,a31,a32,a33,a34
      real*8 a42,a43,a44,a,rprob,si,sq,su,sv
      real*8 costp,sintp,phip
      real*8 bmu,b,ri1,ri3,cosi3,sini3,cosb2,sinbt,sini2,bott,cosdph
      real*8 cosi2,sin2i3,sin2i2,cos2i3,cos2i2,sin2,cos2,sin2cos1
      real*8 cos2sin1,cosi1,sini1,sin2i1,cos2i1

      real*8 ran2

      wght=fi
      fi=fi/wght
      fq=fq/wght
      fu=fu/wght
      fv=fv/wght

      ! uncommented this section as isotropic scattering is assumed
c***** isotropic scattering ******************************
      cost=2.d0*ran2(iseed)-1.d0
      sint=(1.d0-cost*cost)
      if(sint.le.0.d0)then
        sint=0.d0
      else
        sint=dsqrt(sint)
      endif

      phi=twopi*ran2(iseed)
      sinp=dsin(phi)
      cosp=dcos(phi)

      nxp=sint*cosp
      nyp=sint*sinp
      nzp=cost

      goto 100
c***********************************************************

      costp=cost
      sintp=sint
      phip=phi

c***** electron scattering ********************************
10    continue
c      bmu=1.d0-2.d0*ran2(iseed)
c      cosb2=bmu*bmu
c      b=cosb2-1.d0
c      p1=1.d0+cosb2
c      p2=b
c      p3=2.d0*bmu
c      p4=0.d0
c***********************************************************

c***** dust scattering ********************************
      bmu=((1.d0+g2)-((1.d0-g2)/(1.d0-hgg+2.d0*hgg*ran2(iseed)))**2)
      bmu = bmu / (2.d0 * hgg)  ! added because above line got too long
      cosb2=bmu**2
      b=cosb2-1.d0
      call dustmat(p1,p2,p3,p4,bmu,cosb2,pl,pc,sc,hgg,g2,pi)
      a=p1
c***********************************************************

      if(abs(bmu).gt.1.d0) then
         if(bmu.gt.1.d0) then
            bmu=1.d0
            cosb2=1.d0
            b=0.d0
         else
            bmu=-1.d0
            cosb2=1.d0
            b=0.d0
         end if
      end if
      sinbt=dsqrt(1.d0-cosb2)
      ri1=twopi*ran2(iseed)

      if(ri1.gt.pi) then
         ri3=twopi-ri1
         cosi3=dcos(ri3)
         sini3=dsin(ri3)
         sin2i3=2.d0*sini3*cosi3
         cos2i3=2.d0*cosi3*cosi3-1.d0
         a11=p1
         a12=p2*cos2i3
         a13=p2*sin2i3

c****** for electron scattering **********
c         rprob=a11*fi+a12*fq+a13*fu
c         if((2.d0*ran2(iseed)).gt.rprob) goto 10
c         a=rprob
c******************************************

         if(bmu.eq.1.d0) then
            goto 100
         else
            if(bmu.eq.-1.d0) then
               fu=-fu
               goto 100
            end if
         end if

         cost=costp*bmu+sintp*sinbt*cosi3
         if(dabs(cost).lt.1.d0) then
            sint=dabs(dsqrt(1.d0-cost*cost))
            sini2=sini3*sintp/sint
            bott=sint*sinbt
            cosi2=costp/bott-cost*bmu/bott
         else
            sint=0.d0
            sini2=0.d0
            if(cost.ge.1.d0)  cosi2=-1.d0
            if(cost.le.-1.d0) cosi2=1.d0
         end if

         cosdph=-cosi2*cosi3+sini2*sini3*bmu
         if(abs(cosdph).gt.1.d0) then
            if(cosdph.gt.1.d0) then
               cosdph=1.d0
            else
               cosdph=-1.d0
            end if
         end if
         phi=phip+dacos(cosdph)
         if(phi.gt.twopi) phi=phi-twopi
         if(phi.lt.0.d0)    phi=phi+twopi

         sin2i2=2.d0*sini2*cosi2
         cos2i2=2.d0*cosi2*cosi2-1.d0
         sin2=sin2i2*sin2i3
         cos2=cos2i2*cos2i3
         sin2cos1=sin2i2*cos2i3
         cos2sin1=cos2i2*sin2i3

         a21=p2*cos2i2
         a22=p1*cos2-p3*sin2
         a23=p1*cos2sin1+p3*sin2cos1
         a24=-p4*sin2i2
         a31=-p2*sin2i2
         a32=-p1*sin2cos1-p3*cos2sin1
         a33=-p1*sin2+p3*cos2
         a34=-p4*cos2i2
         a42=-p4*sin2i3
         a43=p4*cos2i3
         a44=p3

c      elseif(ri1.le.pi) then
      else  

         cosi1=dcos(ri1)
         sini1=dsin(ri1)
         sin2i1=2.d0*sini1*cosi1
         cos2i1=2.d0*cosi1*cosi1-1.d0
         a11=p1
         a12=p2*cos2i1
         a13=-p2*sin2i1

c************* for electron scattering ****************
c         rprob=a11*fi+a12*fq+a13*fu
c         if((2.d0*ran2(iseed)).gt.rprob) goto 10
c         a=rprob
c*******************************************************

         if(bmu.eq.1.d0) then
            goto 100
         else
            if(bmu.eq.-1.d0) then
               fu=-fu
               goto 100
            end if
         end if

         cost=costp*bmu+sintp*sinbt*cosi1
         if(dabs(cost).lt.1.d0) then
            sint=dabs(dsqrt(1.d0-cost*cost))
            sini2=sini1*sintp/sint
            bott=sint*sinbt
            cosi2=costp/bott-cost*bmu/bott
         else
            sint=0.d0
            sini2=0.d0
            if(cost.ge.1.d0)  cosi2=-1.d0
            if(cost.le.-1.d0) cosi2=1.d0
         end if

         cosdph=-cosi1*cosi2+sini1*sini2*bmu
         if(dabs(cosdph).gt.1.d0) then
            if(cosdph.gt.1.d0) then
               cosdph=1.d0
            else
               cosdph=-1.d0
            end if
         end if
         phi=phip-acos(cosdph)
         if(phi.gt.twopi) phi=phi-twopi
         if(phi.lt.0.d0)    phi=phi+twopi

         sin2i2=2.d0*sini2*cosi2
         cos2i2=2.d0*cosi2*cosi2-1.d0
         sin2=sin2i2*sin2i1
         cos2=cos2i2*cos2i1
         sin2cos1=sin2i2*cos2i1
         cos2sin1=cos2i2*sin2i1

         a21=p2*cos2i2
         a22=p1*cos2-p3*sin2
         a23=-p1*cos2sin1-p3*sin2cos1
         a24=p4*sin2i2
         a31=p2*sin2i2
         a32=p1*sin2cos1+p3*cos2sin1
         a33=-p1*sin2+p3*cos2
         a34=-p4*cos2i2
         a42=p4*sin2i1
         a43=p4*cos2i1
         a44=p3

      end if

      si=(a11*fi+a12*fq+a13*fu)/a
      sq=(a21*fi+a22*fq+a23*fu+a24*fv)/a
      su=(a31*fi+a32*fq+a33*fu+a34*fv)/a
      sv=(a42*fq+a43*fu+a44*fv)/a

      fi=si*wght
      fq=sq*wght
      fu=su*wght
      fv=sv*wght

      cosp=dcos(phi)
      sinp=dsin(phi)

      nxp=sint*cosp
      nyp=sint*sinp
      nzp=cost

100   continue

      return
      end


	subroutine dustmat(p1,p2,p3,p4,cost,cost2,pl,pc,sc,hgg,g2,pi)

c				a.d. code october 25, 1989
c				revised baw apr 10, 1990
c	**********************************************************************
c
c	this program calculates the elements of the phase matrix for a
c	simple representation of the mrn dust mixture using the algorithms
c	for the ultraviolet region due to richard l. white ap.j. 229, 954,
c	1979.
c
c	***********************************************************************
c
c	    cost = cos(angle) of scattering (i.e. angle between incident 
c	        photon and scattered photon)
c	    g = mean value of cosine of scattering angle (henyey-greenstein)
c	    pl = peak linear polarization
c	    pc = peak value of linear to circular conversion
c	    sc = asymmetry of the circular polarization.
c	    p1 = intensity phase function
c	    p2 = polarization function
c	    p3 = skew polarization
c	    p4 = circular polarization
c
c	    the scattering matrix for (i,q,u,v) is of the form
c
c		 p1    p2    0	  0
c		 p2    p1    0	  0
c		 0     0     p3  -p4
c		 0     0     p4   p3
c
c	**********************************************************************

      implicit none

      real*8 p1,p2,p3,p4,cost,cost2,pl,pc,sc,hgg,g2,pi
      real*8 phi,f,f2,c

      p1 = (1.d0 - g2)/(1.d0+g2-2.d0*hgg*cost)**1.5
      p2 = -pl*p1*(1.d0-cost2)/(1.d0+cost2)
      p3 = p1*2.d0*cost/(1.d0+cost2)
      if(dabs(cost).gt.1.d0) then 
        print *, 'in dustmat, cost.gt.1',cost
        if(cost.gt.1.d0) then     
           cost=1.d0
        else    
           cost=-1.d0
        end if 
      end if
c     angle in degrees!
      phi=dacos(cost)*180.d0/pi
      f=3.13d0*phi*dexp(-7.d0*phi/180.d0)
c     now convert to radians
      f2=(phi+sc*f)*pi/180.d0
c     fphi= (1.d0+3.d13*sc*dexp(-7.d0*phi/pi))*phi
      c=(dcos(f2))**2
      p4 = -pc*p1*(1.d0-c)/(1.d0+c)

c******* Isotropic scattering
c      p1 = 1.0
c      p3 = 1.0
c      p2 = 0.0
c      p4 = 0.0  

	return
	end


