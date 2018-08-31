! adapted from code from Kenneth Wood:
! http://www-star.st-and.ac.uk/~kw25/research/montecarlo/points/points.html
! accessed 04.07.2018

      subroutine scattp(sint,cost,sinp,cosp,phi,
     +      fio,fqo,fuo,fvo,sinto,costo,sinpo,cospo,phio,calpha,
     +      pl,pc,sc,hgg,g2,pi,twopi)

      implicit none

      real*8 sint,cost,sinp,cosp,phi,sinto,costo,sinpo,cospo,phio
      real*8 fio,fqo,fuo,fvo

      real*8 pi,twopi
      real*8 p1,p2,p3,p4,a11,a12,a13,a21,a22,a23,a24,a31,a32,a33,a34
      real*8 a42,a43,a44,a,rprob,si,sq,su,sv
      real*8 costp,sintp,phip
      real*8 bmu,b,ri1,ri3,cosi3,sini3,cosb2,sinbt,sini2,bott,cosdph
      real*8 cosi2,sin2i3,sin2i2,cos2i3,cos2i2,sin2,cos2,sin2cos1
      real*8 cos2sin1,cosi1,sini1,sin2i1,cos2i1
      real*8 calpha,sinmu,wght,pl,pc,sc,hgg,g2

      wght=fio
      fio=fio/wght
      fqo=fqo/wght
      fuo=fuo/wght
      fvo=fvo/wght

c***** isotropic scattering ******************************
c      p.a.cost=2.*ran2(iseed)-1.
c      p.a.sint=sqrt(1.-p.a.cost*p.a.cost)
c      p.a.phi=2.*3.14159265*ran2(iseed)
c      p.a.sinp=sin(p.a.phi)
c      p.a.cosp=cos(p.a.phi)
c
c      p.n.x=p.a.sint*p.a.cosp
c      p.n.y=p.a.sint*p.a.sinp
c      p.n.z=p.a.cost
c      goto 100
c***********************************************************

      costp=cost
      sintp=sint
      phip=phi

c***** electron scattering ********************************
10    continue
c      bmu=1.-2.*ran2(iseed)
c      cosb2=bmu*bmu
c      b=cosb2-1.
c      p1=1.+cosb2
c      p2=b
c      p3=2.*bmu
c      p4=0.
c***********************************************************

c***** dust scattering ********************************
c      bmu=((1.+hg2)-((1.-hg2)/(1.-hgg+2.*hgg*ran2(iseed)))**2)/(2.*hgg)
c      cosb2=bmu**2
c      b=cosb2-1.
c      call dustmat(p1,p2,p3,p4,bmu,cosb2)
c      a=p1
c***********************************************************

c      muscatt=bmu
      bmu=calpha
      sinmu=sqrt(1.-bmu*bmu)
      cosb2=bmu**2
      b=cosb2-1.
      call dustmat(p1,p2,p3,p4,bmu,cosb2,pl,pc,sc,hgg,g2,pi)
      a=p1

      if(abs(bmu).gt.1.) then
         if(bmu.gt.1.) then
            bmu=1.
            cosb2=1.
            b=0.
         else
            bmu=-1.
            cosb2=1.
            b=0.
         end if
      end if
      sinbt=sqrt(1.-cosb2)
c      ri1=twopi*ran2(iseed)

c      cosi1=(p.a.cost-po.a.cost*bmu)/(po.a.sint*sinmu)
c      sini1=sin(po.a.phi-p.a.phi)*p.a.sint/sinmu
      cosi1=(costo-cost*bmu)/(sint*sinmu)
      sini1=sin(phi-phio-pi)*sinto/sinmu

      ri1=atan2(sini1,cosi1)+pi

      if(ri1.gt.pi) then
         ri3=twopi-ri1
         cosi3=cos(ri3)
         sini3=sin(ri3)
         sin2i3=2.*sini3*cosi3
         cos2i3=2.*cosi3*cosi3-1.
         a11=p1
         a12=p2*cos2i3
         a13=p2*sin2i3

c****** for electron scattering **********
c         rprob=a11*p.s.i+a12*p.s.q+a13*p.s.u
c         if((2.*ran2(iseed)).gt.rprob) goto 10
c         a=rprob
c******************************************

         if(bmu.eq.1.) then
            goto 100
         else
            if(bmu.eq.-1.) then
               fuo=-fuo
               goto 100
            end if
         end if

c         po.a.cost=costp*bmu+sintp*sinbt*cosi3
         if(abs(costo).lt.1.) then
c            po.a.sint=abs(sqrt(1.-po.a.cost*po.a.cost))
            sini2=sini3*sintp/sinto
            bott=sinto*sinbt
            cosi2=costp/bott-costo*bmu/bott
         else
c            po.a.sint=0.
            sini2=0.
            if(costo.ge.1.)  cosi2=-1.
            if(costo.le.-1.) cosi2=1.
         end if

c         cosdph=-cosi2*cosi3+sini2*sini3*bmu
c         if(abs(cosdph).gt.1.) then
c            if(cosdph.gt.1.) then
c               cosdph=1.
c            else
c               cosdph=-1.
c            end if
c         end if
c         po.a.phi=phip+acos(cosdph)
c         if((po.a.phi).gt.twopi) po.a.phi=po.a.phi-twopi
c         if((po.a.phi).lt.0.)    po.a.phi=po.a.phi+twopi

         sin2i2=2.*sini2*cosi2
         cos2i2=2.*cosi2*cosi2-1.
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

         cosi1=cos(ri1)
         sini1=sin(ri1)
         sin2i1=2.*sini1*cosi1
         cos2i1=2.*cosi1*cosi1-1.
         a11=p1
         a12=p2*cos2i1
         a13=-p2*sin2i1

c************* for electron scattering ****************
c         rprob=a11*p.s.i+a12*p.s.q+a13*p.s.u
c         if((2.*ran2(iseed)).gt.rprob) goto 10
c         a=rprob
c*******************************************************

         if(bmu.eq.1.) then
            goto 100
         else
            if(bmu.eq.-1.) then
               fuo=-fuo
               goto 100
            end if
         end if

c         po.a.cost=costp*bmu+sintp*sinbt*cosi1
         if(abs(costo).lt.1.) then
c            po.a.sint=abs(sqrt(1.-po.a.cost*po.a.cost))
            sini2=sini1*sintp/sinto
            bott=sinto*sinbt
            cosi2=costp/bott-costo*bmu/bott
         else
c            po.a.sint=0.
            sini2=0.
            if(costo.ge.1.)  cosi2=-1.
            if(costo.le.-1.) cosi2=1.
         end if

c         cosdph=-cosi1*cosi2+sini1*sini2*bmu
c         if(abs(cosdph).gt.1.) then
c            if(cosdph.gt.1.) then
c               cosdph=1.
c            else
c               cosdph=-1.
c            end if
c         end if
c         po.a.phi=phip-acos(cosdph)
c         if((po.a.phi).gt.twopi) po.a.phi=po.a.phi-twopi
c         if((po.a.phi).lt.0.)    po.a.phi=po.a.phi+twopi

         sin2i2=2.*sini2*cosi2
         cos2i2=2.*cosi2*cosi2-1.
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

      si=(a11*fio+a12*fqo+a13*fuo)/a
      sq=(a21*fio+a22*fqo+a23*fuo+a24*fvo)/a
      su=(a31*fio+a32*fqo+a33*fuo+a34*fvo)/a
      sv=(a42*fqo+a43*fuo+a44*fvo)/a

      if(abs(sq).gt.abs(si)) then
        print *,'scattp: Q > I'
        print *,calpha,ri1
        print *,a11,a12,a13
        print *,a21,a22,a23
        print *,costo,sinto,phio
        print *
c        sq=0.
      endif
c      if(abs(su).gt.abs(si)) then
c        print *,'scattp'
c        su=0.
c      endif

      fio=si*wght
      fqo=sq*wght
      fuo=su*wght
      fvo=sv*wght

100   continue
      return
      end





















