! adapted from code from Kenneth Wood:
! http://www-star.st-and.ac.uk/~kw25/research/montecarlo/points/points.html
! accessed 04.07.2018

      subroutine scatt1(fi,fq,fu,fv,fio,fqo,fuo,fvo,nxp,nyp,nzp,
     +       sint,cost,sinp,cosp,phi,
     +       obsx,obsy,obsz,sinto,costo,sinpo,cospo,phio,hgg,g2,
     +       pl,pc,sc,hgfac,
     +       pi,twopi,fourpi)

      implicit none

      real fi,fq,fu,fv,fio,fqo,fuo,fvo,nxp,nyp,nzp
      real obsx,obsy,obsz,sinto,costo,sinpo,cospo,phio
      real hgg,g2,pl,pc,sc,hgfac,calpha
      real sint,cost,sinp,cosp,phi,pi,twopi,fourpi

c***** calculate calpha=cos(alpha), where alpha is angle between incident
c***** and outgoing (i.e., observed) photon direction
      calpha=nxp*obsx+nyp*obsy+nzp*obsz

c***** weight photon for isotropic, Thomson, or HG scattering
c       hgfac=0.5*(1.+calpha*calpha)                       ! Thomson
c       hgfac=1./fourpi                                    ! isotropic
       hgfac=(1.-g2)/(1.+g2-2.*hgg*calpha)**1.5/fourpi ! HG

      fio=fi
      fqo=fq
      fuo=fu
      fvo=fv
      call scattp(sint,cost,sinp,cosp,phi,
     +      fio,fqo,fuo,fvo,sinto,costo,sinpo,cospo,phio,calpha,
     +      pl,pc,sc,hgg,g2,pi,twopi)

      return
      end

