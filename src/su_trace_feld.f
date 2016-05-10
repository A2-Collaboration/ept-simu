c
c      subroutine    s u _  t r a c e _ f e l d . f
c
c      Routine verfolgt ein Elektron durch die Feldkarte
c      Parameter der Routine
c        l: steuert bei erstem Aufruf der Routine (l=0.)
c           das Einlesen der Feldkarte.
c      Danach ist es die Schrittlaenge beim Tracing in Index-Einheiten
c      1 entspricht 1.1055 cm
c      it ist die Zahl der Iterationsschritte beim Suchen des Feldwertes.
c      Bei einem Fehler in der Routine wird in=0 gesetzt.
c      e ist die Energie des Elektrons in MeV.
c        x0,z0,a0: Startwerte fuer das Elektron am Rand der Feldkarte (bei z=1).
c        x,z,a Zielpunkt an einem Rand der Feldkarte (Ergebnis)
c
c       
      subroutine trace_feld(xl,it,x0,z0,a0,e,x,z,a)
c
      implicit none
      real feld(200,200),xl,incr,e,x0,z0,a0,x,z,a,x1,z1,
     &  x2,z2,alpha,gamma,pi,pi2,feld0,find,r,xm,zm,find0
      integer ix,iz,ixm,izm,in,ischritt,i,j,it,ie,ifind
      character*20 mfile
      save
c
c      Initialisieren der Subroutine - Einlesen der Feldkarte
c      Die Karte ist vom cmap**.map Typ, sie ist eine Liste Indexz,Indexx,Feld/Gauss
c      die x-Achse laeuft entlang der Eintrittskante und zaehlt in Ablenkrichtung.
c      Die z-Achse laeuft _|_ dazu grob in Strahlrichtung.
c      Der Orts-Bezugspunkt, Kante am Polschuh, ist (0.06,10.66).
c
      pi=3.1415926
      pi2=pi/2
      if(xl.eq.0.)then
        print*,'map-file (z.B. cmap14.map) [a20]:'
        read '(a)', mfile
        open(1,file=mfile,status='old')
        read(1,*)ixm,izm,incr
        do i=1,ixm
          do j=1,izm
            read(1,*,err=2,end=2)ix,iz,feld(ix,iz)
            if(i.ne.ix.or.j.ne.iz)goto 1
          end do
        end do
        print*,'i,j,incr',i,j,incr
        goto 2
1       print*,'Fehler mit Anzahl der Werte im Mapfile'
        print*,'i soll, ist  ',i,ix,'  j soll ist  ',j,iz
        it=0
2       return
      end if
c      Hier beginnt das Tracing
c      Suche des Anfangsfeldes feld0        
c      Punkt im Feld?
      x=x0
      z=z0
      ischritt=0
10    if(x.lt.1.or.x.gt.80.or.z.lt.1.or.z.gt.54)then
        print*,'Punkt ausserhalb des Feldes  ',x,z
        return
      end if
c      Feld suchen
      call such_feld(x,z,feld,find0)
      r=e/(3.E-4*find0*incr)  ! r in Index-Koordinate
      alpha=atan(a0)
4     ischritt=ischritt+1
      print*,'ischritt',ischritt
      x1=x+xl*cos(alpha)
      z1=z+xl*sin(alpha)
      print*,'x1,z1,feld',x1,z1,find0
      call such_feld(x1,z1,feld,find)
      print*,'find',find
      do i=1,it
        find0=(find0+find)/2.
        r=e/(3.E-4*find0)   ! r/cm
        xm=x+r*cos(alpha-pi2)
        zm=z+r*sin(alpha-pi2)
        gamma=atan(xl/r)
        print*,'find0,r,xm,zm,gamma',find0,r,xm,zm,gamma
        x2=r*cos(pi2-alpha+gamma)
        z2=r*sin(pi2-alpha+gamma)
        print*,'Zwischenwert x2,z2',x2,z2
        x2=xm-x2
        z2=zm+z2
        print*,'x2,z2,find0',x2,z2,find0
        find=find0
        call such_feld(x2,z2,feld,find0)
        print*,'find0',find0
        read*,xl
      end do
      x=x2
      z=z2
      ie=e
      ifind=find0
      write(10,*)x,z,alpha,r,ifind,ie
      alpha=alpha-gamma
      a=tan(alpha)
      print*,'alpha',alpha/pi*180,' Grad'
      if(x.lt.1.or.x.gt.80.or.z.lt.1.or.z.gt.54)return
      goto 4
      end
