c
c      subroutine      s u _ t r a c e _ t r a c k . f
c
c      Routine verfolgt ein Elektron durch die Feldkarte
c      Die gefundene Trajektorie wird ueber das Feld
c      track(2,300) uebergeben. Teilweise wird
c      doppelt genau gerechnet.
c      Parameter der Routine
c       xl: steuert bei erstem Aufruf der Routine (l=0.)
c           das Einlesen der Feldkarte durch die
c           Unterroutine su_such_feld.
c      Danach ist es die Schrittlaenge beim Tracing in Index-Einheiten
c      1 entspricht incr in cm (hier 1.1055 cm).
c      it ist die Zahl der Iterationsschritte beim Suchen des Feldwertes.
c      Bei einem Fehler in der Routine wird it=0 gesetzt.
c      e ist die Energie des Elektrons in MeV.
c        x0,z0,a0: Startwerte fuer das Elektron am Radiator
c        dabei ist a0 der Winkel/rad (nicht der tan).
c        x,z,a Zielpunkt, auch hier a/rad
c
c      Das Ergebnis wird ins Feld track(2,300) uebergeben
c       
      subroutine su_trace_track(xl,it,x0,z0,a0,e,x,z,a,track)
c
      implicit none
      real*8 ald,al1d,rd,gamd,xd,zd,x1d,z1d,stepd,xld
      real xl,incr,e,x0,z0,a0,x,z,a,b,x1,z1,xg,zg,
     &  x2,z2,alpha,gamma,pi,pi2,feld0,find,r,xm,zm,find0,
     &  alpha1,aa,track(2,300),find1,h,bl
      integer ix,iz,ixm,izm,in,ischritt,i,j,it,ie,ifind,ja,
     &  ired,isch
      save
c
c      Initialisieren der Subroutine - Einlesen der Feldkarte
c      Die Karte ist eine Liste Indexx,Indexz,Feld/Gauss
c      die x-Achse laeuft entlang der Eintrittskante und zaehlt
c      in Ablenkrichtung.
c      Die z-Achse laeuft _|_ dazu grob in Strahlrichtung.
c
      pi=3.1415926
      bl=0.   ! Bahnlaenge entlang der Trajektorie
      if(xl.eq.0.)then  !Routine wird initialisiert
        find=1.E8       !su_such_feld wird initialisiert
        call su_such_feld(x,z,find)
        if(x.eq.0.)then   !Fehler aus su_such_feld
          it=0
          return
        end if
        xg=x    ! Grenzindices fuer die Feldkarte
        zg=50   ! Uebernahme aus su_such_feld.f
        incr=z
        x0=incr ! Uebergabe des Inkrements an das Hauptprogramm
       return
      end if
c      Hier beginnt das Tracing
c      Suche des Anfangsfeldes feld0        
c      Punkt im Feld?
      x=x0
      z=z0
      xd=x0
      zd=z0
      ired=.5/xl+1   ! Reduktionszahl fuer Anzahl Ergebnis-Punkte
      do i=1,2       ! Nullen des Ergebnisfeldes
        do j=1,300
          track(i,j)=0.
        end do
      end do
      track(1,1)=x
      track(2,1)=z
      ischritt=2
      isch=0
      alpha=a0
      ald=a0
      xld=xl
      stepd=xl
c      Feld suchen
6     call su_such_feld(x,z,find0)
      if(x.eq.0)then
        it=0
        return
      end if
      find=find0
4     isch=isch+1
      do i=1,it
        rd=e/(2.9979E-4*(find0+find)/2.*incr)   ! r/Index
        gamd=datan(stepd/rd)
        al1d=ald-gamd
        x1=xd+rd*(dsin(ald)-dsin(al1d))
        z1=zd+rd*(dcos(al1d)-dcos(ald))
        call su_such_feld(x1,z1,find1)
        if(x1.eq.0.) goto 7
        ifind=find1
        if(abs(ifind).lt.5) then
          stepd=stepd+xld
          goto 6
        end if
        find=find1
      end do
      xd=x1
      zd=z1
      x=x1
      z=z1
      ald=al1d
      ie=e
      ifind=find
      find0=find
      if(isch.ge.ired)then
        track(1,ischritt)=xd
        track(2,ischritt)=zd
c        Bestimmung des Abstandes zwischen Punkten
        h=sqrt((track(1,ischritt-1)-xd)**2+(track(2,ischritt-1)-zd)**2)
        bl=bl+h    ! Bahnlaenge
c        write(22,*)xd,zd,bl,find0  ! Feldverlauf f(Ort)
        isch=0
        ischritt=ischritt+1
        if(ischritt.gt.300) then
          print*,'Schrittzahl in su_trace_track zu hoch!'
          print*,'ischritt,e,xg,zg',ischritt,e,xg,zg
          stop
        end if
      end if
      aa=dtan(ald)
      a=ald
c     Abfrage, ob x,z Werte noch innerhalb der Feldkarte sind,
c      wenn ja, weiterrechnen
      if(x.ge.1.and.x.le.xg.and.z.ge.1.and.z.le.zg) goto 4
7     if(x.le.xg.and.z.lt.2) then
          b=z-aa*x
          z=-10.
          x=(z-b)/aa
          if(ischritt.gt.300)print*,'su_tr_track 2',ischritt
          track(1,ischritt)=x
          track(2,ischritt)=z
        end if
        return
      end
