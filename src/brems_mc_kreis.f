c
c      b r e m s _ m c _ k r e i s . f
c

c      Programm erwuerfelt Brems-Ereignisse. Fragt nach den
c      Partner-Elektronen von Photonen, die durch einen Kollimator
c      passen, verfolgt diese Elektronen durch ein Feld (Geraden,
c      Kreise), sucht nach dem Schnittpunkt der Elektronen-Bahn
c      mit der Bildlinie (0,0) und erstellt ein Histogramm,
c      normiert auf die Zahl der Zuege und eine Hist-Breite von .002.
c      das in brems_mc_kreis.dat ausgegeben wird.
c
c      Das Programm ruft die Routine:
c        su_brems_mc.f : Erwuerfelt das Brems-Ereignis, bestimmt
c        den Startpunkt und -winkel des PBE in Ablenkrichtung und
c        den Winkel des Brems-Photons. Letzterer wird benoetigt,
c        um ds/(dk dO) (br2bs) als Gewicht des Ereignisses
c        zu bestimmen.
c        Diese Routine steht in diesem File!
c        Der Abstand von Index_x = 0 in Index-Einheiten wird bestimmt.
c        Letzterer wird fuer die Erstellung des Histogramms benoetigt.
c
c        In fort.16 stehen die Parameter der ermittelten Response als
c        Energie des PBE/MeV, Schwerpunkt, sigma, Orte min,max,diff/cm
c        Laden als:
c        g77 -fbounds-check -o brems_mc_kreis.exe brems_mc_kreis.f brems.f
c           
      implicit none
      real e0,e,k,pi,me,hist(50000),incr,x0,z0,a0,a00,zeff,
     &  x,z,a,tb,xr,xr0,txbe,xl,zz,ab,bb,d0,histep,alpha0,
     &  xrr,ar,zr,g,d,track(2,300),summ,ssum,sig,h1,h2,bm,
     &  emin,emax,estep,xxr,dxr,sumtr,hizt(31),xyr(2),txybe(2),
     &  dx1,dx2,dx3,radius
      integer i,j,l,ll,ihi,it,ihi0,izug,ja,ihimi,ihima,il,izi
c
      pi=3.1415926
      me=.511006
      incr=1.1055
c
c      Initialisieren der routine su_brems_mc.f
      call su_brems_mc(0.,0.,tb,xyr,txybe)
c
      e0=1600.  ! Einschuss-Energie des Elektrons vor BP
      bm=14000.  ! konstantes Feld/Gauss im Magneten
      x0=5.41   ! Sollstrahl-Eintritt und Richting an Polkante, omap**.map
      z0=10.66  ! in Index-Einheiten
      a00=78.    ! Startwinkel Sollstrahl in Grad!
      xr0=10.    ! Lage des Radiators vor der Polkante in cm
      zz=28.    ! Ladungszahl Radiator
      ab=0.     ! Bildgerade, in Index
      bb=0.
      print*,'mittleres Feld/Gauss im Tagger [',bm,']:'
      read*,bm
      print*,'effektive Kante [bei Polkante',z0,']:'
      zeff=z0
      read*,zeff
      print*,'Vorlauf Random Generator [1]:'
      i=1
      read*,i
      do j=1,i
        x=rand()
      end do
      print*,'Ladungszahl des Radiators (screening) [',zz,']:'
      read*,zz
c      print*,'Bildlinie a,b in Index bei',ab,bb,
c      Gesamt-Laenge der Bildlinie zwischen index_x 0 und 90 fuer Hist
 10   d0=sqrt(90**2+(90.*ab)**2)
      print*,'Gesamtlaenge der Bildlinie in Index_Einheiten'
      print*,'zwischen index_x = 0 bis 90',d0
      i=d0
      histep=1./49000*i
      ihi0=d0/histep  ! Zahl der Histo-Schritte auf Bildlinie
      print*,'Schrittweite fuer Histogramm [',histep,']:'
      read*,histep
      ihi0=d0/histep
      if(ihi0.gt.49990.) goto 10
      print*,'das ergibt ',ihi0,' Histogramm-Schritte'
      print*,'bei 50000 vorgesehenen'
      print*,'Ort/Index Eintritt und Richtung (Grad rel. x-Achse)'
      print*,'x0,z0,a0/Grad [',x0,z0,a00,']:'
      read*,x0,z0,a00
      print*,'Radiator [',xr0,' cm] vor der Polkante'
      read*,xr0
      alpha0=a00/180.*pi ! Sollstart-Richtung am Radiator
      x0=x0-xr0/incr*cos(alpha0) ! Radiator Sollort fuer xr0/cm
      z0=z0-xr0/incr*sin(alpha0) ! vor Polkante
      print*,'Zahl der Zuege fuer MC [20000]:'
      izug=20000
      read*,izug
      print*,'einl. Elektronenenergie/MeV  [',e0,']:'
      read*,e0
      print*,'ausl. Elektronenenergie/MeV min,max,step'
      emin=10.
      emax=180.
      estep=10.
      print*,'[',emin,emax,estep,']:'
      read*,emin,emax,estep
      open(1,file="brems_mc_kreis.dat")
      do e=emin,emax,estep
        k=e0-e
        radius=e/(3.E-4*bm*incr)  ! Radius im Feld bm in Index
        do i=1,50000
          hist(i)=0.
        end do
        summ=0.       ! Groessen nullen fuer response, Schwerpunkt, sigma ...
        ssum=0.      ! Tracklaenge
        ihimi=50000
        ihima=0
        sig=0.
        do i=1,izug
          call su_brems_mc(e0,k,tb,xyr,txybe)
c   e0,k Energien, tb,xrr,txbe Laborwinkel Photon, Ort-x auf Radiator,
c   Laborwinkel Elektron, PBE
          xr=x0+xyr(1)/incr*sin(alpha0)  ! xz-Koordinaten 
          zr=z0+xyr(1)/incr*cos(alpha0)  ! auf dem Radiator
          ar=alpha0+txybe(1)   ! Eintritts-Winkel Polkante
          call br2bs(g,zz,tb,k/me,e0/me,0)
c   Jetzt Rechnung mit Kries und Geraden
          dx1=(zeff-zr)/tan(ar)
          dx2=radius*sin(ar)
          dx3=zeff/tan(ar)
          d=xr+dx1+2.*dx2+dx3
          ihi=d/histep+1.
          if(ihi.lt.1)ihi=1.
          if(ihi.gt.50000)ihi=50000
          hist(ihi)=hist(ihi)+g
          summ=summ+g
          ssum=ssum+g*d
          if(ihi.gt.ihima)ihima=ihi
          if(ihi.lt.ihimi)ihimi=ihi
        end do
        ssum=ssum/summ  ! Schwerpunkt der Verteilung in Index!
        il=0
        do i=1,50000
          if(hist(i).ne.0.)then
            write(1,*) histep*(-.5+i)*incr,
     &        hist(i)/izug/histep*.002
            sig=sig+hist(i)*(ssum-histep*(-.5+i))**2
            il=il+1
          end if
        end do
         sig=sqrt(sig/summ)*incr
        h1=histep*(-.5+ihimi)*incr
        h2=histep*(-.5+ihima)*incr
        write(1,*) '       *      *   '
        write(16,*)e,ssum*incr,sig,h1,h2,h2-h1
        print*,'Ergebnis Apfu in brems_mc_kreis.dat fuer e0 und e/MeV',
     &    e0,e
        print*,'Schwerpunkt ',ssum*incr
        print*,'sigma       ',sig
        print*,'Orte min,max,diff',h1,h2,h2-h1
        print*,'Diese obigen Daten in fort.16 in cm als'
        print*,'E/MeV,Schwerpunkt,sigma,min.,max.Ort,deren Differenz'
      end do
      close(1)
      stop
      end
c      
c     ********************** Routinen *****************************
c
c
c      subroutine s u _ b r e m s _ m c . f
c
c      Subroutine zur Nutzung beim Erwuerfeln von Bremsstrahlungs(B)-
c      Ereignissen im Zusammenhang mit einem Tagger.
c
c      Elektronen treffen auf den Radiator (r) nach Massgabe
c      der Emittanz (s) des Elektronenstrahls. Dort werden also
c      Ort und Winkel fuer horizontal (x=Ablenkrichtung)
c      und vertikal (y) erwuerfelt.
c      Die entsprechenden Parameter sind xr,txr,yr,tyr
c
c      Durch Vielfachstreuung kann das Elektron vor dem B-Prozess
c      eine Aenderung des Winkels erfahren. Die entsprechenden
c      Winkelkomponenten sind txv und tyv.
c      
c      Von dem so gegebenen Ort auf dem Radiator startet
c      ein Bremsstrahlungsphoton, das durch einen Kollimator (k) teffen
c      muss. Der Ort im Kollimator wird gleichverteilt als (xk,yk) erwuerfelt.
c
c      Somit ist die Bahn des Photons bestimmt. Daraus folgen
c      die Winkelkomponenten txrk und tyrk. Die Komponenten
c      fuer den Bremsstrahlungswinkel ergeben sich dann als
c      txb=txrk-txr-txv und tyb=tyrk-tyr-tyv. Hier spielen die Vorzeichen
c      eigentlich keine Rolle, da alle Groessen erwuerfelt sind.
c
c      Der Bremsstrahlungswinkel fuer die Errechnung des WQ ist dann
c      tb=sqrt(txb**2+tyb**2).
c
c      Der Winkel des PBE in Ablenkrichtung(x) im Labor ergibt sich als
c      txbe=-txb*k/e-txr-txv, wieder sind die Vorzeichen ohne Belang.
c
c      Dieser Winkel muss nach Massgabe der Vielfachstreuung des
c      weichen PBE veraendert werden. Hierzu wird die Restlaenge im Radiator
c      verwendet.
c      
c      Das Gewicht des Ereignisses folgt aus dem Brems-Querschnitt
c      dsigma/(dk dO) mit den Parametern k,e0,tb (br2bs) in barn/sr
c      (wird ausserhalb dieser Routine bestimmt).
c
c      Die Parameter der Anordnung stehen im File su_brems_mc.geo,
c      der beim Initialisieren der Routine aufgerufen wird (siehe unten).
c      Hierzu werden die Energien k und e0 mit 0 angegeben.
c      Tritt beim Einlesen dieses Files ein Fehler auf, folg g=0.
c
c      In dem File su_brems_mc.geo steht:
c      Emittanz horizontal, vertikal, als sx,stx,sy,sty,
c      Radius des Kollimators
c      und dessen Entfernung vom Radiator, als rk,lk,
c      fuer die Vielfachstreuung die Strahlungslaenge
c      t0 und die Dicke des Radiators dr in cm.
c
c      
      subroutine su_brems_mc(e0,k,tb,xyr,txybe)
c
c      Dabei bedeuten als Eingabe:
c      e0,k die Energien/MeV des einlaufenden Elektrons und des Photons,
c      (bei e0=k=0. wird die Routine initialisiert)
c      als Ausgabe:
c      tb,xyr(2),txybe(2) der Labor-Winkel des B-Photons, der Radiator-Ort
c      des Elektrons und die Winkel-Komponente in Ablenkrichtung des PBE.
c      Dabei steht der Index (1),(2) fuer die x- und die y-Richtung
c      
      implicit none
      real e0,e,k,sx,sy,stx,sty,rk,lk,xr,yr,txr,tyr,rax,
     &  txv,tyv,xk,yk,txrk,tyrk,txb,tyb,tb,pi,me,xyr(2),txybe(2),
     &  r,r1,r2,z,txbe,g,t0,dr,fi,tv,teb,ddr,ddrv,tvp,tvpx,tvpy
      integer i,j,l
c
      save
      pi=3.1415926
      me=.511006
c      
      if(k+e0.eq.0.)then
        g=1.
        i=0
        open(1,file="su_brems_mc.geo",status="old")
        read(1,*,err=20)sx,stx,sy,sty
        read(1,*,err=20)rk,lk,t0,dr
        close(1)
        rax=1.
        print*,'Radiator-Breite in Ablenkrichtung einschraenken?'
        print*,'halbe Breite/cm [',rax,']:'
        read*,rax
        return
 20     print*,'Fehler beim Einlesen des su_brems_mc.geo'
        return
      else
        i=i+1   ! laufende Nummer der Zuege
        e=e0-k  ! PBE
c        Radiator Reaktionspunkt (xr,yr) Komponenten Polarwinkel (txr,tyr)
 1      xr=sx*sqrt(-2.*log(rand()))*cos(2.*pi*rand())
        if(abs(xr).gt.rax) goto 1
        yr=sy*sqrt(-2.*log(rand()))*cos(2.*pi*rand())
        xyr(1)=xr
        xyr(2)=yr
        txr=stx*sqrt(-2.*log(rand()))*cos(2.*pi*rand())
        tyr=sty*sqrt(-2.*log(rand()))*cos(2.*pi*rand())
c        Vielfachstreuung im Radiator
        ddr=dr*rand()   ! Erwuerfeln des Reaktionsortes im Radiator
        tv=20./e0*sqrt(ddr/t0) ! mittlerer Vf-Streuwinkel
        tv=tv*sqrt(-2*log(rand())) ! obiger verwuerfelt
        ddrv=dr-ddr     ! Restlaenge fuer weiches PBE
        fi=2.*pi*rand()
        txv=tv*cos(fi)   ! x-y Komponenten des Vf-Streuwinkels
        tyv=tv*sin(fi)
c        Ort des Photons im Kollimator
        r=rk*sqrt(rand())
        fi=2.*pi*rand()   ! r und fi fuer Gleichverteilung im Kollimator
        xk=r*cos(fi)
        yk=r*sin(fi)
c        Photon fliegt von (xr,yr) nach (xk,yk) bei Abstand lk
c        Winkelkomponenten hierzu
        txrk=atan((xk-xr)/lk)
        tyrk=atan((yk-yr)/lk)
c        dann Winkel-Komponenten fuer den Brems-Prozess und
c        Winkel des Photons
        txb=txrk-txr-txv ! Winkel x
        tyb=tyrk-tyr-tyv ! Winkel y
        tb=sqrt(txb**2+tyb**2)   ! Photonenwinkel fuer Brems
c        Elektron-Winkel aus Brems
        teb=tb*k/e
c        Elektronenwinkel in Ablenkrichtung im Labor 
        txybe(1)=-txb*k/e-txr-txv
        txybe(2)=tyb*k/e-tyr-tyv
        tvp=20./e*sqrt(ddrv/t0) ! Vielfachstreuwinkel fuer weiches PBE
        tvp=tvp*sqrt(-2*log(rand()))
        fi=2.*pi*rand()
        tvpx=tvp*cos(fi)
        tvpy=tvp*sin(fi)
        txybe(1)=txybe(1)+tvpx
        txybe(2)=txybe(2)+tvpy
        return
        end if
        end
        
c
