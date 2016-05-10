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
