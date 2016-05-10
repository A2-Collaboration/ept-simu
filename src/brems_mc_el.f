c
c      b r e m s _ m c _ e l . f
c
c      Programm erwuerfelt Brems-Ereignisse fuer zufaellige
c      Elektronenenergien. Fragt nach den dazugehoerigen
c      Photonen, die durch einen Kollimator
c      passen, verfolgt dann die Elektronen durch eine Feldkarte,
c      omap**.map, sucht nach dem Schnittpunkt der Elektronen-Bahn
c      mit der Detektorlinie und erstellt Histogramme fuer die getroffenen
c      Detektoren als normierte Ereignisse ueber der El-Energie.
c      Die Hisogramme werden in brems_mc_el.dat ausgegeben.
c      Aus den Histogrammen werden mittlere Energie und sigma 
c      bestimmt und im File fort.16 ausgegeben
c
c      In dem Programm werden ausschliesslich Index-Einheiten benutzt,
c      ausser in der Routine su_brems_mc.f .
c      
c      Das Programm ruft Routinen:
c        su_trace_track.f: Bestimmt die Bahn des Elektrons track(2,300)
c        in der Felkarte.
c        Es ruft su_such_feld.f, das den Feldwert an einem Ort
c        interpoliert. Das Ergebnis steht in track(2,300).
c        su_brems_mc.f  : Erwuerfelt das Brems-Ereignis, bestimmt
c        den Startpunkt und -winkel des PBE in Ablenkrichtung und
c        den Winkel des Brems-Photons. Letzterer wird benoetigt,
c        um ds/(dk dO) (br2bs) als Gewicht des Ereignisses
c        zu bestimmen.
c        su_cut_bildlin.f: Bestimmt den Schnittpunkt mit der
c        Bildlinie und den Abstand von Index_x = 0 in Index-Einheiten.
c        Diese Index-Einheiten werden nach Massgabe eines Files,
c        z.B. detektor_eichung.dat, in den Detektoren entsprechende
c        Energie-Werte umgesetzt.
c        Der File z.B. "deichung_140.dat" wird nach Massgabe eines Resultats
c        aus brems_mc.f (fort.17 -> eich140.dat)in eich_leiter.f erstellt.
c        Er enthaelt:  Det.Nr.,Em/MeV,dE,Ortm/Index_x,dOrt
c        Die zugeordneten El.Energien werden
c        dann zur Einordnung in das Histogramm benutzt.
c        Die Dimension der Histogramme wird ebenfalls aus dem Eich-file
c        bestimmt.
c
c        g77 -fbounds-check -o brems_mc_el.exe brems_mc_el.f
c            su_trace_track.f su_such_feld.f brems.f
c            su_cut_bildlin.f su_brems_mc.f
c           
      implicit none
      real e0,e,k,pi,me,hist(47,75),nhist(47),incr,x0,z0,a0,a00,
     &  x,z,a,tb,xr,xr0,txbe,xl,zz,ab,bb,d0,histep,alpha0,
     &  xrr,ar,zr,g,d,track(2,300),summ,ssum,sig,h1,h2,
     &  emin,emax,estep,xxr,dxr,sumtr,hizt(31),sum3,en(47,75),
     &  xyr(2),txybe(2),em(47),dem(47),om(47),dom(47),sum4,
     &  himax,fwhm,sum0(47),sum1,sum2,sulueck,result(47,4),znr
      integer i,j,l,ll,ihi,it,ihi0,izug,ja,ihimi,ihima,il,
     &  izi,iz,istop
      character*25 eichfile
c
      pi=3.1415926
      me=.511006
c
      open(1,file="stop.dat",status='old')
      write(1,*)0,0
      close(1)
      print*,'raytrace Feldkarte omap**#.map'
c      Initialisieren der Feldkarte in su_trace_feld
      call su_trace_track(0.,it,incr,z0,a0,e,x,z,a,track)
c      Initialisieren der routine su_brems_mc.f
      call su_brems_mc(0.,0.,tb,xyr,txybe)
c
      print*,'Eichfile-Name (z.B. deichung_140.dat) [a25]:'
      read '(a)', eichfile
      open(2,file=eichfile,status='old')
      do i=1,47
        read(2,*,err=1,end=1) znr,em(i),dem(i),om(i),dom(i)
      end do
      goto 2
1     i=i-1
      print*,'in Eich-File: nur ',i,'Werte gefunden'
      print*,'letzter Wertesatz  D-Nr, Emit, dE, Ortmit, dOrt:'
      print*,znr,em(i),dem(i),om(i),dom(i)
c      Vorgabe von Werten
2     xl=.11     ! Runku Inkrement in Index
      it=5       ! Runku Iterationen
      e0=1600.  
      x0=5.45    ! Sollstrahl-Eintritt und Richting an Polkante, omap**.map
      z0=11.6    ! in Index-Einheiten
      a00=78.    ! Startwinkel Sollstrahl in Grad!
      xr0=10.4   ! Abstand des Radiators vor der Polkante in index_y
c      Quellort Radiator
      alpha0=a00/180.*pi ! Sollstart-Richtung am Radiator
      x0=x0-xr0/tan(alpha0) ! Radiator Sollort fuer xr0
      z0=z0-xr0             ! vor Polkante
      print*,'Radiator bei x0,z0',x0,z0
c      
      zz=28.     ! Ladungszahl Radiator
      ab=0.      ! Detektor-Linie Steigung
      bb=2.82    ! Achsabschnitt in Index_y
c      
      print*,'Vorlauf Random Generator [1]:'
      i=1
      read*,i
      do j=1,i
        x=rand()
      end do
      print*,'Ladungszahl des Radiators (screening) [',zz,']:'
      read*,zz
      print*,'Zahl der Zuege fuer MC [20000]:'
      izug=20000
      read*,izug
      print*,'einl. Elektronenenergie/MeV  [',e0,']:'
      read*,e0
      print*,'Elektronenenergie/MeV verwuerfeln zwischen min,max'
      emin=23.
      emax=180.
      print*,'[',emin,emax,']:'
      read*,emin,emax
      do i=1,47
        nhist(i)=0
        do j=1,75
          hist(i,j)=0.
        end do
      end do
c
      sulueck=0.
c      
        do iz=1,izug
          e=emin+(emax-emin)*rand()   ! erwuerfelte El.-Energie
c          do e=emin,emax,(emax-emin)/75
          k=e0-e
c   e0,k Energien, tb,xrr,txbe Laborwinkel Photon, Ort-x auf Radiator,
c   Laborwinkel Elektron post-Brems ! Achtung ! hier Orte in cm
          xr=x0        ! xz-Koordinaten in Index
          zr=z0        ! auf dem Radiator
          call su_brems_mc(e0,k,tb,xyr,txybe)
c          xyr sind hier die Abweichungen vom Sollwert in cm !!!!
          ar=alpha0+txybe(1)     ! Eintritts-Winkel Polkante
          xr=xr+xyr(1)/incr
          call br2bs(g,zz,tb,k/me,e0/me,0)  ! Gewicht des Ereignisses
c
          call su_trace_track(xl,it,xr,zr,ar,e,x,z,a,track)
          call su_cut_bildlin(ab,bb,x,z,d,track)
c          die hier benoetigte Groesse ist x (muss hier =d)
          do i=1,47
            if(x.gt.om(i)-dom(i).and.x.lt.om(i)+dom(i))goto 10
          end do
          sulueck=sulueck+1.
          goto 12
10        ihi=(e-em(i))*51./(2.*dem(i))+37.5
          if(ihi.gt.75)ihi=75
          if(ihi.lt.1) ihi=1
          hist(i,ihi)=hist(i,ihi)+g
          nhist(i)=nhist(i)+1.
12      continue
        if(iz/5000*5000.eq.iz)then
          print*,'iz,sulueck',iz,sulueck
          open(3,file='stop.dat',status='old')
          read(3,*)istop
          close(3)
          if(istop.ne.0)goto 30
          end if
c        end do
        end do
        goto 31
c
 30   print*,'Wuerfeln gestoppt bei iz =',iz
c
 31   do i=1,47
        write(17,*)i,nhist(i)
        do j=1,4
          result(i,j)=0.
        end do
      end do
c
      do i=1,47
        sum0(i)=0.
        sum1=0.
        sum2=0.
        sum3=0.
        sum4=0.
        if(nhist(i).ne.0.)then
          do j=1,75
            en(i,j)=em(i)+2.*dem(i)*(j-37.5)/51.
            sum0(i)=sum0(i)+hist(i,j)
          end do
          do j=1,75
            sum1=sum1+hist(i,j)/sum0(i)*en(i,j)
            sum2=sum2+hist(i,j)/sum0(i)*en(i,j)**2
        end do
        result(i,1)=sum1
        result(i,2)=sqrt(sum2-sum1**2)
        result(i,3)=result(i,2)*2.35
        end if
      end do
      do i=1,47
        write(16,*)i,em(i),(result(i,j),j=1,3),sum0(i)/nhist(i)
      end do
c
      open(1,file="brems_mc_el.dat")
      do i=1,47
        if(nhist(i).ne.0.)then
          do j=1,75
          write(1,*)en(i,j),hist(i,j)/sum0(i)
        end do
        write(1,100)
 100    format('     *       *   ')
        end if
      end do
      close(1)
      stop
      end
c      
