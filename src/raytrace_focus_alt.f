c
c      r a y t r a c e _ f o c u s . f
c
c     Programm sucht 2 Elektronenbahnen in einer Feldkarte, vorgegeben durch
c     omap**.map.
c      die beiden Bahnen werden bestimmt fuer Sollwinkel + und - einen
c      kleinen Winkelschritt. Danach wird der Schnittpunkt
c      beider Bahnen als Focus gesucht.
c      Die Rouitine liefert eine Bahn als Liste (x,z) fuer einen
c      Startwinkel (Index id). Diese Bahnen werden fuer das
c      Plotprogramm raytrace_focus.gle in das Feld bahn(2,2,300)
c      kopiert. Dabei ist
c      1.Index id d.h. bezeichnet den Startwinkel der Bahn
c      2.Index x,z
c      3.Index  Laufindex
c      Der inhalt von bahn(2,2,300) wird in fort.20 abgelgt.
c      Ergebnis fuer die Schnittpunkte im File fort.14
c
c     Laden als:
c     g77 -fbounds-check -o raytrace_focus.exe raytrace_focus.f
c      su_trace_track.f su_such_feld.f su_find_focus.f
c      
      implicit none
      real xl,x0,z0,a0,a00,x,z,a,aa(2),e,pi,alpha0,alphae,dalpha,
     &  emin,emax,estep,ee,bahn(2,2,300),wink,al,dw,win,xr0,incr,
     &  x1,x2,z1,z2,abild,bbild,track(2,300)
      integer it,id,i,j,ie,l
c
      pi=3.1415926
      print*,'raytrace Feldkarte omap**.map'
      print*,'Suche nach Punkt-Punkt Fokus'
      call su_trace_track(0.,it,x0,z0,a0,e,x,z,a,track)  ! initialisieren
      xl=.1
      it=10
      x0=5.41   ! Koordinaten Sollstrahl an Polkante
      z0=10.66
      a00=78.     ! Richtung Sollstrahl/Grad
      xr0=0.    ! Radiator Abstand/cm vor Polkante
      incr=1.1055  ! (1 Index) in Feldkarte entspricht incr cm
      print*,'Eingabe Inkrement l/Index, Iterationen it'
      print*,'[',xl,it,']:'
      read*,xl,it
      print*,'Eintritt Sollstrahl an Polkante'
      print*,'Eingabe x0,z0,a0/Grad [',x0,z0,a00,']:'
      read*,x0,z0,a00
      a0=a00/180.*pi
      alpha0=a0
      print*,'Radiator in Strahlrichtung [',xr0,'cm] vor der Polkante'
      read*,xr0
      x0=x0-xr0/incr*cos(alpha0)  ! Radiator-Ort (index)  fuer xr0/cm
      z0=z0-xr0/incr*sin(alpha0)
      print*,'Winkelschritt, um den die Sollrichtung + und -'
      print*,'geaendert wird [1 Grad]:'
      wink=1
      read*,wink
      print*,'Energien/MeV emin,emax,estep [*,*,*]:'
      read*,emin,emax,estep
      dw=wink/180.*pi
      aa(1)=a0+dw
      aa(2)=a0-dw
      ie=0
      do e=emin,emax,estep
        ie=ie+1
        ee=e
        do i=1,300
          do j=1,2
            do id=1,2
              bahn(id,j,i)=0.0
            end do
          end do
        end do
        do id=1,2
          call su_trace_track(xl,it,x0,z0,aa(id),ee,x,z,a,track)
          if(it.eq.0.)stop
          do i=1,300
            bahn(id,1,i)=track(1,i)
            bahn(id,2,i)=track(2,i)
          end do
        end do
        do l=1,300
          write(21,*)bahn(1,1,l),bahn(1,2,l),
     &      bahn(2,1,l),bahn(2,2,l)
        end do
        call su_find_focus(bahn,x,z)
        write(14,*) x,z,e
        if(ie.eq.2) then
          x1=x
          z1=z
        end if
        print*,'e-Schnittpunkt  e,x,z ',e,x,z
        x2=x
        z2=z
        do l=1,300
          if(bahn(1,1,l)*bahn(1,2,l)*
     &      bahn(2,1,l)*bahn(2,2,l).ne.0.)then
          write(20,*)bahn(1,1,l),bahn(1,2,l),
     &      bahn(2,1,l),bahn(2,2,l)
        else
        if(bahn(1,1,l)*bahn(1,2,l).ne.0.) 
     &    write(20,102) bahn(1,1,l),bahn(1,2,l)
        if(bahn(2,1,l)*bahn(2,2,l).ne.0.)
     &    write(20,102) bahn(2,1,l),bahn(2,2,l)
        end if
      end do
      abild=(z2-z1)/(x2-x1)
      bbild=z1-abild*x1
      print*,'Bildlinie zwischen 2. und letzten Schnittpunkt'
      print*,'a,b,',abild,bbild
      write(20,100)
 100  format( '   *     *     *     * ')
 101  format( '   *     *   ',2f10.3)
 102  format(2f10.3,'   *    *  ')
      end do
      print*,'Ergebnis in fort.14, in z.B. focus14.dat umbenennen'
      stop
      end
