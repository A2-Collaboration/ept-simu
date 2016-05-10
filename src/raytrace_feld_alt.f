c
c      r a y t r a c e _ f e l d . f
c
c     Programm sucht Elektronenbahn in einer Feldkarte, vorgegeben durch
c     omap**.map.
c     Ergebnis steht in File raytrace.dat
c
c     Laden als:
c     g77 -fbounds-check -o raytrace_feld.exe raytrace_feld.f
c      su_trace_track.f su_such_feld.f
c      
      implicit none
      real xl,x0,z0,a0,a00,x,z,a,e,pi,alpha0,alphae,dalpha,
     &  emin,emax,estep,alpha,track(2,300)
      integer it,i
c
      pi=3.1415926
      print*,'raytrace Feldkarte omap**.map'
      call su_trace_track(0.,it,x0,z0,a0,e,x,z,a,track)
      if(it.eq.0)stop
      xl=.1
      it=10
      e=1500.
      x0=4.83
      z0=7.932
      a00=78.
      emin=10.
      emax=160.
      estep=10.
1     print*,'Eingabe Inkrement l/Index, Iterationen it'
      print*,'[',xl,it,']:'
      read*,xl,it
      print*,'Eingabe x0,z0,a0/Grad [',x0,z0,a00,']:'
      read*,x0,z0,a00
      a0=a00/180*pi
      print*,'emin,emax,estep/MeV [',emin,emax,estep,']:'
      read*,emin,emax,estep
      do e=emin,emax,estep
        call su_trace_track(xl,it,x0,z0,a0,e,x,z,a,track)
        if(it.eq.0)stop
        alpha0=a0
        alphae=a
        dalpha=alpha0-alphae
        print*,'Winkel alpha_0,alpha_end,dalpha/rad',
     &    alpha0,alphae,dalpha
        do i=1,300
          if((track(1,i)+track(2,i)).eq.0.) then
            write(30,100)
            write(22,101)
            goto 10
          else
            write(30,*)track(1,i),track(2,i)
          end if
        end do
 10     continue
 100    format(1x,'    *    *')
 101    format(1x,'    *    *    *    *')
      end do
      print*,'fuer Verwendung in GLE :'
      print*,'Gesamtergebnis fuer Trajektorien in fort.30 fuer'
      print*,'fuer Feldwerte entlang Weg in fort.22 als'
      print*,'als Ortx,Ortz,Weg,Feldwert'
      goto 1
      end
