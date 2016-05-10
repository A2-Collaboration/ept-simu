c
c      s u c h _ f e l d . f
c
c      Programm nutzt die Routine su_such_feld, um Feldwerte
c      fuer vorgegebene Orte zu finden
c
c      Nach Initialisierung der Routine wird nach Eingabe eines
c      Punktes in der Feldkarte der Feldwert durch lineare
c      Interpolatin bestimmt.
c
c      Laden mit:
c      g77 -fbounds-check -o such_feld.exe such_feld.f su_such_feld.f
c
      implicit none
      real xl,x,z,find,xg,zg,a,b,xmax,xmin,xstep,xhall,zhall,
     &     fmax
      integer it,ja
c
      a=1.578
      b=-29.63
      xhall=47.37
      zhall=45.11
      print*,'Nutzung der Routine su_such_feld.f'
      find=1.E30   ! Routine wird initialisiert
      call su_such_feld(x,z,find)
      fmax=find
      if(x.eq.0.)then
        print*,'Fehler beim Initialisieren der Routine'
        stop
      end if
1     open(1,file='such_feld.dat',status='old',access='append')
      print*,'Index-Werte (x,z), fuer die das Feld gesucht wird'
      print*,'wird (0.,0.) eingegeben, fragt das Programm nach '
      print*,'einer Geraden entlang derer die Feldwerte bestimmt'      
      print*,'werden'
      print*,'[x<=80,z<=50] Ort-Hall-Sonde [',xhall,zhall,']:'
      xg=xhall
      zg=zhall
      read*,xg,zg
      if(xg+zg.eq.0.) goto 10
      x=xg
      z=zg
      call su_such_feld(x,z,find)
      if(x.eq.0.)then
        print*,'Fehler beim Suchen des Feldes'
        stop
      end if
      print*,'Feldwert ist  ',find
      print*,'Abweichung vom Maximal-Feld in %o',
     &  (fmax-find)/fmax*1000.
      print*,'Ergebnis Max-Feld, Hall-Feld, Abweichung'
      print*,'in File such_feld.dat schreiben ? [ja=1]:'
      ja=1
      read*,ja
      if(ja.eq.1)write(1,100) fmax,find,(fmax-find)/fmax*1000.
 100  format(2f12.2,f8.2)
      goto 1
10    print*,'Eingabe von a,b in Index-Einheiten'
      print*,'Standard Messgerade fuer Hallsonde [',a,b,']:'
      read*,a,b
      print*,'Feldwerte bestimmen fuer Index_x, min,max,step'
      xmin=40.
      xmax=50.
      xstep=.1
      print*,'[',xmin,xmax,xstep,']:'
      read*,xmin,xmax,xstep
      do x=xmin,xmax,xstep
        z=x*a+b
        print*,'x,z  ',x,z
        call su_such_feld(x,z,find)
        print*,'x,z,find  ',x,z,find
        if(x.eq.0.)then
          print*,'Fehler beim Suchen des Feldes'
          stop
        end if
        print*,x,z,find
        write(25,*)x,z,find
      end do
      print*,'Ergebnis steht in fort.25 als x,z,Feldwert'
      goto 1
      end
