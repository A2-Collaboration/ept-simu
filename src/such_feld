c
c       t e s t _ s u _ s u c h _ f e l d . f
c
c      Programm testet oder nutzt die Routine su_such_feld
c
c      Nach Initialisierung der Routine wird nach Eingabe eines
c      Punktes in der Feldkarte der Feldwert durch lineare
c      Interpolatin bestimmt.
c
c      Laden mit:
c      g77 -fbounds-check -o test_su_such_feld.exe
c        test_su_such_feld.f su_such_feld.f
c
      implicit none
      real xl,x,z,find,xg,zg
      integer it
c
      print*,'test der Routine su_such_feld.f'
      find=1.E30   ! Routine wird initialisiert
      call su_such_feld(x,z,find)
      if(x.eq.0.)then
        print*,'Fehler beim Initialisieren der Routine'
        stop
      end if
1     print*,'Index-Werte (x,z), fuer die das Feld gesucht wird'
      print*,'[x<=80,z<=50]:'
      read*,xg,zg
      x=xg
      z=zg
      call su_such_feld(x,z,find)
      if(x.eq.0.)then
        print*,'Fehler beim Suchen des Feldes'
        stop
      end if
      print*,'Feldwert ist  ',find
      goto 1
      end
