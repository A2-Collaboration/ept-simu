c
c      i n _ s u _ b r e m s _ m c . f
c
c      Programm erstellt file su_brems_mc.geo, der bei der Initialisierung
c      der Routine su_brems_mc.f eingelsen wird.
c
      implicit none
      real emx,emy,sx,stx,sy,sty,rk,lk,t0,dr,z,v
c
      emx=9.E-3
      emy=4.E-4
      v=1
      rk=.25
      lk=250.
      dr=1.E-4
      t0=1.44
      z=28.
      print*,'Emittanzen in x (hor.) und y (vert.) in mm*mrad'
      print*,'[',emx,emy,']:'
      read*, emx,emy
      print*,'Aufteilen in Verhaeltnis Ort/Winkel [1]:'
      read*,v
      stx=sqrt(emx/v)
      sx=stx*v
      sty=sqrt(emy/v)
      sy=sty*v
      print*,'Radius Photonen-Kolli, Abstand Radiator-Kollimator'
      print*,'in cm [',rk,lk,']:'
      read*,rk,lk
      print*,'fuer die Vielfachstreuung im Radiator'
      print*,'Radiatordicke, Strahlungslaenge in cm'
      print*,'fuer Ni [',dr,t0,']:'
      read*,dr,t0
      print*,'Dicke entspricht ',dr/t0,' Strahlungslaengen'
      print*,'Z des Radiators [',z,']:'
      read*,z
      open(1,file='su_brems_mc.geo',status='new')
      write(1,*)sx*.1,stx*.001,sy*.1,sty*.001
      write(1,*)rk,lk,t0,dr,z
      close(1)
      stop
      end
