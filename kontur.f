c
      character file*20
      print*,'Skalieren und schieben einer Kontur'
      print*,'z.B. polschuh.geo oder kartenrand.geo ...'
      print*,'File-Name a20 [*];'
      read '(a)', file
      read file '(a)'
      open(1,file=file,status='old')
3     print*,'xadd,yadd [*,*]:'
      read*,xadd,yadd
      open(2,file='kontur.geo')
1     read(1,*,end=2,err=2)x,y
      x=x+xadd
      y=y+yadd
      write(2,*)x,y
      goto 1
2     close(2)
      rewind(1)
      goto 3
      end
