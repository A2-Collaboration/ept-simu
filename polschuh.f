c
      open(1,file='polschuh.geo',status='old')
3     print*,'xadd,yadd [*,*]:'
      read*,xadd,yadd
      open(2,file='ps.geo')
1     read(1,*,end=2,err=2)x,y
      x=x+xadd
      y=y+yadd
      write(2,*)x,y
      goto 1
2     close(2)
      rewind(1)
      goto 3
      end
