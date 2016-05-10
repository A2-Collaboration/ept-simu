c
c      Erstellt einen file raster.dat mit Koordinaten, die
c      die ein Raster mit vorgegebenem Inkrement geben
c
      implicit none
      real incr,xmin,xmax,ymin,ymax,x,y
c
      print*,'Feldgroesse in x und y'
      print*,'xmin,xmax [*,*]:'
      read*,xmin,xmax
      print*,'ymin,ymax [*,*]:'
      read*,ymin,ymax
      print*,'Schrittweite (Inkrement)*]:'
      read*,incr
      open (1,file='raster.dat')
      do x=xmin,xmax,incr
        do y=ymin,ymax,incr
          write(1,*)x,y
        end do
      end do
      stop
      end
      
