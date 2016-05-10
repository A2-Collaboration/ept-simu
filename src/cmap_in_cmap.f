c
c       c m a p _ i n _ c m a p . f
c
c      Das Programm kann in einem File cmap**.map die Indices und den Feldwert
c      verschieben oder umeichen. Es liest also Index-x, Index_y, Feld/Gauss
c      ein und kann die einzelnen Parameter mit einem Faktor mulitiplizieren
c      oder etwas addieren. Damit kann z.B. durch Multiplikation mit einem
c      Imkrement von Index auf Ort umgeeicht werden.
c      
c      Startfile ist immer ein cmap**.map mit einer Comment-Zeile am Anfang
c      in der steht: Zahl der Spalten, Zahl der Zeilen im Feld(i,j) und das
c      Inkrement Index in Ort/cm
c      
c      Das Ergebnis steht im File fort.15, es muss dann einen sinnvollen
c      anderen Namen erhalten (siehe auch endpoint.txt). Alle Elemente
c      im Ergebnisfile sind jetzt real!
c
c      Laden alss  g77 -o cmap_in_cmap.exe cmap_in_cmap.f
c
      implicit none
      real feld(200,200),inccm,sum(3),fakt(3),fneu,x,z,
     &  xmin,xmax,zmin,zmax,fmin,fmax
      integer i,j,inx,inz,ix,iz
      character file*20
c
      print*,'Startfile z.B. cmap10.map, a20 [**]:'
      read '(a)', file
      open(1,file=file,status='old')
      read(1,*)inx,inz,inccm
      fmin=1.E30
      fmax=-1.E30
      do i=1,inx
        do j=1,inz
          read(1,*,err=2,end=2)ix,iz,feld(ix,iz)
          if(i.ne.ix.or.j.ne.iz) goto 2
          if(feld(ix,iz).lt.fmin)fmin=feld(ix,iz)
          if(feld(ix,iz).gt.fmax)fmax=feld(ix,iz)
        end do
      end do
      print*,'Feld_min, Feld_max  ',fmin,fmax
      print*,'Das Inkrement Index in Ort/cm ist  ',inccm
      print*,'Die Elemente der Liste koennen geaendert werden durch'
      print*,'Element_neu = (Element_alt + sum) * Faktor'
      do i=1,3
        sum(i)=0.
        fakt(i)=1.
      end do
      print*,'3 Summanden  [',sum,' ]:'
      read*,sum
      print*,'3 Faktoren   [',fakt,']:'
      read*,fakt
      fmin=1.E30
      fmax=-1.E30
      xmin=1.E30
      xmax=-1.E30
      zmin=1.E30
      zmax=-1.E30
      do ix=1,inx
        do iz=1,inz
          x=(ix+sum(1))*fakt(1)
          z=(iz+sum(2))*fakt(2)
          fneu=(feld(ix,z)+sum(3))*fakt(3)
          if(x.lt.xmin)xmin=x
          if(x.gt.xmax)xmax=x
          if(z.lt.zmin)zmin=z
          if(z.gt.zmax)zmax=z
          if(fneu.lt.fmin)fmin=fneu
          if(fneu.gt.fmax)fmax=fneu
          write(15,*)x,z,fneu
        end do
      end do
      print*,'Minimal- und Maximal-Werte'
      print*,'x     ',xmin,xmax
      print*,'z     ',zmin,zmax
      print*,'feld  ',fmin,fmax
      print*,'Ergebnisfile ist fort.15, alle Werte jetzt real!'
      stop
2     print*,'Etwas falsch mit dem Indexing von Feld(i,j)'
      stop
      end
