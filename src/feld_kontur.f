c
c      f e l d _ k o n t u r . f
c
c      Programm fasst die Programme
c          feld_kontur_punkte.f  und  feld_kontur_linie.f     zusammen.
c      Es liefert ebenso die Ergebnisfiles fort.11 und fort.12
c
c                      feld_kontur_punkte.f      
c      Programm sucht in mit convert_map.f aus map.map erstellten
c      File omap**.map (umgedrehte cmap)
c      eine Daten-Liste fuer Punkte mit gleichem Feldwert.
c      Hierzu werden solche Koordinaten aufsteigend und absteigend
c      entlang fester x- und z-Werte gesucht und im File fort.11
c      abgelegt. Damit sind diese Werte dann  !!nicht!! entlang
c      Linien geordnet sondern in der Reihenfolge der Findung.
c
c                      feld_kontur_linie.f
c      Setzt die ungeordnete Punktefolge aus feld_kontur.f, fort.11,
c      in Linien mit festem Feld um. Das Ergebnis steht dann in fort.12
c      
c      Laden als   g77 -o feld_kontur.exe feld_kontur.f
c      


      implicit none
      real feld(200,200),xkontur(400),zkontur(400),wert,xs,zs,
     &  xinc,zk(500),xk(500),dimi,dist
      integer ix,iz,i,j,inz,inx,l,ii,imi
      character file*20
c      
      print*,'Start-File (z.B. omap14.map [a20]:'
      read '(a)', file
      open(1,file=file,status='old')
      print*,'Suche Kontur fuer Wert/Gauss [10000]:'
      wert=10000
      read*,wert
      read(1,*)inx,inz,xinc
      do i=1,inx
        do j=1,inz
          read(1,*)ix,iz,feld(ix,iz)
          if(i-ix+j-iz.ne.0)print*,'i,j,ix,iz',i,j,ix,iz
        end do
      end do
      l=0
      do iz=1,inz
        do ix=1,inx-1
          if(feld(ix,iz).le.wert.and.feld(ix,iz+1).ge.wert)then
            zs=iz+(feld(ix,iz)-wert)/(feld(ix,iz)-feld(ix,iz+1))
            l=l+1
            if(l.gt.400)stop
            zkontur(l)=zs
            xkontur(l)=ix
            write(11,*)xkontur(l),zkontur(l)
          end if
          if(feld(ix,iz).ge.wert.and.feld(ix,iz+1).le.wert)then
            zs=iz+(feld(ix,iz)-wert)/(feld(ix,iz)-feld(ix,iz+1))
            l=l+1
            if(l.gt.400)stop
            zkontur(l)=zs
            xkontur(l)=ix
            write(11,*)xkontur(l),zkontur(l)
          end if
        end do
      end do
      do ix=1,inx
        do iz=1,inz-1
          if(feld(ix,iz).le.wert.and.feld(ix+1,iz).ge.wert)then
            xs=ix+(feld(ix,iz)-wert)/(feld(ix,iz)-feld(ix+1,iz))
            l=l+1
            if(l.gt.400)stop
            xkontur(l)=xs
            zkontur(l)=iz
            write(11,*)xkontur(l),zkontur(l)
          end if
          if(feld(ix,iz).ge.wert.and.feld(ix+1,iz).le.wert)then
            xs=ix+(feld(ix,iz)-wert)/(feld(ix,iz)-feld(ix+1,iz))
            l=l+1
            if(l.gt.400)stop
            xkontur(l)=xs
            zkontur(l)=iz
            write(11,*)xkontur(l),zkontur(l)
          end if
        end do
      end do
2     ii=l
      print*,ii,'  Punkte gefunden'
c      Jetzt Punkte zur Linie ordnen
      xk(1)=xkontur(1)
      zk(1)=zkontur(1)
      dimi=1.E30
      xkontur(1)=0.
      zkontur(1)=0.
      j=1
4     l=0
      do i=1,ii
        if(zkontur(i).ne.0.)then
          l=l+1
          dist=sqrt((zk(j)-zkontur(i))**2+(xk(j)-xkontur(i))**2)
          if(dist.lt.dimi)then
            dimi=dist
            imi=i
          end if
        end if
      end do
      if(l.eq.0)goto 3
      j=j+1
      if(j.gt.499)then
        print*,'Ziel-Zahl zu hoch  ',j
        stop
      end if
      zk(j)=zkontur(imi)
      xk(j)=xkontur(imi)
      write(12,*)xk(j),zk(j)
      xkontur(imi)=0.
      zkontur(imi)=0.
      dimi=1.E30
      goto 4
3     write(12,*)xk(1),zk(1)
      write(12,*)xk(2),zk(2)
      print*,'die Ergebnisse stehen'
      print*,'als Punkte in fort.11, als Linien in fort.12'
      stop
      end
