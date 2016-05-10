c
c      a s y m a p s . f
c
c      Das Programm bestimmt die "Asymmetrie" zwischen zwei Feldkarten.
c      D.h., die Feldkarten werden eingelesen, dann wird Punkt fuer Punkt
c      die Asymmetrie gebildet und es werden verschiedene Werte in den
c      File asym.dat geschrieben.
c
c      Laden als:  g77 -o asymaps.exe asymaps.f
      implicit none
      real feld1(80,54),feld2(80,54),feldr(80,54),asmin,asmax,asym,
     &  h1,h2,h
      integer ix,iz,i,j,l,imin,imax,jmin,jmax,im,jm,ix1,iz1,ix2,iz2
      character*20 map1,map2
c
      print*,'Fledkarte 1 (z.B. omap14.map) [20A1]:'
      read '(a)', map1
      print*,'Feldkarte 2 [20A1]:'
      read '(a)', map2
      open(1,file=map1,status='old')
      read(1,*)ix1,iz1,h
      do i=1,ix1
        do j=1,iz1
          read(1,*)ix,iz,feld1(ix,iz)
          if(i.ne.ix.or.j.ne.iz) then
            print*,'Feldkarte 1'
            print*,'etwas falsch mit den Indices i,ix,j,iz',i,ix,j,iz
            stop
          end if
        end do
      end do
      close(1)
      open(1,file=map2,status='old')
      read(1,*)ix2,iz2,h
      if(ix1.ne.ix2.or.iz1.ne.iz2)then
        print*,'Felder passen nicht zusammen ix1,ix2,iz1,iz2'
        print*,ix1,ix2,iz1,iz2
        stop
      end if
      do i=1,ix2
        do j=1,iz2
          read(1,*)ix,iz,feld2(ix,iz)
          if(i.ne.ix.or.j.ne.iz) then
            print*,'Feldkarte 2'
            print*,'etwas falsch mit den Indices i,ix,j,iz',i,ix,j,iz
            stop
          end if
        end do
      end do
      close(1)
      im=ix1/2
      jm=iz1/2
      print*,'zentraler Wert bei [',im,jm,'] ?:'
      read*,im,jm
      h1=feld1(im,jm)-feld2(im,jm)
      h2=(feld1(im,jm)+feld2(im,jm))/2.
      open(1,file='asym.dat')
      asmin=1.E30
      asmax=-1.E30
      do i=1,ix1
        do j=1,iz1
          asym=(feld1(i,j)-feld2(i,j)-h1)/h2
          feldr(i,j)=asym
          write(1,*) i,j,abs(asym)
          if (asym.lt.asmin)then
            asmin=asym
            imin=i
            jmin=j
          end if
          if (asym.gt.asmax)then
            asmax=asym
            imax=i
            jmax=j
          end if
        end do
      end do
      print*,'Extremwerte :'
      print*,'Minimum  ',asmin,' bei ix,iz ',imin,jmin
      print*,'Maximun  ',asmax,' bei ix,iz ',imax,jmax
      print*,'zentrale Werte bei',im,jm,feld1(im,jm),feld2(im,jm)
      print*,'dort Differenz und Mittelw.',h1,h2
      print*,'Ergebnis in asym.dat als ix,iz,asym'
1     print*,'betrachten Einzelwerte bei ix,iz [*,*]:'
      read*,ix,iz
      print*,feld1(ix,iz),feld2(ix,iz),feld1(ix,iz)-feld2(ix,iz),
     &  feldr(ix,iz)
      goto 1
      stop
      end
      
