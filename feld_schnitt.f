c
c      f e l d _ s c h n i t t . f
c
c      Programm benutzt Feldkarte cmap**.map und bestimmt
c      Schnitte Feld=F(Ort)
c      cmap benutzt Laborkoordinaten x(horizontal) und z(Strahlrichtung)
c      dort steht iz,ix,Wert/Gauss.
c      Ergebnis steht in fort.15
c      
c      Laden als:   g77 -o feld_schnitt.exe feld_schnitt.f
c
      implicit none
      real feld(300,300),felds(3,300),incr,x,z,ibx,ibz
      integer ix,iz,iicm,ixmax,izmax,izstart,izstop,
     &  ixstart,ixstop
      character mapfile*20
c
      ibx=-.06
      ibz=10.66
      print*,'map-file (z.B. cmap14.map, a20) [*** ***]:'
      read '(a)', mapfile
      open(1,file=mapfile,status='old')
      read(1,*)ixmax,izmax,incr
 1    read(1,*,err=2,end=2)ix,iz,feld(ix,iz)
      goto 1
 2    if(iz.ne.izmax.or.ix.ne.ixmax)then
        print*,'Fehler mit Anzahl der Werte im Mapfile'
        print*,'izmax soll, ist  ',izmax,iz,
     &    '  ixmax soll ist  ',ixmax,ix
        stop
      end if
c
      print*,'Ausgabe in Index = 0 oder cm =/= 0'
      print*,'dann relativ zum Bezugspunkt (Ibx,Ibz)=(.06,10.66) [0]:'
      iicm=0
      read*,iicm
      print*,'index z, min, max (step ist 1) [1,54]:'
      izstart=1
      izstop =54
      read*,izstart,izstop
      print*,'index x, min, max (step ist 1) [1,80]:'
      ixstart=1
      ixstop =80
      read*,ixstart,ixstop
c      
      do iz=izstart,izstop,1
        do ix=ixstart,ixstop,1
          x=ix
          z=iz
          if(iicm.ne.0)then
            x=(x-ibx)*incr
            z=(z-ibz)*incr
          end if
          write(15,*)x,z,feld(ix,iz)
        end do
      end do
      print*,'Das Ergebnis steht in fort.15'
      stop
      end
      
