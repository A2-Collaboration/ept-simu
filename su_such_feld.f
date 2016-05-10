c
c      subroutine s u _ s u c h _ f e l d . f
c
c      Routine sucht Feldwerte in einer Feldkarte ****.map
c      mit Hilfe linearer Interpolation.
c
c      Die Routine wird initialisiert durch Eingabe
c      eines Wertes fuer find > 1.E5. Dann wird nach dem
c      Namen des Files mit der Feldkarte gefragt (mapfile).
c      Der File wird in das Feld feld(200,200) eingelsen.
c      In der ersten Zeile steht ixm,izm,incr.
c      ixm und izm besagen, dass die Indices des Feldes
c      von 1 bis ixm und 1 bis izm laufen. Beide muessen also
c      <=200 sein. incr ist das Inkrement von Index-Einheiten
c      zu cm. 1 Index =^= incr cm. Nach der Initialisierung
c      wird x=ixm, z=izm und find=Maximal-Feld gesetzt.
c      
c      Eingabe-Parameter:
c      (x,z) Ort/Index-Einheiten, fuer den das Feld gesucht wird
c      Ausgabe-Parameter:
c      find, der ermittelte Feldwert in Gauss
c
c      Treten Fehler auf, werden diese gemeldet und x=0. gesetzt.
c
      subroutine su_such_feld(x,z,find)
c
      implicit none
      real x,z,feld(400,400),find,find1,find2,feld1,feld2,
     &  fmin,fmax,incr,h
      integer ix,iz,ixm,izm,i,j
      character*20 mapfile
      save
c
      if(find.gt.1.E5) then
        fmin=1.E30
        fmax=-1.E30
        x=1.
        print*,'Name des Files (Feldkarte) ****.map [a20]:'
        read '(a)',mapfile
        open(1,file=mapfile,status='old')
        read(1,*)ixm,izm,incr
        print*,'Kennzeile mapfile:',ixm,izm,incr
        do i=1,ixm
          do j=1,izm
            read(1,*,err=2,end=2)ix,iz,feld(ix,iz)
            if(i.ne.ix.or.j.ne.iz) then
              print*,'Fehler bei Einlesen der Feldkarte'
              print*,'i,ix,j,iz',i,ix,j,iz
              x=0.
              return
            else
              if(feld(ix,iz).lt.fmin)fmin=feld(ix,iz)
              if(feld(ix,iz).gt.fmax)fmax=feld(ix,iz)
            end if
          end do
        end do
        do i=ixm+1,ixm+5
          do j=1,izm
            feld(i,j)=feld(ixm,j)
          end do
        end do
        ixm=ixm+5
        print*,'Feldkarte erweitert auf Index-x = ',ixm
        print*,'Min-, Max-Feldwerte/Gauss', fmin,fmax
        x=ixm
        z=incr
        find=fmax
        return
2       print*,'Fehler beim Einlesen des Feldes'
        x=0.
        return
c        Ab hier Interpolation des Feldwertes
        else
          ix=x
          iz=z
          if(ix.lt.1.or.ix.gt.ixm.or.iz.lt.1.or.iz.gt.izm) then
c            print*,'aus su_such_feld'
c            print*,'Punkt ausserhalb des Feldes  '
c            print*,'ix,iz,x,z',ix,iz,x,z
            x=0.
            return
          else
            feld1=feld(ix,iz)+(x-ix)*(feld(ix+1,iz)-feld(ix,iz))
            feld2=feld(ix,iz+1)+(x-ix)*(feld(ix+1,iz+1)-feld(ix,iz+1))
            find1=feld1+(z-iz)*(feld2-feld1)
            feld1=feld(ix,iz)+(z-iz)*(feld(ix,iz+1)-feld(ix,iz))
            feld2=feld(ix+1,iz)+(z-iz)*(feld(ix+1,iz+1)-feld(ix+1,iz))
            find2=feld1+(x-ix)*(feld2-feld1)
            find=(find1+find2)/2.
            return
          end if
        end if
        end
        
