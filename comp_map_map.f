c
c      c o m p _ m a p _ m a p . f
c
c      Programm vergleicht Feldkarten, bestimmt Mittlewerte aus 
c      waehlbaren Segmenten (x,z)_min (x,z)_max, Streuung der
c      Elemente im Segment, Eichfaktor zwischen beiden Feldkarten
c      und dann die Streuung der Unterschiede zwischen den Feldkarten,
c      vofuer groessere Segmente gewaehlt werden koennen.
c
c      Laden als  g77 -o comp_map_map.exe comp_map_map.f
c      
      implicit none
      real x,z,feld(2,400,400),sum(2),isum(2),incr,sig(2),
     &  mwert(2),eich,ssum,hist(101),dmin,dmax,d,dhi,his,dhi0,
     &  diff(400,400)
      integer i,j,l,ixm,izm,isxmin,isxmax,iszmin,iszmax,ix(2),iz(2),
     &  ipxmin,ipxmax,ipzmin,ipzmax,jx,jz,ihi,
     &  ixmax,ixmin,izmax,izmin
      character*20 map1,map2
c
      print*,'Vergleich zwischen Feldkarten'
      print*,'Die Eichung gleicht der Bezugsfeldkarte an'
      print*,'Name der Bezugs-Feldkarte ****.map [a20]:'
      read '(a)',map1
      print*,'Name der Vergleichs-Feldkarte      [a20]:'
      read '(a)',map2
      open(1,file=map1,status='old')
      open(2,file=map2,status='old')
      read(1,*)ix(1),iz(1),incr
      read(2,*)ix(2),iz(2),incr
      if(ix(1).ne.ix(2).or.iz(1).ne.iz(2))then
        print*,'Feldkarten haben verschiedene Dimensionen'
        stop
      end if
      print*,'Kennzeile mapfiles:',ix(1),iz(1),incr
      do l=1,2
        do i=1,ix(l)
          do j=1,iz(l)
            read(l,*)jx,jz,feld(l,i,j)
            if(jx.ne.i.or.jz.ne.j)then
              print*,'etwas falsch mit den Feldern',l,i,jx,j,jz
              stop
            end if
          end do
        end do
        sum(l)=0.
        isum(l)=0.
      end do
      isxmin=ix(1)/2-1
      isxmax=ix(1)/2+1
      iszmin=iz(1)/2-1
      iszmax=iz(1)/2+1
      print*,'Segment, ueber das gemittelt werden soll'
      print*,'xmin,xmax [',isxmin,isxmax,']:'
      read*,isxmin,isxmax
      print*,'zmin,zmax [',iszmin,iszmax,']:'
      read*,iszmin,iszmax
      do l=1,2
        do i=isxmin,isxmax
          do j=iszmin,iszmax
            sum(l)=sum(l)+feld(l,i,j)
            isum(l)=isum(l)+1.
          end do
        end do
        mwert(l)=sum(l)/isum(l)
        sum(l)=0.
      end do
      do l=1,2
        do i=isxmin,isxmax
          do j=iszmin,iszmax
            sum(l)=sum(l)+(mwert(l)-feld(l,i,j))**2
            ssum=ssum+feld(l,i,j)
          end do
        end do
        sig(l)=sqrt(sum(l)/isum(l))
      end do
      print*,'Feld-Mittelwerte, sigma-Werte',mwert,sig
      eich=mwert(1)/mwert(2)
      ipxmin=2
      ipxmax=ix(1)-1
      ipzmin=2
      ipzmax=iz(1)-1      
      print*,'Differenz, Mittelwert ',mwert(1)-mwert(2),
     &  (mwert(1)+mwert(2))/2.
      print*,'Eichfaktor [',eich,']:'
      read*,eich
      print*,'Segment, ueber das geprueft werden soll'
      print*,'xpmin,xpmax [',ipxmin,ipxmax,']:'
      read*,ipxmin,ipxmax
      print*,'zpmin,zpmax [',ipzmin,ipzmax,']:'
      read*,ipzmin,ipzmax
      sum(1)=0.
      isum(1)=0.
      ssum=0.
      dmin=1.E30
      dmax=-1.E30
      ixmin=100
      ixmax=-10
      izmin=100
      izmax=-10
      do i=ipxmin,ipxmax
        do j=ipzmin,ipzmax
          d=feld(1,i,j)-eich*feld(2,i,j)
          diff(i,j)=d/(feld(1,i,j)+feld(2,i,j))/2.
          if(d.le.dmin)then
            ixmin=i
            izmin=j
            dmin=d
          end if
          if(d.ge.dmax)then
            ixmax=i
            izmax=j
            dmax=d
          end if
          sum(1)=sum(1)+d
          isum(1)=isum(1)+1.
          ssum=ssum+feld(1,i,j)+eich*feld(2,i,j)
        end do
      end do
      mwert(1)=sum(1)/isum(1)
      mwert(2)=ssum/2./isum(1)
      print*,'Mittelwert, mittlere Abweichung, relativ '
      print*,mwert(2),mwert(1),mwert(1)/mwert(2)
      print*,'Minimal-Wert',dmin,' bei',ixmin,izmin
      print*,'Maximal-Wert',dmax,' bei',ixmax,izmax
      dhi=abs(dmin)
      if(d.lt.abs(dmax))dhi=abs(dmax)
      dhi0=dhi/50.     ! Schrittweite Histogramm
      dhi=5.
      print*,'Vorschlag Schrittweite Histogramm ',dhi0
      print*,'Vorgabe   Schrittweite Histogramm [',dhi,']:'
      read*,dhi
      sum(2)=0.
      sum(1)=0.
      isum(2)=0.
      do i=1,101
        hist(i)=0.
      end do
      do i=ipxmin,ipxmax
        do j=ipzmin,ipzmax
          sum(2)=sum(2)+(mwert(1)-feld(1,i,j)+eich*feld(2,i,j))**2
          isum(2)=isum(2)+1.
          sum(1)=sum(1)+(mwert(2)-feld(1,i,j)+eich*feld(2,i,j))**2
          ihi=(feld(1,i,j)-eich*feld(2,i,j))/dhi+50.51
          if(ihi.lt.1)ihi=1
          if(ihi.gt.101)ihi=101
          hist(ihi)=hist(ihi)+1.
        end do
      end do
      do i=1,101
        write(31,*) dhi*(-51.+i),hist(i)
      end do
      sig(2)=sqrt(sum(2)/isum(2))
      sig(1)=sqrt(sum(1)/isum(2))
      print*,'sigma - Verteilung Differenzen ',sig(2)
      print*,'sigma-Verteilung Werte         ',sig(1)
      print*,'Histogramm in fort.31'
      do i=ipxmin,ipxmax
        do j=ipzmin,ipzmax
          write(32,*)i,j,diff(i,j)
        end do
      end do
      print*,'Feld mit rel. Differenzen in fort.32'
      stop
      end
      
