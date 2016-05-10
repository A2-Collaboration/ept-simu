c
c      m a p _ o m a p . f
c
c      Die gemessenen Feldkarten werden in Feldern map**.map
c      abgelegt. Dieses Programm liest Karten map**.map ein, die durch
c      entsprechendes Umkopieren von den Original-Karten z.B. map.map
c      oder, wie sie aus der Messung hervorgehen,
c      erhalten werden. Die Feldkarte wird
c      umgeschrieben in praktischere Koordinaten.
c      Die Feldwerte werden in ein Feld(Indexx,Indexz)
c      eingeordnet. Der Ergebnisfile ist z.B. omap**.map
c      (der Name kann frei gewaehlt werden), in dem die Achsen
C      entsprechend dem Labor liegen.
c      In map.map steht urspruenglich:
c      3 Kommentar-Zeilen: Filename bei Messung
c                          Datum: Start der Messung
c                          Zeit : Ende der Messung
c      Dann folgen Zeilen mit Index1, Index2, Feldwert/Gauss
c      wobei die Indices die Inkrementschritte der Koordinatenmaschine
c      ab einem maschineneigenen Nullpunkt zaehlen.
c      Ein Inkrementschritt der Koordinatenmaschine entspricht 0.2211 mm
c      (hier in data abgelegt).
c      Die erste Achse liegt im Labor _|_ zur Eintrittskante des
c      Endpoint-Taggers, entspricht dort also der z-Achse,
c      sie zaehlt entgegengesetzt der Strahlrichtung.
c      Die zweite Achse entspricht im Labor der x-Achse, || Eintrittskante,
c      sie zaehlt in Richtung der Ablenkung. Diese Richtung ist 
c      entgegengesetzt zur Laborrichtung (Strahleintritt in den ET rechts)
c      wird aber aus Darstellungsgruenden so beibehalten,
c      d.h. der Enpointtagger lenkt in den Rechnungen nach rechts ab,
c      in der praktischen Installation aber nach liks, wir der
c      Main Tagger, MT.
c      Es treten in manchen Original-Listen Anfangs- und End-Feld-Werte
c      ohne Bedeutung auf, oder es fehlt der letzte Wert.
c      Diese Werte sind zwar unwichtig, werden aber durch lineare
c      Interpolation bestimmt.
c
c      Laden als  g77 -fbounds-check -o map_omap.exe map_omap.f
c      
      implicit none
      real increment,pi,x0,y0,amp,x1,y1,f,dxy,fmin,fmax,xmin,xmax,
     &  ymin,ymax,xn,yn,feld(200,200),inccm,xe,ye 
      integer l,i,j,iwert,inx,iny,ix,iy
      character*25 comment,file_ein,file_aus
      data increment/0.02211/,pi/3.1415926/
      print*,'*********** Programm   convert_map.f ***********'
      print*,'Abmessungen in cm, Energie in MeV, Winkel Grad'
      print*,'Programm liest einen original Mess-File map.map'
      print*,'setzt ihn nach Wahl eines Anfangspunktes als'
      print*,'Index (1,1)'
      print*,'in Index-x, -z, Feld/Gauss um. In der Kommentarzeile'
      print*,'des Ergebnisses stehen die Maximal-Werte'
      print*,'fuer die Indices und die Konversion Index in Laenge/cm'
      print*,'Der Ergenisfile ist fort.10'
c
      do i=1,200
        do j=1,200
          feld(i,j)=0.
        end do
      end do
c
      print*,'Start(original)-File z.B. map14.map [A25]:'
      read '(a)',file_ein
      open(1,file=file_ein,status='old')
      iwert=2
      do l=1,3
        read(1,'(a)') comment
        print '(1x,a25)',comment
      end do
      read(1,*) y0,x0,f
      print*,'Anfangspunkt  ',x0,y0
      xmin=x0
      xmax=x0
      ymin=y0
      ymax=y0
c
1     read(1,*,err=2,end=2)y1,x1,f
      if(iwert.eq.2) then
        xe=x1
        ye=y1
        dxy=x1-x0
        inccm=dxy*increment
        print*,'Schrittweite/cm  ',inccm
        fmin=f
        fmax=f
      end if
      if(x1.lt.xmin)xmin=x1
      if(x1.gt.xmax)xmax=x1
      if(y1.lt.ymin)ymin=y1
      if(y1.gt.ymax)ymax=y1
      if(f.lt.fmin)fmin=f
      if(f.gt.fmax)fmax=f
      iny=(ymax-ymin)/dxy+1.5
      inx=(xmax-xmin)/dxy+1.5
      iwert=iwert+1
      goto 1
2     print*, iwert  ,'  Werte eingelesen'
      print*,'Endpunkt bei  ',x1,y1
      print*,'Minimal- und Maximal-Werte'
      inx=(xmax-xmin)/dxy+1.1
      iny=(ymax-ymin)/dxy+1.1
      print*,'x     ',xmin,xmax, '   also ',inx  ,' Werte'
      print*,'y     ',ymin,ymax, '   also ',iny  ,' Werte'
      print*,'feld  ',fmin,fmax
      close(1)
c
      print*,'Ausgabe-File z.B. omap14.map [A25]:'
      read '(a)', file_aus
      open(2,file=file_aus)
      yn=ymin
      xn=xmax
      print*,'Punkt (1,1) als xmax,ymin zum Anfangspunkt im Labor ?'
      print*,'[',xn,yn,']:'
      read*,xn,yn
      print*,'Punkt (1,1) jetzt bei  ',xn,yn
c      
      open(1,file=file_ein,status='old')
      do l=1,3
        read(1,'(a)')comment
      end do
      iwert=1
3     read(1,*,err=4,end=4)y1,x1,f
      ix=(xn-x1)/dxy+1.5    ! da xn > x1
      iy=(y1-yn)/dxy+1.5    ! da yn < y1
c      ix=inx-ix+1   ! Aenderung Laufrichtung der Achsen
      iy=iny-iy+1    ! Aenderung Laufrichtung der Achse
c      print*,'x1,xn,ix,y1,yn,iy',x1,xn,ix,y1,yn,iy
      feld(ix,iy)=f
      iwert=iwert+1
      goto 3
4     write(2,*) inx,iny,inccm
      print*,'inx,iny,inccm',inx,iny,inccm
      feld(inx,1)=((2.*feld(inx,2)-feld(inx,3))+
     &  (2.*feld(inx-1,1)-feld(inx-2,1)))/2.
      feld(inx,iny)=((2.*feld(inx-1,iny)-feld(inx-2,iny))+
     &  (2.*feld(inx,iny-1)-feld(inx,iny-2)))/2.
      do ix=1,inx
        do iy=1,iny
          write(2,*)ix,iy,feld(ix,iy)
        end do
      end do
      stop
      end
