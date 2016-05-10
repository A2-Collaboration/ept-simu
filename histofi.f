c      h i s t o f i . f
c
c      Programm macht aus einer Spalte einer bis
c      zu 6-spaltigen Zahlenliste
c      ein Histogramm, indem es die Liste einliest,
c      minimalen und maximalen Wert bestimmt, dann 
c      Rahmenwerte fuer ein Histogramm vorschlaegt,
c      auf deren Basis man dann das Histogramm
c      erstellt. Das Resultat steht in histofi.dat
c
c      laden als:
c      g77 -o histofi.exe histifi.f
c
      implicit none
      real hist(101),x(6),xmin,xmax,xm,dx,ddx,xlow,xhigh,xx
      integer i,j,l,ihi,isp,iaus
      character*20 file
c
      print*,'Erstellen eines Histograms aus einer Zahlenspalte'
      print*,'Startfile A20 [**]:'
      read '(a)',file
      open(1,file=file,status='old')
      print*,'File hat wieviele Splaten? (maximal 6)[1]:'
      isp=1
      read*,isp
      print*,'welche Spalte auswerten? [1]:'
      iaus=1
      read*,iaus
      xmin=1.E30
      xmax=-1.E30
      i=0
1     read(1,*,end=2)(x(j),j=1,isp)
      i=i+1
      if(x(iaus).le.xmin)xmin=x(iaus)
      if(x(iaus).ge.xmax)xmax=x(iaus)
      goto 1
2     close(1)
      print*,i,'Werte gefunden und eingelsen'
      print*,'Minimal- und Maximal-Wert',xmin,xmax
      xm=(xmin+xmax)/2.
      dx=(xmax-xmin)*1.05
      ddx=dx/101.
      xlow=xm-50.*ddx
      xhigh=xm+50.*ddx
      print*,'Vorschlag, low,high,step :',xlow,xhigh,ddx
      print*,'waehle Werte'
      read*,xlow,xhigh,ddx
      do i=1,101
        hist(i)=0.
      end do
      open(1,file='fort.10',status='old')
3     read(1,*,end=4)(x(j),j=1,isp)
      ihi=(x(iaus)-ddx/2.-xlow)/ddx+.5
      if(ihi.lt.1)ihi=1
      if(ihi.gt.101)ihi=101
      hist(ihi)=hist(ihi)+1.
      goto 3
4     close(1)
      open(1,file='histofi.dat')
      do i=1,101
        xx=(i-.5)*ddx+ddx/2.+xlow
        write(1,*)xx,hist(i)
      end do
      stop
      end
