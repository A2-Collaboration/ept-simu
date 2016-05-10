c
c      b l d f i . f
c
c      Programm erstellt einen 3-spaltigen File
c      aus einem Start-File mit 1 bis mehreren Spalten
c      dabei kann man Konstanten an die einzelnen Elemente
c      multiplizieren, dividieren, addieren
c      es koennen Funktionen auf die Elemente angewendet werden
c      Ziel-File ist bldfi.dat
c
c      Laden mit
c      g77 -fbounds-check -o bldfi.exe bldfi.f
c
      implicit none
      real x(20),y(10),fv(10),tv(10),sv(10),fn(10),tn(10),
     &  sn(10),a(10),hilf,pot,xmin,xmax,ymin,ymax,
     &  xanf,xstep,h1
      character*30 plfi
      integer i,isp,nri,ind(10),jf(10),jafe,jsf
c
      print*,'erstellen von Plotfiles  bldfi.dat mit x,y,dy'
      print*, 'aus Files mit einer oder mehr Spalten'
      print*, 'fuer die einezelnen Elemente x,y,dy koennen'
      print*, 'Faktoren, Teiler, Summanden, Funktionen'
      print*, 'angegeben werden.'
      print*, 'Start-Filename (A30):'
      read '(a)', plfi
      print 101, plfi
 101  format(1x,30a1)
      open(1,file=plfi,status='old')
      do i=1,4
        read(1,102)x
        write(0,103)x
      end do
 102  format(20a4)
 103  format(1x,20a4)
      rewind 1
      print*,'wie viele Spalten ? [3]:'
      isp=3
      read*,isp
      print*,'wieviele Spalten in den File bldfi.dat ? '
      print*,'normal 3 Spalten als x, y, dy'
      nri=3
      read*,nri
      print*,'welche Spalten in file bldfi.dat?'
      print*,'bei 0 fuer x nach Massgabe '
      print*,'bei 0 fuer dy wird 1. geschrieben'
      print*,'bei gleichem Index fuer y und dy wird abgefragt,'
      print*,'ob sqrt(y) geschrieben werden soll'
      do i=1,10
        ind(i)=i
        fv(i)=1.
        tv(i)=1.
        sv(i)=0.
        fn(i)=1.
        tn(i)=1.
        sn(i)=0.
        jf(i)=0
        a(i)=1.
      end do
      print*,'Indices [',(ind(i),i=1,nri),']:'
      read*,(ind(i),i=1,nri)
      if(ind(2).ne.ind(3))goto 70
      print*,'Fehler sei sqrt(y) ? [ja=1]:'
      jafe=1
      read*,jafe
 70   if(ind(1).ne.0)goto 31
      print*,'xstart,xstep [1,1]:'
      xanf=1.
      xstep=1.
      read*,xanf,xstep
      hilf=xanf
 31   print*,'3 Funktionen  nein=0,log=1,exp=2,**a=3,abs=4 [0,0,0]:'
      read*,(jf(i),i=1,nri)
      jsf=0
      do 13 i=1,nri
        jsf=jsf+jf(i)
        if(jf(I).ne.3)goto 13
        print*,'fuer Element',i,'Exponent [1.]:'
        read*, a(i)
 13   continue
      if(jsf.ne.0)print*,'vor Funktionsbildung'
      print*,' Faktoren ? [',(fv(i),i=1,nri),']:'
      read*,(fv(i),i=1,nri)
      print*,'3 Teiler ? [',(tv(i),i=1,nri),']:'
      read*,(tv(i),i=1,nri)
      print*,'3 Summanden ? [',(sv(i),i=1,nri),']:'
      read*,(sv(i),i=1,nri)
      if(jsf.eq.0)goto 14
      print*,'nach Funktionsbildung'
      print*,'3 Faktoren ? [',(fn(i),i=1,nri),']:'
      read*,(fn(i),i=1,nri)
      print*,'3 Teiler ? [',(tn(i),i=1,nri),']:'
      read*,(tn(i),i=1,nri)
      print*,'3 Summanden ? [',(sn(i),i=1,nri),']:'
      read*,(sn(i),i=1,nri)
 14   print*,'mit Potenz der 1. Zielspalte multiplizieren [0] :'
      pot=0.
      read*,pot
      xmin=1.E30
      xmax=-1.E30
      ymin=1.E30
      ymax=-1.E30
      open(2,name='bldfi.dat')
 1    read(1,*,end=10)(x(i),i=1,isp)
      do i=1,nri
        if(ind(i).ne.0)y(i)=x(ind(i))
      end do
      if(ind(1).eq.0.)y(1)=hilf
      hilf=hilf+xstep
      if(y(1).eq.99999.)goto 5
      do 16 i=1,nri
             y(i)=y(i)*fv(i)/tv(i)+sv(i)
             goto(21,22,23,24),jf(i)
17           y(i)=y(i)*fn(i)/tn(i)+sn(i)
             goto 16
 21          y(i)=log(y(i))
             goto 17
 22          y(i)=exp(y(i))
             goto 17
 23          h1=1.
             if(a(i).ne.0.)h1=y(i)**a(i)
             y(i)=h1
             goto 17
 24          y(i)=abs(y(i))
             goto 17
 16        continue
           if(ind(3).eq.0)goto 3
           if(ind(2).eq.ind(3).and.jafe.eq.1)y(3)=sqrt(y(2))
           goto 2
 3         y(3)=1.
 2         if(pot.eq.0.)goto 30
           y(2)=y(2)*y(1)**pot
           y(3)=y(3)*y(1)**pot
 30        if(y(1).lt.xmin)xmin=y(1)
           if(y(1).gt.xmax)xmax=y(1)
           if(y(2).lt.ymin)ymin=y(2)
           if(y(2).gt.ymax)ymax=y(2)
 5         write(2,*)(y(i),i=1,nri)
           goto 1
 10        close(1)
           close(2)
           print*,'X-Werte zwischen',xmin,' und',xmax
           print*,'Y-Werte zwischen',ymin,' und',ymax
           stop
         end
         
