c
c      e t _ k m _ m t . f
c
c      Programm prueft, wie der Strahl durch
c      Endpoit-Tagger- ET, Korrektur-Magnet- KM, Main-Tagget- MT,
c      zum Strahlfaenger gelangt.
c      Rechnung mit Kreisen und Geraden.
c      MT hat den Eintritt bei (x,z)=(0,0), er lenkt in x-Richtung,
c      wie auch MT, ab.
c      Die z-Richtung ist die Richtung des einlaufenden Strahls.
c      ET und KM werden duch jeweils zwei z-Werte gegeben,
c      die den Eintritt und Austritt angeben. Die Ablenkung wird durch
c      Radien angegeben, die alle positiv anzusetzen sind, obwohl
c      KM zuruecklenkt!
c      Der Soll-Eintritt in ET ist bei x=0, in MT bei (z,x)=(0,0)
c
c      Die Geoemetrieangaben stehen in einem File etkmmt.geo
c      alle Laengen in cm
c      als ET: z-Eintritt und Austritt
c          KM: z-Eintritt und Austritt
c          MT: Eintrittskante: xe,ze,re  (Mittelpunkt und Radius)
c              Austrittskante: xa,za,ra (Mittelpunkt und Radius)
c              Austritt ist regulaer bei xmta0=225
c          Ablenkradien: ret,rkm,rmt  (alle positiv angeben)
c
c      Das Ergebnis steht in et_km_mt.dat      
c
c      Laden als
c      g77 -o et_km_mt.exe et_km_mt.f sub_g_k.f
c      
      implicit none
      real xete,xeta,xkme,xkma,xmte,xmta,xe,ze,re,xa,za,ra,
     &  zete,zeta,zkme,zkma,ret,rkm,rmt,wet,wi,xet,zet,
     &  dz1,dz2,dx1,dx2,fi,xkm,zkm,a,b,vz,t,zm,xm,bahn,
     &  alp,zmta,zmte,pi,al0,a0,b0,xmta0,rmt0,zk,xk,al,d,
     &  db,zmta0,xs,zs,xg,zg(3),del,h1,a0s,b0s,z0,x0,z,x,
     &  zmte0,xmte0
c
c      der Sollstrahl hat die Koordinaten
c      fuer z<0  stets x=0.
c      fuer z>0 Kreis mit zm,xm,rm  0.,280.,280. bis x=225.
c      Ablenkwinkel alpha=78.672 Grad
c      dann Gerade x=z*a+b mit a,b  tan(alpha)=4.992, -572.48
c      bahn ist die Laenge der Bahn des Sollstrahls.
c      Gegen diese Groesse sind Werte bezogen auf den
c      Sollstrahl aufzutragen (z.B. Winkel, Abstand vom Sollstr.)
c      
      pi=3.1415926
      rmt0=280.         ! Parameter Sollstrahl
      zmte0=0.
      xmte0=0.
      xmta0=225.
      al0=acos((rmt0-xmta0)/rmt0)
      zmta0=rmt0*sin(al0)
      a0=tan(al0)
      b0=xmta0-a0*zmta0
      print*,'alpha0/Grad,a0,b0,zmta0',al0/pi*180.,a0,b0,zmta0
c      
      open(1,file='etkmmt.geo',status='old')
      read(1,*)zete,zeta
      read(1,*)zkme,zkma
      read(1,*)ze,xe,re
      read(1,*)za,xa,ra
      read(1,*)ret,rkm,rmt
      close(1)
      print*,'Laengen in cm, Winkel in Grad'
      print*,'Radien in ET, KM, MT'
      print*,'unendlich, also kein Feld ist r>1.E30'
      print*,'alle Radien postiv angeben! '
      print*,'ret,rkm,rmt [',ret,rkm,rmt,']:'
      read*,ret,rkm,rmt
c
      open(1,file='et_km_mt.dat')
c      Strahl durch ET
      xete=0.
      wet=asin((zeta-zete)/ret)
      xeta=ret*(1.-cos(wet))
      xkme=xeta+(zkme-zeta)*tan(wet)
c      Strahlverlauf notieren
      write(1,100)-200.,0.,0.,-200.,0.
      write(1,100)zete,xete,0.,zete,xete
      do wi=0,wet,wet/30.
        xet=ret*(1.-cos(wi))
        zet=zete+ret*sin(wi)
        write(1,100)zet,xet,wi/pi*180.,zet,xet
      end do
      write(1,100)zeta,xeta,wet/pi*180.,zeta,xeta
c      Strahl durch KM
      dz1=rkm*sin(wet)
      dx1=rkm*cos(wet)-xkme
      dz2=zkma-zkme-dz1
      print*,'dz1,dx1,dx2',dz1,dx1,dx2
      fi=asin(dz2/rkm)
      xkma=rkm*cos(fi)-dx1
      xmte=xkma+zkma*tan(fi)  ! dies gilt nur fuer Sollstrahl
c      Strahlverlauf notieren
      write(1,100)zkme,xkme,wet/pi*180.,zkme,xkme
      do wi=wet,-fi,-(wet+fi)/30.
        xkm=rkm*cos(wi)-dx1
        zkm=zkme+dz1-rkm*sin(wi)
        write(1,100)zkm,xkm,wi/pi*180.,zkm,xkm
      end do
      write(1,100)zkma,xkma,-fi/pi*180.,zkma,xkma
c      Strahl durch MT
      a=tan(-fi)
      b=xmte
      vz=-1.
      call cut_g_k(a,b,ze,xe,re,vz,zmte,xmte,t)
      print*,'zmte,xmte',zmte,xmte
      write(1,100)zmte,xmte,-fi/pi*180.,zmte,xmte
      zm=zmte+rmt*sin(fi)
      xm=xmte+rmt*cos(fi)
      xmta=xmta0
      del=asin((xm-xmta)/rmt)
      alp=pi/2-del+fi
      zmta=rmt*cos(del)+zm
      print*,'zm,xm,zmta,xmta',zm,xm,zmta,xmta
      bahn=0.
c      Bahnverfolgung durch MT
      do wi=-fi,alp-fi,alp/75.
        zk=rmt*sin(wi)+zm
        xk=-rmt*cos(wi)+xm
        if(xk.lt.xmta0)then
          h1=sqrt(zk**2+(rmt0-xk)**2) 
          d=rmt0-h1
          zs=zk*rmt0/h1
          bahn=rmt0*asin(zs/rmt0)
          write(1,100)zk,xk,wi/pi*180.,bahn,d
        end if
      end do
      bahn=bahn+sqrt((zmta-zk)**2+(xmta0-xk)**2)
      a=tan(alp-fi)
      b=xmta-zmta*a
      a0s=-1./a0
      do x0=xmta0,700.,(700.-xmta)/10.
        z0=(x0-b0)/a0
        b0s=z0*(a0-a0s)+b0   ! Linie senkrecht zu Sollstrahl
        z=(b0s-b)/(a-a0s)
        x=z*a+b
        d=sqrt((z0-z)**2+(x0-x)**2)
        if(z.gt.z0)d=-d
        db=sqrt((zmta0-z0)**2+(xmta0-x0)**2)
        write(1,100)z,x,(alp-fi)/pi*180.,bahn+db,d
      end do
      stop
100   format(5f12.5)
      end



