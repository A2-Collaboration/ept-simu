c      sub_g_k.f
c
c      Subroutinen zu Geraden und Kreisen und Punkten
c      Geraden sind gegeben durch  y=ai*x+bi
c      Kreise durch (x-xi)^2 + (y-yi)^2 = ri^2
c      Schnittpunkte:
c      werden ausgegeben (x,y), ebenso, dort, wo Kreise
c      im Spiel sind die entsprechenden Winkel der Radien
c      zum Schnittpunkt in rad (ti)
c      Die Subroutinen heissen dann sinngemaess:
c      cut_g_g, cut_g_k, cut_k_k
c      sie enthalten in ihren Prameterlisten gegebenenfalls eine
c      Zahl vz, die das Vorzeichen einen Wurzel angibt.
c      Man muss die vz bestimmen um sinnvolle L"osungen zu erhalten.
c      Falls es keine Loesung gibt, wird vz=0 gesetzt
c
c      Gerade durch zwei Punkte, Kreis durch drei Punkte
c      die Subroutinen heissen dann sinngemaess:
c      gerade, kreis und haben die Parameterlisten
c      (x1,y1,x1,y2,vz,a,b) und (x1,y1,x2,y2,x3,y3,vz,xk,yk,r)
c      wobei vz die obigen Bedeutungen annehmen kann
c      
c      Rotation eines Koordinatensystems
c
c      Folgende Routinen sind im Paket enthalten      
c      subroutine cut_g_g(a1,b1,a2,b2,vz,x,y)
c      subroutine cut_g_k(a,b,x1,y1,r,vz,x,y,t)
c      subroutine cut_k_k(x1,y1,r1,x2,y2,r2,vz1,vz2,x,y,t1,t2)
c      subroutine cut_k_k_trig(x1,y1,r1,x2,y2,r2,vz1,x,y)
c      subroutine gerade(x1,y1,x2,y2,vz,a,b)
c      subroutine kreis(x1,y1,x2,y2,x3,y3,vz,xk,yk,r)
c      subroutine rot(x0,y0,x,y,alpha)
c
c
c
      subroutine cut_g_g(a1,b1,a2,b2,vz,x,y)
c
c      Schnitt Gerade-Gerade
c
      if(a1-a2.eq.0.)then
        vz=0.
        return
      endif
      x=(b2-b1)/(a1-a2)
      y=a1*x+b1
      return
      end
c      
c      
       subroutine cut_g_k(a,b,x1,y1,r,vz,x,y,t)
c
c      Schnitt Gerade-Kreis
c       
       e=(a*(b-y1)-x1)/(a**2+1.)
       f=(r**2-(b-y1)**2-x1**2)/(a**2+1.)
       if(f+e**2.lt.0.)then
         vz=0.
         return
       else
         x=-e+vz*sqrt(f+e**2)
         y=a*x+b
       endif
       t=atan2((y-y1),(x-x1))
       return
       end
c     
       subroutine cut_k_k(x1,y1,r1,x2,y2,r2,vz1,vz2,x,y,t1,t2)
c
c      Schnitt Kreis-Kreis
c       
       c=((r1**2-r2**2)-(y2-y1)**2+(x2**2-x1**2))/(2.*(x1-x2))
       d=(y2-y1)/(x1-x2)
       e=(c-x2*d**2)/(1.+d**2)
       f=(r2**2*d**2-c**2-x2**2*d**2)/(1.+d**2)
       if(f+e**2.lt.0.)then
         vz1=0.
         return
       endif
       x=-e+vz1*sqrt(f+e**2)
       if(r2**2-(x-x2)**2.lt.0.)then
         vz2=0.
         return
       endif
       y=y2+vz2*sqrt(r2**2-(x-x2)**2)
       t1=atan((y-y1)/(x-x1))
       t2=atan((y-y2)/(x-x2))
       return
      end
c      
c     Schnitt Kreis-Kreis mit Trigonometrie
c
      subroutine cut_k_k_trig(x1,y1,r1,x2,y2,r2,vz1,x,y)
c
c      trigonometrischer Ansatz zum Schnittpunkt zweier Kreise
c      Rechnung am Dreieck mit Bestimmung der Basislinie aus
c      den xi,yi und der darauf stehenden Hoehe, dann Winkel usw.
c
      implicit none
      real x1,x2,y1,y2,r1,r2,x,y,vz1,vz,d,d2,bet,alph,gam,
     &  x0,y0,h
c      
      d=sqrt((x1-x2)**2+(y1-y2)**2)
      vz=1.
c      if(x2.lt.x1.or.y2.lt.y1)vz=-1
      d=vz*d
      bet=atan2(y2-y1,x2-x1)
      d2=(d**2-r1**2+r2**2)/(2.*d)
      h=vz1*sqrt(r2**2-d2**2)
      alph=atan2(h,d2)
      gam=alph-bet
      x0=r2*cos(gam)
      y0=r2*sin(gam)
c      print*,'d,bet,d2,h,alph,gam,x0,y0'
c      print*,d,bet,d2,h,alph,gam,x0,y0
      x=x2-x0
      y=y2+y0
        return
        end
c        
      subroutine gerade(x1,y1,x2,y2,vz,a,b)
c
c      Zweipunkt-Geradengleichung
c
      if(x2-x1.eq.0.)then
        vz=0
        return
      endif
      a=(y2-y1)/(x2-x1)
      b=y1-a*x1
      return
      end
c
      subroutine kreis(x1,y1,x2,y2,x3,y3,vz,xk,yk,r)
c
c     Kreis durch drei Punkte
c      
      if(y2-y1.eq.0..or.y3-y1.eq.0.)then
        vz=0.
        return
      endif
      a=(x1-x2)/(y2-y1)-(x1-x3)/(y3-y1)
      b=.5*((y3**2+x3**2-y1**2-x1**2)/(y3-y1)+
     &      (y2**2+x2**2-y1**2-x1**2)/(y2-y1))
      xk=b/a
      yk=xk*(x1-x2)/(y2-y1)+.5*(y2**2+x2**2-y1**2-x1**2)/(y2-y1)
      r=sqrt(xk**2+yk**2)
      return
      end
c
c
      subroutine rot(x0,y0,x,y,alpha)
c
c     Rotation des  !Bezugssystems! im Uhrzeigersinn um alpha
c     im Bogenmass
c
      x=x0*cos(alpha)-y0*sin(alpha)
      y=x0*sin(alpha)+y0*cos(alpha)
      return
      end
c      
