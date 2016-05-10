c      b r e m s . f
c
c      Routinen zur Berechnung von Brems-Querschnitten
c      Formeln nach   Koch and Motz, Rev.Mod.Phys., 31(1959)920
c      tragen die dortige Formelbezeichnung
c      hinter den Buchstaben  'br'. z.B. Formel 3bs(e) wird in der
c      routine  'br3bse'  berechnet. die Wirkungsquerschnitte sind
c      jeweils in barn angegeben, die Eingangsenergien in mec**2,
c      Winkel in rad.
c
c      Bezeichnungen
c      e0 Energie des einlaufenden Elektrons
c      xk die des auslaufenden Photons
c      und  e  des auslaufenden  Elektrons
c       e0=xk+e
c      z Ordnungszahl
c      t Emissionswinkel des Photons und
c      b der des Elektrons
c      i Index  =0 ds/dOmrga  =1 ds/dt(oder db fuer das Elektron)
c
c
c	br2bn
c
c	ds/(dk*dt) noscreen,Born
c	entspricht brglhu fuer den unpolarisierten Teil
c	(as=99999.)
c
	subroutine br2bn(as,az,at,axk,ae0,i)
        implicit none
        real as,az,at,axk,ae0
        real*8 z,t,xk,e0,st,e,e02,e2,xk2,p02,p2,p0,p,ee0,
     &    pp0,ct,d0,d02,d04,eps,q2,q,xl,epsq,vf,s,s1
        integer i
c
        
	z=az
	t=at
	xk=axk
	e0=ae0
	st=sin(t)
	e=e0-xk
	e02=e0**2
	e2=e**2
	xk2=xk**2
	p02=e02-1.
	p2=e2-1.
	p0=sqrt(p02)
	p=sqrt(p2)
	ee0=e*e0
	pp0=p*p0
	ct=cos(t)
	d0=e0-p0*ct
	d02=d0**2
	d04=d02**2
	eps=log((e+p)/(e-p))
	q2=p02+xk2-2.*p0*xk*ct
	q=sqrt(q2)
	xl=log((ee0-1.+pp0)/(ee0-1.-pp0))
	epsq=log((q+p)/(q-p))
	vf=2.3058e-5*z**2*p/p0/xk
	if(i.eq.1)vf=vf*6.2832*st
	s=8.*st**2*(2.*e02+1.)/(p02*d04)-2.*(5.*e02+2.*ee0+3.)/
     *   (p02*d02)-2.*(p02-xk2)/(q2*d02)+4.*e/(p02*d0)
	s1=4.*e0*st**2*(3.*xk-p02*e)/(p02*d04)+
     *   (4.*e02*(e02+e2)+2.-2.*(7.*e02-3.*ee0+e2))/(p02*d02)+
     *   2.*xk*(e02+ee0-1.)/(p02*d0)
	s=s+xl/pp0*s1-4.*eps/(p*d0)+epsq/(p*q)*(4./d02-6.*xk/d0-
     *     2.*xk*(p02-xk2)/(q2*d0))
	as=s*vf
	return
	end
c
c
c      br2bna
c
c      ds/(dk*dt)  kleine Winkel,Born,noscreen,extr.rel.
c
       subroutine br2bna(s,z,t,xk,e0,i)
c
       implicit none
       real s,z,t,xk,e0,e,y,y2,y1,vf
       integer i
c       
       e=e0-xk
       y=e0*t
       y2=y*y
       y1=(1.+y2)**2
       vf=3.6893E-4*z**2*e/(xk*y1)
       if(i.ne.0)vf=vf*t*6.2832
       s=vf*(16.*y2*e0/y1-(e0+e)**2/e+
     *  ((e*e+e0*e0)/e-4.*y2*e0/y1)*2.*log(2.*e*e0/xk))
       return
       end
c
c
c       br2bnb
c
c       ds/(dk*dt)
c       differential in photon energy and angle
c       Born,no screen,large angles !!! extr.rel.
c
        subroutine br2bnb(s,z,t,xk,e0,i)
c
        implicit none
        real s,z,t,xk,e0,st,ct,ct1,e,e02,e2,e0e,q2,q,
     &    vf,x
        integer i
c        
        st=sin(t)
        ct=cos(t)
        ct1=1.-ct
        e=e0-xk
        e02=e0**2
        e2=e**2
        e0e=e0*e
        q2=e2+2.*xk*e0*ct1
        q=sqrt(q2)
        vf=z**2*4.612e-5*e/(e02*e0*xk*ct1**2)
        if(i.ne.0)vf=vf*6.2833*st
        x=(e02+e2)/e0e*st**2*log(2.*e0e/xk)-(5.*e0+2.*e)/e0+
     *    2.*(e2-2.*e02*log(xk/e0))/e0e*ct1-e*(q2+e*(xk-e0*ct))
     *    /(e0*q2*ct1)-e0*xk*ct1/(e*q2*q)*(3.*q2+e*(e0+xk))*
     *    log((q+e)/(q-e))
       s=vf*x
       return
       end
c
c
c      br4bs
c      Firad, Born,compl.screen.,extr.rel.
c
       subroutine br4bs(s,z)
c
       implicit none
       real s,z
c       
       s=2.318E-3*z**2*(log(183/(z**.333333))+1./18.)
       return
       end
c
c
c      br4bn
c      step wird fuer eine interne Integration gebraucht
c      step=0 setzt step=100
c      Firad, Born,noscreen
c
       subroutine br4bn(s,z,e0,step)
c
       implicit none
       real s,z,e0,step,p0,ep,sep,epl,x,dx,ss,y
c       
       if(step.eq.0.)step=100.
       p0=sqrt(e0**2-1.)
       ep=e0*p0
       sep=e0+p0
       epl=log(sep)
       x=2.*p0*sep
       dx=x/step
       ss=0.
       do 10 y=dx*.5,x,dx
 10    ss=ss+log(y+1.)/y
       ss=ss*dx
       s=5.7952e-4*z**2*((12.*e0**2+4.)/(3.*ep)*epl-
     *   (8.*e0+6.*p0)/(3.*ep*p0)*epl**2-
     *   1.33333+2./ep*ss)
       return
       end
c
c
c      br4bna
c      Firad, Born,noscreen,non rel.
c
       subroutine br4bna(s,z)
c
       implicit none
       real s,z
c       
       s=3.091E-3*z**2
       return
       end
c
c
c      br4bnb
c      Firad, Born,noscreen, extr.rel.
c
       subroutine br4bnb(s,z,e0)
c
       implicit none
       real s,z,e0
c       
       s=2.318E-3*z**2*(log(2.*e0)-.3333333)
       return
       end
c
c	br3bn
c
c	ds/dk Born, no screening
c
	subroutine br3bn(ys,yz,yxk,ye0)
c
        implicit none
        real ys,yz,yxk,ye0
        real*8 z,xk,e0,pi,e,e2,xk2,e02,p02,p2,p0,p,vf,e0e,
     &    p0p,xl,eps0,eps,s1,s2,s
c        
	z=yz
	xk=yxk
	e0=ye0
	pi=3.14159
	e=e0-xk
	e2=e*e
	xk2=xk*xk
	e02=e0*e0
	p02=e02-1.
	p2=e2-1.
	p0=sqrt(p02)
	p=sqrt(p2)
	vf=.9977/4./pi/137*z**2/xk*p/p0
	e0e=e0*e
	p0p=p0*p
	xl=2*log((e0e+p0p-1.)/xk)
	eps0=log((e0+p0)/(e0-p0))
 	eps=log((e+p)/(e-p))
	s1=8.*e0e/(3.*p0p)+xk2*(e0e**2+p0p**2)/p0p**3+.5*xk/p0p*
     *    (eps0*(e0e+p02)/p0**3-eps*(e0e+p2)/p**3+2.*xk*e0e/p0p**2)
	s2=4./3.-2.*e0e*(p02+p2)/p0p**2+eps0*e/p0**3+eps*e0/p**3-
     *    eps0*eps/p0p
	s=vf*(s2+xl*s1)
	ys=s
	return
	end
c
c      br3bsa
c      ds/dk  Born,approx.screen,extr.rel.
c
       subroutine br3bsa(s,z,xk,e0)
c
       implicit none
       real s,z,xk,e0,ee0
c       
       ee0=(e0-xk)/e0
       s=2.318e-3*z**2/xk
       s=((1.+ee0**2-.666667*ee0)*log(183.*z**(-.33333))+ee0/9.)*s
       return
       end
c
c      br3bse
c      ds/dk  Born,approx.screen,extr.rel
c
       subroutine br3bse(s,z,xk,e0)
c
       implicit none
       real s,z,xk,e0,e,h1,h2,b,b2,atb,h,ee0,vf
c       
       e=e0-xk
       h1=xk/(2.*e*e0)
       h2=z**.333333/111.
       b=h2/h1
       b2=b*b
       atb=atan(b)/b
       h=log(1./(h1*h1+h2*h2))
       ee0=e/e0
       vf=2.*z**2/(xk*1725.48)
       s=(1.+ee0**2.-.66667*ee0)*(h+1.-2.*atb)+
     * (2./b2*log(1.+b2)+1.33333/b2*(2.-b2)*atb-2.66667/b2+.2222)*ee0
       s=s*vf
       return
       end
c
c
c      br2bs
c      ds/(dk*domega) Born,approx.screen,small angles,extr.rel.
c
       subroutine br2bs(as,az,ate,axk,ae0,i)
c
       implicit none
       real as,az,ate,axk,ae0
       real*8 z,te,xk,e0,e,y,e02,y2,y1,x,z1,s
       integer i
c       
	z=az
	te=ate
	xk=axk
	e0=ae0
       e=e0-xk
       y=e0*te
       e02=e0*e0
       y2=y*y
       y1=(1.+y2)**2
       x=1./((xk/(2.*e0*e))**2+(z**.33333/111.)**2/y1)
       x=log(x)
       z1=z
       if(z.eq.0.)z1=1.
       s=3.6893e-4*e02*z1**2/(xk*y1)*
     *   (16.*y2*e/(e0*y1)-(e0+e)**2/e02+
     *   ((e02+e*e)/e02-4.*y2*e/(e0*y1))*x)
       if(i.ne.0)s=s*te*6.2833
	as=s
       return
       end
c
c      br2bsd
c      ds/(dk*domega) Born,approx.screen,small angles,extr.rel.
c      Version voll doppelt genau
c       
       subroutine br2bsd(as,az,ate,axk,ae0,i)
c
       implicit none
       real*8 as,az,ate,axk,ae0
       real*8 z,te,xk,e0,e,y,e02,y2,y1,x,z1,s
       integer i
c
       z=az
       te=ate
       xk=axk
       e0=ae0
       e=e0-xk
       y=e0*te
       e02=e0*e0
       y2=y*y
       y1=(1.+y2)**2
       x=1./((xk/(2.*e0*e))**2+(z**.33333/111.)**2/y1)
       x=log(x)
       z1=z
       if(z.eq.0.)z1=1.
       s=3.6893e-4*e02*z1**2/(xk*y1)*
     &      (16.*y2*e/(e0*y1)-(e0+e)**2/e02+
     &      ((e02+e*e)/e02-4.*y2*e/(e0*y1))*x)
       if(i.ne.0)s=s*te*6.2833
       as=s
       return
       end
c      
c
c      br2cn
c      ds/(dk*dt) extr.rel.,small angles
c
       subroutine br2cn(s,z,t,xk,e0)
c
       implicit none
       real s,z,t,xk,e0,e,e02,p0,u,u2,zeta,zet2,del,
     &   gam
c       
       e=e0-xk
       e02=e0**2
       p0=sqrt(e02-1.)
       u=p0*t
       u2=u*u
       zeta=1./(1.+u2)
       zet2=zeta**2
       del=xk/(2.*e0*e)
       gam=log(1./del)-2.-1.2021*(z/137)**2
       s=2.3234e-3*z**2/xk*p0*u*zet2/e02*
     *   (-2.*e0*e*(1.+4.*u2*zet2*gam)+(e02+e**2)*(3.+2.*gam))
       return
      end
      
c
c      z i r p
c
c      function  zirp(e0,xk,z) bestimmt den mittleren Zirkular-
c      polarisationsgrad von Bremsquanten nach Integration
c      ueber Elektronen- und Photonen-Richtungen
c      bei 100% Polarisation der einl. Elektronen
c      nach Fronsdal,Ueberall Phys.Rev.111(1958)580
c
c      e0 - Energie der einlaufenden Elektronen in mec**2
c      xk - Energie der Photonen in mec**2
c      z - Ordnugszahl Bremstarget
c
       function zirp(e0,xk,z)
c
       implicit none
       real zirp,e0,xk,z,x,be,q,q2,h,h1,a,b,h2
c       
       x=xk/e0
       be=111./z**.333333
       q=.5*be/e0*x/(1.-x)
       q2=q**2
       h=log(be**2/(1.+q2))
       h1=2*q*atan(1./q)
       a=1.+h-h1
       b=.666667+h+4.*q2-3.*q2*log(1.+1./q2)-2.*q2*h1
       h2=.666667*(1.-x)*b
       zirp=((2.-x)*a-h2)*x/((1.+(1.-x)**2)*a-h2)
       return
       end
c
c
c      z i r p t
c
c     Function, die die Zirkularpolarisation von Bremsquanten
c     erzeugt von 100% longitudinal polarisierten Elektronen
c     bestimmt.
c     e0 - Energie der einl. Elektronen in mec**2
c     xk - Energie der auslaufenden Photonen in mec**2
c     z - Ordnungszahl des Bremstargets
c     t - Emissionswinkel der Photonen in rad
c
       function zirpt(e0,xk,z,t)
c
       implicit none
       real zirpt,e0,xk,z,t,x,u,u2,xa,h,h1
c       
       x=xk/e0
       u=e0*t
       u2=(1.+u**2)**2
       xa=(z**(1./3.)/111.)**2+(x*.5/(e0*(1.-x)))**2*u2
       h=log(u2/xa)-1.
       h1=2.*(1.-x)*(2.*u**2/u2*(h-3.)+1.)
       zirpt=((2.-x)*h-h1)*x/((1.+(1.-x)**2)*h-h1)
       return
       end
c
c      z i r p i t
c
c     Function, die die Zirkularpolarisation von Bremsquanten
c     erzeugt von 100% longitudinal polarisierten Elektronen
c     bestimmt.
c     Dabei wird bis zu dem angegebenen Winkel integriert !!!! (it)
c     e0 - Energie der einl. Elektronen in mec**2
c     xk - Energie der auslaufenden Photonen in mec**2
c     z - Ordnungszahl des Bremstargets
c     t1,t2 - Emissionswinkel der Photonen in rad,
c     zwischen denen integriert wird.	
c
       function zirpit(e0,xk,z,t1,t2)
c
       implicit none
       real e0,xk,z,t1,t2,ssbh,ssbl,x,t,u,u2,xa,h,h1,
     &   sbh,sbl,zirpit
c       
	ssbh=0.
	ssbl=0.
       x=xk/e0
	do 10 t=t1,t2,(t2-t1)/500.
       u=e0*t
       u2=(1.+u**2)**2
       xa=(z**(1./3.)/111.)**2+(x*.5/(e0*(1.-x)))**2*u2
       h=log(u2/xa)-1.
       h1=2.*(1.-x)*(2.*u**2/u2*(h-3.)+1.)
	sbh=((1.+(1.-x)**2)*h-h1)/x*sin(t)/u2
	sbl=((2.-x)*h-h1)*sin(t)/u2
	ssbh=ssbh+sbh
10	ssbl=ssbl+sbl
       zirpit=ssbl/ssbh
       return
       end
c
c      briksc
c
c      Screening-Korrektur zum Integral ueber alle Photonen-(k)-
c      Richtungen. Die Korrektur ist additiv zum ungescreenten
c      Integral, das in der Function  brik  steht.
c      Zu Grunde liegt der Atomformfaktor nach Thomas-Fermi.
c      Integral geloest von Maximon.
c      Die Korrektur ist somit differentiell in  Elektronenenergie
c      und Elektronen(raum)winkel
c      dscreen/(de*dt(o)).
c      Laeuft nur fuer extreme Energien und Winkel
c
       function briksc(e0,xk,z,b,i)
c
       implicit none
       real a(3),bb(3),be(3),vf,e,ee0,d,qm,zzqm,su,s,
     &   dqm,briksc,z,b,e0,xk
       integer i,j
       data a/.1,.55,.35/,be/6.,1.2,.3/
c
       vf=9.223e-5*z**2
       if(i.ne.0)vf=vf*6.2832*sin(b)
       e=e0-xk
       ee0=e*e0
       d=xk*.5/ee0
       qm=d+(sin(.5*b))**2/d
       zzqm=(z**.33333/(121.*qm))**2
       do 10 i=1,3
 10    bb(i)=zzqm*be(i)**2
       su=0.
       do 11 i=1,3
       s=0.
       do 12 j=1,3
       if(i.eq.j)goto 12
       s=s+a(i)*a(j)*((1.+bb(j))/(bb(i)-bb(j))*log(1.+bb(j))+.5)
 12    continue
 11    su=su+2.*s-a(i)**2*log(1.+bb(i))
       dqm=d/qm
       briksc=vf/(e0**2*xk)*(e0**2+e**2-4.*ee0*dqm*(1.-dqm))*su
       return
       end
c
c      brglhu
c
c       Bremsquerschnitt integriert ueber alle Richtungen
c       p (auslaufendes Elektron)
c       nach  R.L.Gluckstern und M.H.Hull, Phys.Rev. 90(1953)1030
c       fuer as=99999. wird nur au bestimmt, sonst nur as,ap und
c	au=as+ap
c
c       Parameter:
c       Eingang - az =^= z
c                 ate=^=te, ak=^=xk, ae0=^=e0
c       Ausgang - au gesamter Querschnitt (au=as+ap)
c                 as,ap senkrecht und paralell pol. W.Q.
c       **** fuer as=99999. wird nur au bestimmt ****
c
       subroutine brglhu(au,as,ap,az,ate,ak,ae0,i)
c
       implicit none
       real au,as,ap,az,ate,ak,ae0
       real*8 e0,xk,teta0,z,e,p0,p,st0,ct0,del0,t,r,eps,
     &   epst,st02,del02,del04,e02,e2,xk2,p02,t2,vf,x,x1,x2
       integer i
c       
       e0=ae0
       xk=ak
       teta0=ate
       z=az
        e=e0-xk
        p0=sqrt(e0*e0-1.)
        p= sqrt(e *e -1.)
        st0=sin(teta0)
        ct0=cos(teta0)
        del0=e0-p0*ct0
       t=sqrt(p0*p0+xk*xk-2.*xk*p0*ct0)
        r=log((e*e0-1.+p*p0)/(e*e0-1.-p*p0))
        eps=log((e+p)/(e-p))
        epst=log((t+p)/(t-p))
       st02=st0*st0
       del02=del0*del0
       del04=del02*del02
       e02=e0*e0
       e2=e*e
       xk2=xk*xk
       p02=p0*p0
       t2=t*t
       vf=z**2*2.306e-5*p/(p0*xk)
       if(i.ne.0)vf=vf*st0*6.283
	if(as.ne.99999.)goto 20
        au=vf*(8.*st02*(2.*e02+1.)/(p02*del04)
     *   -2.*(5.*e02+2.*e0*e+3.)/(p02*del02)
     *   -2.*(p02-xk2)/(t2*del02)+4.*e/(p02*del0)
     *   +(r/(p*p0))*(4.*e0*st02*(3.*xk-p02*e)/(p02*del04)
     *   +(4.*e02*(e02+e2)-2.*(7.*e02-3.*e*e0+e2)+2.)/(p02*del02)
     *   +2.*xk*(e02+e*e0-1.)/(p02*del0))
     *   +(epst/(p*t))*(4./del02-6.*xk/del0-2.*xk*(p02-xk2)/(t2*del0))
     *   -4.*eps/(p*del0))
       if(as.eq.99999.)return
20       x=8.*st02*(2.*e02+1.)/(p02*del04)
     *  -(5.*e02+2.*e*e0+5.)/(p02*del02)
     *  -(p02-xk2)/(t2*del02)+2.*(e+e0)/(p02*del0)
       x1=r/(p*p0)*(4.*e0*st02*(3.*xk-p02*e)/(p02*del04)+(2.*e02*
     *  (e02+e2)-(9.*e02-4.*e0*e+e2)+2.)/(p02*del02)
     *  +xk*(e02+e*e0)/(p02*del0))
       x2=epst/(p*t)*(4./del02-
     *  7.*xk/del0-xk*(p02-xk2)/(t2*del0)-4.)-4.*eps/(p*del0)
     *  +1./(p02*st02)*((2.*r/(p*p0))*(2.*e02-e0*e-1-xk/del0)
     *  -4.*epst*(del0-e)**2/(p*t)-2.*eps*(del0-e)/p)
       ap=vf*(x+x1+x2)
       as=vf*(-1.*(5.*e02+2.*e*e0+1.)/(p02*del02)
     *  -(p02-xk2)/(t2*del02)-2.*xk/(p02*del0)+r/(p*p0)*
     *  ((2.*e02*(e02+e2)-(5.*e02-2.*e*e0+e2))/(p02*del02)
     *  +xk*(e02+e*e0-2.)/(p02*del0))+epst/(p*t)*
     *  (xk/del0-xk*(p02-xk2)/(t2*del0)+4.)-1./(p02*st02)*
     *  (2.*r/(p*p0)*(2.*e02-e*e0-1.-xk/del0)-
     *  4.*epst*(del0-e)**2/(p*t)-2.*eps*(del0-e)/p))
	au=as+ap
       return
       end
c
c      brik
c
c      Bremsquerschnitt integriert ueber Richtungen und
c      Polarisationen des auslaufenden Photons
c      nach Maximon und Isabelle  Phys.Rev 133B(1964)1344
c      die Parameter sind : a0=^=e0, ak=^=xk, az=^=z,
c                           ab=^=b(elektronenwinkel)
c                           io=0 ds/do, io=/=0 ds/db
c
       function brik(a0,ak,az,ab,io)
c
       implicit none
       real a0,ak,az,ab,brik
       real*8 e0,xk,z,teta,e,p0,p,st,ct,st2,pp0,ee0,
     &   xl,xl2,x2l,e0e,e0p,p0p,p0e,xk2,x,y,c,d,f
       integer io
c
       e0=a0
       xk=ak
       z=az
       teta=ab
       e=e0-xk
       p0=sqrt(e0**2-1.)
       p =sqrt(e **2-1.)
       st=sin(teta)
       ct=cos(teta)
       st2=st**2
       pp0=p*p0
       ee0=e*e0
       xl=ee0-pp0*ct-1.
       xl2=xl**2
       x2l=sqrt(xl*(xl+2.))
       e0e=e0+e
       p0p=p0**2+p**2
       e0p=2.*e0*p**2
       p0e=2.*e*p0**2
       xk2=xk**2
       x=(xk2*xl+(2.*ee0-xl)*(xl+1.))/(xl2*x2l)*log(xl+1.+x2l)*2.
       y=xk/(xl2*p**3)*(-2.*xk*e*p**2+3.*e*e0e-p0p+(2.*ee0-xl)*
     *  (ee0+p**2)+p0e*st2*(e0p-3.*xk*e**2)/xl2)*log(e+p)
       c=-xk/(xl2*p0**3)*(2.*xk*e0*p0**2+3.*e0*e0e-p0p+(2.*ee0-xl)*
     *  (ee0+p0**2)+e0p*st2*(p0e+3.*xk*e0**2)/xl2)*log(e0+p0)
       d=2.*xk2*st2/(xl2*pp0)**2*(2.*pp0**2*(e0**2+e**2-ee0)
     *  +3.*xk2*e0e**2)
       f=xk2/(xl2*pp0**2)*(p0p*(ee0+1.)-(2.*ee0-xl)
     *  *(p0p+ee0+1.)-3.*e0e**2)-2./xl2*(2.*ee0-xl)
       brik=(x+y+c+d+f)*p/(p0*xk)*9.2234E-5*z**2
       if(io.ne.0)brik=brik*st*6.2832
       return
       end
c
c
c
c       subroutine moliere(c,z,q2)
c
c       Bestimmt den Screening-Korrektur-Faktor (1-f(q))/q^2, der
c       quadriert an Stelle der generellen 1/q^4 Abhaengikeit
c       zu setzten ist. Nach Moliere.
c       c ist der KorrekturFaktor (1-f(q))/q^2
c       z die Ladung des Kerns
c       q2 das Quadrat des Impulsuebertrags.
c
       subroutine moliere(c,z,q2)
       implicit none
       real q2,c,z,bet(3),beta(3),alpha(3)
       integer i
c
       data alpha/.1,.55,.35/,bet/6.,1.2,.3/
c
       c=0.
       do i=1,3
         beta(i)=z**(1./3.)/121.*bet(i)
         c=c+alpha(i)/(q2+beta(i)**2)
       end do
       return
      end
c
c
      
