c      
c     ********************** Routinen *****************************
c
c
c      subroutine s u _ c u t _ b i l d l i n . f
c
c      bestimmt Schnittpunkt zwischen Elektronenbahn aus
c      track(2,300) (erstellt wie in raytrace_feld.f)
c      und einer vorgegebenen Bildlinie (gegeben durch die Parameter
c      ab, bb).
c      Das Feld track(2,300) wird im Hauptprogramm erstellt und hier
c      wieder verwendet. Dann wird auf der Elektronenbahn rueckwaerts
c      nach dem Schnittpunkt gesucht, der dann als x, z
c      als Ergebnis ausgegeben wird.
c      Zur Nutzung in einem folgenden Histogramm wird noch der Abstand
c      auf der Bildlinie von Index_x = 0 in Index-Einheiten bestimmt.
c
      subroutine su_cut_bildlin(ab,bb,x,z,d,track)
c
      implicit none
      real ab,bb,track(2,300),x,z,a,b,d
      integer i,it
c
      do i=1,300
        if(track(1,i)+track(2,i).eq.0.)goto 1
      end do
 1    it=i-1
      do i=it,2,-1
        a=(track(2,i)-track(2,i-1))/(track(1,i)-track(1,i-1))
        b=track(2,i)-a*track(1,i)
        x=(bb-b)/(a-ab)
        z=ab*x+bb
        if(z.le.track(2,i)) then
          d=sqrt(x**2+(z-bb)**2)
          close(2)
          return
        end if
        if(z.le.track(2,i-1).and.z.ge.track(2,i)) then
         d=sqrt(x**2+(z-bb)**2)
          close(2)
          return
        end if
      end do
      print*,'aus su_cut_bildlin.  keine Loesung gefunden'
      x=0.
      z=0.
      return
      end
c      
