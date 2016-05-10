c
c      s u _ f i n d _ f o c u s . f
c
c      Subroutine bestimmt den Schnittpunkt aus zwei Bahnen, die mit
c      Winkelabweichungen bestimmt wurden. Der Schnittpunkt entspricht
c      dann dem Ort der Punkt-Punkt Abbildung. Die Bahnen stehen im
c      Feld bahn(2,2,300), dabei steht der 1. Index fuer die beiden Winkel,
c      der 2. fuer die Ortskoordinaten und der 3. ist der Laufindex.
c
c      Die Routine beginnt mit den jeweils groessten x-Werten und sucht
c      dann das Punktepaar mit dem kleinsten Abstand.
c      Dann nimmt es die davor liegenden Punkte dazu und bestimmt den
c      gesuchten Schnittpunkt. Dieser wird als Ergebis zurueckgegeben.
c
      subroutine su_find_focus(bahn,x,z)
c
      implicit none
      real bahn(2,2,300),x,z,x1,z1,x2,z2,a(2,3),b(2,3),
     &  za(2),d,dist,ax(3),az(3)
      integer i1,i2,la(2),ls(2),l,id,ll(2)
c
      za(1)=0.
      za(2)=0.
      do id=1,2
        do l=1,300
          if(bahn(id,2,l).gt.za(id))then
            za(id)=bahn(id,2,l)
            la(id)=l
          end if
          if(bahn(id,2,l).eq.0.)then
            ls(id)=l-1
            goto 1
          end if
        continue
      end do
 1    continue
      end do
      print*,'la(1),ls(1),la(2),ls(2)',la(1),ls(1),la(2),ls(2)
      dist=1.E30
      do i1=la(1),ls(1),1
        x1=bahn(1,1,i1)
        z1=bahn(1,2,i1)
        do i2=la(2),ls(2),1
          x2=bahn(2,1,i2)
          z2=bahn(2,2,i2)
          d=(x1-x2)**2+(z1-z2)**2
          if(d.lt.dist)then
            dist=d
            ll(1)=i1
            ll(2)=i2
          end if
        end do
      end do
      do l=1,3
        do id=1,2
          a(id,l)=(bahn(id,2,ll(id))-bahn(id,2,ll(id)-l))/
     &      (bahn(id,1,ll(id))-bahn(id,1,ll(id)-l))
          b(id,l)=bahn(id,2,ll(id))-a(id,l)*bahn(id,1,ll(id))
        end do
        ax(l)=(b(2,l)-b(1,l))/(a(1,l)-a(2,l))
        az(l)=a(1,l)*ax(l)+b(1,l)
      end do
      x=(ax(1)+ax(2)+ax(3))/3.
      z=(az(1)+az(2)+az(3))/3.
      print*,'ax',ax
      print*,'az',az
      return
      end
c      
