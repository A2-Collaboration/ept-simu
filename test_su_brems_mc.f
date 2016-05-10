c
c      t e s t _ s u _ b r e m s _ m c . f
c
c      Programm tetstet die routine su_brems_mc.f
c
c      Laden als:
c g77 -o test_su_brems_mc.exe test_su_brems_mc.f su_brems_mc.f brems.f
c      
      implicit none
      real e0,k,me,pi,tb,xyr(2),txybe(2)
      integer i,j
c
      pi=3.1415926
      me=.511006
      print*,'test_su_brems_mc.f'
c
      call su_brems_mc(0.,0.,tb,xyr,txybe)
c
      print*,'Energie/MeV des einl. Elektrons '
      print*,'und des B-Photons [*,*]:'
      read*,e0,k
      print*,'i,tb,xyr(1),xyr(2),txybe(1),txybe(2)'
      do i=1,100
        call su_brems_mc(e0,k,tb,xyr,txybe)
        print 100,i,tb,xyr,txybe
        write(10,100)i,tb,xyr,txybe
 100    format(i5,1p5e14.4)
      end do
      stop
      end
      
