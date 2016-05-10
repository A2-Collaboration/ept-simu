c
c      t e s t _ s u _ t r a c e _ t r a c k . f
c
c      Test-Programm  fuer Routine su_trace_track.f
c
c      Laden als:
c      g77 -fbounds-check -o test_su_trace_track.exe
c        test_su_trace_track.f su_trace_track.f su_such_feld.f
c      
      implicit none
      real xl,x0,z0,a0,x,z,a,e,track(2,300),pi
      integer it
c
      pi=3.1415926
      print*,'test Routine su_trace_track.f'
      xl=0.
      it=1
      call su_trace_track(xl,it,x0,z0,a0,e,x,z,a,track)
      if(it.eq.0)stop
1     print*,'Eingabe Schrittweite l, Iterationen it und e [*,*,*]:'
      read*,xl,it,e
      print*,'Eingabe x0,z0/Index,a0/Grad [*,*,*]:'
      read*,x0,z0,a0
      a0=a0/180.*pi
      call su_trace_track(xl,it,x0,z0,a0,e,x,z,a,track)
      if(x.eq.0)stop
      print*,'x,z,a  ',x,z,a
      do i=1,300
        if(track(1,i)+track(2,i).eq.0.)stop
        if(track(1,i).eq.0)then
          write(10,*)'     *   ',track(2,i)
        else
          if(track(2,i).eq.0.)then
            write(10,*) track(1,i),'   *  '
          else
            write(10,*)track(1,i),track(2,i)
          endif
        endif
      end do  
      goto 1
      end
c
