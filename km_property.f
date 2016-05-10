c      k m _ p r o p e r t y . f
c
c      Programm gibt Schaetzwerte fuer den Ablenkwinkel
c      im KM, wenn man hierfuer den Dipole DS4 verwendet,
c      der eine Polschuhlaenge von 40 cm hat.
c      Das Ergebnis steht in km_property.dat
c
      implicit none
      real r,alpha,d,t,l
c
      l=40.
      open(1,file='km_property.dat')
      do e=1000,1600,200
        do t=1,1.5,.05
          r=e/t/3.
          alpha=asin(l/r)
          write(1,*)t,alpha/3.1415926*180.,e
        end do
        write(1,*)'    *     *     *'
      end do
      stop
      end
