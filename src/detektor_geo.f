c
c      d e t e k t o r _ g e o . f
c
c      Programm erstellt eine Datenliste zum Zeichnen der
c      Detektoren in einem GLE-Plot. Basis sind die Dimensionen
c      gegeben in Zeichnung endpoint_d_20.dwg
c
c      Laden als:
c      g77 -fbounds-check -o detektor_geo.exe detektor_geo.f
c
c      Das Ergebnis (Plot-File) steht in detektor_geo.geo
c
      implicit none
      real xz(5,2),x0,z0,step,xx
      integer i,j,k
c
c      Grunddaten:
      x0=79.47
      z0=2.81      ! absolute Index-Koordinaten Mittelpunkt Z.# 47
      xz(1,1)=-.63
      xz(1,2)=.14
      xz(2,1)=.52
      xz(2,2)=.39
      xz(3,1)=.63
      xz(3,2)=-.14   ! relative (zum Mittelpunkt) Koordinaten der Ecken
      xz(4,1)=-.52
      xz(4,2)=-.39
      xz(5,1)=-.63
      xz(5,2)=.14
      step=2.68/2.
c
      open(1,file='detektor_geo.geo')
      do i=1,47
        xx=x0-step*(47-i)
        do j=1,5
          write(1,100) xx+xz(j,1),z0+xz(j,2),48-i
        end do
        write(1,101)
      end do
 100  format(2f14.3,i5)
 101  format('       *         *          * ')
      stop
      end
      
