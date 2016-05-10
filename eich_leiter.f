c
c      e i c h _ l e i t e r . f
c
c      Programm benutzt ein Rechenergebnis aus brems_mc
c      z.B. gerechnet mit omap140.map fuer E_el 10 bis 170 step 2,
c      das in einer Liste  E/Mev, Lage_Max/Index_x, sigma/Index_x
c      steht.
c      Es ordnet es den Inedx_x-Orten der Detektoren zu (reine Gometrie):
c      1. Detektor bei Index_x 18.79, letzter (47.) bei 80.47,
c      alle Detektoren haben die Breite +- 0.62 Index_x, und den
c      Abstand 1.34 Index_x .
c      es werden bestimmt die Rand- und der Mitte-Energiewert.
c      Als Ergebnis im File eich_leiter.dat steht
c      # des Detektors (von 1-47), der Mittelwert/MeV, die Randwerte/MeV
c      und deren Differenz/MeV als energetische Detektorbreite.
c
c      Benoetigt wird also der Ergebisfile aus brems_mc (!in Index_x!)
c      dort fort.17 umbenannt z.B. in max_sig_14.dat.
c
c      Laden als:
c      g77 -fbounds-check -o eich_leiter.exe eich_leiter.f
c
      implicit none
      real max(100),sig(100),ebrems(100),emin,emax,mmin,mmax,
     &  ort,lage(3,47),elage(3,47)
      integer i,j,l,liste
      character*20 infile
c
      print*,'Bestimmung der Energiezuordnung fuer die 47 Zaehler'
      print*,'aus deren Geometrie und einer Eichrechnung mit'
      print*,'brems_mc.exe (Bestimmung der Verteilung der'
      print*,'PB-Elektronen fester Energien am Radiator),'
      print*,'die im einem File abgelegt ist als'
      print*,'Energie/MeV,Maximum und sigma des Resultats/Index_x'
      print*,'Name des Files aus brems_mc (z.B. max_eich_14.dat)'
      print*,'[a20:]'
      read '(a)',infile
      open(1,file=infile,status='old')
      i=1
1     read(1,*,err=2,end=2)ebrems(i),max(i),sig(i)
      if(i.eq.1)then
        emin=ebrems(i)
        mmin=max(i)
      end if
      i=i+1
      goto 1
2     liste=i-1
      close(1)
      emax=ebrems(liste)
      mmax=max(liste)
      print*,'File eingelesen, emin,emax,mmin,mmax'
      print*,emin,emax,mmin,mmax
      do i=1,47
        ort=18.79+1.34*(i-1)
        lage(1,i)=ort-0.62
        lage(2,i)=ort
        lage(3,i)=ort+0.62
        do j=1,3
          elage(j,i)=0.
        end do
      end do
      do j=1,47
        do l=1,3
          ort=lage(l,j)
          do i=1,liste-1
            if(max(i).lt.ort.and.max(i+1).gt.ort)then
              elage(l,j)=ebrems(i)+(ebrems(i+1)-ebrems(i))*
     &          (ort-max(i))/(max(1+i)-max(i))
            end if
          end do
        end do
      end do
      open(1,file='e_eich_leiter.dat')
      open(2,file='detektor_eichung.dat')
      do j=1,47
        write(1,*) elage(2,j),lage(2,j),(lage(3,j)-lage(1,j))/2,
     &    (elage(3,j)-elage(1,j))/2.
        write(2,101)j,elage(2,j),(elage(3,j)-elage(1,j))/2.,
     &    lage(2,j),(lage(3,j)-lage(1,j))/2.
101     format(I4,4f9.2)
      end do
      print*,'Das Ergebnis steht in zwei files:'
      print*,'in e_eich_leiter.dat als'
      print*,'              E/MeV,Ort/Index_x,dOrt,dE'
      print*,'in detektor_eichung.dat als'
      print*,'              Det #,E/MeV,dE,Ort/Index_x,dOrt '
      stop
      end
     
