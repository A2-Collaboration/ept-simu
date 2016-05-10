c      h e x a _ d e z . f
c
c      Programm setzt hexadezimale Zahlen in dezimale um.
c      Es werden die hexa-Stellen in steigenden Potenzen eingegeben
c      in dezimal umgesetzt und aufaddiert
c
c      Laden als  g77 -o hexa_dez.exe hexa_dez.f
c
      implicit none
      real hz,dz
      integer n
c
1     dz=0
      print*,'wird als hexa-Stelle ein Wert <0 eingegeben,'
      print*,'bedeutet das das obere Ende der Zahl'
      do n=0,7
        print*,' 16**',n
        read*,hz
        if(hz.lt.0)goto 2
        dz=dz+hz*16.**n
        print*,' bisher dz= ',dz
      end do
2     print*,'Wert der Dezimalzahl',dz
      goto 1
      end
      
