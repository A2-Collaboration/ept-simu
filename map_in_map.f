c
c       m a p _ i n _ m a p . f
c
c      Das Programm convert_map.f ordnet die Messwerte aus der Koordinaten-
c      Maschine in ein Feld mit Wert(1,1) dem ersten Wert der Messung.
c      Der erste Index gibt die Richtung parallel zur geraden Kante der
c      (Eintrittskante) des Endpoint-Taggers (ET), der zweite gibt die Richtung
c      des einlaufenden Strahls (abgesehen vom Drehwinkel des Taggers).
c      Der erste Index laeuft in Richtung der Ablenkung am Eintrittspunkt,
c      der zweite entgegen der Strahlrichtung.
c
c      Das hier vorliegende Programm gibt zu den Indices Laengen, die von einem
c      Nullpunkt (z.B. Strahleintritt an Magnetkante) aus gemessen werden.
c      Es werden also zu diesen Indices je Felder mit Laengen bestimmt, wie sie
c      sinnvoll entsprechend der Geometrie in der Halle sind.
c      Dabei steht der erste Index fuer die X-Achse in der Halle,
c      horizontal, senkrecht zur Strahlachse, vom Ursprung nach links
c      (Ablenkrichtung, des Haupttaggers).
c      Der zweite Index entspricht der z-Achse in der Halle (Strahlrichtung),
c      er zaehlt in Richtung des Strahls.
c
c      Das Programm liest also den Ergebnisfile aus convert_map.f ein
c      und setzt ihn entsprechend um. Hierfuer werden der Faktor
c      Winkel-Inkrement (Koordinaten-Maschine) in Laenge/mm und
c      der Drehwinkel des Taggers benoetigt.
c
c      Man erhaelt dann als Ergebnis ein Feld (x-Wert/cm, z-Wert/cm, Feld/Gauss)
