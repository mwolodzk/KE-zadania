#!/usr/bin/python
# -*- coding: utf-8 -*-

# Treść zadania
#Korzystając z wiadomości przekazanych podczas
#wykładu KE04 obliczyć napięcie indukowane w obciążeniach (R NE = R FE = 50 ) 2-przewodowej linii
#TEM zbudowanej z przewodów oddalonych od siebie o s = 1.27 mm i ułożonych na płaszczyźnie XY
#zgodnie z Rys. 2. Przyjąć, że linia jest oświetlana falą płaską i natężenie pola E wynosi 10 V/m. Przyjąć
#harmoniczne pobudzenie. Długość pojedynczego odcinka przewodu B = 20 cm.
#Wskazówka: podzielić obwód na 3 sekcje o długości B, dla każdej sekcji określić wartość pojemności
#linii (taka sama dla każdej sekcji) oraz wydajności źródeł prądowych i napięciowych (ze względu na
#różną orientację jednej z sekcji względem padającej fali własności źródeł tej sekcji będą inne). Na
#podstawie długości jednej sekcji oszacować maksymalną częstotliwość analizy i przeprowadzić
#obliczenia tylko do tej częstotliwości.

from math import log,pi
import matplotlib.pyplot as plt
import numpy as np

class Constants:
    """Stałe potrzebne do jednoznacznego rozwiązania zadania."""
    c0 = 299792458            # [m/s] prędkość światła w próżni
    c = c0                     # [m/s] prędkość światła w przyjętym ośrodku rozchodzenia się fali
    E = 10.0                # [V/m] natężenie pola elektrycznego oświetlającej fali płaskiej
    Z0 = 120.0*pi            # przenikalność próżni
    H = E/Z0                # natężenie pola magnetycznego oświetlającej fali płaskiej
    Eps0 = 1.0E-9/(36.0*pi)    # [F/m] przenikalność elektryczna próżni
    Eps = Eps0                # [F/m] przenikalność elektryczna ośrodka rozchodzenia się fali
    Mi0 = pi*4.0E-7            # [H/m] przenikalność magnetyczna próźni
    Mi = Mi0                 # [H/m] przenikalność magnetyczna ośrodka rozchodzenia się fali
    
    def __init__(self):
        return
        
class UkladKartezjanski:
    """
    Struktura modelująca układ współrzędnych kartezjańskich.
    do zastosowania w rozwiązaniu zadania drugiego.
    """
        
    def __init__(self,x,y,z):
        'Inicjalizacja położenia w układzie kartezjańskim.'
        self.x = x        # [m] 
        self.y = y        # [m]
        self.z = z        # [m]
        
#def ortogonalnosc(v1,v2):
    #'Sprawdzenie ortogalnosci wektorow w ukladzie kartezjanskim.'
    #if ()

class FalaPlaska():
    """
    Model fali płaskiej TEM rozchodzącej się w ośrodku izotropowym. Częstotliwość fali może być zewnętrznie zmieniana.
    """
    
    def __init__(self,k,E,H,f=1.0):
        'Ustalenie właściwości fali płaskiej.'
        self.k = k                # kierunek propagacji fali w układzie kartezjańskim
        self.E = E                 # wektor zmian komponentu elektrycznego fali
        self.H = H                 # wektor zmian komponentu magnetycznego fali
        self.f = f                # [Hz] częstotliwość fali płaskiej
        
        # sprawdzenie czy fala jest TEM
        if (self.k.x*self.E.x + self.k.y*self.E.y + self.k.z*self.E.z != 0.0):
            raise Exception("Wektory k i E nie są ortogonalne! To nie jest fala TEM!")
        if (self.k.x*self.H.x + self.k.y*self.H.y + self.k.z*self.H.z != 0.0):
            raise Exception("Wektory k i H nie są ortogonalne! To nie jest fala TEM!")
        if (self.E.x*self.H.x + self.E.y*self.H.y + self.E.z*self.H.z != 0.0):
            raise Exception("Wektory E i H nie są ortogonalne! To nie jest fala TEM!")

class LiniaTEM():
    """
    Model linii 2-przewodowej z zadania drugiego. Składa się z czterech fragmetów modelowanych
    jako linie krótkie z źródłami prądowymi i napięciowymi. 
    """
    
    Rne = 50.0                    # [Ohm] obciążenie linii TEM, na wrotach wejściowych pierwszego segmentu
    Rfe = 50.0                     # [Ohm] drugie obciążenie linii TEM, na wrotach wyjściowych ostatniego segmentu
    
    def __init__(self):
        'Utworzenie linii 2-przewodowej z segmentów krótkich linii TEM.'
        
        self.segmenty = (\
            self.SegmentTEM(UkladKartezjanski(1.0,0.0,0.0)),\
            self.SegmentTEM(UkladKartezjanski(1.0,0.0,0.0)),\
            self.SegmentTEM(UkladKartezjanski(1.0,0.0,0.0)),\
            self.SegmentTEM(UkladKartezjanski(0.0,-1.0,0.0))\
            )
        
    def Vne(self, falaPadajaca):
        """
        Obliczenie napięcia indukowanego przez falę padającą na linię TEM
        na obciążeniu Rne (na wrotach wejściowych pierwszego segmentu).
        Z uwagi na uproszczony, czyli pozbawiony elementów biernych, model 
        linii, mogę prosto zsumować wszystkie źródła napięciowe i prądowe w obwodzie.
        """
        V = 0.0                # [V] napięcie wypadkowe indukowany przez pole magnetyczne na linii TEM
        I = 0.0                # [I] prąd wypadkowy powodowany przez pole elektryczne na linii TEM
        Rne = self.Rne
        Rfe = self.Rfe
        for segment in self.segmenty:
            V += segment.V(falaPadajaca)
            I += segment.I(falaPadajaca)
        Vne = (Rne)/(Rne+Rfe)*V - (Rne*Rfe)/(Rne+Rfe)*I
        return Vne
        
    def Vfe(self, falaPadajaca):
        """
        Obliczenie napięcia indukowanego przez falę padającą na linię TEM
        na obciążeniu Rfe (na wrotach wyjściowych ostatniego segmentu).
        Z uwagi na uproszczony, czyli pozbawiony elementów biernych, model 
        linii, mogę prosto zsumować wszystkie źródła napięciowe i prądowe w obwodzie.
        """
        V = 0.0                # [V] napięcie wypadkowe indukowany przez pole magnetyczne na linii TEM
        I = 0.0                # [I] prąd wypadkowy powodowany przez pole elektryczne na linii TEM
        Rne = self.Rne
        Rfe = self.Rfe
        for segment in self.segmenty:
            V += segment.V(falaPadajaca)
            I += segment.I(falaPadajaca)
        Vfe = -(Rfe)/(Rne+Rfe)*V - (Rne*Rfe)/(Rne+Rfe)*I
        return Vfe
    
    class SegmentTEM():
        """
        Segment krótkiej linii TEM, modelowany jako źródło napięciowe i prądowe.
        """
        
        B = 20.0E-2                # [m] długość segmentu linii TEM
        s = 1.27E-3                # [m] odległość między przewodami linii TEM
        d = 1.00E-3                # [m] przekrój poprzeczny przewodu linii TEM
        
        def __init__(self, kierunekSegmentu):
            """
            Konstruktor segmentu linii TEM, w którym obliczane są wszystkie jej właściwości.
            Kierunek segmentu ustalany jest wedle punktu odniesienia, który znajduje się na 
            wrotach segmentu. Wektor kierunku wskazuje kierunek zwrócenia segmentu.
            Na użytek zadania przyjmuję, że kierunek segmentu o składowej x != 0.0 oznacza
            ułożenie segmentu na płaszczyźnie OXZ, a o składowej y != 0.0 na płaszczyźnie OYZ.
            """
            self.C = (pi*Constants.Eps) / (log((2.0*self.s) / (self.d)))    # [F] pojemność własna segmentu
            self.kierunek = kierunekSegmentu                                # kierunek segmentu w układzie kartezjańskim
            self.H = None                                                    # pole magnetyczne indukujące napięcie w segmencie
            self.E = None                                                    # pole elektryczne powodujące przepływ prądu zmiennego w segmencie
            self.f = None                                                    # częstotliwość fali oświetlającej segment
            
        def V(self,falaPadajaca):
            'Obliczenie wartości źródła napięciowego spowodowanego zmiennym polem magnetycznym (wedle prawa Faradaya).'
            if self.kierunek.x != 0.0:
                self.H = falaPadajaca.H.y * Constants.H
            elif self.kierunek.y != 0.0:
                self.H = falaPadajaca.H.x * Constants.H
            else:
                self.H = 0.0
            f = falaPadajaca.f
            V = 1j * 2.0*pi*f * Constants.Mi * self.H * self.B*self.s
            return V
            
        def I(self,falaPadajaca):
            'Oblicza wartość źródła prądowego spowodowanego zmiennym polem elektrycznym.'
            if self.kierunek.x != 0.0:
                self.E = falaPadajaca.E.z * Constants.E
            elif self.kierunek.y != 0.0:
                self.E = falaPadajaca.E.z * Constants.E
            else:
                self.E = 0.0
            f = falaPadajaca.f
            I = 1j * 2.0*pi*f * self.C * self.E * self.B*self.s
            return I
            
# -----------------------------------------------------------
# Rozwiązanie zadania i rysowanie wyników na wykresie
# -----------------------------------------------------------

class Rozwiazanie:
    """
    Klasa reprezentująca rozwiązanie zadania drugiego.
    
    Treść zadania:
    Korzystając z wiadomości przekazanych podczas
    wykładu KE04 obliczyć napięcie indukowane w obciążeniach (R NE = R FE = 50 Ohm) 2-przewodowej linii
    TEM zbudowanej z przewodów oddalonych od siebie o s = 1.27 mm i ułożonych na płaszczyźnie XY
    zgodnie z Rys. 2. Przyjąć, że linia jest oświetlana falą płaską i natężenie pola E wynosi 10 V/m. Przyjąć
    harmoniczne pobudzenie. Długość pojedynczego odcinka przewodu B = 20 cm.
    """
    
    punktyNaDekadeF = 10            # [Hz] ilość punktów na dekadę częstotliwości w których obliczane będą napięcia na linii TEM
    fmin = 1.0                        # [Hz] częstotliwość początkowa obliczeń
    
    def __init__(self):
        """
        Przygotuj potrzebne obiekty i zmienne do obliczeń. 
        Do obliczenia częstotliwości maksymalnej symulacji zakładam warunek, że 
        minimalna długość fali padającej musi być większa od ośmiu długości pojedynczego
        segmentu, czyli lambda > 8 * B.
        """
        self.linia = LiniaTEM()
        self.fmax = Constants.c / (8.0*LiniaTEM.SegmentTEM.B)
        self.fala = FalaPlaska(\
            k = UkladKartezjanski(x = 1.0, y = 1.0, z = 0.0),\
            E = UkladKartezjanski(x = -1.0, y = 1.0, z = 0.0),\
            H = UkladKartezjanski(x = 0.0, y = 0.0, z = 1.0),\
            f = self.fmin
        )
        
    def Vne(self):
        'Wyznaczenie zależności napięcia na obciążeniu Rne od częstotliwości podanej fali padającej.'
        fstart = round(log(self.fmin,10))        # [Hz] wykładnik częstotliwości początkowej obliczeń
        fstop = round(log(self.fmax,10))        # [Hz] wykładnik częstotliwości końcowej obliczeń
        self.f = np.logspace(\
            start=fstart,\
            stop=fstop,\
            num=self.punktyNaDekadeF * fstop
            )
        Vne = []                                # [V] lista obliczonych napięć na obciążeniu Rne
        for f in self.f:
            self.fala.f = f
            Vne.append(self.linia.Vne(self.fala))
        return Vne
        
    def Vfe(self):
        'Wyznaczenie zależności napięcia na obciążeniu Rfe od częstotliwości podanej fali padającej.'
        fstart = round(log(self.fmin,10))        # [Hz] wykładnik częstotliwości początkowej obliczeń
        fstop = round(log(self.fmax,10))        # [Hz] wykładnik częstotliwości końcowej obliczeń
        self.f = np.logspace(\
            start=fstart,\
            stop=fstop,\
            num=self.punktyNaDekadeF * fstop
            )
        Vfe = []                                # [V] lista obliczonych napięć na obciążeniu Rne
        for f in self.f:
            self.fala.f = f
            Vfe.append(self.linia.Vfe(self.fala))
        return Vfe

# Przydatne narzędzia
def amplitudy(listaWartosciZespolonych):
    'Zwraca listę amplitud liczb zespolonych przesłanych w postaci obiektu wyliczeniowego.'
    absolutne = []
    for wartosc in listaWartosciZespolonych:
        absolutne.append(abs(wartosc))
    return absolutne

if __name__ == '__main__':
    
    rozwiazanie = Rozwiazanie()
    
    # Z uwagi na niejasność orientacji fali padającej i obwodu oświetlanego.
    fala = FalaPlaska(\
            k = UkladKartezjanski(x = 1.0, y = 1.0, z = 1.0),\
            E = UkladKartezjanski(x = -0.5, y = -0.5, z = 1.0),\
            H = UkladKartezjanski(x = -1.0, y = 1.0, z = 0.0),\
            f = rozwiazanie.fmin
    )
    rozwiazanie.fala = fala
    rozwiazanie.Rne = 573.0
    # ---------------------------------------------
    
    Vne = rozwiazanie.Vne()
    Vfe = rozwiazanie.Vfe()
    f = rozwiazanie.f
    
    Vne = amplitudy(Vne)
    Vfe = amplitudy(Vfe)
    
    Vne = np.asarray(Vne)
    Vfe = np.asarray(Vfe)
    f = np.asarray(f)
    
    plt.figure()
    
    plt.subplot(211)
    plt.plot(f,Vne, label='Napiecie na Rne')
    plt.ylabel('Vne [V]')
    plt.xlabel('f [Hz]')
    plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    plt.legend(loc='upper center')
    
    plt.title('Napiecia indukowane na obciazeniach')
    
    plt.subplot(212)
    plt.plot(f,Vfe,label='Napiecie na Rfe')
    plt.ylabel('Vfe [V]')
    plt.xlabel('f [Hz]')
    plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    plt.legend(loc='upper center')
    
    plt.show()
    
