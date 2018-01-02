#!/usr/bin/python
# -*- coding: utf-8 -*-

# Treść zadania:
#Zadanie 1 (termin oddania najpóźniej 8-I-2018). Korzystając z wiadomości zawartych podczas
#wykładu KE02 oraz KE03 napisać krótki skrypt w środowisku Matlab (Octave lub innym podobnym),
#który obliczy natężenie pola E w strefie dalekiej (np. dla r = 10 km) wytworzone przez ścieżki
#zasilające umieszczone na przykładowej płytce PCB. Kształt przewodnika przyjąć za Rys. 1. Należy
#założyć, że cały przewód jest krótki w porównaniu do długości fali, a więc amplituda prądu płynącego
#w każdym odcinku jest stała. Założyć, że amplituda prądu wynosi A = 20 mA, a długość pojedynczego
#odcinka przewodu B = 10 cm. Częstotliwość harmonicznego sygnału pobudzającego to 500 kHz, 50
#MHz, 500 MHz. Dlaczego przybliżenie oparte na prostych wzorach z wykładu załamuje się dla
#większych częstotliwości?

# Koncepcja rozwiązania:
#Przyjmując założenia, że
#1. cały przewód jest krótki w porównaniu do długości fali, a więc amplituda prądu płynącego
#w każdym odcinku jest stała,
#2. każdy z odcinków przewodu należy potraktować jak dipol Hertza,
#3. pole EM w dowolnym punkcie P w przestrzeni w strefie dalekiej jest sumą fal emitowanych przez wszystkie odcinki przewodu, 
#4. przyjmujemy UPROSZCZENIE, że poszczególne dipole nie oddziaływają na siebie,
# zadanie sprowadza się do obliczenia pola elektrycznego od każdego fragmentu przewodu, traktując taki fragment jako dipol Hertza, w zadanej odległości rx, gdzie x to numer fragmentu, a następnie zsumowanie wszystkich składowych pola elektrycznego.

from cmath import pi,sin,exp,sqrt
import matplotlib.pyplot as plt
import numpy as np

t = 0.0        # [s] moment w propagacji fali

class Constants:
    """Stałe potrzebne do jednoznacznego rozwiązania zadania."""
    c0 = 299792458    # [m/s] prędkość światła w próżni
    c = c0             # [m/s] prędkość światła w przyjętym ośrodku rozchodzenia się fali
    Eps0 = 1.0E-9/(36.0*pi)        # [F/m] przenikalność elektryczna próżni
    Eps = Eps0        # [F/m] przenikalność elektryczna ośrodka rozchodzenia się fali
    q0 = 1.602176E-19            # [C] ładunek elektryczny elementarny
    
    def __init__(self):
        return

# Funkcje globalne, powszechnie przydatne    
def beta(f):
    'Obliczenie liczby falowej (Bety) dla zadanej częstotliwości.'
    return 2.0*pi*f/Constants.c0

def I0(f):
    'Obliczenie prądu odniesienia (I0) dla zadanej częstotliwości.'
    return Constants.q0 * 2.0 * pi * f


class UkladPolarny:
    """
    Struktura modelująca układ współrzędnych polarnych.
    do zastosowania w rozwiązaniu zadania pierwszego.
    """
    
    def __init__(self,theta,phi,r):
        'Inicjalizacja położenia w układzie polarnym.'
        self.theta = theta        # [rad] 
        self.phi = phi            # [rad]
        self.r = r                # [m]

# Zauważam, ze pozostając w zgodzie z założeniami 
# zadaną pętlę z prądem można rozłożyć na pary 
# przewodów z prądami różnicowymi. 

class PetlaZPradem:
    'Model pętli z prądem złożonej z dipoli Hertza.'
    
    def __init__(self,dlugosc,polozenieWUkladziePolarnym):
        'Konstruktor pętli z prądem o podanej długości i położeniu w układzie polarnym.'
        self.dlugosc = dlugosc
        self.theta = polozenieWUkladziePolarnym.theta
        self.phi = polozenieWUkladziePolarnym.phi
        self.r = polozenieWUkladziePolarnym.r

class DipolHertza:
    'Model dipolu Hertza.'
    
    def __init__(self,dlugosc,polozenieWUkladziePolarnym):
        'Konstruktor dipola Hertrza o podanej długości i położeniu w układzie polarnym.'
        self.dlugosc = dlugosc                                    # [m]
        self.theta = polozenieWUkladziePolarnym.theta            # [rad]
        self.phi = polozenieWUkladziePolarnym.phi                # [rad]
        self.r = polozenieWUkladziePolarnym.r                    # [m]
        
    def E(self,f,punktWUkladziePolarnym,prad):
        """
        Zwraca wartość pola elektrycznego od dipola Hertza w punkcie zadanym w położeniu polarnym
        i prądzie płynącym w zadanym kierunku.
        Wartość pola elektrycznego jest obliczona jako absolutna względem polożenia dipola.
        """
        r = punktWUkladziePolarnym.r - self.r                # [m]
        theta = punktWUkladziePolarnym.theta - self.theta    # [rad]
        phi = punktWUkladziePolarnym.phi - self.phi            # [rad]
        I = prad                                            # [A] prąd dipola
        
        stalaWzoru = (I0(f) * self.dlugosc * beta(f)**2)\
            / (4.0 * pi * Constants.Eps * 2.0 * pi * f * r)\
            * sin(theta)
            
        wykladnik = exp(1j*(2.0*pi*f * t - beta(f)*r))

        wartoscPolaElektrycznego = I * stalaWzoru * wykladnik
        
        return wartoscPolaElektrycznego
        
        
# Przedstawienie badanego obwodu elektrycznego

class EmitujacyObwod:
    """
    Model analizowanego obwodu jako suma dipoli Hertza.
    
    Każdy z odcinków przewodu należy potraktować jak dipol Hertza. Pole EM w dowolnym
    punkcie P w przestrzeni w strefie dalekiej jest sumą fal emitowanych przez wszystkie odcinki
    przewodu (przyjmujemy UPROSZCZENIE, że poszczególne dipole nie oddziaływają na siebie).
    Ponieważ każdy z odcinków jest w nieco innym miejscu, więc docierające do punktu P fale są
    przesunięte względem siebie w fazie. Zmieniają się także ich amplitudy, ponieważ każdy z odcinków
    przewodu (dipoli Hertza) może być inaczej zorientowany w przestrzeni, a dipol Hertza nie jest
    promiennikiem izotropowym.
    
    Przyjmuję za biegun układu współrzędnych polarnych lewy środkowy róg obwodu.
    Zakładam, że układ leży w kącie radialnym phi = 0. Kąt theta równy 0 jest skierowany
    ,,do góry'', czyli od bieguna do dłuższego boku obwodu.
    Na poniższym szkicy oznaczony jako 'o'
    
        <--  3B -->
        ___.___.___
    (8)|    (1)    |(2)  ^__ B
       |___     ___|     v
       o(7)|   |(3)
        (6)|___|(4)
            <B>
            (5)
            
    Odcinki numeruję począwszy od górnego lewego rogu, zaczynając od 1. Zatem pierwszy odcinek
    będzie miał długość 3B, drugi prostopadły do niego B, pozostałe również B.
    """
    
    B = 10.0E-2        # [m] długość pojedynczego odcinka
    A = 20.0E-3        # [A] natężenie prądu w obwodzie
    
    f = (\
        500.0E3,    # pierwsza badana częstotliwość fali \
        50.0E6,        # druga badana częstotliwość fali \
        500.0E6     # trzecia badana częstotliwość fali \
        )
    
    def __init__(self):
        'Konstruktor analizowanego obwodu.'
        # Obliczam odległości polarne 'r' środków poszczególnych odcinków obwodu
        odleglosc_srodka_1 = sqrt(self.B*self.B + (1.5*self.B)**2).real
        odleglosc_srodka_2 = sqrt(3.0*self.B + (0.5*self.B)**2).real
        odleglosc_srodka_3 = 2.5*self.B
        odleglosc_srodka_4 = sqrt((2.0*self.B)**2 + (-0.5*self.B)**2).real
        odleglosc_srodka_5 = sqrt((1.5*self.B)**2 + (-1.0*self.B)**2).real
        odleglosc_srodka_6 = sqrt(self.B**2 + (-0.5*self.B)**2).real
        odleglosc_srodka_7 = 0.5*self.B
        odleglosc_srodka_8 = 0.5*self.B
        
        self.dipol_1 = DipolHertza( \
            dlugosc = 3.0*self.B, \
            polozenieWUkladziePolarnym = UkladPolarny(phi = 0,theta=0.5*pi,r=odleglosc_srodka_1) )
        self.dipol_2 = DipolHertza( \
            dlugosc = self.B, \
            polozenieWUkladziePolarnym = UkladPolarny(phi = 0, theta=pi, r=odleglosc_srodka_2) )
        self.dipol_3 = DipolHertza( \
            dlugosc = self.B, \
            polozenieWUkladziePolarnym = UkladPolarny(phi = 0, theta = 1.5*pi, r=odleglosc_srodka_3) )
        self.dipol_4 = DipolHertza( \
            dlugosc = self.B, \
            polozenieWUkladziePolarnym = UkladPolarny(phi = 0, theta = pi, r=odleglosc_srodka_4) )
        self.dipol_5 = DipolHertza( \
            dlugosc = self.B, \
            polozenieWUkladziePolarnym = UkladPolarny(phi = 0, theta = 1.5*pi, r=odleglosc_srodka_5) )
        self.dipol_6 = DipolHertza( \
            dlugosc = self.B, \
            polozenieWUkladziePolarnym = UkladPolarny(phi = 0, theta = 0, r=odleglosc_srodka_6) )
        self.dipol_7 = DipolHertza( \
            dlugosc = self.B, \
            polozenieWUkladziePolarnym = UkladPolarny(phi = 0, theta = 1.5*pi, r=odleglosc_srodka_7) )
        self.dipol_8 = DipolHertza( \
            dlugosc = self.B, \
            polozenieWUkladziePolarnym = UkladPolarny(phi = 0, theta = 0, r=odleglosc_srodka_8) )
            
        self.obwod = (\
            self.dipol_1,\
            self.dipol_2,\
            self.dipol_3,\
            self.dipol_4,\
            self.dipol_5,\
            self.dipol_6,\
            self.dipol_7,\
            self.dipol_8,\
            )
            
    def E(self,punktWUkladziePolarnym):
        """
        Zwraca wartość amplitudy pola elektrycznego, promieniowanego od obwodu,
        w zadanym w układzie polarnym punkcie.
        """
        wartoscPolaEWPunkcie = 0.0        # [V/m] wartość natężenia pola, inicjalizacja
        rozwiazania = []                # lista rozwiązań dla wszystkich podanych częstotliwości fali
        for fx in self.f:
            for dipol in self.obwod:
                wartoscPolaEWPunkcie += dipol.E(fx,punktWUkladziePolarnym,self.A)
            amplitudaPolaEWPunkcie = abs(wartoscPolaEWPunkcie)
            rozwiazania.append(amplitudaPolaEWPunkcie)
            
        return rozwiazania

            
# -----------------------------------------------------------
# Rozwiązanie zadania i rysowanie wyników na wykresie
# -----------------------------------------------------------

class Rozwiazanie:
    """
    Klasa reprezentująca rozwiązanie zadania pierwszego.
    
    Treść zadania:
    Zadanie 1 (termin oddania najpóźniej 8-I-2018). Korzystając z wiadomości zawartych podczas
    wykładu KE02 oraz KE03 napisać krótki skrypt w środowisku Matlab (Octave lub innym podobnym),
    który obliczy natężenie pola E w strefie dalekiej (np. dla r = 10 km) wytworzone przez ścieżki
    zasilające umieszczone na przykładowej płytce PCB. Kształt przewodnika przyjąć za Rys. 1. Należy
    założyć, że cały przewód jest krótki w porównaniu do długości fali, a więc amplituda prądu płynącego
    w każdym odcinku jest stała. Założyć, że amplituda prądu wynosi A = 20 mA, a długość pojedynczego
    odcinka przewodu B = 10 cm. Częstotliwość harmonicznego sygnału pobudzającego to 500 kHz, 50
    MHz, 500 MHz. 
    Dlaczego przybliżenie oparte na prostych wzorach z wykładu załamuje się dla
    większych częstotliwości?
    """

    rozdzielczoscTheta = 1            # [stopien] rozdzielczość za jaką będzie rozwiązywane zadanie, musi być typu integer
    
    def __init__(self,r,phi=0.0):
        'Przygotuj potrzebne obiekty i zmienne do obliczeń.'
        self.obwod = EmitujacyObwod()
        self.r = r                    # [m] odległość obliczanego rozwiązania od bieguna
        self.phi = phi                # [rad] kąt obrotowy obliczanego rozwiązania względem osi wiodącej
        self.listaKatow = [(pi/180.0)*kat for kat in range(0, 361, self.rozdzielczoscTheta)]    # [rad]
        
    def poleE(self):
        'Obliczenie natężenia pola E (składowa E(theta)) wokół układu przewodów.'
        self.listaKatow = [(pi/180.0)*kat for kat in range(0, 361, self.rozdzielczoscTheta)]    # [rad]
        self.listaPunktow = []         # lista punktów w układzie polarnym, dla których wykonuję obliczenia
        poleE = []                    # obliczone pole elektryczne wokół badanego obwodu w zadanej odległości i nachyleniu
        for x in range(len(self.listaKatow)):
            self.listaPunktow.append(\
                UkladPolarny(r=self.r,phi=self.phi,theta=self.listaKatow[x]) )
        for punkt in self.listaPunktow:
            poleE.append(self.obwod.E(punkt))
        return poleE
        
    #def rysujRozwiazanie(self):
        #'Rysuje rozwiązanie zadania pierwszego na wykresie polarnym.'
        #return 0.0
        
# Przydatne narzędzia
def transponuj(macierz):
    'Transponuje macierz.'
    return [list(i) for i in zip(*macierz)]

if __name__ == '__main__':
    
    print("W jakiej odległości od promieniującego obwodu obliczyć pole E? (np. 30.0E4)")
    odleglosc = float(input("Podaj wartość w metrach [m]: "))
    print("")
    
    rozwiazanie = Rozwiazanie(r=odleglosc)
    
    dane = rozwiazanie.poleE()
    dane = transponuj(dane)
    theta = rozwiazanie.listaKatow                
    
    indeks = 0
    
    for f in rozwiazanie.obwod.f:
        plt.figure()
        
        theta = np.asarray(theta)            # zamiana na typ array, oś theta
        r = np.asarray(dane[indeks])        # zamiana na typ array, oś r
        
        tytul = 'Czestotliwosc ' + '{:.2e}'.format(f) + '[Hz]'
        plt.title(tytul)
        plt.ylabel('Pole E [V/m]')
        plt.ticklabel_format(style='sci')
        plt.polar(theta,r)                # rysuj wykres polarny
        indeks += 1
        
    print("Przybliżenie oparte na prostych wzorach załamuje się dla większych częstotliwości,\n \
        ponieważ dla większych częstotliwości nie jest spełniony warunek na znacznie większą\n \
        długość fali od rozmiarów promieniującego obwodu. W takim wypadku proste przybliżenia nie\n \
        oddają dosyć dobrze rozkładu pól w strefie bliskiej obowdu. Na przykład dla częstotliwości\n \
        500 MHz długość fali w próżni to {:.3e} [m], co jest porównywalnym wymiarem z rozmiarami\n \
        badanego obwodu.".format(Constants.c0/500.0E6))
    
    plt.show()
