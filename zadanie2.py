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

from cmath import pi,log

class Constants:
	"""Stałe potrzebne do jednoznacznego rozwiązania zadania."""
	E = 10.0				# [V/m] natężenie pola elektrycznego oświetlającej fali płaskiej
	Z0 = 120.0*pi			# przenikalność próżni
	H = self.E/self.Z0		# natężenie pola magnetycznego oświetlającej fali płaskiej
	Eps0 = 1.0E-9/(36.0*pi)	# [F/m] przenikalność elektryczna próżni
	Eps = Eps0				# [F/m] przenikalność elektryczna ośrodka rozchodzenia się fali
	Mi0 = pi*4.0E-7			# [H/m] przenikalność magnetyczna próźni
	Mi = Mi0 				# [H/m] przenikalność magnetyczna ośrodka rozchodzenia się fali
	fmax = 					# [Hz] maksymalna częstotliwość analizy
	
	def __init__(self):
		return
		
class UkladKartezjanski:
	"""
	Struktura modelująca układ współrzędnych kartezjańskich.
	do zastosowania w rozwiązaniu zadania drugiego.
	"""
		
	def __init__(self,x,y,z):
		'Inicjalizacja położenia w układzie kartezjańskim.'
		self.x = x		# [m] 
		self.y = y		# [m]
		self.z = z		# [m]
		
def ortogonalnosc(v1,v2):
	'Sprawdzenie ortogalnosci wektorow w ukladzie kartezjanskim.'
	if ()

class FalaPlaska():
	"""
	Model fali płaskiej TEM rozchodzącej się w ośrodku izotropowym. Częstotliwość fali może być zewnętrznie zmieniana.
	"""
	
	def __init__(self,k,E,H,f=0.0):
		'Ustalenie właściwości fali płaskiej.'
		self.k = k				# kierunek propagacji fali w układzie kartezjańskim
		self.E = E 				# wektor zmian komponentu elektrycznego fali
		self.H = H 				# wektor zmian komponentu magnetycznego fali
		self.f = f				# [Hz] częstotliwość fali płaskiej
		
		# sprawdzenie czy fala jest TEM
		if (self.k.x*self.E.x + self.k.y*self.E.y + self.k.z*self.E.z != 0.0):
			raise self.Exception("Weself.ktory self.k i self.E nie są ortogonalne! To nie jest fala Tself.EM!")
		if (self.k.x*self.H.x + self.k.y*self.H.y + self.k.z*self.H.z != 0.0):
			raise self.Exception("Weself.ktory self.k i self.H nie są ortogonalne! To nie jest fala Tself.EM!")
		if (self.E.x*self.H.x + self.E.y*self.H.y + self.E.z*self.H.z != 0.0):
			raise Exception("Wektory E i H nie są ortogonalne! To nie jest fala TEM!")

class LiniaTEM():
	"""
	Model linii 2-przewodowej z zadania drugiego. Składa się z czterech fragmetów modelowanych
	jako linie krótkie z źródłami prądowymi i napięciowymi. 
	"""
	
	Rne = 50.0					# [Ohm] obciążenie linii TEM
	Rfe = 50.0 					# [Ohm] drugie obciążenie linii TEM
	
	def __init__(self):
		'Utworzenie linii 2-przewodowej z segmentów krótkich linii TEM.'
		
		self.segmenty = (\
			SegmentTEM(UkladKartezjanski(1.0,0.0,0.0)),\
			SegmentTEM(UkladKartezjanski(1.0,0.0,0.0)),\
			SegmentTEM(UkladKartezjanski(1.0,0.0,0.0)),\
			SegmentTEM(UkladKartezjanski(0.0,-1.0,0.0))\
			)
		
	
	
	class SegmentTEM():
		"""
		Segment krótkiej linii TEM, modelowany jako źródło napięciowe i prądowe.
		"""
		
		B = 20.0E-2				# [m] długość segmentu linii TEM
		s = 1.27E-3				# [m] odległość między przewodami linii TEM
		d = 1.00E-3				# [m] przekrój poprzeczny przewodu linii TEM
		
		def __init__(self, kierunekSegmentu):
			"""
			Konstruktor segmentu linii TEM, w którym obliczane są wszystkie jej właściwości.
			Kierunek segmentu ustalany jest wedle punktu odniesienia, który znajduje się na 
			wrotach segmentu. Wektor kierunku wskazuje kierunek zwrócenia segmentu.
			Na użytek zadania przyjmuję, że kierunek segmentu o składowej x != 0.0 oznacza
			ułożenie segmentu na płaszczyźnie OXZ, a o składowej y != 0.0 na płaszczyźnie OYZ.
			"""
			self.C = (pi*Constants.Eps) / (log((2.0*self.s) / (self.d)))	# [F] pojemność własna segmentu
			self.kierunek = kierunekSegmentu								# kierunek segmentu w układzie kartezjańskim
			self.H = 0.0													# pole magnetyczne indukujące napięcie w segmencie
			self.E = 0.0													# pole elektryczne powodujące przepływ prądu zmiennego w segmencie
			
		def V(self,falaPadajaca):
			'Obliczenie wartości źródła napięciowego spowodowanego zmiennym polem magnetycznym (wedle prawa Faradaya).'
			if self.kierunek.x != 0.0:
				self.H = falaPadajaca.H.y
			elif self.kierunek.y != 0.0:
				self.H = falaPadajaca.H.x
			else:
				self.H = 0.0
			V = 1j * 2.0*pi*f * Constants.Mi * self.H * self.B*self.s
			return V
			
		def I(self,falaPadajaca)
			'Oblicza wartość źródła prądowego spowodowanego zmiennym polem elektrycznym.'
			if self.kierunek.x != 0.0:
				self.E = falaPadajaca.E.z
			elif self.kierunek.y != 0.0:
				self.E = falaPadajaca.E.z
			else:
				self.E = 0.0
			I = 1j * 2.0*pi*f * self.C * self.E * self.B*self.s
			return I
			
		
