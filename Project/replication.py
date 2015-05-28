import random 
import numpy.random as nr
import processes as proc
import molecules as mol
from processes import Process
from molecules import BioMoleculeCount
from molecules import BioMolecule

'''
Klassen, mit Objekten (DNA, Helicase, Polymerase), Prozess (Replikation). Elternklasse (z.B. BioMolecule, BioMoleculeCount). 
Elternklasse vererbt
'''
class DNA(BioMolecule):
    def __init__(self, mid, name, length, nucleotides, mass=0):
        super(DNA, self).__init__(mid, name, mass)
        self._length=length/2
        self._nucleotides = nucleotides

    @property
    def nucleotides(self):
        return self._nucleotides

    @nucleotides.setter
    def nucleotides(self, value):
        self._nucleotides = value

    @property
    def length(self):
        return self._length

    @length.setter
    def length(self, value):
        self._length = value

class PolymeraseIII(BioMoleculeCount):
	#position=int; bound=boolean
    def __init__(self, mid, name, count=0, position=0, bound=False):
        super(PolymeraseIII, self).__init__(mid, name, count)
        self._position=position
        self._bound=bound

    @property
    def bound(self):
        return self._bound

    @bound.setter
    def bound(self, value):
        self._bound = value

    @property
    def position(self):
        return self._position

    @position.setter
    def position(self, value):
        self._position = value

class Helicase(BioMoleculeCount):
	#position=int; bound=boolean
    def __init__(self, mid, name, count=0, position=0, bound=False):
        super(Helicase, self).__init__(mid, name, count)
        self._position=position
        self._bound = bound

    @property
    def bound(self):
        return self._bound

    @bound.setter
    def bound(self, value):
        self._bound = value

    @property
    def position(self):
        return self._position

    @position.setter
    def position(self, value):
        self._position = value


class Replication(Process):

    def __init__(self, id, name, ATP, NT, double=False):
        super(Replication, self).__init__(id, name)
        self._double = double
        self.ATP_molecules = ATP
        self.Nucleotide = NT

    @property
    def double(self):
        return self._double

    @double.setter
    def double(self, value):
        self._double = value


    def update(self, model):
    	#return position helicase und polymeraseIII
        #Uebergeben der Objekte, vom Ori aus wird in zwei Richtungen repliziert (anticlock und clock)
        self.PolymeraseIII_anticlock = model.states['PolymeraseIII_0']
        self.PolymeraseIII_clock = model.states['PolymeraseIII_1']
        self.Helicase_anticlock = model.states['Helicase_0']
        self.Helicase_clock = model.states['Helicase_1']
        self.DNA = model.states['DNA']

        if self.double ==False: # Replikation darf nur stattfinden, wenn DNA noch nicht veerdoppelt ist
        #in anticlock_Richtung
            if (self.PolymeraseIII_anticlock.bound == False) and (self.PolymeraseIII_anticlock.position <self.DNA.length):
                self.Helicase_anticlock, self.PolymeraseIII_anticlock=self.initiate(DNA, self.PolymeraseIII_anticlock, self.Helicase_anticlock)
            elif self.Helicase_anticlock.bound == True:
                self.Helicase_anticlock, self.PolymeraseIII_anticlock=self.elongate(DNA, self.PolymeraseIII_anticlock, self.Helicase_anticlock)   
            else:
                self.Helicase_anticlock, self.PolymeraseIII_anticlock=self.terminate(DNA, self.PolymeraseIII_anticlock, self.Helicase_anticlock)    
        #in clock_Richtung
            if (self.PolymeraseIII_clock.bound == False) and (self.PolymeraseIII_anticlock.position <self.DNA.length):
                self.Helicase_clock, self.PolymeraseIII_clock = self.initiate(DNA, self.PolymeraseIII_clock, self.Helicase_clock)
            elif self.Helicase_clock.bound == True:
                self.Helicase_clock, self.PolymeraseIII_clock = self.elongate(DNA, self.PolymeraseIII_clock, self.Helicase_clock)
            else:
                self.Helicase_clock, self.PolymeraseIII_clock = self.terminate(DNA, self.PolymeraseIII_clock, self.Helicase_clock)
        #wenn beide Richtungen komplett sind, wird double auf True gesetzt, um weitere Replikation zu verhindern
            if (self.PolymeraseIII_anticlock.position >= self.DNA.length) and (self.PolymeraseIII_clock.position >= self.DNA.length):
                self.double = True

        self.DNA.nucleotides = 4*self.DNA.length + 2*self.PolymeraseIII_anticlock.position + 2*self.PolymeraseIII_clock.position
        #print self.DNA.nucleotides
        #print self.Helicase_anticlock.position, self.PolymeraseIII_anticlock.position, self.Helicase_clock.position, self.PolymeraseIII_clock.position, self.DNA.nucleotides

    def initiate(self, DNA, Pol, Hel):
    	"""
        wird aufgerufen, wenn Polymerase noch nicht gebunden ist, Helikase bindet mit def. Bindungswahrscheinlichkeit und startet 
        Strangauftrennung (und damit ATP-Verbrauch). Falls Abstand Helikase-Polymerase Mindestabstand erreicht hat, bindet Polymerase mit def. Bindungswahrscheinlichkeit.
        """
        Helicase = Hel
        PolymeraseIII =Pol
        if Helicase.bound == False:
            x_number =nr.randint(1,10)#Bindungswahrscheinlichkeit Helikase
            if x_number ==1:
                Helicase.bound =True
        elif self.ATP_molecules >= 100 and (Helicase.position - PolymeraseIII.position) < 3000:
            Helicase.position += 100
            self.ATP_molecules -= 100
        elif self.ATP_molecules > 0 and (Helicase.position - PolymeraseIII.position) < 3000:
            Helicase.position += self.ATP_molecules
            self.ATP_molecules -= self.ATP_molecules
        if Helicase.position >=1500: # 1500 ist Mindestabstand zwischen Helicase und PolyIII
            y_number =nr.randint(1,5)# Bindungswahrscheinlichkeit Polymerase
            if y_number ==1:
                PolymeraseIII.bound = True

        if Helicase.position > self.DNA.length:
            self.ATP_molecules=self.ATP_molecules+(Helicase.position -self.DNA.length)
            Helicase.position = self.DNA.length
        #print ('ATP:',self.ATP_molecules,'NT:',self.Nucleotide)
        return Helicase, PolymeraseIII

    def elongate(self,DNA, Pol, Hel):
        """
        Wird aufgerufen, wenn Polymerase und Helicase gebunden sind. Testet, ob genug ATP Molekuele und Nukelotide vorhanden sind.
        Verlaengert pro Step um 100 Nukelotide oder die maximal moegliche Anzahl bei ATP/Nukleotid Begrenzung.
        Der maximale Abstand zwischen Helikase und Polymerase ist 3000, der minimale 1500.
        """
        Helicase = Hel
        PolymeraseIII = Pol
        if self.ATP_molecules >= 100 and (Helicase.position - PolymeraseIII.position) < 3000: #genug ATP, Abstand klein genug
            Helicase.position += 100 
            self.ATP_molecules -= 100
            if self.Nucleotide >= 200 and (Helicase.position - PolymeraseIII.position) > 1500: #genug  Nucleotide (>=200)
                PolymeraseIII.position += 100
                self.Nucleotide -= 200
            elif self.Nucleotide > 1 and (Helicase.position - PolymeraseIII.position) > 1500: #nicht genug Nucleotide (1-199)
                PolymeraseIII.position += self.Nucleotide/2
                Helicase.position = Helicase.position -100 +self.Nucleotide/2
                self.ATP_molecules =self.ATP_molecules+100-self.Nucleotide/2
                self.Nucleotide -= 2*(self.Nucleotide/2)
        
        elif self.ATP_molecules >= 0 and (Helicase.position - PolymeraseIII.position) < 3000: #nicht genug ATP, Abstand klein genug
            Helicase.position += self.ATP_molecules
            if self.Nucleotide >= 200 and (Helicase.position - PolymeraseIII.position) > 1500:  #genug Nucleotide
                PolymeraseIII.position += 100
                self.Nucleotide -= 200
            elif self.Nucleotide > 1 and (Helicase.position - PolymeraseIII.position) > 1500: #nicht genug Nucleotide
                PolymeraseIII.position += self.Nucleotide/2
                Helicase.position = Helicase.position -self.ATP_molecules +self.Nucleotide/2
                self.ATP_molecules -=self.Nucleotide/2
                self.Nucleotide -= 2*(self.Nucleotide/2)
            self.ATP_molecules -= self.ATP_molecules

        if Helicase.position > self.DNA.length:
            self.ATP_molecules=self.ATP_molecules+(Helicase.position -self.DNA.length)
            Helicase.position = self.DNA.length

        if Helicase.position >= self.DNA.length:
            Helicase.bound =False
        #print ('ATP:',self.ATP_molecules,'NT:',self.Nucleotide)
        return Helicase, PolymeraseIII



    def terminate(self,DNA, Pol, Hel):
        """
        Beendet die Replikation. Wird aufgerufen, wenn die Helicase nicht mehr gebunden ist.
        Wenn genug Nucleotide vorhanden sind, wird die Polymerase pro Step um 100 Nukelotide verschoben,
        sonst um die  maximal moegliche Anzahl. Ist die Mitte der DNA erreicht, wird die Polymerase abgeloest
        """
        Helicase = Hel
        PolymeraseIII = Pol
        #print self.DNA.length, PolymeraseIII.position
    	#wenn pos= length helicase + polyIII fallen ab
        if PolymeraseIII.bound == True:
            if self.Nucleotide >= 200 :
                PolymeraseIII.position += 100
                self.Nucleotide -= 200
            elif self.Nucleotide>1 :
                PolymeraseIII.position += self.Nucleotide/2
                self.Nucleotide -= 2*self.Nucleotide/2

        if PolymeraseIII.position >= self.DNA.length:
            self.Nucleotide += (PolymeraseIII.position-self.DNA.length)*2
            PolymeraseIII.position=self.DNA.length
            PolymeraseIII.bound = False

        return Helicase, PolymeraseIII

    def gene_check(self,DNA,Pol_ac,Pol_c,gene_begin,gene_end): 
        """
        Testet, ob ein Gen schon doppelt oder einfach vorliegt. Uebergeben werden muessen die DNA-Polymerasen,
        sowie Start und Endpunkt des Gens auf dem Strang.
        return=2 fuer Gene, die schon repliziert wurden
        return=1 fuer noch nicht repliziert wurden
        """
        PolymeraseIII_ac = Pol_ac
        PolymeraseIII_c = Pol_c
        if (gene_end < PolymeraseIII_c.position) or (gene_begin > (2*self.DNA.length-PolymeraseIII_ac.position)):
            return 2
        else:
            return 1
