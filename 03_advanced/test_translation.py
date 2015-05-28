import random
import sys

class BioMolecule(object):
    """
    A generic molecule that has basic attributes like id, name and
    mass.

    @type id: int
    @type name: str
    @type mass: float
    """
    def __init__(self, id, name, mass=None):
        self._id = id
        self.name = name
        self.mass = mass

    @property
    def _id(self):
        return self.__id
    
    @_id.setter
    def _id(self, value):
        if not isinstance(value, int):
            raise TypeError("Id must be Integer.")
        self.__id = value

    @property
    def name(self):
        return self.__name
    
    @name.setter
    def name(self, value):
        if not isinstance(value, str):
            raise TypeError("Name must be String.")
        self.__name = value

    @property
    def mass(self):
        return self.__mass
    
    @mass.setter
    def mass(self, value):
        if not isinstance(value, float):
            raise TypeError("Mass must be Float.")
        self.__mass = value
    
class Polymer(BioMolecule):
    """
    A polymer molecule that has a sequence attribute which is
    accessible via indexing the object. 

    @type id: int
    @type name: str
    @type sequence: str
    @type mass: float
    """
    
    def __init__(self, id, name, sequence, mass=None):
        
       super(Polymer, self).__init__(id, name, mass)
       self._sequence = sequence

    @property
    def _sequence(self):
        return self.__sequence
    
    @_sequence.setter
    def _sequence(self, value):
        if not isinstance(value, str):
            raise TypeError("Mass must be Float.")
        self.__sequence = value    
    
    def __getitem__(self, value):
        """
        Makes the sequence accessible via the indexing operators:
<        p[10] returns the tenth character in the sequence.
        """
        return self.sequence[value]

    def __setitem__(self, key, value):
        """
         Enables changing of sequence characters via the indexing operators.       
        """
        self.sequence[key] = value

class MRNA(Polymer):
    def __init__(self, id, name, sequence, mass=None):
        super(MRNA, self).__init__(id, name, sequence, mass)
        
        # 7. Create a list that stores if a ribosome is bound for each
        # codon (triplet).
        self.binding = [0] * (len(sequence)/3)# use this attribute for 7.

    def calculate_mass(self):
        NA_mass = {'A': 1.0, 'U': 2.2, 'G':2.1, 'C':1.3}
        codon_dict = {'A':0, 'B'}
        while i in len(sequence):


protein = BioMolecule(111, 'protein', 100.0)
print protein._id
protein2 = Polymer(123, 'protein2', 'AGAG', 100.0)
print protein2._sequence