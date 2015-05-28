import random
import molecules_mod1 as mol
import numpy.random as npr
import time
from Input.KnowledgeBase import KnowledgeBase as kb
import model


class Process(object):
    """
    Parent for all cellular processes.
    """
    def __init__(self, id, name):
        self.id = id
        self.name = name


        self.enzyme_ids = []
        self.substrate_ids = []

    def set_states(self, substrate_ids, enzyme_ids):
        self.enzyme_ids = enzyme_ids
        self.substrate_ids = substrate_ids

    def update(self, model):
        """
        Has to be implemented by child class.
        """
        pass


class Translation(Process):
    """
    Translation is instantiated in the Cell to produce proteins.

    Defines Translation process. It iterates over all ribosomes and decides what
    they should do. They either bind to mRNA or elongate/terminate a protein if
    they are already bound.

    """
    code = dict([('UCA', 'S'), ('UCG', 'S'), ('UCC', 'S'), ('UCU', 'S'),
                 ('UUU', 'F'), ('UUC', 'F'), ('UUA', 'L'), ('UUG', 'L'),
                 ('UAU', 'Y'), ('UAC', 'Y'), ('UAA', '*'), ('UAG', '*'),
                 ('UGU', 'C'), ('UGC', 'C'), ('UGA', '*'), ('UGG', 'W'),
                 ('CUA', 'L'), ('CUG', 'L'), ('CUC', 'L'), ('CUU', 'L'),
                 ('CCA', 'P'), ('CCG', 'P'), ('CCC', 'P'), ('CCU', 'P'),
                 ('CAU', 'H'), ('CAC', 'H'), ('CAA', 'Q'), ('CAG', 'Q'),
                 ('CGA', 'R'), ('CGG', 'R'), ('CGC', 'R'), ('CGU', 'R'),
                 ('AUU', 'I'), ('AUC', 'I'), ('AUA', 'I'), ('AUG', 'M'),
                 ('ACA', 'T'), ('ACG', 'T'), ('ACC', 'T'), ('ACU', 'T'),
                 ('AAU', 'N'), ('AAC', 'N'), ('AAA', 'K'), ('AAG', 'K'),
                 ('AGU', 'S'), ('AGC', 'S'), ('AGA', 'R'), ('AGG', 'R'),
                 ('GUA', 'V'), ('GUG', 'V'), ('GUC', 'V'), ('GUU', 'V'),
                 ('GCA', 'A'), ('GCG', 'A'), ('GCC', 'A'), ('GCU', 'A'),
                 ('GAU', 'D'), ('GAC', 'D'), ('GAA', 'E'), ('GAG', 'E'),
                 ('GGA', 'G'), ('GGG', 'G'), ('GGC', 'G'), ('GGU', 'G')])

    def __init__(self, id, name):
        super(Translation, self).__init__(id, name)

        # declare attributes
        self.__ribsomes = []



    def update(self, model):
        """
        Update all mrnas and translate proteins.
        """
        self.__ribosomes = model.states[self.enzyme_ids[0]]
        for mrna_id in self.substrate_ids:
            prot = None
            mrna = model.states[mrna_id]
            if not mrna.binding[0]:
                self.initiate(mrna)
            else:
                prot = self.elongate(mrna)
            if isinstance(prot, molecules.Protein):
                if prot.id in model.states:
                    model.states[prot.id].append(prot)
                else:
                    model.states[prot.id] = [prot]

    def initiate(self, mrna):
        """
        Try to bind to a given MRNA. Binding probability corresponds to the ribosome count.

        @type mrna: MRNA
        """
        if not mrna.binding[0]:  #  no mrna bound yet and target mrna still free at pos 0
            # bind a nascent protein to the 0 codon
            if npr.poisson(self.__ribosomes.count) > 1: # at least one binding event happens in time step
                # if a ribosome binds the position a new Protein is created and stored on the
                # position as if it were bound to the ribosome
                mrna.binding[0] =  molecules.Protein("Protein_{0}".format(mrna.id),
                                                     "Protein_{0}".format(mrna.id),
                                                     self.code[mrna[0:3]])
                self.__ribosomes.count -= 1

    def elongate(self, mrna):
        """
        Elongate the new protein by the correct amino acid. Check if an
        MRNA is bound and if ribosome can move to next codon.
        Terminate if the ribosome reaches a STOP codon.

        @type return: Protein or False
        """

        # TODO: this needs to update in a random order
        for i, ribosome in enumerate(mrna.binding):
            if isinstance(ribosome, molecules.Protein):
                codon = mrna[i*3:i*3+3]
                aa = self.code[codon]
                if aa == "*":  # terminate at stop codon
                    return self.terminate(mrna, i)

                if not mrna.binding[i + 1]:  # if the next rna position is free
                    mrna.binding[i] + aa
                    mrna.binding[i + 1] = mrna.binding[i]
                    mrna.binding[i] = 0
        return 0

    def terminate(self, mrna, i):
        """
        Splits the ribosome/MRNA complex and returns a protein.

        @type mrna: MRNA
        """

        protein = mrna.binding[i]  # bound mRNA
        mrna.binding[i] = 0
        self.__ribosomes.count += 1
        return protein








class Transcription(Process):
    """
    Implements mRNA transcription from genes on the chromosome.
    """

    def __init__(self, id, name):
        super(Transcription, self).__init__(id, name)

        self.__mRNAs = []




        
    def transcribe(self):

        #global MRNA_Liste

        StopCodonDict = ['TAG', 'TAA' ,'TGA']

        StartCodon = 'ATG'

        DNA = ''


        while DNA[0:3]!=StartCodon:
            GeneNumber = npr.randint(1,10)
            SetInd = 3-len(str(GeneNumber))
            GeneInd = 'MG_'+'0'*SetInd+str(GeneNumber)
            GetSeq = kb()
            try:

                DNA = GetSeq.get_sequence(GeneInd)
            except:
                pass

        mRNA = mol.MRNA(0,'','',0,0)

        #mRNA.mid = int(GeneNumber)

        mRNA.name = 'mRNA_'+str(GeneNumber)

        TransTime = float(len(DNA))*0.075



        StopCodonDict = ['TAG', 'TAA' ,'TGA']

        StartCodon = 'ATG'

        j = 0


        if str(DNA[j*3:j*3+3]) == StartCodon:
            while (3*j)<(len(DNA)) and str(DNA[j*3:j*3+3]) not in StopCodonDict:
                mRNA.sequence += str(DNA[j*3:j*3+3])
                j=j+1
            mRNA.sequence += str(DNA[j*3:j*3+3])
            mRNA.sequence = mRNA.sequence.replace('T','U')
        else:
            while str(DNA[j*3:j*3+3]) not in StartCodon:
                j = j+1
            while (3*j)<(len(DNA)):
                if str(DNA[j*3:j*3+3]) not in StopCodonDict:
                    mRNA.sequence += str(DNA[j*3:j*3+3])
                    j=j+1
            mRNA.sequence += str(DNA[j*3:j*3+3])
            mRNA.sequence = mRNA.sequence.replace('T','U')


        #if  mRNA.sequence in MRNA_Liste:
            #SEQ_Pos = MRNA_Liste.index(mRNA.sequence)
            #MRNA_Liste[SEQ_Pos+2] = MRNA_Liste[SEQ_Pos+2]+1
        #else:
            #MRNA_Liste += [mRNA.mid, mRNA.name, mRNA.sequence, TransTime, mRNA.count]

        
        #print MRNA_Liste
        #print str(TransTime) + "s"

        return mRNA



    def update(self, model):
        prot = None

        j=0

        while j < 20:
            prot = self.transcribe()
            if prot.name in model.states:
                model.states[prot.name].append(prot)
            else:
                model.states[prot.name] = [prot]
            j = j+1









if __name__ == "__main__":
    c = model.Model()
    MRNA_Liste = []
    TR = Transcription(1,'digga')
    TR.update(c)
    #TR.update_trans(c)




    




    # TODO: implement transcription