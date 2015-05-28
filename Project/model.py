import processes_mod2 as proc
import molecules as mol

class Model(object):
    """
    Initializes the states and processes for the model and lets the processes update their corresponding states.
    """
    def __init__(self):
        self.states = {}
        self.processes = {}

        # initiate states
        self.ribosomes = {'Ribosomes': mol.Ribosome('Ribosomes', 'Ribosomes', 10)}
        self.mrnas = {}

       
        #self.mrnas = {'MRNA_1':mol.MRNA(1,'MRNA_1',"AGTGGTATTATATAG")}
#                    {MRNA_1, 2, 3 usw...          : ID, Name  (MRNA_1, 2, 3...) , sequence }     von 0 bis 50
#  strings zusammenbauen "{0} {1}".format("a","b")
       


        self.states.update(self.ribosomes) # updates 
        self.states.update(self.mrnas) #

        # initiate processes
        translation = proc.Translation(1, "Translation")
        transcription = proc.Transcription(1, "Transcription")
        

        translation.set_states(self.mrnas.keys(), self.ribosomes.keys())
        

        self.processes = {"Translation":translation, "Transcription":transcription}

    def step(self):
        """
        Do one update step for each process.

        """
        for p in self.processes:
            self.processes[p].update(self)

    def simulate(self, steps, log=True):
        """
        Simulate the model for some time.

        """
        for s in xrange(steps):
            self.step()
            if log:
             # This could be an entry point for further logging

                # print count of each protein to the screen
                print '\r{}'.format([len(self.states[x]) for x in self.states.keys() if "mRNA_" in x]),
            #        string aus: Laenge von Array (states)     x laeuft durch states.keys, BDG: Protein_ kommt in states.key-Eintrag vor
            #                                   Anzahl der Proteine

if __name__ == "__main__":

    #print 'sds'

    c = Model()
    c.simulate(30, log=True)


