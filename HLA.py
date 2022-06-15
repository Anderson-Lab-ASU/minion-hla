from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

class HLA:                                                                      
    def __init__(self, gene, field_1, field_2, field_3, field_4, hla_sequence, imgt_id):
        self.gene = gene                                                        
        self.field_1 = field_1                                                  
        self.field_2 = field_2                                                  
        self.field_3 = field_3                                                  
        self.field_4 = field_4                                                  
        self.hla_sequence = hla_sequence                                         
        self.imgt_id = imgt_id                                                  
                                                                                
    def printHLA(self):                                                             
        return ("HLA-"+self.gene+"*"+self.field_1+":"+self.field_2+":"+self.field_3+":"+self.field_4)

    def is01last(self):
        if(self.field_4 == "01" or self.field_4 == ""):
            return True
        else: 
            return False 
    
    def getGene(self):
        return self.gene
    
    def getAllele(self):
        return self.field_1

    def toSeqRecord(self):
        return SeqRecord(Seq(self.hla_sequence),id=self.imgt_id, description=self.printHLA())

