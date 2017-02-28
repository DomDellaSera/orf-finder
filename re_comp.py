import pdb
import sys
import logging










 	
def getKey(custom): #For sequence class
    return custom.length


class nonspecifc_Sequence(object): # This will be used to create a tuple later indexed by len(u)
    def __init__(self, sequence):
        self.sequence = sequence
        self.length = len(sequence)
        self.ORF = False
        
    def __repr__(self): # This forces the sequence to be printed, abstracting the attributes
        return('{}'.format(self.sequence))
        

#logging.disable() = True
logger = logging.getLogger()
#print("Arguements \n\n", str(sys.argv), "\n\n")
if "-debug"  in str(sys.argv):
    logging.basicConfig(level=logging.DEBUG)
    logging.basicConfig(format='%(asctime)s %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p')
    
def translation(rna_3_letter_code):
    #input:
        #3 letter dna code
    #output: protein
    RNA_codon_table = { #
		'UUU': 'Phe', 'UCU': 'Ser', 'UAU': 'Tyr', 'UGU': 'Cys', 
		'UUC': 'Phe', 'UCC': 'Ser', 'UAC': 'Tyr', 'UGC': 'Cys', 
		'UUA': 'Leu', 'UCA': 'Ser', 'UAA': '---', 'UGA': '---', 
		'UUG': 'Leu', 'UCG': 'Ser', 'UAG': '---', 'UGG': 'Trp', 
		'CUU': 'Leu', 'CCU': 'Pro', 'CAU': 'His', 'CGU': 'Arg', 
		'CUC': 'Leu', 'CCC': 'Pro', 'CAC': 'His', 'CGC': 'Arg', 
		'CUA': 'Leu', 'CCA': 'Pro', 'CAA': 'Gln', 'CGA': 'Arg', 
		'CUG': 'Leu', 'CCG': 'Pro', 'CAG': 'Gln', 'CGG': 'Arg', 
		'AUU': 'Ile', 'ACU': 'Thr', 'AAU': 'Asn', 'AGU': 'Ser', 
		'AUC': 'Ile', 'ACC': 'Thr', 'AAC': 'Asn', 'AGC': 'Ser', 
		'AUA': 'Ile', 'ACA': 'Thr', 'AAA': 'Lys', 'AGA': 'Arg', 
		'AUG': 'Met', 'ACG': 'Thr', 'AAG': 'Lys', 'AGG': 'Arg', 
		'GUU': 'Val', 'GCU': 'Ala', 'GAU': 'Asp', 'GGU': 'Gly', 
		'GUC': 'Val', 'GCC': 'Ala', 'GAC': 'Asp', 'GGC': 'Gly', 
		'GUA': 'Val', 'GCA': 'Ala', 'GAA': 'Glu', 'GGA': 'Gly', 
		'GUG': 'Val', 'GCG': 'Ala', 'GAG': 'Glu', 'GGG': 'Gly'
	}
    singleletter = {
		'Cys': 'C', 'Asp': 'D', 'Ser': 'S', 'Gln': 'Q', 'Lys': 'K', 'Trp': 'W', 'Asn': 'N', 'Pro': 'P',
        'Thr': 'T', 'Phe': 'F', 'Ala': 'A', 'Gly': 'G', 'Ile': 'I', 'Leu': 'L', 'His': 'H', 'Arg': 'R',
        'Met': 'M', 'Val': 'V', 'Glu': 'E', 'Tyr': 'Y', '---': '*', "NNN" : 'X'
	}
    def codon_error_tolerator(tetracodon):
        inferred_amino_acid = {
            "UUN" : "NNN",
            "CUN" : "Leu",
            "GUN" : "Val",
            "AUN" : "NNN",
            "UCN" : "Ser",
            "CCN" : "Pro",
            "ACN" : "Thr",
            "GCN" : "Ala",
            "UAN" : "NNN",
            "CAN" : "NNN",
            "AAN" : "NNN",
            "GAN" : "NNN",
            "UGN" : "NNN",
            "CGN" : "Arg",
            "AGN" : "NNN", 
            "GGN" : "Gly"}
        if tetracodon[2] is "N" and tetracodon.count("N") == 1:
            correction = inferred_amino_acid[tetracodon]
            return(correction)
        elif "N" in tetracodon:
            return("NNN")
        else: return(tetracodon)
	
    
    rna_nuc = ""
    for i in rna_3_letter_code:# We can't work with T
        
        if i is "T": rna_nuc += "U"
        else: rna_nuc += i
    if "N" in rna_nuc:
        aa = codon_error_tolerator(rna_nuc)
    elif "N" not in rna_nuc:
        aa = RNA_codon_table[rna_nuc]
            
    
    return(singleletter[aa])
	

    
def reverse_complement(DNA): #HELPER FUNCTION ONLY - don't call this directly
    #Input:
            #Single DNA sequence
    #Output:
            # Sequence with of biological bonding partners in reverse order


    def complement(n): #Not sure if this shoud be its own function but I put it inside reverse_compliment to 
                       # clarify the flow of function calls
        u = None
        if n == "A": u = "T"
        elif n == "T": u = "A"
        elif n == "G": u = "C"
        elif n == "C": u = "G"
        elif n == "N" : u = n
        #else: Exc
        
        
        return(u)
        
        
    transcribed_rna = ""
    logging.debug("""
    ##############################
    #Reverse Compliment Debugging#
    ##############################
    """)
    #logging.debug("input: "+DNA[0])
    
    
    for i in DNA:
        try:
            transcribed_rna += complement(i)
        except:
            print(i)
    return(transcribed_rna[::-1])#This [::-1] alone is doing the reversal here-
                                    #Useful for strings
logging.debug('output: '+ reverse_complement("T"))
                                   


def reading_frame_generator(sequence):
    # Calls Reverse complement, builds list data structure of length two and in each adds 3 strings, each at a different
    # amino acid starting point as such all combinations of reading frames are produced
    #input:
        #DNA sequence
    #output:
        # 6 DNA sequences, 3 for the sequence and 3 for the reverse complement
        
    sequence_rc = reverse_complement(sequence)
	
    HexRF = []
    for i in [sequence, sequence_rc]:
        #print(i)
    
        HexRF.append(i)
        HexRF.append(i[1:])
        HexRF.append(i[2:])
    return(HexRF)

def reading_frame_translator(dna): #HELPER FUNCTION FOR reading_frame_generator

    #input:
        #linear DNA sequence
    #output:
        #Single Peptide peptide
    #function calls:
        #translation
        
    peptide = ""
    frame = ""
    for i in dna:
        frame += i
        #print(frame)
        if len(frame) == 3:
            peptide +=translation(str(frame))
            frame = ""
        else: continue
    return(peptide)


def peptide_frame_chooser(DNA_FRAMES):
    logging.debug("peptide_frame_chooser function input: ")
    logging.debug(DNA_FRAMES)
    #input: a list of 6 strings, each a different possible string
    #output: the longest open reading frame(s)
    
    maximum_length = len(DNA_FRAMES[0]) #If the whole sequence is translatable, it is probably the ORF
    peptide_frames = []
    max_segment=0
    
    for i in DNA_FRAMES:
        peptide_frames.append(reading_frame_translator(i)) #Translation function called here
        #print(peptide_frames)
        
    open_reading_frames= [] #This is the list that will be restarted as the maximum peptide length grows
    #Depreciated as part of spagetti code# max_orf_length = 0 # This is just to keep track of our largest ORF without having to index
    
    for i in peptide_frames: 
        # Each split at stop codons, and reiterated to find the largest split
        max_segment_length = 0
        possible_proteins = i.split("*")
        
        for j in possible_proteins:
            open_reading_frames.append(j)
    ORFs_Indexed = []
    
    def ORF_Indexer(sequence): #I'm depreciating this for a class instead for now
        indexed_orf = (sequence, len(sequence))
        return(indexed_orf)
        
    orf_tmp = []
    for i in open_reading_frames:
        f = nonspecifc_Sequence(i) #Object Initilized here; converge objs in the future
        orf_tmp.append(f)
   
            
    
    sorted_orf =open_reading_frames.sort
    
    logging.debug("Final open reading frames:")
    logging.debug("Sorted By greatest length")
    logging.debug(sorted(orf_tmp, key=getKey, reverse=True))
    
    logging.debug("""
    #############################
    # Open Reading Frame Finder #
    #############################
    """)
    
    
    
    
    orfs_Len = sorted(orf_tmp, key=getKey, reverse=True)
    max_orf_length = orfs_Len[0].length #This shouldn't be mutable now
    logging.info("Protein Sequence   "+ " "* (max_orf_length-len('Protein Sequence '))+"Length")
    logging.info('Top 5 Open Reading Frames')
    for i in orfs_Len[:6]: #Describing the sequences as ORFs# This is just to create an ORF Table
        spacer = " "*((3+ max_orf_length)-i.length)
        
        logging.info(i.sequence+spacer+ str(i.length))
        if i.length == max_orf_length:
            i.ORF = True

        logging.info("\nThe Max ORF is %i" % max_orf_length)
    logging.info([x for x in orfs_Len if x.ORF is True])

    return([x.sequence for x in orfs_Len if x.ORF is True])
       
dna = "TTATAGCGCATGTTAGCGCATGCATGCTATGCGCGATGTGTATAGTGACTGATATACAATGCCATGCATGCAAAATATAGCGCATGTTAGCGCATGCATGCTATGCGCGATGTGTATAGTGACTGATATACAATGCCATGCATGCAAAAATAGCGCATGTTAGCGCATGCATGCTATGCGCGATGTGTATAGTGACTGATATACAATGCCATGCATGCAAAATATAGCGCATGTTAGCGCATGCATGCTATGCGCGATGTGTATAGTGACTGATATACAATGCCATGCATGCAAAA"  

rf = reading_frame_generator(dna)
returned_obj = peptide_frame_chooser(rf)
logging.debug("Variable stored in assigned returner: ")
logging.debug(returned_obj)

#print(rf)
	
	
	
	
	
	
	
	
	
	
