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
        'Met': 'M', 'Val': 'V', 'Glu': 'E', 'Tyr': 'Y', '---': '*'	
	}
	
    rna_nuc = ""
    for i in rna_3_letter_code:
        
        if i is "T": rna_nuc += "U"
        else: rna_nuc += i
    
    
    
    single_letter_aa = singleletter[\
		RNA_codon_table[rna_nuc]\
		]
    
    return(single_letter_aa)
	
def complement(nucleotide):
    if nucleotide == "A": return("T")
    elif nucleotide == "T": return("A")
    elif nucleotide == "G": return("C")
    elif nucleotide == "C": return("G")
    
def reverse_complement(DNA):
    transcribed_rna = ""
    for i in DNA:
        transcribed_rna += complement(i)
    return(transcribed_rna[::-1])

def reading_frame_generator(sequence):
    # Calls Reverse complement, builds list data structure of length two and in each adds 3 strings, each at a different
    # amino acid starting point as such all combinations of reading frames are produced
    #input:
        #DNA sequence
    #output:
        # 6 DNA sequences, 3 for the sequence and 3 for the reverse compliment
        
    sequence_rc = reverse_complement(sequence)
	
    HexRF = []
    for i in [sequence, sequence_rc]:
        #print(i)
    
        HexRF.append(i)
        HexRF.append(i[1:])
        HexRF.append(i[2:])
    return(HexRF)

def reading_frame_translator(dna):

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
        #print(ORF_Indexer(i))
        f = nonspecifc_Sequence(i)
        orf_tmp.append(f)
        
    
        
        
    spagetti_code= """
        #This is producing too much spagetti code for something that isn't done recursively
        
        for j in enumerate(possible_proteins): #This section finds the largest continuous sequence in each frame
            # I'm using enumerate to keep track of some sort of index
            #segment_Length = len(j[1])
            if len(j) > max_segment_length:
                #print("Our old segment length is", max_segment_length, "and our new is", len(j))
                max_segment_length = len(j)
                max_segment = j
                #print(j, "is our max segment at ", max_segment_length, "residues")
                
                
        #print(max_segment_length, max_orf_length)        
            if max_segment_length > max_orf_length:
                #If we find a larger continuous stretch we want to nullify our old data
                max_orf_length = max_segment_length
                #print("Larger ORF found")
                open_reading_frames = None
                open_reading_frames = [] #Erase old ORFs because we found a longer one here
                
            
            #message(max_segment)
            elif max_segment_length == max_orf_length:
                #print(open_reading_frames, "There exists an ORF of a similar length")
                open_reading_frames.append(max_segment)# There could be two possible ORFs of equal length
                #print(open_reading_frames)
                if max_segment_length > 40:
                    pdb.set_trace()
                """
            
    #print("Final open reading frames")
    #print(open_reading_frames)
    sorted_orf =open_reading_frames.sort
    #logging.basicConfig(format='%(asctime)s %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p')
    #print(open_reading_frames)
    #handler = logging.StreamHandler()
    

    #logging.basicConfig(format='%(asctime)s %(name)-12s %(levelname)-8s %(message)s')
    #logging.debug("1###########################################")

    #handler.setFormatter(formatter)
    #logging.debug("2###########################################")

    #logger.addHandler(handler)
    #logger.setLevel(logging.DEBUG)
    
    
    #logging.debug("3###########################################")
    logging.debug("Final open reading frames:")
    #logging.debug("Sequences and their respective lengths")
    logging.debug([x for x in orf_tmp])
    orfs_List = [x for x in orf_tmp]
    logging.debug("Sorted By greatest length")
    logging.debug(sorted(orfs_List, key=getKey, reverse=True))
    
    logging.debug("""
    #############################
    # Open Reading Frame Finder #
    #############################
    """)
    
    
    
    
    orfs_Len = sorted(orfs_List, key=getKey, reverse=True)
    max_orf_length = orfs_Len[0].length #This shouldn't be mutable now
    logging.debug("Protein Sequence   "+ " "* (max_orf_length-len('Protein Sequence '))+"Length")
    for i in orfs_Len: #Describing the sequences as ORFs# This is just to create an ORF Table
        spacer = " "*((3+ max_orf_length)-i.length)
        logging.debug(i.sequence+spacer+ str(i.length))
        if i.length == max_orf_length:
            i.ORF = True
            #print(i,"\nClass ORF Status: ", i.ORF)
    logging.debug("\n\nThe Max ORF is %i" % max_orf_length)
    logging.debug([x for x in orfs_Len if x.ORF is True])
    #logging.debug([x.length for x in orf_tmp])
    #logger.warning(orf_tmp)
    #logging.debug()
   
    #for i in open_reading_frames:
     #   print(i)
        #print(i[1])
        
       
dna = "TTATAGCGCATGTTAGCGCATGCATGCTATGCGCGATGTGTATAGTGACTGATATACAATGCCATGCATGCAAAATATAGCGCATGTTAGCGCATGCATGCTATGCGCGATGTGTATAGTGACTGATATACAATGCCATGCATGCAAAAATAGCGCATGTTAGCGCATGCATGCTATGCGCGATGTGTATAGTGACTGATATACAATGCCATGCATGCAAAATATAGCGCATGTTAGCGCATGCATGCTATGCGCGATGTGTATAGTGACTGATATACAATGCCATGCATGCAAAA"  

rf = reading_frame_generator(dna)
peptide_frame_chooser(rf)
#print(rf)
	
	
	
	
	
	
	
	
	
	