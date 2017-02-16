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
    open_reading_frames= []
    max_orf_length = 0
    
    for i in peptide_frames: 
        max_segment_length = 0
        possible_proteins = i.split("*")
        
        for j in enumerate(possible_proteins):
            segment_Length = len(j[1])
            if segment_Length > max_segment_length:
                max_segment_length = segment_Length
                max_segment = j
                #print(j)
                
        #print(max_segment_length, max_orf_length)        
        if max_segment_length > max_orf_length:
            max_orf_length = max_segment_length
            open_reading_frames = None
            open_reading_frames = [] #Erase old ORFs because we found a longer one here
            #message(max_segment)
        elif max_segment_length == max_orf_length:
            #print(open_reading_frames)
            open_reading_frames.append(max_segment)# There could be two possible ORFs of equal length
            
    try:
        
        orf = [x[1] for x in open_reading_frames]
    except:
        orf = "fuck"
    
    return(orf)
	
	
	
	
	
	
	
	
	
	
	
	