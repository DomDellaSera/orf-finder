# orf-finder

Finds Open Reading Frames in a Fasta File of DNA Nucleotides. This is the first program I've made whose complexity could be nicely handled with object-oriented abstractions, and thus is still a work in progress. 

Open Reading Frames are found by assuming UAA/UAG/UGA sequences are stop codons and taking the longest continuous sequence from the possible reading frames from the sequence and its reverse compliment.
