# -*- coding: utf-8 -*-
"""
YOUR HEADER COMMENT HERE

@author: Paige Pfenninger

"""

import random
from amino_acids import aa, codons, aa_table   # you may find these useful
from load import load_seq


def shuffle_string(s):
    """Shuffles the characters in the input string
        NOTE: this is a helper function, you do not
        have to modify this in any way """
    return ''.join(random.sample(s, len(s)))

# YOU WILL START YOUR IMPLEMENTATION FROM HERE DOWN ###


def get_complement(nucleotide):
    """ Returns the complementary nucleotide

        nucleotide: a nucleotide (A, C, G, or T) represented as a string
        returns: the complementary nucleotide
    >>> get_complement('A')
    'T'
    >>> get_complement('C')
    'G'
    >>> get_complement('G')
    'C'
    >>> get_complement('T')
    'A'
    """
    if nucleotide == 'A':
        return 'T'
    elif nucleotide == 'T':
        return 'A'
    elif nucleotide == 'C':
        return 'G'
    else: # nucleotide == 'G'
        return 'C'


def get_reverse_complement(dna):
    """ Computes the reverse complementary sequence of DNA for the specfied DNA
        sequence

        dna: a DNA sequence represented as a string
        returns: the reverse complementary DNA sequence represented as a string
    >>> get_reverse_complement("ATGCCCGCTTT")
    'AAAGCGGGCAT'
    >>> get_reverse_complement("CCGCGTTCA")
    'TGAACGCGG'
    >>> get_reverse_complement("ATTCGAGT")
    'ACTCGAAT'
    """
    dna2 = ''
    for i in range(len(dna)):
        dna2 = dna2 + get_complement(dna[i])

    index = len(dna2) - 1
    reverse_dna = ''
    while index >= 0:
        reverse_dna = reverse_dna + dna2[index]
        index = index - 1
    return reverse_dna


def rest_of_ORF(dna):
    """ Takes a DNA sequence that is assumed to begin with a start
        codon and returns the sequence up to but not including the
        first in frame stop codon.  If there is no in frame stop codon,
        returns the whole string.

        dna: a DNA sequence
        returns: the open reading frame represented as a string
    >>> rest_of_ORF("ATGTGAA")
    'ATG'
    >>> rest_of_ORF("ATGAGATAGG")
    'ATGAGA'
    >>> rest_of_ORF("ATGGGCCCCAAATAA")
    'ATGGGCCCCAAA'
    """
    new_dna = ''

    index = 0

    while index < len(dna) - 2:
        codon = dna[index] + dna[index + 1] + dna[index + 2]

        if codon == 'TAA' or codon == 'TAG' or codon == 'TGA': 
            return new_dna
        else: #codon is not a stop codon
            new_dna = new_dna + codon

        index = index + 3
    
    return dna #only happens if there is no stop codon 


def find_all_ORFs_oneframe(dna):
    """ Finds all non-nested open reading frames in the given DNA
        sequence and returns them as a list.  This function should
        only find ORFs that are in the default frame of the sequence
        (i.e. they start on indices that are multiples of 3).
        By non-nested we mean that if an ORF occurs entirely within
        another ORF, it should not be included in the returned list of ORFs.

        dna: a DNA sequence
        returns: a list of non-nested ORFs
    >>> find_all_ORFs_oneframe("ATGCATGAATGTAGATAGATGTGCCC")
    ['ATGCATGAATGTAGA', 'ATGTGCCC']
    >>> find_all_ORFs_oneframe("ATGCATGAATGTAGATAGAAAAAAATGTGCCC")
    ['ATGCATGAATGTAGA', 'ATGTGCCC']
    """
    index = 0
    new_dna = ''
    dna_list = []
    while index < len(dna) - 2:
        codon = dna[index: index + 3]

        if codon == 'ATG':
            new_dna = dna[index: len(dna)] #gets a string of dna that begins with a start codon and extends to the end of the dna string
            cut_dna = rest_of_ORF(new_dna) #cuts dna at stop codon

            dna_list.append(cut_dna) #adds dna to the list of ORF

            index = index + len(cut_dna) + 3
        else:
            index = index + 3

    return dna_list


def find_all_ORFs(dna):
    """ Finds all non-nested open reading frames in the given DNA sequence in
        all 3 possible frames and returns them as a list.  By non-nested we
        mean that if an ORF occurs entirely within another ORF and they are
        both in the same frame, it should not be included in the returned list
        of ORFs.

        dna: a DNA sequence
        returns: a list of non-nested ORFs

    >>> find_all_ORFs("ATGCATGAATGTAG")
    ['ATGCATGAATGTAG', 'ATGAATGTAG', 'ATG']
    >>> find_all_ORFs("")
    []
    """
    frame1 = dna
    frame2 = dna[1: len(dna)]
    frame3 = dna[2: len(dna)]

    orf = []

    orf.extend(find_all_ORFs_oneframe(frame1)) #adds the ORF in frame one to the complete list of ORF
    orf.extend(find_all_ORFs_oneframe(frame2)) #adds the ORF in frame two to the complete list of ORF
    orf.extend(find_all_ORFs_oneframe(frame3)) #adds the ORF in frame three to the complete list of ORF

    return orf


def find_all_ORFs_both_strands(dna):
    """ Finds all non-nested open reading frames in the given DNA sequence on both
        strands.

        dna: a DNA sequence
        returns: a list of non-nested ORFs
    >>> find_all_ORFs_both_strands("ATGCGAATGTAGCATCAAA")
    ['ATGCGAATG', 'ATGCTACATTCGCAT']
    >>> find_all_ORFs_both_strands("")
    []
    """
    strand1 = dna 
    strand2 = get_reverse_complement(dna)

    orf = []
    orf.extend(find_all_ORFs(strand1)) # adds the ORF in strand one to the complete list of ORF
    orf.extend(find_all_ORFs(strand2)) # adds the ORF in strand two to the complete list of ORF

    return orf


def longest_ORF(dna):
    """ Finds the longest ORF on both strands of the specified DNA and returns it
        as a string
    >>> longest_ORF("ATGCGAATGTAGCATCAAA")
    'ATGCTACATTCGCAT'
    """
    length = 0
    longest = ''
    orf = find_all_ORFs_both_strands(dna)
    for element in orf:
        if len(element) > length:
            longest = element
            length = len(element)

    return longest


def longest_ORF_noncoding(dna, num_trials):
    """ Computes the maximum length of the longest ORF over num_trials shuffles
        of the specfied DNA sequence

        dna: a DNA sequence
        num_trials: the number of random shuffles
        returns: the maximum length longest ORF """
    i = 1
    max_length = 0
    while i <= num_trials:
        new_dna = shuffle_string(dna) #shuffles dna
        l = longest_ORF(new_dna) #l is the longest orf in the shuffeled dna 

        if len(l) > max_length:
            max_length = len(l)

        i = i + 1

    return max_length



def coding_strand_to_AA(dna):
    """ Computes the Protein encoded by a sequence of DNA.  This function
        does not check for start and stop codons (it assumes that the input
        DNA sequence represents an protein coding region).

        dna: a DNA sequence represented as a string
        returns: a string containing the sequence of amino acids encoded by the
                 the input DNA fragment

        >>> coding_strand_to_AA("ATGCGA")
        'MR'
        >>> coding_strand_to_AA("ATGCCCGCTTT")
        'MPA'
    """
    index = 0
    protein = ''
    aF = ['TTT', 'TTC']
    aL = ['TTA', 'TTG', 'CTT', 'CTC', 'CTA', 'CTG']
    aI = ['ATT', 'ATC', 'ATA']
    aM = ['ATG']
    aV = ['GTT', 'GTC', 'GTA', 'GTG']
    aS = ['TCT', 'TCC', 'TCA', 'TCG', 'AGT', 'AGC']
    aP = ['CCT', 'CCC', 'CCA', 'CCG']
    aT = ['ACT', 'ACC', 'ACA', 'ACG']
    aA = ['GCT', 'GCC', 'GCA', 'GCG']
    aY = ['TAT', 'TAC']
    aH = ['CAT', 'CAC']
    aQ = ['CAA', 'CAG']
    aN = ['AAT', 'AAC']
    aK = ['AAA', 'AAG']
    aD = ['GAT', 'GAC']
    aE = ['GAA', 'GAG']
    aC = ['TGT', 'TGC']
    aW = ['TGG']
    aR = ['CGT', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG']
    aG = ['GGT', 'GGC', 'GGA', 'GGG']

    while index < len(dna) - 2:
        codon = dna[index: index + 3]

        if codon in aF:
            protein = protein + 'F' # Phenylalanine
        elif codon in aL:
            protein = protein + 'L' # Leucine
        elif codon in aI:
            protein = protein + 'I' # Isoleucine
        elif codon in aV:
            protein = protein + 'V' # Valine
        elif codon in aM: 
            protein = protein + 'M' # Methionine (start codon)
        elif codon in aC:
            protein = protein + 'C' # Cysteine
        elif codon in aA:
            protein = protein + 'A' # Alanine
        elif codon in aG:
            protein = protein + 'G' # Glycine
        elif codon in aP:
            protein = protein + 'P' # Proline
        elif codon in aT:
            protein = protein + 'T' # Threonine
        elif codon in aS:
            protein = protein + 'S' # Serine
        elif codon in aY:
            protein = protein + 'Y' # Tyrosine
        elif codon in aW:
            protein = protein + 'W' # Tryptophan
        elif codon in aQ:
            protein = protein + 'Q' # Glutamine
        elif codon in aN:
            protein = protein + 'N' # Asparagine
        elif codon in aH:
            protein = protein + 'H' # Histidine
        elif codon in aE:
            protein = protein + 'E' # Glutamic acid
        elif codon in aD:
            protein = protein + 'D' # Aspartic acid
        elif codon in aK:
            protein = protein + 'K' # Lysine
        elif codon in aR:
            protein = protein + 'R' # Arginine 
        else:
            protein = protein + '-'

        index = index + 3

    return protein



def gene_finder(dna):
    """ Returns the amino acid sequences that are likely coded by the specified dna

        dna: a DNA sequence
        returns: a list of all amino acid sequences coded by the sequence dna.
    """
    
    threshold = longest_ORF_noncoding(dna, 1500)
    print threshold

    orf_list = find_all_ORFs_both_strands(dna)
    new_orf_list = []
    protein = []

    for i in range(0, len(orf_list)):
        s = orf_list[i]
        if len(s) >= threshold:
            new_orf_list.append(s)
            protein.append(coding_strand_to_AA(s))

    
    return protein


from load import load_seq
dna = load_seq("./data/X73525.fa")
print gene_finder(dna)

# if __name__ == "__main__":
#     import doctest
#     doctest.run_docstring_examples(longest_ORF, globals(), verbose = True)
#     doctest.run_docstring_examples(coding_strand_to_AA, globals(), verbose = True)
