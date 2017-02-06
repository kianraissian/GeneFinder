# -*- coding: utf-8 -*-
"""
YOUR HEADER COMMENT HERE

@author: YOUR NAME HERE

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
    """
    if nucleotide=='A':
        return 'T'
    if nucleotide=='T':
        return 'A'
    if nucleotide=='G':
        return 'C'
    if nucleotide=='C':
        return 'G'



def get_reverse_complement(dna):
    """ Computes the reverse complementary sequence of DNA for the specfied DNA
        sequence

        dna: a DNA sequence represented as a string
        returns: the reverse complementary DNA sequence represented as a string
    >>> get_reverse_complement("ATGCCCGCTTT")
    'AAAGCGGGCAT'
    >>> get_reverse_complement("CCGCGTTCA")
    'TGAACGCGG'
    """

    dna=dna[::-1]
    reversecomplement=''
    for i in dna:
        complement=get_complement(i)
        reversecomplement+=complement
    return reversecomplement


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
    """

    for i in range(0, len(dna), 3):
        b=dna[i:i+3]
        if b=='TAG':
            return dna[0:i]
        if b=='TAA':
            return dna[0:i]
        if b=='TGA':
            return dna[0:i]
    return dna



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
    """

    ORFs=[]
    i=0
    while i<len(dna):
        b=dna[i:i+3]
        if b=='ATG':
            new_dna=rest_of_ORF(dna[i:])
            ORFs.append(new_dna)
            i+=len(new_dna)
        else:
            i+=3
    return ORFs


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
    """
    ORFlist=[]
    for i in [0, 1, 2]:
        ORFlist += find_all_ORFs_oneframe(dna[i:])

    return ORFlist

def find_all_ORFs_both_strands(dna):
    """ Finds all non-nested open reading frames in the given DNA sequence on both
        strands.

        dna: a DNA sequence
        returns: a list of non-nested ORFs
    >>> find_all_ORFs_both_strands("ATGCGAATGTAGCATCAAA")
    ['ATGCGAATG', 'ATGCTACATTCGCAT']
    """
    ORFlist=[]
    ORF=get_reverse_complement(dna)
    newORF=find_all_ORFs(ORF)
    ORFlist.append(newORF)
    otherORF=find_all_ORFs(dna)
    ORFlist.append(otherORF)
    ORFlist.reverse()
    return ORFlist


def longest_ORF(dna):
    """ Finds the longest ORF on both strands of the specified DNA and returns it
        as a string
    >>> longest_ORF("ATGCGAATGTAGCATCAAA")
    'ATGCTACATTCGCAT'
    """

    ORFlist=find_all_ORFs_both_strands(dna)
    ORFlongest=sorted(ORFlist, key=len)

    return ORFlongest[-1]



def longest_ORF_noncoding(dna, num_trials):
    """ Computes the maximum length of the longest ORF over num_trials shuffles
        of the specfied DNA sequence

        dna: a DNA sequence
        num_trials: the number of random shuffles
        returns: the maximum length longest ORF """
    while i<num_trials:
        list(dna)
        list.shuffle(dna)
        newdna=join(dna)
        i+=1
    return newdna




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
    aa=''
    while t<len(dna):
        for t in range(0, len(dna), 3):
            d=dna[i:i+3]
            aa_table[d]
        aa+=d
    return aa



def gene_finder(dna):
    """ Returns the amino acid sequences that are likely coded by the specified dna

        dna: a DNA sequence
        returns: a list of all amino acid sequences coded by the sequence dna.
    """
    scrambled=longest_ORF_noncoding(dna, 1500)
    scrambledbothstrands=find_all_ORFs_both_strands(scrambled)
    aminoacids=coding_strand_to_AA(scrambledbothstrands)
    return aminoacids



if __name__ == "__main__":
    import doctest
    doctest.testmod()
