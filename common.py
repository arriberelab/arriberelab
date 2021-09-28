#!/usr/bin/env python3
"""
Joshua Arribere Aug 10, 2013
Converted to python 3: Mar 24, 2020

Several commonly used scripts
"""
from pyx import color
from pickle import load, dump


def reverseTranslate(seq):
    """Given AA sequence, will reverse translate into optimal codons
    for worm"""
    AAs = {'F': 'TTT', 'L': 'CTT', 'I': 'ATT', 'V': 'GTT', 'C': 'TGT',
           'M': 'ATG', 'A': 'GCT', 'G': 'GGA', 'T': 'ACA', 'W': 'TGG',
           'S': 'TCA', 'Y': 'TAT', 'P': 'CCA', 'H': 'CAT', 'E': 'GAA',
           'Q': 'CAA', 'D': 'GAT', 'N': 'AAT', 'K': 'AAA', 'R': 'AGA',
           '*': 'TAA'}
    aa = ''
    seq = seq.upper()
    for letter in seq:
        aa += AAs[letter]
    bb = ''
    for ii in range(0, len(aa), 3):
        bb += aa[ii:ii + 3] + ' '
    return bb.strip()


def scrambler(seq):
    """Will shuffle a sequence until a stop-codon free version
    is generated"""
    import random
    seq = ''.join(random.sample(seq, len(seq))).upper()
    for ii in range(10000):
        seq = ''.join(random.sample(seq, len(seq)))
        if not translate(seq)[1].endswith('*'):
            return seq


def reCode(seq):
    """Will take a sequence and translate it to AA space, then
    reverse translate it back to codon space, changing all
    but the Met/Trp codons"""
    import random, copy
    seq = seq.upper()
    geneticCode = {'GGG': 'G', 'GGA': 'G', 'GGC': 'G', 'GGT': 'G',
                   'GAG': 'E', 'GAA': 'E', 'GAC': 'D', 'GAT': 'D',
                   'GCG': 'A', 'GCA': 'A', 'GCC': 'A', 'GCT': 'A',
                   'GTG': 'V', 'GTA': 'V', 'GTC': 'V', 'GTT': 'V',
                   'AGG': 'R', 'AGA': 'R', 'AGC': 'S', 'AGT': 'S',
                   'AAG': 'K', 'AAA': 'K', 'AAC': 'N', 'AAT': 'N',
                   'ACG': 'T', 'ACA': 'T', 'ACC': 'T', 'ACT': 'T',
                   'ATG': 'M', 'ATA': 'I', 'ATC': 'I', 'ATT': 'I',
                   'CGG': 'R', 'CGA': 'R', 'CGC': 'R', 'CGT': 'R',
                   'CAG': 'Q', 'CAA': 'Q', 'CAC': 'H', 'CAT': 'H',
                   'CCG': 'P', 'CCA': 'P', 'CCC': 'P', 'CCT': 'P',
                   'CTG': 'L', 'CTA': 'L', 'CTC': 'L', 'CTT': 'L',
                   'TGG': 'W', 'TGA': '*', 'TGC': 'C', 'TGT': 'C',
                   'TAG': '*', 'TAA': '*', 'TAC': 'Y', 'TAT': 'Y',
                   'TCG': 'S', 'TCA': 'S', 'TCC': 'S', 'TCT': 'S',
                   'TTG': 'L', 'TTA': 'L', 'TTC': 'F', 'TTT': 'F'}
    reverseGeneticCode = {}
    for codon in geneticCode:
        AA = geneticCode[codon]
        if AA not in reverseGeneticCode:
            reverseGeneticCode[AA] = []
        reverseGeneticCode[AA].append(codon)
    # now reverse translate
    revSeq = ''
    for ii in range(len(seq)):
        if ii % 3 == 0:
            codon = seq[ii:ii + 3]
            AA = geneticCode[codon]
            if AA not in ['M', 'W']:
                temp = copy.copy(reverseGeneticCode[AA])
                # temp.remove(codon)
                revSeq += random.choice(temp)
            else:
                revSeq += reverseGeneticCode[AA][0]
    return revSeq


def translate(seq):
    geneticCode = {'GGG': 'G', 'GGA': 'G', 'GGC': 'G', 'GGT': 'G',
                   'GAG': 'E', 'GAA': 'E', 'GAC': 'D', 'GAT': 'D',
                   'GCG': 'A', 'GCA': 'A', 'GCC': 'A', 'GCT': 'A',
                   'GTG': 'V', 'GTA': 'V', 'GTC': 'V', 'GTT': 'V',
                   'AGG': 'R', 'AGA': 'R', 'AGC': 'S', 'AGT': 'S',
                   'AAG': 'K', 'AAA': 'K', 'AAC': 'N', 'AAT': 'N',
                   'ACG': 'T', 'ACA': 'T', 'ACC': 'T', 'ACT': 'T',
                   'ATG': 'M', 'ATA': 'I', 'ATC': 'I', 'ATT': 'I',
                   'CGG': 'R', 'CGA': 'R', 'CGC': 'R', 'CGT': 'R',
                   'CAG': 'Q', 'CAA': 'Q', 'CAC': 'H', 'CAT': 'H',
                   'CCG': 'P', 'CCA': 'P', 'CCC': 'P', 'CCT': 'P',
                   'CTG': 'L', 'CTA': 'L', 'CTC': 'L', 'CTT': 'L',
                   'TGG': 'W', 'TGA': '*', 'TGC': 'C', 'TGT': 'C',
                   'TAG': '*', 'TAA': '*', 'TAC': 'Y', 'TAT': 'Y',
                   'TCG': 'S', 'TCA': 'S', 'TCC': 'S', 'TCT': 'S',
                   'TTG': 'L', 'TTA': 'L', 'TTC': 'F', 'TTT': 'F'}
    pep = ''
    for ii in range(0, len(seq) - 2, 3):
        codon = geneticCode[seq[ii:ii + 3]]
        if codon != '*':
            pep += geneticCode[seq[ii:ii + 3]]
        else:
            pep += '*'
            return seq[:ii + 3], pep
    return seq, pep


def writeFasta(dict1, outFile):
    """Given a dict of {name:seq} and an outFile, will write a fasta-formatted file"""
    with open(outFile, 'w') as f:
        for (name, seq) in dict1.iteritems():
            f.write(f'>{name}\n{seq}\n')


def digest(list1, length):
    list1.sort()
    bb = []
    for ii in range(len(list1) - 1):
        bb.append(list1[ii + 1] - list1[ii])
    bb.append(length - list1[-1] + list1[0])
    bb.sort()
    return bb


def partialDigest(list1, length):
    list1.sort()
    bb = []
    for ii in range(len(list1) - 2):
        bb.append(list1[ii + 2] - list1[ii])
    bb.append(length - list1[-2] + list1[0])
    bb.append(length - list1[-1] + list1[1])
    bb.sort()
    return bb


def unPickle(file):
    with open(file, 'r') as f:
        return load(f)


def rePickle(obJect, file):
    with open(file, 'wb') as f:
        dump(obJect, f, protocol=2)


def parseGeneList(file):
    """Will parse a line-delimited list of entries to a dict"""
    with open(file, 'r') as f:
        return dict((line.strip().split()[0], f) for line in f)


def colors(i):
    a = [color.cmyk(1, 0.5, 0, 0),  # blue
         color.cmyk(0, 0.8, 1, 0),  # vermillion
         color.cmyk(0.97, 0, 0.75, 0),  # blue green
         color.cmyk(0.1, 0.7, 0, 0),  # reddish purple
         color.cmyk(0.1, 0.05, 0.9, 0),  # yellow
         color.cmyk(0.8, 0, 0, 0),  # sky blue
         color.cmyk(0, 0.5, 1, 0),  # orange
         color.cmyk(0, 0, 0, 1),  # black
         color.cmyk(1, 0.5, 0, 0),  # blue
         color.cmyk(0.1, 0.7, 0, 0),  # reddish purple
         color.cmyk(0, 0.8, 1, 0),  # vermillion
         color.cmyk(1, 0.5, 0, 0),  # blue
         color.cmyk(0.97, 0, 0.75, 0),  # blue green
         color.cmyk(0.1, 0.7, 0, 0),  # reddish purple
         color.cmyk(0, 0, 0, 1),  # black
         color.cmyk(0, 0.8, 1, 0),  # vermillion
         color.cmyk(1, 0.5, 0, 0),  # blue
         color.cmyk(0.97, 0, 0.75, 0),  # blue green
         color.cmyk(0.1, 0.7, 0, 0),  # reddish purple
         color.cmyk(0, 0.8, 1, 0),  # vermillion
         color.cmyk(1, 0.5, 0, 0),  # blue
         color.cmyk(0, 0, 0, 1),  # black
         color.cmyk(0.1, 0.05, 0.9, 0),  # yellow
         color.cmyk(0, 0.5, 1, 0),  # orange
         color.cmyk(0, 0.8, 1, 0),  # vermillion
         color.cmyk(0.1, 0.7, 0, 0),  # reddish purple
         color.cmyk(1, 0.5, 0, 0),  # blue
         color.cmyk(0.97, 0, 0.75, 0),  # blue green
         color.cmyk(0.8, 0, 0, 0),  # sky blue
         color.cmyk(0.5, 1, 0, 0),  # purple
         color.cmyk(0, 0, 0, 1),  # black
         color.cmyk(0.97, 0, 0.75, 0),  # blue green$
         color.cmyk(0, 0.8, 1, 0),  # vermillion
         color.cmyk(1, 0.5, 0, 0),  # blue
         color.cmyk(0, 0, 0, 1),  # black
         color.cmyk(0.97, 0, 0.75, 0),  # blue green$
         color.cmyk(0.97, 0, 0.75, 0),  # blue green$
         color.cmyk(0, 0.5, 1, 0),  # orange
         color.cmyk(0, 0.5, 1, 0),  # orange
         color.cmyk(0, 0.5, 1, 0),  # orange
         color.cmyk(0.1, 0.7, 0, 0),  # reddish purple
         color.cmyk(0.1, 0.7, 0, 0),  # reddish purple
         color.cmyk(0.1, 0.7, 0, 0),  # reddish purple
         color.cmyk(1, 0.5, 0, 0),  # blue
         color.cmyk(1, 0.5, 0, 0),  # blue
         color.cmyk(1, 0.5, 0, 0),  # blue
         color.cmyk(0.8, 0, 0, 0),  # sky blue
         color.cmyk(0, 0.8, 1, 0),  # vermillion
         color.cmyk(0.1, 0.05, 0.9, 0),  # yellow
         color.cmyk(0.5, 1, 0, 0),  # purple
         color.cmyk(0, 0, 0, 1),  # black
         ]
    """
    a=[    
            color.cmyk(1,0.5,0,0),#blue
            color.cmyk(0.8,0,0,0),
            color.cmyk(0.17,0,0.29,0.51),#green
            color.cmyk(0.17,0,0.33,0.34),
            color.cmyk(0,0.36,0.71,0),#poppy
            color.cmyk(0,0.16,0.71,0.02),
            color.cmyk(0,0.8,1,0),#vermillion
            color.cmyk(0,0.8,0.8,0.1),
            color.cmyk(0.1,0.7,0,0),
            color.cmyk(0.1,0.05,0.9,0),
            color.cmyk(0.97,0,0.75,0),
            color.cmyk(0,0,0,1)
            ]
    """
    return a[i % len(a)]


def mkTuples(dict1):
    """Given a dict1={x,y}, will return [(x,y)] s.t. the x are increasing"""
    a = dict1.keys()
    a.sort()
    return [(key, dict1[key]) for key in a]


def parseFasta(file):
    """Will parse a fasta file to a dict of {name:seq}"""
    with open(file, 'r') as f:
        a = {}
        for line in f:
            line = line.strip()
            if line[0] == '>':
                name = line[1:]  # .strip().split()[0]
                a[name] = ''
            else:
                a[name] = ''.join([a[name], line])
    return a


def revCompl(seq):
    """Will return the reverse complement of a sequence"""
    a = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G',
         'a': 't', 't': 'a', 'g': 'c', 'c': 'g',
         'N': 'N', 'n': 'n', '.': '.', '-': '-'}
    return ''.join([a[seq[-i]] for i in range(1, len(seq) + 1)])


def flipStrand(strand):
    """Will return the opposite strand"""
    if strand == '+':
        return '-'
    elif strand == '-':
        return '+'
    else:
        return 'na'
