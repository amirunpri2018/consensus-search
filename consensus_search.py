#!/usr/bin/env python

"""
A very simple consensus sequence search with required letters.
------------------------------------------------------------------------------

Information
===========

This script searches a fasta genome for instances of a consensus sequence while
also requiring a particular 3' sequence.

Under the hood, we use pygr (for loading the genome) and motility (for
performing the search).


Download
========

If you're on our `github page <https://github.com/uci-cbcl/consensus_search/>`_, 
you can click on the "Download Zip" button.  Otherwise, try this link::

   https://github.com/uci-cbcl/consensus_search/archive/master.zip


Usage
=====

To perform the search used in the *genesis* paper, first download the xenopus genome::

   wget ftp://ftp.xenbase.org/pub/Genomics/JGI/Xentr7.1/xenopus_tropicalis_v7.1.tar.gz
   tar xfz xenopus_tropicalis_v7.1.tar.gz
   mv 20100930/sequences/Xenopus_tropicalis.main_genome.scaffolds.fasta .
   rm -r xenopus_tropicalis_v7.1.tar.gz 20100930

Then install `motility <https://github.com/ctb/motility>`_ and
pygr (`easy_install pygr`).

Finally, you can reproduce the results from the paper by running::

    python consensus_search.py --genome Xenopus_tropicalis.main_genome.scaffolds.fasta \
            --consensus GGAACTGGCCCCTGCAAACA --required_3p_seq NGG --mismatches 5 \
            --outfile results.bed

This search will allow up to 5 mismatches to the tyrosinase site above enforcing  a degenerate
NGG PAM sequence at the 3' end. To search for mismatch sites to your own target site of
interest, substitute your sequence in place of the above consensus after --consensus.

Note that any of the IUPAC letters can be used in the 
consensus and required sequences. Specifically,

==========   ==============
IUPAC code   Allowed letter
==========   ==============
A            Adenine
C            Cytosine
G            Guanine
T            Thymine
U            Uracil (converted to T for DNA search)
R            Purine (A or G)
Y            Pyrimidine (C or T)
M            C or A
K            T, or G
W            T, or A
S            C or G
B            C, T, or G (not A)
D            A, T, or G (not C)
H            A, T, or C (not G)
V            A, C, or G (not T)
N            Any base (A, C, G, or T)
==========   ==============



For additional usage instructions, run::

    python consensus_search.py --help

Or shoot me an email at jake.biesinger@gmail.com

"""

__author__ = "Jake Biesinger, jake.biesinger@gmail.com"
__license__ = "Apache"

import argparse
import sys

try:
    import motility
except ImportError:
    print "Couldn't import motility! please download and install it from https://github.com/ctb/motility"
    sys.exit()

try:
    from pygr.seqdb import SequenceFileDB
except ImportError:
    print "Couldn't import pygr! Please install it by typing `sudo easy_install pygr`"
    sys.exit()

DNA_LETTERS = 'ACGT'
IUPAC_LETTERS = dict(A='A', C='C', G='G', T='T', U='T', R='AG', Y='CT', S='GC',
                     W='AT', K='GT', M='AC', B='CGT',
                     D='AGT', H='ACT', V='ACG', N='ACGT')
IUPAC_SCORES = {k: [0 if l in v else -1 for l in DNA_LETTERS]
                for k, v in IUPAC_LETTERS.items()}

# required letters are just PWM entries with -1000000 in them instead of -1
REQUIRED_SCORES = {k: [elem * 1000000 for elem in v]
                   for k, v in IUPAC_SCORES.items()}


def make_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument("--genome", required=True,
                        help="fasta file with genome to search")
    parser.add_argument("--mismatches", required=True, type=int,
                        help="The number of mismatches to allow")
    parser.add_argument("--consensus", required=True,
                        help="The consensus sequence to search for")
    parser.add_argument("--required_3p_seq", required=True,
                        help="the sequence *required* on the 3' end")
    parser.add_argument("--outfile", required=True,
                        help="where to save the results")
    return parser


def main(argv=None):
    parser = make_parser()
    args = parser.parse_args(argv)
    genome = SequenceFileDB(args.genome)


    pwm = [IUPAC_SCORES[l] for l in args.consensus]
    pwm.extend([REQUIRED_SCORES[l] for l in args.required_3p_seq])
    pwm = motility.PWM(pwm)

    # find all matches
    with open(args.outfile, 'w') as outfile:
        for chrom in genome.keys():
            chromseq = str(genome[chrom])
            print "searching ", chrom, "of length", len(chromseq)
            matches = pwm.find(chromseq, -args.mismatches)
            for start, stop, strand, seq in matches:
                score = pwm.calc_score(seq)
                outfile.write('\t'.join(
                    [chrom, str(start), str(stop), seq, str(score),
                     '+' if strand == 1 else '-']) + '\n')


if __name__ == '__main__':
    main()
