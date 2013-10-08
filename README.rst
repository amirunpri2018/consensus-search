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

Which will search for the tyrosine site #2 (consensus sequence) with the
N[A/C]GG sequence at the 3' end.  

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
Y            Pyrimidine (C, T, or U)
M            C or A
K            T, U, or G
W            T, U, or A
S            C or G
B            C, T, U, or G (not A)
D            A, T, U, or G (not C)
H            A, T, U, or C (not G)
V            A, C, or G (not T, not U)
N            Any base (A, C, G, T, or U)
==========   ==============


For additional usage instructions, run::

    python consensus_search.py --help

Or shoot me an email at jake.biesinger@gmail.com
