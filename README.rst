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
            --consensus GCACAAAAGCCGGAGCTTC --required_3p_seq NMGG --mismatches 5 \
            --outfile results.bed

Which will search for the tyrosine site #2 (consensus sequence) with the
N[A/C]GG sequence at the 3' end.

For additional usage instructions, run::

    python consensus_search.py --help

Or shoot me an email at jake.biesinger@gmail.com
