Automatic Multi-Locus Species Tree (autoMLST) Overview
====================================================
autoMLST is a webserver and analysis pipeline for constructing species phylogenies from genomic data. The pipeline scans for conserved single copy housekeeping genes, selects comparable genomes in the database, and builds a phylogeny using either a concatenated gene matrix or coalescent tree method.

For full instructions see `autoMLST on bitbucket`_.

Simplified wrapper instructions
-------------------------------
The ``simplified_wrapper.py`` script constructs an MLSA-based tree from a directory of Genbank files. Conserved single copy housekeeping genes are identified and a tree built from a concatenated gene matrix. However, no comparable genomes are added, nor are existing genomes removed.
Usage: ``simplified_wrapper.py directory``

NOTE: ``reducedcore.zip`` needs to be unpacked before running. 

.. _https://bitbucket.org/ziemertlab/automlst/src/master/


