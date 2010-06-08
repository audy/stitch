# Stitch

Austin Glenn Davis-Richardson
<harekrishna@gmail.com>

An experiment:

Illumina sequencing has the ability to generate paired-end reads.  If the sum of the DNA fragments is roughly less than that of 2x an illumina read (~100bp),  they can be said to "overlap."  This script also preserves the quality scores (the 5'-most overrides the 3'-most)

This script uses BLAST to find how a set of paired end reads overlaps, and then writes out a contig based on the best BLAST hit.

## Installation

1. Make sure you have NCBI's BLAST version 2.2.23 or greater.
   Get it here: <ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/>

2. Once in the stitch directory do:
   > $ sudo python setup.py install

3. Stitch should have installed successfully meaning that you can now do this:
   > $ stitch -h
   to display help information

## Usage

Do this:

> $ python test.py fastqfile1 fastqfile2 contigfile contig_failed_file

Where fastqfile1 is the 5' most and fastqfile2 is the 3' most.

(Not to self: Switch to argparse or something!)

## Bugs

 - Quality scores of the left (5') most read override those of the right (3') read.

## Know your rights

Stitch is free and open-source.
See LICENSE file (GNU GPL v3)