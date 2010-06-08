# Stitch

Austin Glenn Davis-Richardson
<harekrishna@gmail.com>

An experiment:

Illumina sequencing has the ability to generate paired-end reads.  If the sum of the DNA fragments is roughly less than that of 2x an illumina read (~100bp),  they can be said to "overlap."  This script also preserves the quality scores (the 5'-most overrides the 3'-most)

This script uses BLAST to find how a set of paired end reads overlaps, and then writes out a contig based on the best BLAST hit.

## Usage

python test.py fastqfile1 fastqfile2 contigfile contig_failed_file

Where fastqfile1 is the 5' most and fastqfile2 is the 3' most.

(Not to self: Switch to argparse or something!)

## Files

 - dnaobj.py - object representing a FASTA or FASTQ record
 - doubleblast.py - Uses NCBI+ bl2seq to perform a blast between two sequences.
 - fastitr.py - Iterates through a FASTA/Q file, generating dnaobj objects
 - gcblast.py - This was to be the main script
 - test.py - This IS the main script right now.

## Pre-reqs

Uses blast+ 2.2.23 from NCBI
Get it here: <ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/>

## Bugs

 - Quality scores of the left (5') most read override those of the right (3') read.

## Know your rights

Stitch is free and open-source.
See LICENSE file (GNU GPL v3)