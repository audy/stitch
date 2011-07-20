# Stitch

Austin Glenn Davis-Richardson  
<heyaudy@gmail.com>

## Description:

Stitch assembles overlapping paired-end reads into a single contig for each pair. This increases the read length and hopefully the quality of a _de novo_ or reference assembly. Stitch is multi-threaded and will automatically use all cores on your system unless told otherwise. Stitch currently only reads FASTQ format. QSEQ and FASTA formats to come. Reads that are not found to overlap are dumped in a file called `<prefix>-singletons` and are in FASTQ format. These can then be trimmed and combined with contigs to do a _de novo_ assembly.

## License

Stitch is licensed under the GNU GPL v3.

Stitch is a __pre-publication__ release meaning (a) It's still buggy (b) Be nice and cite me (c) [send a message](mailto:adavisr@ufl.edu) if you need to find out how.

### Score Calculation:

Scores are calculated with the following equation:

`matches - mismatches / (length-of-overlap)`

The default minimum score is 0.6 (60%).

Stitch calculates all possible overlaps and chooses the one with
the highest score.  Any matches/mismatches containing an `N` do not
contribute to the score.
		
## Requirements

- Python 2.6, 2.7
- Mac OSX or Linux (Windows might work but hasn't been tested)

## Usage

Invoke _comme ca_:

    $ stitch -i <fastq file 1> -j <fastq file 2> -o <output prefix>

Where `<fastqfile1>` is the 5' most and `<fastqfile2>` is the 3' most.

`<outputprefix>` will result in three files: contigs, and two files for rejects

Other options:

 - `-h` prints help
 - `-t [THREADS]` specify number of cores to use (Default=all)
 - `-s [SCORE]` minimum score (Default=0.6)
 - `-p pretty output` prints a user-friendly output.  Useful for adjusting
   score.

## Bugs/Feature requests

 - If you have any, [let me know](https://github.com/audy/stitch/issues). Thx!


