# Stitch

Austin G. Davis-Richardson  
<harekrishna@gmail.com>

Stitch assembles overlapping paired-end reads into a single contig for each pair. This increases the read length and hopefully the quality of a _de novo_ or reference assembly. Stitch is multi-threaded and will automatically use all cores on your system unless told otherwise. Stitch currently only reads FASTQ format. QSEQ and FASTA formats to come. Reads that are not found to overlap are dumped in a file called `<prefix>-singletons` and are in FASTQ format. These can then be trimmed and combined with contigs to do a _de novo_ assembly.

Stitch is licensed under the GNU GPL v3.

Stitch is a __pre-publication__ release meaning (a) It's still buggy (b) Be nice and cite me (c) [send a message](mailto:adavisr@ufl.edu) if you need to find out how.
		
## Requirements

- Python 2.6, 2.7
- Mac OS X or Linux (Windows might work but hasn't been tested)

## Alignment Algorithm

Stitch aligns overlapping paired end reads by counting the number of matching nucleotides in an overlapping window. The window that provides the highest number of matching nucleotides wins.

The consensus sequence is generated thusly,

    A ================= ===================================
    B                   =================================== ============
    C ================= =================================== ============
                        <--------- "the middle" ----------->
                        
    In the region dubbed "the middle", the nucleotide with the highest
    corresponding quality score is used (if there is a mismatch). If there is
    a match, then that nucleotide is used (and the highest quality score is 
    given). In the case of a tie, a 'N' is used and the quality score is 
    unchanged.
    
Stitch assumes that read A is 3'-5' and read B is 5'-3'. So read B is automatically reverse-complemented in the alignment procedure. (I know I should make this an option but I haven't yet).

## Usage

**NOTE** - Stitch expects reads to be of the same length!

Invoke, _comme ca_

    Usage: stitch.py -i <fastq file 1> -j <fastq file 2> -o <output prefix>

More options,

    -h, --help            show this help message and exit
    -i FILEA, --first=FILEA
                        first fastq file
    -j FILEB, --second=FILEB
                        second fastq file
    -o PREFIX, --output=PREFIX
                        output prefix (omit to print to stdout)
    -t THREADS, --threads=THREADS
                        number of threads (default = all available)
    -p, --pretty_output   displays overlapping contigs in a nice way.
    -s SCORE, --score=SCORE
                        minimum percent identity (default = 25)
    -b TABLE, --table=TABLE
                        output overlap length to a text file

## Bugs/Feature requests

 - If you have any, [let me know](https://github.com/audy/stitch/issues). Thx!


