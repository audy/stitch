# Stitch

**No longer maintained**. Please use [Pandsaseq](https://github.com/neufeld/pandaseq).

Austin G. Davis-Richardson
<harekrishna@gmail.com>

Stitch assembles overlapping paired-end reads into a single contig for each
pair. This increases the read length and hopefully the quality of a _de novo_
or reference assembly. Stitch is multi-threaded and will automatically use all
cores on your system unless told otherwise. Stitch currently only reads FASTQ
format. QSEQ and FASTA formats to come. Reads that are not found to overlap are
dumped in a file called `<prefix>-singletons` and are in FASTQ format. These
can then be trimmed and combined with contigs to do a _de novo_ assembly.

Stitch is licensed under the MIT open source license. See `LICENSE.md` for
details.

# Citing Stitch

Please us the following citation:

Brown, Christopher T., et al. "Gut microbiome metagenomics analysis suggests a
functional model for the development of autoimmunity for type 1 diabetes." PloS
one 6.10 (2011): e25792.

Or, in BibTex:

```
@article{brown2011gut,
  title={Gut microbiome metagenomics analysis suggests a functional model for the development of autoimmunity for type 1 diabetes},
  author={Brown, Christopher T and Davis-Richardson, Austin G and Giongo, Adriana and Gano, Kelsey A and Crabb, David B and Mukherjee, Nabanita and Casella, George and Drew, Jennifer C and Ilonen, Jorma and Knip, Mikael and others},
  journal={PloS one},
  volume={6},
  number={10},
  pages={e25792},
  year={2011},
  publisher={Public Library of Science}
}
```

## Requirements

- Python 2.6, 2.7
- Mac OS X or Linux (Windows might work but hasn't been tested)

## Alignment Algorithm

Stitch aligns overlapping paired end reads by counting the number of matching
nucleotides in an overlapping window. The window that provides the highest
number of matching nucleotides wins.

The consensus sequence is generated thusly,

    A 5' =============== =================================> 3'
    B                3' <================================= ===========  5'
    C 5' =============== ================================= ===========> 3'
                        <--------- "the middle" ---------->

In the region dubbed "the middle", the nucleotide with the highest
corresponding quality score is used (if there is a mismatch). If there is a
match, then that nucleotide is used (and the highest quality score is  given).
In the case of a tie, a 'N' is used and the quality score is unchanged.

Stitch assumes that read A is 3'-5' and read B is 5'-3'. So read B is
automatically reverse-complemented in the alignment procedure. (I know I should
make this an option but I haven't yet).

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
