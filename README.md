# Stitch

Austin Glenn Davis-Richardson  
<adavisr@ufl.edu>  
[@audyyy](http://www.twitter.com/audyyy)

Working under Dr. Eric Triplett & Dr. Adriana Giongo

## Description:

Stitch creates contigs from overlapping paired-end Illumina reads using a simple (but fast) algorithm.

Stitch is a __pre-publication__ release meaning (a) It's buggy (b) You should cite us (c) [send a message](mailto:adavisr@ufl.edu) if you need to cite us.

### Score Calculation:

Scores are calculated with the following equation:

`matches - mismatches / (length-of-overlap)`

The default minimum score is 0.6 (60%).

Stitch calculates all possible overlaps and chooses the one with
the highest score.  Any matches/mismatches containing an `N` do not
contribute to the score.
		
## Installation

1. Once in the stitch directory do:
   > `$ python setup.py install`

	(you may have to sudo if you do not have write privileges to `/usr/local/bin`)

2. Stitch should have installed successfully meaning that you can now do this:
   > `$ stitch -h`
   to display help information

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

## Importing

			from stitch import *
			stitch(options={
		  	'filea': 'filenamea',
		  	'fileb': 'filenameb',
		  	'prefix': 'data/stitched.fastq',
		  	'pretty': False,
		  	'threads': None,
		  	'score': None })

## Bugs/Feature requests

 - If you have any, [let me know](https://github.com/audy/stitch/issues). Thx!

## Features

 - Automatically joins overlapping paired ends!
 - **Multithreaded** - Automatically uses all cores on your system!
 - Generates 3 files: `contigs`, and two files for unassembled pairs.

## Plans

 - Accept interleaved fastq format.
 - Support piped IO.
 - Implement with MapReduce.

## License

Stitch is free and open-source.
See LICENSE (GNU GPL v3)