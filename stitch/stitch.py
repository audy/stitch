#!/usr/bin/env python

# STITCH
#  Constructs contigs from overlapping paired-end Illumina sequencing reads.

# Austin G. Davis-Richardson
#  adavisr@ufl.edu
#  working under Drs. Adriana Giongo & Eric Triplett
#  at The University of Florida, Gainesville, FL, USA

# https://www.github.com/audy/stitch

from fasta import Dna, Fasta
from itertools import izip, imap, dropwhile
from multiprocessing import Pool
import os
import sys
from time import time


def stitch(*args, **kwargs):
    ''' The stitcher '''
    
    filea = kwargs.get('filea')
    fileb = kwargs.get('fileb')
    prefix = kwargs.get('prefix', None)
    score = kwargs.get('score', 0.6)
    pretty = kwargs.get('pretty', None)
    threads = kwargs.get('threads', None)
    table = kwargs.get('table', None)

    
    if not (filea or fileb):
        raise Exception, 'stitch(filea=\'filea\', fileb=\'fileb\')'
    
    if prefix:
        dudsa = open('%s-nh-s1.fastq' % prefix, 'w')
        dudsb = open('%s-nh-s2.fastq' % prefix, 'w')
        outfile = open('%s-contigs.fastq' % prefix , 'w')
    
    if table:
        htable = open(table, 'w')
    
    seqsa = open(filea, 'r')
    seqsb = open(fileb, 'r')
    
    # Ready.. Set..
    numcontigs, numtotes, overlaps = 0, 0, 0
    starttime = time()
    p = Pool(threads)
    
    # Go!
    for i in p.imap(doStitch, izip(Fasta(seqsa), Fasta(seqsb))):
        numtotes += 1
        if i.score > score:
            numcontigs += 1
            overlaps += i.overlap
            
            if prefix:
                print >> outfile, '%s' % i.record
            
            if pretty:
                print >> sys.stdout, '>%s (%.3f)' % (i.reca.header, i.score)
                print >> sys.stdout, i.pretty
            
            if table:
                print >> htable, i.overlap
        else:
            reca, recb = i.originals
            
            if prefix:
                print >> dudsa, reca
                print >> dudsb, recb
    
    # Clean-up & inform the user
    duration = time() - starttime
    print >> sys.stderr, \
        'Made %s contigs out of %s reads in %.2f seconds (%.2f per sec)' % \
        (numcontigs, numtotes, duration, numtotes/duration)
    try:
        print >> sys.stderr, \
            'Average overlap was %.2f\n' % (float(overlaps)/numcontigs)
    except ZeroDivisionError:
        print >> sys.stderr, 'no contigs :(\n'
    

class Stitch:
    ''' Stitches together two overlapping Illumina reads using Doubleblast '''
    def __init__(self, reca, recb):
        self.reca = reca
        self.recb = recb
        self.record = False
        self.overlap = 0
        self.pretty = ''
        self.score = 0.0
        self.find_overlaps()
    
    @property
    def originals(self):
        return (self.reca, self.recb)
    
    def find_overlaps(self):
        ''' Alignment algorithm, returns new DNA object of contig '''
        
        # reverse complement second sequence
        # this should be made into an option
        seqa, seqb = self.reca.seq, self.recb.revcomp
        seqb = self.recb.revcomp
        
        if len(seqa) != len(seqb):
            raise Exception, 'Stitch requires reads to be of the same length!'
        
        # get the minimum length for score generation
        length = min([len(i) for i in (seqa, seqb)])
        
        scores, hits = [], {}
        
        # Find overlap
        for i in range(len(seqa)):
            score = 0
            for na, nb in zip(seqa[i-1:], seqb[:-i+1]):
                if ('N' in (na, nb)):
                    continue
                if (na == nb):
                    score += 1
                if (na != nb):
                    score -= 1
            hits[score] = i
        
        # BUG: this scoring scheme sometimes favors small alignments over larger ones
        # with slightly worse scores. We prefer the larger ones. How to score differently
        # so that length is taken into account? Normalize by overlap lenght?
        
        score = max(hits.keys())
        i = hits[score] # i stands for index of best overlap
        
        self.overlap = len(seqa) - i
        self.score = float(score)/(len(seqa)-i)
        
        # fuggered code to generate the contig
        beg = seqa[0:i-1]
        end = seqb[-i+1:]
        qbeg = self.reca.qual[0:i-1]
        qend = self.recb.qual[::-1][-i+1:]
        smida = seqa[i-1:]
        smidb = seqb[0:len(seqb)-i+1]
        qmida = self.reca.qual[i-1:]
        qmidb = self.recb.qual[::-1][0:len(seqb)-i+1]
        mid, midq = [], []
        
        for (na, qa), (nb, qb) in \
                        zip(zip(smida, qmida), zip(smidb, qmidb)):
            if qa > qb:
                mid += na
                midq += qa
            elif qa < qb:
                mid += nb
                midq += qb
            else:
                if qa == qb:
                    mid += na
                    midq += qa
                else:
                    mid += 'N'
                    midq += qa
        
        newseq = beg + ''.join(mid) + end
        newqual = qbeg + ''.join(midq) + qend
        
        # generate pretty print view
        self.pretty = '1:%s\n2:%s\nC:%s\n' % \
                        (seqa + '-'*(i-1),
                        '-'*(i-1) + seqb, newseq)
        
        # create a new record sequence for the contig
        self.record = Dna(self.reca.header, newseq, newqual)

def get_args():
    from optparse import OptionParser
    ''' Parses command-line arguments, returns parser object '''
    
    parser = OptionParser(
        description="""Stitch - Tool for creating contigs from overlapping
        paried-end illumina reads.""",
        usage='-i <fastq file 1> -j <fastq file 2> -o <output prefix>')
    
    # TODO: add option to specify orientation of reads.
    parser.add_option('-i', '--first', dest='filea',
        help='first fastq file')
    
    parser.add_option('-j', '--second', dest='fileb',
        help='second fastq file')
    
    parser.add_option('-o', '--output', dest='prefix',
        help='output prefix (omit to print to stdout)')
    
    parser.add_option('-t', '--threads', dest='threads', default=None,
        type=int, help='number of threads (default = all available)')
    
    parser.add_option('-p', '--pretty_output', dest='pretty', default=False,
        action='store_true',
        help='displays overlapping contigs in a nice way.')
    
    parser.add_option('-s', '--score', dest='score', default=0.6,
        help='minimum percent identity (default = 25)', type=float)
    
    parser.add_option('-b', '--table', dest='table', default=None,
        help='output overlap length to a text file')
    
    return parser

def doStitch(recs):
    ''' Used by Pool.imap to create stitch jobs '''
    
    try:
        reca, recb = recs
        return Stitch(reca, recb)
    
    except KeyboardInterrupt:
        # BUG Pool() has a habit of not exiting when a CTRL-C is passed
        # So far, I haven't found a way around this except by killing
        # the process by hand
        # Note: os.abort() may work
        print 'Ouch!'
        quit()

def main():
    ''' run from command-line '''
    
    # parse arguments
    parser = get_args()
    (options, args) = parser.parse_args()
    
    # verify two files specified
    if not (options.filea and options.fileb):
        print >> sys.stderr, 'Usage: %s %s' % \
            (parser.get_prog_name(), parser.usage)
        sys.exit()
    
    # output everything to stdout if no output prefix is specified
    if not (options.prefix):
        print >> sys.stderr, 'Warning: no outputfile'
        dudsa, dudsb, outfile = sys.stdout, sys.stdout, sys.stdout
    
    # initiate the stiching!
    try:
        stitch( filea   = options.filea,
                fileb   = options.fileb,
                prefix  = options.prefix,
                threads = options.threads,
                pretty  = options.pretty,
                score   = options.score,
                table   = options.table)
    
    except KeyboardInterrupt:
        print >> sys.stderr, 'Ouch!'
        quit()

if __name__ == '__main__':
    main()