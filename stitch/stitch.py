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
    score = kwargs.get('score', 35)
    pretty = kwargs.get('pretty', None)
    threads = kwargs.get('threads', None)
    table = kwargs.get('table', None)

    score = 20

    
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
    #p = Pool(threads)
    
    # Go!
    for i in imap(doStitch, izip(Fasta(seqsa), Fasta(seqsb))):
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
        a, b = self.reca.seq, self.recb.revcomp
        
        # convert quality score to integers
        qa, qb = [ [ ord(i) for i in j ] for j in [ self.reca.qual, self.recb.qual[::-1] ] ]
        
        scores = {
            'eq': +1,
            'un': -1,
            'N': 0,
        }
        
        alignments = {}
        
        for n in range(len(a)):
            
            # get overlapping region
            ta = a[n:]
            
            # because 'string'[:-0] == ''
            if n > 0:
                tb = b[:-n]
            else:
                tb = b
            
            # score overlap
            score = 0
            for i, j in zip(ta, tb):
                if i == j:
                    score += scores['eq']
                elif i != j:
                    score += scores['un']
                if 'N' in [i, j]:
                    score += scores['N']
        
            alignments[score] = n
        
        best_score = max(alignments.keys())
        best_index = alignments[best_score]
        
        # GENERATE CONTIG
        
        # beginning
        if best_index == 0:
            beginning = a
            qual_beg  = qa
        else:
            beginning = a[:best_index]
            qual_beg  = qa[:best_index]
        
        # middle
        middle = []
        qual_middle = []
        for (i, qi), (j, qj) in zip(zip(a[best_index:], qa[best_index:]), \
                    zip(b[:-best_index], qb[:-best_index])):
            if i == j:
                middle.append(i)
                qual_middle.append(max([qi, qj]))
            elif i != j:
                # take best quality
                if qi > qj: # i wins
                    middle.append(i)
                    qual_middle.append(qi)
                elif qi < qj: # j wins
                    middle.append(j)
                    qual_middle.append(qj)
                elif qi == qj: # tie
                    middle.append('N')
                    qual_middle.append(qi)
                else:
                    raise Exception
            else:
                raise Exception
                
        middle = ''.join(middle)
        qual_middle = qual_middle
        
        assert len(middle) == len(qual_middle)
        
        # end
        if best_index == 0:
            end = ''
            qual_end = []
        else:
            end = b[best_index:]
            qual_end = qb[best_index:]
        
        # concatenate
        newseq  = beginning + middle + end
        newqual = ''.join(chr(i) for i in qual_beg + qual_middle + qual_end)
        
        # double-check
        assert len(newseq) == len(newqual)

        # print >> sys.stderr, " b: %s \n m: %s \n e: %s\n\n"  % ( beginning, middle, end )
        # print >> sys.stderr, newseq 

        # generate pretty print view
        self.pretty = '1:%s\n2:%s\nC:%s\n' % \
                        (a + '-'*(best_index-1),
                        '-'*(best_index-1) + b, newseq)
        
        # create a new record sequence for the contig
        self.record = Dna(self.reca.header, newseq, newqual)
        self.score = best_score

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
    
    parser.add_option('-s', '--score', dest='score', default=35,
        help='minimum percent identity (default = 25)', type=float)
    
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
                score   = options.score)
    
    except KeyboardInterrupt:
        print >> sys.stderr, 'Ouch!'
        quit()

if __name__ == '__main__':
    main()
