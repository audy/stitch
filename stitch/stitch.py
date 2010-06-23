from fasta import *
from dna import *
from itertools import izip, imap, dropwhile
from optparse import *
from multiprocessing import Pool
import os
import sys
import time
    
def main():
    ''' Spaghetti & Meatballs '''
    parser = OptionParser(
        description="""Stitch - Tool for creating contigs from overlapping
paried-end illumina reads.""",
        usage='-i <fastq file 1> -j <fastq file 2> -o <output prefix>')
    parser.add_option('-i', '--first', dest='filea')
    parser.add_option('-j', '--second', dest='fileb')
    parser.add_option('-o', '--output', dest='prefix')

    (options, args) = parser.parse_args()
    
    if not (options.filea and options.fileb and options.prefix):
        print >> sys.stderr, 'Usage: %s %s' % \
            (parser.get_prog_name(), parser.usage)
        sys.exit()
    
    seqsa = open(options.filea, 'r')
    seqsb = open(options.fileb, 'r')
    
    numcontigs, numtotes = 0, 0
    
    dudsa = open('%s-nh-s1.fastq' % options.prefix, 'w')
    dudsb = open('%s-nh-s2.fastq' % options.prefix, 'w')
    outfile = open('%s-contigs.fastq' % options.prefix , 'w')
    
    p = Pool()
     
    for i in imap(doStitch, izip(Fasta(seqsa), Fasta(seqsb))):
        if i.hits:
            print >> outfile, i.record
            numcontigs += 1
        else:
            pass
        numtotes += 1
    dudsa.close()
    dudsb.close()
    outfile.close()
    
    print 'Made %s contigs, out of %s reads' % \
        (numcontigs, numcontigs + numtotes)
    
    
class Stitch:
    ''' Stitches together two overlapping Illumina reads using Doubleblast '''
    def __init__(self, reca, recb):
        self.reca = reca
        self.recb = recb
        self.hits = False
        self.contig = []
        self.quality = []
        result = self.find_overlaps()
        
    def find_overlaps(self):
        ''' Alignment algorithm, returns new DNA object of contig '''
        print self.reca.seq, self.recb.seq
        raise StopIteration

def doStitch(recs):
    ''' Used by Pool.imap to create stitch jobs '''
    try:
        reca, recb = recs
        return Stitch(reca, recb)
    except KeyboardInterrupt:
        return        
        
if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        print >> sys.stderr, 'Ouch!'
        quit()