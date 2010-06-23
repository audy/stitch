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
    
    if not (options.filea and options.fileb):
        print >> sys.stderr, 'Usage: %s %s' % \
            (parser.get_prog_name(), parser.usage)
        sys.exit()
    if not (options.prefix):
        print >> sys.stderr, 'Warning: no outputfile'
        dudsa, dudsb, outfile = sys.stdout, sys.stdout, sys.stdout
    else:
        dudsa = open('%s-nh-s1.fastq' % options.prefix, 'w')
        dudsb = open('%s-nh-s2.fastq' % options.prefix, 'w')
        outfile = open('%s-contigs.fastq' % options.prefix , 'w')

    seqsa = open(options.filea, 'r')
    seqsb = open(options.fileb, 'r')
    
    numcontigs, numtotes = 0, 0
    
    p = Pool()
     
    for i in p.imap(doStitch, izip(Fasta(seqsa), Fasta(seqsb))):
        numtotes += 1
        if i.record:
            numcontigs += 1
        else:
            print 'loser!'
        
    
    print 'Made %s contigs, out of %s reads' % \
        (numcontigs, numcontigs + numtotes)
        
    
class Stitch:
    ''' Stitches together two overlapping Illumina reads using Doubleblast '''
    def __init__(self, reca, recb):
        self.reca = reca
        self.recb = recb
        self.record = False
        self.contig = []
        self.quality = []
        self.find_overlaps()
        
    def find_overlaps(self):
        ''' Alignment algorithm, returns new DNA object of contig '''
        seqa, seqb = self.reca.seq, self.recb.revcomp
        seqb = self.recb.revcomp
        
        length = min([len(i) for i in (seqa, seqb)])
        scores = []
        
        hits = {}
        
        for i in range(len(seqa)):
            score = 0
            for na, nb in zip(seqa[i-1:], seqb[:-i+1]):
                if na == nb:
                    score += 1
            hits[score] = i
        
        
        score = max(hits.keys())
        i = hits[score]
        
        if ((score > 20)):
            print 'winner! %s' % (i)
            print seqa + '-'*(i-1)
            print '-'*(i-1) + seqb
            self.record = self.reca



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