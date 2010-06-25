from fasta import *
from dna import *
from itertools import izip, imap, dropwhile
from optparse import *
from multiprocessing import Pool
import os
import sys
from time import time

    
def main():
    ''' Spaghetti & Meatballs '''
    parser = OptionParser(
        description="""Stitch - Tool for creating contigs from overlapping
paried-end illumina reads.""",
        usage='-i <fastq file 1> -j <fastq file 2> -o <output prefix>')
    parser.add_option('-i', '--first', dest='filea'
        help='first fastq file')
    parser.add_option('-j', '--second', dest='fileb'
        help='second fastq file')
    parser.add_option('-o', '--output', dest='prefix'
        help='output prefix (omit to print to stdout)')
    parser.add_option('-t', '--threads', dest='threads', default=None,
        type=int, help='number of threads (default = all available)')
    parser.add_option('-p', '--pretty_output', dest='pretty',
        help='displays overlapping contigs in a nice way.')
    parser.add_option('-m', '--overlap', dest='overlap', default=30,
        help='minimum overlap (default = 30)', type = int)
    parser.add_option('-s', '--score', dest='score', default=90,
        help='minimum percent identity (default = 90)' type=int)

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
    
    p = Pool(options.threads)
    
    starttime = time()
    overlaps = 0
    
    for i in p.imap(doStitch, izip(Fasta(seqsa), Fasta(seqsb))):
        numtotes += 1
        if i.record:
            numcontigs += 1
            overlaps += i.overlap
            print >> outfile, '%s' % i.record,
        else:
            reca, recb = i.originals
            print >> dudsa, reca,
            print >> dudsb, recb,
            
    duration = time() - starttime
    
    print >> sys.stderr, \
        'Made %s contigs out of %s reads in %.2f seconds (%.2f per sec)' % \
        (numcontigs, numtotes, duration, numtotes/duration)
    print >> sys.stderr, \
        'Average overlap was %.2f' % (float(overlaps)/numcontigs)
        
    
class Stitch:
    ''' Stitches together two overlapping Illumina reads using Doubleblast '''
    def __init__(self, reca, recb):
        self.reca = reca
        self.recb = recb
        self.record = False
        self.overlap = 0
        self.find_overlaps()
        self.pretty = ''
        
    @property
    def originals(self):
        return (self.reca, self.recb)
        
    def find_overlaps(self):
        ''' Alignment algorithm, returns new DNA object of contig '''
        seqa, seqb = self.reca.seq, self.recb.revcomp
        seqb = self.recb.revcomp
        length = min([len(i) for i in (seqa, seqb)])
        scores, hits = [], {}
        
        for i in range(len(seqa)):
            score = 0
            for na, nb in zip(seqa[i-1:], seqb[:-i+1]):
                if na == nb:
                    score += 1
            hits[score] = i
        
        score = max(hits.keys())
        i = hits[score]
        
        self.overlap = score
        
        if ((score > 25)):
            
            beg = seqa[0:i-1]
            end = seqb[-i+1:]
            qbeg = self.reca.qual[0:i-1]
            qend = self.recb.qual[-i+1:][::-1]
            
            smida = seqa[i-1:]
            smidb = seqb[0:len(seqb)-i+1]
            qmida = self.reca.qual[i-1:]
            qmidb = self.recb.qual[0:len(seqb)-i+1][::-1]
            
            # Vote!
            mid, midq = [], []
            for (na, qa), (nb, qb) in zip(zip(smida, qmida), zip(smidb, qmidb)):
                if qa>qb:
                    mid+=na
                    midq+=qa
                elif qa<qb:
                    mid+=nb
                    midq+=qb
                else:
                    if qa==qb:
                        mid+=na
                        midq+=qa
                    else:
                        mid+='N'
                        midq+=qa
            
            newseq = beg + ''.join(mid) + end
            newqual = qbeg + ''.join(midq) + qend
            # Would it be faster to recycle an object?
            self.record = Dna(self.reca.header, newseq, newqual)


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