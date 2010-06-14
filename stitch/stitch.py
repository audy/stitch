from fastitr import *
from doubleblast import *
from itertools import izip, imap
from optparse import *
from multiprocessing import *
import os
import sys
import time
    
def main():
    parser = OptionParser(
        description="""Stitch - Tool for creating contigs from overlapping
paried-end illumina reads.""",
        usage='-i fastqfile1 -j fastqfile2')
    parser.add_option('-i', '--first', dest='filea')
    parser.add_option('-j', '--second', dest='fileb')

    (options, args) = parser.parse_args()
    
    if not (options.filea or options.fileb):
        print >> sys.stderr, 'Usage: %s %s' % \
            (parser.get_prog_name(), parser.usage)
        quit()
    
    seqsa = open(options.filea, 'r')
    seqsb = open(options.fileb, 'r')
    
    p = Pool()
        
    for i in imap(doStitch, izip(Fastitr(seqsa), \
            Fastitr(seqsb))):
        if i.hits:
            print i.reca.seq
            print i.recb.revcomp
            print i.contig
        else:
            # Send to duds
            print '.',

def doStitch(recs):
    reca, recb = recs
    return Stitch(reca, recb)

class Stitch:
    def __init__(self, reca, recb):
        self.reca = reca
        self.recb = recb
        self.hits = False
        result = Doubleblast.query(reca, recb)
        if result:
            for key in result:
                setattr(self, key, result[key])
            self._generate_contig()
            
    def _generate_contig(self): 
        if self.qstart < self.sstart:
            self.hits = True
            subject, subqual = self.recb.seq, self.recb.qual
            query, quequal = self.reca.revcomp, self.reca.qual[::-1]
            
            
            
            for nuc in subject[]
            
            finalsec, finalqual = [], []
            
            contig = { "seq": finalseq, "qual": finalqual}
           
            setattr(self, 'contig', contig)
            
            
            
if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        print >> sys.stderr, 'Ouch!'
        quit()