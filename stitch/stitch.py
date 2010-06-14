from fasta import *
from doubleblast import *
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
     
    for i in imap(doStitch, izip(Fasta(seqsa), Fasta(seqsb))):
        if i.hits:
            print i.record
        else:
            # Send to duds
            pass


def doStitch(recs):
    ''' Used by Pool.imap to create stitch jobs '''
    try:
        reca, recb = recs
        return Stitch(reca, recb)
    except KeyboardInterrupt:
        return
    
class Stitch:
    ''' Stitches together two overlapping Illumina reads using Doubleblast '''
    def __init__(self, reca, recb):
        self.reca = reca
        self.recb = recb
        self.hits = False
        self.contig = []
        self.quality = []
        results = Doubleblast.query(reca.seq, recb.revcomp)
        
        # BUG: See README.rst
        
        # Todo: Get rid of results that don't make sense.
        
        # Grab most e-valued
        if results:
            self.result = sorted(results, key=itemgetter('evalue'))[0]
            for key in self.result:
                setattr(self, key, self.result[key])
            self.record = self._generate_contig()
            
    def _generate_contig(self): 
        ''' Generate le contig '''
        if self.qstart >= self.sstart: return
        if self.qstart != 1: return
            
        self.hits = True

        subject, subqual = self.recb.seq, self.recb.qual
        query, quequal = self.reca.revcomp, self.reca.qual[::-1]    
                
        # The beginning.
        self.contig.append(query[0:self.sstart])
        self.quality.append(quequal[0:self.sstart])
        
        # The middle - TODO implement voting.
        self.contig.append(query[self.sstart:])
        self.quality.append(quequal[self.sstart:])
        # The end
        self.contig.append(subject[self.qend:])
        self.quality.append(subqual[self.qend:])
        
        finalseq, finalqual = [], []
        self.contig = ''.join(self.contig)
        self.quality = ''.join(self.quality)
        
        return Dna(self.reca.header, self.contig, self.quality)
            
        
if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        print >> sys.stderr, 'Ouch!'
        quit()