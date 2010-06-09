from fastitr import *
from doubleblast import *
from itertools import izip, imap
from optparse import *
from multiprocessing import *
import os
import time
    
def main():
    parser = OptionParser(description="stitch.py")
    parser.add_option('-i', '--first', dest='filea')
    parser.add_option('-j', '--second', dest='fileb')

    (options, args) = parser.parse_args()

    seqsa = open(options.filea, 'r')
    seqsb = open(options.fileb, 'r')
    
    p = Pool()
        
    for i in imap(doStitch, izip(Fastitr(seqsa), \
            Fastitr(seqsb))):
        if i.hits:        
            # Send to contigs
            print '%g %s %s' % (i.evalue, i.bitscore, i.gaps)
            print 'q.start=%s, q.end=%s, s.start=%s, s.end=%s' % \
                (i.qstart, i.qend, i.sstart, i.send)
            print '%s\n%s' % (i.reca.seq, i.recb.seq)
            print '%s, %s\n' % (i.contig, len(i.contig))
        else:
            # Send to duds
            pass


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
            self.hits = True
            for key in result:
                setattr(self, key, result[key])
            self.generate_contig()
    def generate_contig(self):
        if self.qstart < self.sstart:
            print 'qstart<sstart'
            self.contig = self.recb.revcomp[:self.sstart-self.qstart] + \
                self.reca.seq
            self.qual = self.recb.rqual[:self.sstart-self.qstart] + \
                    self.reca.qual            
        elif self.qstart > self.sstart:
            print 'qstart>sstart'
            self.contig = self.reca.seq + \
                self.recb.revcomp[self.send-self.qend:]
            self.qual = self.reca.qual + \
                self.recb.rqual[self.send-self.qend:]
        elif self.sstart == self.qstart:
            self.contig = self.reca.seq
            self.qual = self.reca.seq
if __name__ == '__main__':
    main()
    