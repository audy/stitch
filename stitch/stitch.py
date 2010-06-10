from fastitr import *
from doubleblast import *
from itertools import izip, imap
from optparse import *
from multiprocessing import *
import os
import sys
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
            print i.reca.seq
            print i.recb.revcomp
            pass
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
            self.generate_contig()
    def generate_contig(self): 
        # i am the shittiest programmer in the world!
        if self.qstart < self.sstart:
            self.hits = True
            subject, subqual = self.reca.seq, self.reca.qual
            query, quequal = self.recb.revcomp, self.recb.qual
            qs, ss, qe, se = self.qstart, self.sstart, self.qend, self.send
            print qs, ss, qe, se
            
            # first region, comes from query sequence
            a = query[0:(qs-ss)]
            aq = quequal[0:(qs-ss)]

            # vote (overlap) region, all should be the same length
            bas = query[qs-ss:len(query)]
            baq = quequal[qs-ss:len(query)]
            bbs = subject[0:len(query)+ss-qs]
            bbq = subqual[0:len(query)+ss-qs]
            
            for i in (bas, baq, bbs, bbq): print len(i),
            
            # last region, comes from subject sequence
            c = subject[len(query):len(query)-qe]
            cq = subqual[len(query):len(query)-qe]
            setattr(self, 'contig', '')

            
            
            
if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        print >> sys.stderr, 'Ouch!'
        quit()