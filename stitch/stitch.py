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
            a = subject[0:ss]
            aq = subject[0:ss]

            # vote (overlap) region, all should be the same length
            bqs = query[0:len(subject)-ss+qs]
            bqq = quequal[0:len(subject)-ss+qs]
            bss = subject[ss-qs:]
            bsq = subqual[ss-qs:]
            
            b = []
            bq = []
            
            # seriously, someone shoot me.
            for ((qn, qq), (sn, sq)) in izip(izip(bqs, bqq), izip(bss, 0bsq)):
                if qq > sq:
                    b += qn
                    bq += qq
                elif qq < sq:
                    b += sn
                    bq += sq
                elif qq == sq:
                    if qn == sn:
                        b += sn
                        bq += sq
                    elif qn != sn:
                        if qn == 'N':
                            b += sn
                            bq += sq
                        elif sn == 'N':
                            b += qn
                            bq += qq
                        if (qn and sn) in ('A', 'G'):
                            b += 'R'
                            bq += max(qq, sq)
                        elif (qn and sn) in ('T', 'C'):
                            b += 'Y'
                            bq += max(qq, sq)
                        elif (qn and sn) in ('A', 'C'):
                            b += 'M'
                            bq += max(qq, sq)
                        elif (qn and sn) in ('G', 'C'):
                            b += 'S'
                            bq += max(qq, sq)
                        elif (qn and sn) in ('A', 'T'):
                            b += 'W'
                            bq += max(qq, sq)
                        else:
                            b += 'N'
                            bq += max(qq, sq)
                                        
            b = ''.join(b)
            
            # last region, comes from subject sequence
            c = query[len(subject)-qs+ss:]
            cq = quequal[len(subject)-qs+ss:]
            
            setattr(self, 'contig', a+b+c)
            
            
            
if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        print >> sys.stderr, 'Ouch!'
        quit()