from fastitr import *
from doubleblast import *
from time import *
import sys
from itertools import izip

def main(argv):
    filea, fileb = sys.argv[1:]
    seqsa = open(filea, 'r')
    seqsb = open(fileb, 'r')
    
    for reca, recb in izip(Fastitr(seqsa, filetype='fastq'), \
        Fastitr(seqsb, filetype='fastq')):
        
        a = Stitch.stitch(reca, recb)
        print a
        
class Stitch:
    @classmethod
    def stitch(self, reca, recb):
        result = Doubleblast.query(reca, recb)
        if not result: return None
        qstart = result['query_start']
        sstart = result['subject_start']
        qend = result['query_end']
        send = result['subject_end']
        
        print qstart, qend, sstart, send
        print reca.revcomp
        print recb.seq
        
if __name__ == '__main__':
    main(sys.argv)
    
    


def return_contig(self):
    ''' takes blast results -> ('sequence', 'quality')'''
    # What if there is no result?  We must return something that tells
        # main() so he can just write out the regular sequences.
    if self.result_list == []:
        return ()
    # return contig as a tuple, create the object in main, or don't!
    # It might be faster to skip the object completely.
    # But then why did I bother writing it?  For fun.
    best_score = 0.0
    best_item = {}
    for result in self.result_list:
        if result['e-value'] > best_score:
            best_item = result
    qstart = best_item['query_start']
    sstart = best_item['subject_start']
    qend = best_item['query_end']
    send = best_item['subject_end']
    if qstart < sstart:     
        contig_sequence = self.blobj2.revcomplement()[:sstart-qstart] + self.blobj1.sequence
        contig_quality = self.blobj2.quality[::-1][:sstart-qstart] + self.blobj1.quality            
    elif qstart > sstart:            
        contig_sequence = self.blobj1.sequence + self.blobj2.revcomplement()[send-qend:]
        contig_quality = self.blobj1.quality + self.blobj2.quality[::-1][send-qend:]
    elif sstart == qstart:
        return (self.blobj1.sequence, self.blobj1.quality)
    return (contig_sequence, contig_quality)