from fastitr import *
from doubleblast import *
from itertools import izip
from optparse import *
from multiprocessing import *
    
def main():
    parser = OptionParser(description="stitch.py")
    parser.add_option('-i', '--first', dest='filea')
    parser.add_option('-j', '--second', dest='fileb')

    (options, args) = parser.parse_args()

    seqsa = open(options.filea, 'r')
    seqsb = open(options.fileb, 'r')
    
    p = Pool()
    
    for i in p.imap(doStitch, izip(Fastitr(seqsa), \
            Fastitr(seqsb))):
        print dir(i)


def doStitch(recs):
    reca, recb = recs
    return Stitch.stitch(reca, recb)
        
class Stitch:
    @classmethod
    def stitch(self, reca, recb):
        self.reca = reca
        self.recb = recb
        self.result = Doubleblast.query(reca, recb)
        if result: 
            for key in result:
                setattr(self, key, result[key])
        return self
        
if __name__ == '__main__':
    main()
    

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