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
        return Doubleblast.query(reca, recb)
		    
if __name__ == '__main__':
    main(sys.argv)
