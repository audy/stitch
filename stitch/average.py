from fasta import *
from dna import *
import sys

def main(argv):
    infile = sys.argv[1]

    scores = [0]*100
    counter = 0
    
    with open(infile, 'r') as handle:
        for record in Fasta(handle):
            for score, i in zip(record.qual, range(len(scores))):
                score = ord(score)
                scores[i] += score
            counter += 1
            
    print [ i/counter for i in scores ]
    
        

if __name__ == '__main__':
    try:
        main(sys.argv)
    except KeyboardInterrupt:
        print 'Ouch!'
        quit()