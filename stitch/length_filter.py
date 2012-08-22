from fasta import *
import sys

filename = sys.argv[1]

for minimum, maximum in ((50, 59), (60, 69), (70, 79), (80, 89), (90, 99), (100, 129), (130, 200)):
    handle = open(filename)
    output = open('%s-%s-%s' % (filename, minimum, maximum), 'w')
    counter = 0
    for fasta in Fasta(handle):
        if len(fasta) > maximum:
            continue
        if len(fasta) < minimum:
            continue
        print >> output, fasta
        counter += 1
    print minimum, maximum, counter
