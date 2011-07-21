#!/usr/bin/env python

from stitch import *

def test_headers():
    ''' make sure headers are parsed correctly '''
    with open('seq_a.txt') as handle:
        records = [i for i in Fasta(handle)]
        
    assert records.pop(0).header == 'sequence_1'

def test_alignments():        
    with open('seq_a.txt') as ha, open('seq_b.txt') as hb:
        records = [(a, b) for (a, b) in zip(Fasta(ha), Fasta(hb))]

    contigs = [ i.strip() for i in open('contigs.txt').readlines() ]

    for i, j in records:
        c = contigs.pop(0)
        res = Stitch(i, j).record.seq
        
        # TODO, also check quality scores
        
        try:
            assert res == c
        except AssertionError:
            print 'ERROR:    %s' % i.header
            print 'expected: %s' % i.seq
            print 'got:      %s' % res


if __name__ == '__main__':
    test_headers()
    test_alignments()