#!/usr/bin/env python

from stitch import *

def test_stitch():
    ''' make sure stitch algorithm behaves rationally '''
    
    ha = open('seq_a.txt')
    hb = open('seq_b.txt')
    
    records = ((a, b) for (a, b) in zip(Fasta(ha), Fasta(hb)))
    contigs = (i.strip() for i in open('contigs.txt').readlines())
    for i, j in records:
        c = contigs.next()
        res = Stitch(i, j).record.seq
        assert res == c
    
    return True
    
def test_dnaio():
    ''' make sure fastq and fasta parsing/output behave rationally '''
    return True
    
def run_tests():
    assert test_stitch()
    assert test_dnaio()
    
run_tests()