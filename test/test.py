#!/usr/bin/env python

from stitch import *
import unittest


class TestDnaioFunctions(unittest.TestCase):
    ''' test DNA IO '''
    def test_headers(self):
        ''' make sure headers are parsed correctly '''
        with open('seq_a.txt') as handle:
            records = [i for i in Fasta(handle)]
            
        self.assertEqual(records.pop(0).header, 'sequence_1')


class TestAlignments(unittest.TestCase):
    ''' test alignment functions for rationality '''
    
    def test_alignments(self):        
        with open('seq_a.txt') as ha, open('seq_b.txt') as hb:
            records = [(a, b) for (a, b) in zip(Fasta(ha), Fasta(hb))]

        contigs = [ i.strip() for i in open('contigs.txt').readlines() ]

        for i, j in records:
            c = contigs.pop(0)
            res = Stitch(i, j).record.seq
            self.assertEqual(res, c)


if __name__ == '__main__':
    unittest.main()