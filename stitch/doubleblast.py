"""
doubleblast.py

Uses NCBI+ bl2seq to perform a blast between two sequences.
Takes DNAobjects from dnaobj.py

Austin Glenn Davis-Richardson
austingdr@gmail.com
Triplett Lab, University of Florida
"""

from subprocess import Popen, PIPE
from random import randint, seed
from commands import getoutput
from operator import itemgetter


BLAST = 'blastn'
THREADS = '2'

class NoBlast(Exception):
    def __str__(self):
        return 'BLAST (%s) either not found or wrong version!' % BLAST

class Doubleblast:
    ''' Doubleblastdoc '''
    @classmethod
    def query(self, seqa, seqb):
        ''' bl2seq '''

        seed()
        fname = 's.%s' % (hex(randint(0,1048575))[2:])
        results = []

        try:
            with open(fname, 'w') as sfile:
                print >> sfile, '%s' % seqb
            pipe = Popen([BLAST,
            '-task', 'blastn-short',
            '-subject', fname,
            '-num_threads', THREADS,
            '-outfmt', '6',
            '-strand', 'plus'],
            stdin=PIPE, stdout=PIPE)
            pipe.stdin.write('%s' % seqa)
            self.hits = pipe.communicate()[0]
        except OSError:
            raise NoBlast
        finally:
            getoutput('rm %s' % fname)

        if not self.hits:
            return None
                    
        for line in self.hits.split('\n'):
            if not line: continue
            percent_identity, length, mismatches, gap_openings, \
                qstart, qend, sstart, send, evalue, bitscore \
                = line.split('\t')[2:]
            result = {'identity': float(percent_identity),
                      'length': int(length),
                      'mismatches': int(mismatches),
                      'gaps': int(gap_openings),
                      'qstart': int(qstart),
                      'qend': int(qend),
                      'sstart': int(sstart),
                      'send': int(send),
                      'evalue': float(evalue),
                      'bitscore': float(bitscore),
                      'raw_output': self.hits}
            results.append(result)
            
        return results