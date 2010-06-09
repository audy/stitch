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
    def query(self, sub, que):
        ''' bl2seq '''
        seed()
        fname = 's.%s' % (hex(randint(0,1048575))[2:])
        results = []
        
        with open(fname, 'w') as sfile:
            print >> sfile, '%s' % sub.seq

        try:
            pipe = Popen([BLAST,
            '-task', 'blastn-short',
            '-subject', fname,
            '-num_threads', THREADS,
            '-outfmt', '6'],
            stdin=PIPE, stdout=PIPE)
            pipe.stdin.write('%s' % que.revcomp)
            hits = pipe.communicate()[0]
        except OSError:
            raise NoBlast
        finally:
            getoutput('rm %s' % fname)

        if not hits:
            return None
        
        for line in hits.split('\n'):
            if not line: continue
            percent_identity, length, mismatches, gap_openings, \
                qstart, qend, sstart, send, evalue, bitscore \
                = line.split('\t')[2:]
            result = {'percent_identity': float(percent_identity),
                      'alignment_length': int(length),
                      'mismatches': int(mismatches),
                      'gap_openings': int(gap_openings),
                      'query_start': int(qstart),
                      'query_end': int(qend),
                      'subject_start': int(sstart),
                      'subject_end': int(send),
                      'e-value': float(evalue),
                      'bitscore': float(bitscore) }
            results.append(result)
        
            
        results = sorted(results, key=itemgetter('e-value'))

        return results[0]
        