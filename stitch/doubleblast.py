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
        fname = 's.%s' % (hex(randint(0,65535))[2:])
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
            hits = pipe.communicate()
        except OSError:
            raise NoBlast
        finally:
            getoutput('rm %s' % fname)
        if not hits[0]:
            print 'no hits'
        else:
            print hits[0],
        print '%s\n%s\n' % (sub.seq, que.revcomp)
        
        
