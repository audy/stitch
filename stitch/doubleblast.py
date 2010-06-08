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
THREADS = '1'

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
            print >> sfile, '%s' % sub.revcomp
        try:
            pipe = Popen([BLAST,
            '-subject', fname,
            '-num_threads', THREADS,
            '-outfmt', '6'],
            stdin=PIPE, stdout=PIPE)
            pipe.stdin.write('%s' % que.seq)
            hits = pipe.communicate()[0].split('\n')
        except OSError:
            raise NoBlast
        finally:
            getoutput('rm %s' % fname)
        print hits
        
        
