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

BLAST = 'blastn'
THREADS = '1'

class NoBlast(Exception):
	return 'BLAST (%s) either not found or wrong version!' % BLAST

class Doubeblast:
	''' Doubleblastdoc '''
	@classmethod
	def query(self, sub, que):
		''' bl2seq '''
		seed()
		id = hex(randint())[2:]
		with open('s.%s' % id, 'w') as sfile:
			print >> sfile, '%s' % sub.seq
		try:
			pipe = Popen([BLAST,
			'-subject', self.subfile.
			'-num_threads', THREADS,
			'-outfmt', '6'],
			stdin=PIPE, stdout=PIPE
			)
			pipe.stdin.write('%s' % query)
			hits = pipe.communicate()[0].split('\n')
		except OSError:
			raise NoBlast
        
        
