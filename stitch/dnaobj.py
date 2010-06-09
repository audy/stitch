'''
dnaobj.py

An object representing a FASTA or FASTQ record.

Austin Glenn Davis-Richardson
austingdr@gmail.com
Triplett Lab, University of Florida
'''

import string

_complement = string.maketrans('GATCRYgatcry','CTAGYRctagyr')

class Dnaobj:
    ''' An object representing either a FASTA or FASTQ record '''
    def __init__(self, header, sequence, quality = False):
        self.header = header[1:]
        self.seq = sequence
        self.quality = quality
        if quality:
            self.type = 'fastq'
        else:
            self.type = 'fasta'
    def __str__(self):
        ''' returns a FASTA/Q formatted string '''
        if not self.quality:
            return ('>%s\nself.sequence\n') % \
                (self.header, self.seq)
        else:
            return('@%s\n%s\n+%s\n%s\n') % \
                (self.header, self.seq, self.header, self.quality)
    def __repr__(self):
        return '<dnaobj.%s instance: %s>' % (self.type, self.header)

    @property
    def complement(self):
        ''' returns complement of sequence '''
        return self.seq.translate(_complement)

    @property 
    def revcomp(self):
        ''' returns reverse complement of sequence '''
        return self.complement[::-1]
        
