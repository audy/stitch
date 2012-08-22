import string
_complement = string.maketrans('GATCRYgatcry','CTAGYRctagyr')

class Fasta:
    ''' iterates through a fastq file, returning dnaobj objects '''

    def __init__(self, handle, filetype='fastq'):
        self.filetype = filetype
        self.handle = handle
        self.county = 0

    def __iter__(self):
        if self.filetype == 'fastq':
            counter = 0
            rec = { 0: '', 1: '', 2: '', 3: '' }
            for line in self.handle:
                if counter < 3:
                    rec[counter] = line.strip()
                    counter += 1
                elif counter == 3:
                    rec[counter] = line.strip()
                    counter = 0
                    yield Dna(rec[0], rec[1], rec[3])


class Dna:
    ''' An object representing either a FASTA or FASTQ record '''

    def __init__(self, header, sequence, quality = False):
        self.header = header.lstrip('@').rstrip('\n')
        self.seq = sequence
        self.qual = quality
        if quality:
            self.type = 'fastq'
        else:
            self.type = 'fasta'

        if len(self.seq) != len(self.qual):
            raise IOError, \
                'Seq length and qual length do not agree: %s' % (self.header)

    def __str__(self):
        ''' returns a FASTA/Q formatted string '''
        if not self.qual:
            return ('>%s\nself.sequence\n') % \
                (self.header, self.seq)
        else:
            return('@%s\n%s\n+%s\n%s') % \
                (self.header, self.seq, self.header, self.qual)

    def __len__(self):
        return len(self.seq)

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

    @property
    def rqual(self):
        ''' returns reverse quality'''
        return self.qual[::-1]
