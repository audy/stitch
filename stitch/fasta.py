from dnaobj import *

class Fasta:
    ''' iterates through a fasta or fastq file, returning dnaobj objects '''
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
                    yield Dnaobj(rec[0], rec[1], rec[3])
                    
        elif self.filetype == 'fasta':
            # implement this later
            pass
