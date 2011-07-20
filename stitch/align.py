from fasta import *

class Stitch:
    ''' Stitches together two overlapping Illumina reads
        using simple, ungapped alignment '''
        
    def __init__(self, reca, recb):
        self.reca = reca
        self.recb = recb

        self.record = False
        self.overlap = 0
        self.pretty = ''
        self.score = 0.0
        
        # run the alignment
        self.find_overlaps()
    
    @property
    def originals(self):
        return (self.reca, self.recb)
    
    @classmethod
    def find_overlaps(self):
        ''' Alignment algorithm, returns new DNA object of contig '''
        
        # reverse complement second sequence
        # this should be made into an option
        seqa, seqb = self.reca.seq, self.recb.revcomp
        seqb = self.recb.revcomp
        
        if len(seqa) != len(seqb):
            raise Exception, 'Stitch requires reads to be of the same length!'
        
        # get the minimum length for score generation
        length = min([len(i) for i in (seqa, seqb)])
        
        scores, hits = [], {}
        
        # Find overlap
        for i in range(len(seqa)):
            score = 0
            for na, nb in zip(seqa[i-1:], seqb[:-i+1]):
                if ('N' in (na, nb)):
                    continue
                if (na == nb):
                    score += 1
                if (na != nb):
                    score -= 1
            hits[score] = i
        
        # BUG: this scoring scheme sometimes favors small alignments over larger ones
        # with slightly worse scores. We prefer the larger ones. How to score differently
        # so that length is taken into account? Normalize by overlap lenght?
        
        score = max(hits.keys())
        i = hits[score] # i stands for index of best overlap
        
        self.overlap = len(seqa) - i
        self.score = float(score)/(len(seqa)-i)
        
        # fuggered code to generate the contig
        beg = seqa[0:i-1]
        end = seqb[-i+1:]
        qbeg = self.reca.qual[0:i-1]
        qend = self.recb.qual[::-1][-i+1:]
        smida = seqa[i-1:]
        smidb = seqb[0:len(seqb)-i+1]
        qmida = self.reca.qual[i-1:]
        qmidb = self.recb.qual[::-1][0:len(seqb)-i+1]
        mid, midq = [], []
        
        for (na, qa), (nb, qb) in \
                        zip(zip(smida, qmida), zip(smidb, qmidb)):
            if qa > qb:
                mid += na
                midq += qa
            elif qa < qb:
                mid += nb
                midq += qb
            else:
                if qa == qb:
                    mid += na
                    midq += qa
                else:
                    mid += 'N'
                    midq += qa
        
        newseq = beg + ''.join(mid) + end
        newqual = qbeg + ''.join(midq) + qend
        
        # generate pretty print view
        self.pretty = '1:%s\n2:%s\nC:%s\n' % \
                        (seqa + '-'*(i-1),
                        '-'*(i-1) + seqb, newseq)
        
        # create a new record sequence for the contig
        self.record = Dna(self.reca.header, newseq, newqual)