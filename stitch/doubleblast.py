"""
doubleblast.py

Uses NCBI+ bl2seq to perform a blast between two sequences.
Takes DNAobjects from dnaobj.py

Austin Glenn Davis-Richardson
austingdr@gmail.com
Triplett Lab, University of Florida
"""

import commands
import random

blastprogram = 'bl2seq'

class Blast:
    ''' 1 v 1 sequence blaster '''
    def __init__(self, dnaobj1, dnaobj2):
        self.blobj1 = dnaobj1
        self.blobj2 = dnaobj2
        self.query = self.blobj1.sequence
        self.subject = self.blobj2.revcomplement()
        self.result_list = []

    def bl2seq(self):        
        ''' performs bl2seq command, returns raw output '''
        
        suffix = str(random.randint(0,1000))
        
        # XXX This actually writes files.  Find a way to get around this.  BECAUSE IO IS 'SPENSIVE! 
        
        with open('seqa' + suffix, 'w') as queryfile:
            queryfile.write(self.query)
            
        with open('seqb' + suffix, 'w') as subjectfile:
            subjectfile.write(self.subject)
            
        command = '%s -i %s -j %s -p blastn -D 1' % (blastprogram, 'seqa' + suffix, 'seqb' + suffix)
        
        self.output = commands.getoutput(command).split('\n')
        
        delcommand = 'rm %s %s' % ('seqa' + suffix, 'seqb' + suffix)
        
        commands.getoutput(delcommand)
        
    def parse(self):
        ''' parses output '''

        if not self.output:
            raise Exception, 'do blast.bl2seq() first!'

        result_list = []

        for line in self.output[3:]:
            percent_identity, alignment_length, mismatches, gap_openings, qstart, qend, sstart, send, evalue, bitscore = line.split('\t')[2:]
            
            results = {
            'percent_identity': float(percent_identity),
            'alignment_length': int(alignment_length),
            'mismatches': int(mismatches),
            'gap_openings': int(gap_openings),
            'query_start': int(qstart),
            'query_end': int(qend),
            'subject_start': int(sstart),
            'subject_end': int(send),
            'e-value': float(evalue),
            'bitscore': float(bitscore)
            }
            
            self.result_list.append(results)

        return self.result_list

        
    def return_contig(self):
        ''' takes blast results -> ('sequence', 'quality')'''
        # What if there is no result?  We must return something that tells main() so he can just write out the regular sequences.
        
        if self.result_list == []:
            return ()
            
        # return contig as a tuple, create the object in main, or don't!  It might be faster to skip the object completely.  But then why did I bother writing it?  For fun.
        
        best_score = 0.0
        best_item = {}
     
        for result in self.result_list:
            if result['e-value'] > best_score:
                best_item = result

        qstart = best_item['query_start']
        sstart = best_item['subject_start']
        qend = best_item['query_end']
        send = best_item['subject_end']

        if qstart < sstart:     
            contig_sequence = self.blobj2.revcomplement()[:sstart-qstart] + self.blobj1.sequence
            contig_quality = self.blobj2.quality[::-1][:sstart-qstart] + self.blobj1.quality            

        elif qstart > sstart:            
            contig_sequence = self.blobj1.sequence + self.blobj2.revcomplement()[send-qend:]
            contig_quality = self.blobj1.quality + self.blobj2.quality[::-1][send-qend:]

            
        elif sstart == qstart:
            return (self.blobj1.sequence, self.blobj1.quality)
            
        return (contig_sequence, contig_quality)
        
      
        
        
