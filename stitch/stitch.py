from fastitr import *
from doubleblast import *
from time import *
import sys
from itertools import izip

def main(argv):
    fastqfile1, fastqfile2, contigfile, uncontigfile = argv[1:]
    start_time = time()
    handle = open(fastqfile1, 'r')
    handle2 = open(fastqfile2, 'r')
    contig_output = open(contigfile, 'w')
    uncontig1_output = open(uncontigfile + '1', 'w')
    uncontig2_output = open(uncontigfile + '2', 'w')
    for recorda, recordb in izip(Fastitr(handle, filetype='fastq'),\
 		Fastitr(handle2, filetype='fastq')):
        blobj = Blast(recorda, recordb)
        blobj.bl2seq()
        blobj.parse()
        contig = blobj.return_contig()
        print contig
        if not contig == ():
            result = dnaobj(recorda.header, contig[0], contig[1])
            contig_output.write(str(result))
        else:
            uncontig1_output.write(str(recorda))
            uncontig2_output.write(str(recordb))
    handle.close()
    handle2.close()
    stop_time = time()
    print start_time, stop_time, (stop_time - start_time)
    
if __name__ == '__main__':
    main(sys.argv)
