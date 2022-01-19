import linecache
import argparse
import os
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

parser = argparse.ArgumentParser(
     prog='createFasta.py',
     usage='''python createFasta.py --fasta [fasta file given by braker] --path [Path of fasta file]''',
     description='''This program pulls out specific sequences from a fasta file, given the fasta file and a list of sequences saved in a text file''',
     epilog='''It requires numpy and biopython libraries''')

parser.add_argument('--fasta', type=str, help='The name of the fasta file', required=True)
parser.add_argument('--path', type=str, help='The path of the fasta file', required=False)

args=parser.parse_args()
fastapath=args.path
fasta=args.fasta

if fastapath==None:
    fastafile=fasta
else:
    fastafile=os.path.join(fastapath, fasta)

id_dict=SeqIO.to_dict(SeqIO.parse(fastafile, "fasta"))
identifiers=[id_dict[seq_record].id for seq_record in id_dict]

print identifiers[:5]

with open('window1','w') as w1, open('window2','w') as w2, open('window3','w') as w3:
    for record in id_dict:
        s1=str(id_dict[record].seq[0:333])
        s2=str(id_dict[record].seq[333:666])
        s3=str(id_dict[record].seq[666:])
        w1.write(">"+id_dict[record].id+"\n")
        w1.write(s1+"\n")
        w2.write(">"+id_dict[record].id+"\n")
        w2.write(s2+"\n")
        w3.write(">"+id_dict[record].id+"\n")
        w3.write(s3+"\n")



        
