from Bio import AlignIO
import sys

inputFasta=open(sys.argv[1], "r")
print(inputFasta)
inputLines = inputFasta.read()
lines = inputLines.split("\n")

length=lines[0].split(" ")[1]

print(len(lines))
#print(lines[0]) 
#print(lines[-1])
#print(lines[-2])
print(length)
first_interval=int(length)/3
second_interval=first_interval*2
last_interval=first_interval*3

print(first_interval)
print(second_interval)
print(last_interval)

with open('window1','w') as w1, open('window2','w') as w2, open('window3','w') as w3:
    for line in lines[1:]:
        if len(line)!=0:
            cols=line.split(" ")
            seq=cols[1]
            id=cols[0]
            if id[0]==">":
                pass
            else:
                id=">"+id
            s1=str(seq[0:first_interval])
            s2=str(seq[first_interval:second_interval])
            s3=str(seq[second_interval:])
            w1.write(id+"\n")
            w1.write(s1+"\n")
            w2.write(id+"\n")
            w2.write(s2+"\n")
            w3.write(id+"\n")
            w3.write(s3+"\n")
