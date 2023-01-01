#Problem 2 - Overlap Graphs
#Guillermo Carrillo Martín - MSc Bioinformatics 22/12/2022 - Python 3.9.12

import re
import os

k = 5 #Minimum number of overlapping matches between 2 sequences
input = "R2_test.fasta" #Input file

##0.Eliminate \n jumps to have one strig for each record
def fasta_no_jump(imput, out):
    fasta_file = open(imput, "rt")
    output = open(out, "wt")

    lin = fasta_file.readline()
    output.write(lin)

    for line in fasta_file:
        match = re.search(r"^>",line)
        if match:
            output.write("\n" + line)
        else:
            output.write(line.replace("\n",""))

    fasta_file.close()
    output.close()

fasta_no_jump(input, "./NoJump.fasta")
input1 = open("./NoJump.fasta", "rt")

##1.Append the tags and the sequences of the input(NoJump) file into two lists
tag_lst = []
seq_lst = []

while True: #for each record
    tag = input1.readline() #Tag
    if not tag: #If there is no tag, break the loop
        break
    seq = input1.readline() #Sequence
    
    tag_lst.append(tag[1:len(tag)].replace("\n","")) #Append tag without > and \n 
    seq_lst.append(seq.replace("\n","")) #Append sequence without \n

input1.close() #close NoJump file
os.remove("./NoJump.fasta") #Remove NoJump file

##2.Create a file with the overlap graph
out = open("./R2_output.txt", "wt") #create output file
blank = "Yes" #At the beggining, the output file has no lines

for x in range(0,len(seq_lst)): #For each pair of sequences in the list (s1 and s2)
    for y in range(0,len(seq_lst)): 
        s1 = seq_lst[x]
        s2 = seq_lst[y]
        if s1 != s2: #If they are not the same sequence...
            if s1[len(s1)-k:len(s1)] == s2[0:k]: #And the end of s1 overlaps the beggining of s2 in "k" bp
                tag1 = tag_lst[x] #tag of s1 is tag1
                tag2 = tag_lst[y] #tag of s2 is tag2

                if blank == "Yes": #If the output is blank, skip the first \n
                    blank = "No" 
                else:
                    out.write("\n") #If not, writte \n in outpt file
                
                out.write(tag1 + " " + tag2) #Writte tail and head tags (v,w) in output file

out.close() #close output file