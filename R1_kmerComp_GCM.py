#Problem 1 - k-Mer Composition
#Guillermo Carrillo Martín - MSc Bioinformatics 22/12/2022 - Python 3.9.12

import re
import os

imput_file = "rosalind_kmer.txt" #imput file
klen = 4 # kmer length

#Eliminate \n jumps to have one strig for each record
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
            output.write(line[0:len(line)-1])

    fasta_file.close()
    output.close()

#Recursive function to put all kmers of lenth k in a dictionary
def kmer(k, base="", kmer_dic={}): 
    if k == 0:
        kmer_dic[base] = 0
        return kmer_dic
    else:
        for c in "ACGT": 
            kmer(k-1,base + c)
    return kmer_dic

#Calculating the kmer profile of a sequence
def kmer_prof(seq):
    kmer_num = len(seq)-klen+1
    all_kmerD = kmer(klen)

    for i in range(kmer_num):
        novo_kmer = seq[i:i+4]
        all_kmerD[novo_kmer] += 1
    
    return all_kmerD


#Printing the kmer profile in a txt file and removing the NoJump file
def print_profile(input):
    fasta = open(input, "rt")
    output = open("./R1_output.txt", "wt")

    fasta.readline() #Skip tag
    sequence = fasta.readline() #Read the sequence in the FASTA file

    final_dic = kmer_prof(sequence) #Compute the kmer profile of the sequence

    for x in final_dic.values(): #Print each dictionary value in a new file
        output.write(str(x) + " ")

    fasta.close() #Close the imput file
    os.remove("./NoJump.fasta") #Remove NoJump file

#Calling the funtions
fasta_no_jump(imput_file, "./NoJump.fasta")
print_profile("./NoJump.fasta")