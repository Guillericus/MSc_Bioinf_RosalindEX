#Problem 3 - Assembly as Shortest Super String 
#Guillermo Carrillo Martín - MSc Bioinformatics 25/12/2022 - Python 3.9.12

import re
import os

input = "R3_test.fasta" #Input file

##0.Eliminate \n jumps to have one strig for each record

def fasta_no_jump(imput, out):
    fasta_file = open(imput, "rt")
    output = open(out, "wt")

    first_line = fasta_file.readline()
    output.write(first_line)

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


##2.Create the overlapping graph (with the number of overlapping bp) as an array

#Calculate the max overlapping bp between 2 sequences
def max_overlap(seq1, seq2):
    min_l = min(len(seq1), len(seq2)) 

    for k in range(0,min_l+1): #for a sequence length (k) from 0 to the minimum length of the two sequences
        if seq1[len(seq1)-k:len(seq1)] == seq2[0:k]: #if the k last bases of s1 are equal to the first k bases of s2
            overlap_bp = k #k is equal to the overlap bp; if there are no overlapping bp, k = 0

    return overlap_bp

#Calculate the overlapping graph as an Array in the way (seq1, seq2, overlap bp)
def overlap_graph_generator(seql, tagl):
    overlap_graph_list = []

    for x in range(0,len(seql)): #For each pair of sequences in the list (s1 and s2)
        for y in range(0,len(seql)): 
            s1 = seql[x]
            s2 = seql[y]
            if s1 != s2: #If they are not the same sequence...
                overlap_bp = max_overlap(s1,s2) #calculate the max overlap
                tag1 = tagl[x] #tag of s1 is tag1
                tag2 = tagl[y] #tag of s2 is tag2

                new_entry = [tag1, tag2, overlap_bp]
                overlap_graph_list.append(new_entry) #Append tail and head tags (v,w) + number of overlapping bp in the list 
    
    return overlap_graph_list

##3. Assembly all the sequences 

#Assembly the two most overlapped sequences into one
def merge_most_overlap(sequence_list, tag_list, ov_graph):
    overlap_val = []

    for record in ov_graph: #for each record
        overlap_val.append(record[2])
    index = overlap_val.index(max(overlap_val)) #calculate the index of the maximum overlapping record
    #If two or more records have the maximum overlapping bp, it chooses the first index to be found

    overlap_record = ov_graph[index] #Extract, from the most overlapping record:
    s1 = sequence_list[tag_list.index(overlap_record[0])] #Tail sequence
    s2 = sequence_list[tag_list.index(overlap_record[1])] #Head sequence
    n = overlap_record[2] #Overlapping bp
    
    new_sequence = s1[0:len(s1)-n] + s2 #Assembly the two overlap sequences
    return new_sequence, overlap_record

#Recursively assembling the genome
def recursive_assembly(seq_lst, tag_lst, count = 0):
    
    if len(seq_lst) == 1: #If there is only one string in the list, return it as he assembly sequence
        assembly = seq_lst[0]
        return assembly
    
    else: #If there are 2 or more...
        overlap_graph = overlap_graph_generator(seq_lst, tag_lst) #Calculate the overlapping graph
        new_seq, record = merge_most_overlap(seq_lst, tag_lst, overlap_graph) #Calculate the new sequence and return the record info

        count += 1 
        new_tag = "MergedSeq_" + str(count) #Create a new tag for the new seqence

        seq_lst.pop(tag_lst.index(record[0])) #Eliminate s1 and tag1
        tag_lst.pop(tag_lst.index(record[0]))

        seq_lst.pop(tag_lst.index(record[1])) #Eliminate s2 and tag2
        tag_lst.pop(tag_lst.index(record[1]))
        
        seq_lst.append(new_seq) #Append new merged sequence
        tag_lst.append(new_tag) #Append new created tag
    
        assembly = recursive_assembly(seq_lst, tag_lst, count) #Repeat the function with the new lists
    
    return assembly

#Writte the result into a file
out = open("./R3_output.txt", "wt")

shortest_superstring = recursive_assembly(seq_lst, tag_lst) 
out.write(shortest_superstring)

out.close()