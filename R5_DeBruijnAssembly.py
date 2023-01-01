#Problem 5 - Genome Assembly with Perfect Coverage
#Guillermo Carrillo Martín - MSc Bioinformatics 01/01/2023 - Python 3.9.12


file = "./R5_testSC.txt"

##1.Create a list with the sequences
input = open(file, "rt")

seq_list = [] #Create an empty list
for line in input: #For each line in the input file
    seq_list.append(line.replace("\n","")) #Append it to the list without "\n"

input.close() #Close the input file


##2.Constructing the DeBruijn Graph
#Create the substrings pairs [for example: TAGC = (TAG, AGC)] and write them in a list with 2 dimensions
debruijn_graph = []

for seq in seq_list: #For each sequence
    substring1 = seq[0:len(seq)-1] #Calculate substring 1
    substring2 = seq[1:len(seq)] #Calculate substring 2
    
    record = [substring1, substring2] 
    debruijn_graph.append(record) #Append both substrings in the list


##3.Sort DeBruijn Graph
#The idea is to concatenate each record so that the final substring of each record matches the initial substring of the next record 
debruijn_sorted = [debruijn_graph[0]] #The sorted DeBruijn graph starts with a random record (in this case, the first one)

for n in range(0,len(debruijn_graph)-1): #For a number of iterations equal to the number of records in the DeBruijn graph minus one...
    link_kmer = debruijn_sorted[n][1] #link_kmer is equal to the final substring of the last record
    for record in debruijn_graph: #Nested loop: For each record in the DeBruijn graph
        if record[0] == link_kmer: #If the final substring of the last record is equal to the first substring of the present record
            debruijn_sorted.append(record) #Append that record in the sorted list


##4.Solving the DeBruijn graph into the cyclic superstring
#Now, 
cyclic_superstring = debruijn_sorted[0][0] + debruijn_sorted[0][1][-1]

for record in debruijn_sorted[1:len(debruijn_sorted)]: #For each missing record...
    cyclic_superstring += record[1][-1] #Add the last base to the string


#Eliminate last k-1 bases because of the 
kmer_len = len(seq_list[0]) #The length of the first kmer (all kmers have the same lenth)
def_cyclic_superstring = cyclic_superstring[0:-kmer_len+1]

#Write the output in a file
out = open("./R5_output.txt", "wt")
out.write(def_cyclic_superstring)
out.close()