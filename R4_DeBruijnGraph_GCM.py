#Problem 4 - Constructing a De Bruijn Graph
#Guillermo Carrillo Martín - MSc Bioinformatics 27/12/2022 - Python 3.9.12

file = "./R4_test.txt"

##1.Build a sorted list with the estandar and the reverse_complement sequences

#Create a set wih the input file sequences to eliminate duplicated sequences
input = open(file, "rt")

seq_set = set() #Create an empty set
for line in input: 
    seq_set.add(line.replace("\n","")) #Append each line of the input file in a set

input.close()

#Create another set with the reverse complement
def reverse_complement(sequence): #Reverse-complement function
    complement = { #Complement dictionary
    "A": "T", 
    "C": "G", 
    "G": "C", 
    "T": "A"
    }

    rev_seq = sequence[::-1] #Reverse the sequence
    rc_seq = ""
    for i in rev_seq: #For each nucleotide in the reverse sequence...
        rc_seq += complement[i] #Append the complement nucleotide
    return rc_seq

rc_set = set() #Create an empty set
for seq in seq_set: #For each sequence in the initial set...
    rc_seq = reverse_complement(seq) #Calculate the reverse-complement
    rc_set.add(rc_seq) #Append it to the reverse-complement set


#Convert the set to list, in orther to sort them
seq_set.update(rc_set) #Merge two sets (avoiding duplicate sequences)

seq_list = list(seq_set) #Convert the set to a list
seq_list.sort() #Sort it


##2. Constructing the DeBruijn Graph

#Create the substrings pairs: from TAGC = (TAG, AGC) and write them in a text file
out = open("./R4_output.txt", "wt") #Create the output file

blank = "Yes" #At the beggining, the output file has no lines

for seq in seq_list: #For each sequence
    substring1 = seq[0:len(seq)-1] #Calculating substring 1
    substring2 = seq[1:len(seq)] #Calculating substring 2
    
    if blank == "Yes": 
        blank = "No"
    else: #If there is not the first sequence
        out.write("\n") #Write an "\n" 
    out.write("(" + substring1 + ", " + substring2 + ")") #Write both substrings in the file

out.close()