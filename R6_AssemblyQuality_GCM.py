#Problem 6 - Assessing Assembly Quality with N50 and N75
#Guillermo Carrillo Martín - MSc Bioinformatics 01/01/2023 - Python 3.9.12

file = "./R6_test.txt"

##1.Create a list with the sequences
input = open(file, "rt")

seq_list = [] #Create an empty list
for line in input: 
    seq_list.append(line.replace("\n","")) #Append each line of the input file in the list

input.close() #Close the input file


##2.Sort the list by the sequence length
def len_values(val): #Create a function that returns the lenght of the sequences
  return len(val)

seq_list.sort(reverse=True, key=len_values) #Sort the sequence by the length


##3.Calculate N50 and N75
total_bp = 0 #Calculate the total bp of the contigs
for seq in seq_list: 
    total_bp += len(seq) #Sum the length of each contig

def N50_calculator(seq_records): #N50 Calculator
    bp_count = 0
    for seq in seq_records:
        bp_count += len(seq) #Sum the length of the contigs until...
        if bp_count >= total_bp*0.5: #... The length is greater or equal than half of the total bp number
            N50 = len(seq) #Then, return the length of the last contig as N50 value
            return N50

def N75_calculator(seq_records): #N75 Calculator
    bp_count = 0
    for seq in seq_records:
        bp_count += len(seq) #Sum the length of the contigs until...
        if bp_count >= total_bp*0.75: #... The length is greater or equal than 75% of the total bp number
            N75 = len(seq) #Then, return the length of the last contig as N75 value
            return N75

N50 = N50_calculator(seq_list)
N75 = N75_calculator(seq_list)

print(str(N50) + " " + str(N75)) #Print output