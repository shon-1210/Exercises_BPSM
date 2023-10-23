#!/usr/bin/python3


# Program to take a DNA sequence and calculate A+T content
# Written by s2572046 on 20 Oct 2023
# -------------------------------------------- #
# Input



#1. Calculating A+T content
# Have a look at the sequence, do an estimate of the A+T content visually: maybe somewhere around 50% (0.5) or so?
# We need to store the information we need - the A and T count, and the length of the DNA sequence
my_dna = "ACTGATCGATTACGTATAGTATTTGCTATCATACATATATATCGATGCGTTCAT"
at_content = (my_dna.count('A') + my_dna.count('T')) / len(my_dna)
print("A+T content is",str(int(100*at_content)),"percent")



#2. Complementing DNA
# Need a sequential replacement approach as we haven't done loops or regular expressions yet!
# Could use the complement lower case, or numbers, or anything unique


my_dna = "ACTGATCGATTACGTATAGTATTTGCTATCATACATATATATCGATGCGTTCAT"
replacement1 = my_dna.replace('A', 't')
replacement2 = replacement1.replace('T', 'a')
replacement3 = replacement2.replace('C', 'g')
replacement4 = replacement3.replace('G', 'c')
print("The complement sequence is " + str(replacement4.upper()))



# Program to take a DNA sequence and complement it
# Written by s123456 on 20 Oct 2023
#--------------------------------------------#
# Input
my_dna = "ACTGATCGATTACGTATAGTATTTGCTATCATACATATATATCGATGCGTTCAT"
# Process AND output
print("The complement sequence is", my_dna.replace('A', 't').replace('T','a').replace('G','c').replace('C','g').upper())


#3. Restriction fragment lengths


my_dna = "ACTGATCGATTACGTATAGTAGAATTCTATCATACATATATATCGATGCGTTCAT"
frag1_length = my_dna.find("GAATTC") + 1
frag2_length = len(my_dna) - frag1_length
print("Length of fragment one is " + str(frag1_length))

print("Length of fragment two is " + str(frag2_length))

my_dna = "ACTGATCGATTACGTATAGTAGAATTCTATCATACATATATATCGATGCGTTCAT"
site = "GAATTC"
print("Lengths of",site,"cleaved fragments are",my_dna.find(site) + 1,"and",len(my_dna) -(my_dna.find(site)+1),"bases")


#4. Splicing out introns

# First, store the DNA and extract the exons as separate variables:
my_dna ="ATCGATCGATCGATCGACTGACTAGTCATAGCTATGCATGTAGCTACTCGATCGATCGATCGATCGATCGATCGATCGATCGATCATGCTATCATCGATCGATATCGATGCATCGACTACTAT"

# Remember that we start counting from zero, and that positions are inclusive at the start and exclusive at the end
exon1 = my_dna[0:63]
exon2 = my_dna[90:]


# Now we can just join the two exons and print them out
print("### Exons joined ###\n" + exon1 + exon2)

# To calculate the coding percentage, we just
# have to take the length of the exons,
# divide by the length of the sequence, and multiply by 100
coding_length = len(exon1 + exon2)
total_length = len(my_dna)
print("### Coding percentage (rounded) ###\n" + str(int((coding_length / total_length) * 100)))


# To print out the upper/lower case version,
# we need to take the middle bit
# i.e. the intron and convert it to lower case,
# then join the bits back together
# We already specified the exons
exon1 = my_dna[0:63]
exon2 = my_dna[90:]
intron = my_dna[63:90]
print("### Exon intron exon ###\n" + exon1 + intron.lower() + exon2)



























