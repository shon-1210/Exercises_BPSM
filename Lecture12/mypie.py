#!/usr/local/bin/python3


import subprocess 

# programme to trim fixed length adapter sequences from start of each sequence
# and keep the output in a file, outputting to screen as we go
# Code written by s123456 on 27 October 2023
#......................................................
# Open and read the file, trimming off adapter "on the fly"
input_txt_contents = open('input.txt').read().upper().replace('ATTCGATTATAAGC','').split()

# Alternative way to do it
input_txt_contents2 = open('input.txt').read().upper().split()
len(input_txt_contents)
len(input_txt_contents2)



# Contents to a single file and screen
# Open the output pipe for writing, I've called it cleanseqs
cleanseqs = open('Cleaned_sequences.txt','w')


# Now loop through and do the outputs required
# We weren't specifically asked for fasta format...!
# Version 1 (adapter gone already)
for cleanones in input_txt_contents :
    cleanseqs.write(cleanones + '\n')
    cleanones


# Version 2: adapter bit not written to file
# Open an output pipe for writing, I've called it cleanseqs2


cleanseqs2 = open('Cleaned_sequences2.txt','w')
for cleanones in input_txt_contents2 :
    cleanseqs2.write(cleanones[14 :] + '\n')
    cleanones[14 :]


# "Better version" using with
with open('Cleaned_sequences3.txt','w') as cleanseqs_w :
 for cleanones in input_txt_contents :
    cleanseqs_w.write(cleanones + '\n')
    cleanones


cleanseqs.closed
cleanseqs2.closed
cleanseqs_w.closed

cleanseqs.close()
cleanseqs2.close()


print(open('Cleaned_sequences.txt').read().rstrip())


import subprocess
subprocess.call("cat Cleaned_sequences.txt",shell=True)



##### Code for TASK 2 ####
##### Code for TASK 2 ####


# Do a couple of things in Linux first
# Check we have what we need : genomic file
wc genomic_dna2.txt
head genomic_dna2.txt

wc exons.txt
head exons.txt


# Open and read the genomic file as a string
genomic_dna2 = open('genomic_dna2.txt').read().upper()
len(genomic_dna2)


# Version 1 : open and read the exons file as a list directly
exons_v1 = open('exons.txt').read().replace(',','\n').split()
exons_v1

# Could process these taking the first/second, third/fourth etc
# using a list(range(1,10,2)) command or similar
# This approach could rapidly get too complicated, think of a better way!



# Version 2 : open and read the exons file as an rstrip() file,
# this comes in as a string
exons_v2 = open('exons.txt').read().rstrip()
exons_v2

# Quick check
counter=0
for exons in exons_v2 :
    counter += 1
    print(counter,exons)


type(exons_v2)


# Version 3 : open and read the exons file as an rsplit() input,
# which comes in as a list, one line per string element
exons_v3 = open('exons.txt').read().rsplit()
type(exons_v3)

exons_v3

# Quick check
counter=0
for exons in exons_v3 :
    counter += 1
    print(counter,exons)


# This looks correct, let's use this version
counter=0
# Open pipe to file for writing output
with open('Exercise2_coding_seqence.fasta','w') as fullcoding :
    # Fasta header line insertion : write to the pipe
    #Then extract the coding segments, need Python3 index positions
    fullcoding.write('>Lecture12_exercise2_codingseq\n')
    for exons in exons_v3 :
        counter += 1
        startexon = int(exons.split(',')[0])-1
        endexon = int(exons.split(',')[1])
        exonwanted = genomic_dna2[startexon :endexon]
        fullcoding.write(exonwanted)
        regionsummary = 'Exon'+str(counter)+' '+str(exons)+' runs from index position '+str(startexon)+' upto but not including '+str(endexon)+ '.'
        print(regionsummary,'\n\t',exonwanted)



fullcoding.closed

print(open('Exercise2_coding_seqence.fasta').read())





##### Code for TASK 3 ####

import os, shutil
shutil.copy("../Lecture11/remote_exon01.fasta", ".")

sorted(os.listdir())

# What is in the file?
open('remote_exon01.fasta').read().split()

# Read in the sequence as a string, ignoring the first list element as it is the fasta header line
AJ223353_seq = open('remote_exon01.fasta').read().split()[1]
len(AJ223353_seq)

# Check the end of the sequence : should be a stop codon
AJ223353_seq[-3:]


# Double-check we know how indexing works
len(AJ223353_seq[0:30])


# Set the windowing parameters required
windowsize = 30
offset = 3


# Generate a numeric range for positions for the window starts
segment_starts = list(range(0,len(AJ223353_seq)-windowsize+1,offset))


# Alternative (better) range to also get the smaller segments at the end
segment_starts = list(range(0,len(AJ223353_seq),offset))
segment_starts

# Create and then close a file object to take all the segment sequences
basic_fasta_header = '_window' + str(windowsize) + '_offset' + str(offset)




alloutfilename = 'AJ223353_coding_all' + basic_fasta_header + '.fasta'
with open(alloutfilename, 'w') as allsegments :
 allsegments.write('')
# Create a blank list to hold the segments
 segments_made = []
# Create a possible yes/no switching variable, i.e. a conditional
 fileswanted = 1
# The for loop
 for seg_bits in segment_starts :
     tempseq = AJ223353_seq[seg_bits :seg_bits+windowsize].upper()
     segments_made = segments_made + [tempseq]
# percentage GC content, convert float to integer value
     tempGC = int(100 * (tempseq.count('G') + tempseq.count('C'))/len(tempseq) )
# make a fasta header line
     descriptionline = 'AJ223353_coding_'+str(len(segments_made))+'_'+tempseq+basic_fasta_header
     fastaheader = '>'+descriptionline
     print(len(segments_made),'\t',tempseq,'\t',tempGC)
# Question : do we want files written? Answer will be either True or False ('maybe' is NOT an option!)
     if fileswanted == 1 :
# open the segment fastafiles
        with open(descriptionline+'.fasta', 'w') as segmentfile :
# output to files
          segmentfile.write(fastaheader+'\n'+tempseq)
          allsegments.write(fastaheader+'\n'+tempseq+'\n\n')
     else :
        allsegments.close()
       


segmentfile.closed
allsegments.closed




# Outside the loop, print the first 10 segments, one per line
print('\n'.join(segments_made[0:10]))

# Print the first 3 segments, separated by the word 'then'
print(' then '.join(segments_made[0:3]))



# Does the combined fasta file output look correct?
# Using a loop makes this easy
for these in 1,2,3 :
    open(alloutfilename).read().split('>')[these]


# Better still, take advantage of the \n\n between each sequence in the combined fasta file...
print('\n'.join(open(alloutfilename).read().split('\n\n')[0:3]))


# Or just use the file object, this will process the whole file
my_file = open(alloutfilename)
for line in my_file :
    print(line.rstrip('\n\n'))




# What files do I have in my directory?
dir()

# NO, that shows the things I have in my Python session/"bubble", not what I have saved!
# We need to be able to 'talk' to the operating system (os), luckily we have that os module to do that

# Have I already imported the os module?
'os' in dir()


os.getcwd()



# Now we can get current working directory again
os.getcwd()



# How many files/directories are there, visible from where we are?
len(os.listdir())


# Show me the files
os.listdir()


os.listdir()[0]
os.listdir()[1]

# Show me the "dna" files
os.listdir('*.dna')




# Need to import another module if we want to use the * wildcard...
import glob

# What functionality does this module have?
# DOT TAB TAB


glob.glob('*.dna')
len(glob.glob('*.dna'))



glob.glob('*.txt')





##### Code for TASK 4 ####


# In Linux shell
pwd

# What files are there that we might need?
ls *.dna

# How many sequences are there in each/all the files?
wc -l *.dna

# In Linux, make a sub-directory for the sequence files
mkdir dna_files

# Move the dna files to the dna_files directory that we have just made
mv *.dna ./dna_files

# Have a quick look at one file to see what the format is
head -n5 ./dna_files/xaa.dna

# Start Python
python3
# Import the various modules we think we might need
import os, sys, subprocess, shutil



# Get a list of the files we might want to process
# We need the path : directory and the file name.
# This loop would open each file
for file_name in os.listdir('dna_files/') :
    if file_name.endswith('.dna') :
        f'Reading sequences from {file_name}'
        dna_file = open('dna_files/' + file_name)
# This loop then looks at each line and gets the length
# Note the indentation...
        for line in dna_file :
            dna = line.rstrip('\n')
            length = len(dna)
            print(f'\tFound a DNA sequence with length {str(length)}' )


# Now add some code that will try each size bin in turn to see if the sequence fits.
# We'll use a range to check the bounds on each size bin

# Go through each file in the directory
for file_name in  sorted(os.listdir('dna_files')) :
# Check if the file name ends with .dna
  if file_name.endswith('.dna') :
    print('Reading sequences from ' + file_name)
# Open the file and process each line
    dna_file = open('dna_files/' + file_name)
# Calculate the sequence length
    for line in dna_file :
      dna = line.rstrip('\n')
      length = len(dna)
      print('\tsequence length is ' + str(length))
# Go through each size bin and check if the sequence belongs in it
      for bin_lower in list(range(100,1000,100)) :
        bin_upper = bin_lower + 99
        if length >= bin_lower and length <= bin_upper :
          print('\t\tbin is ' + str(bin_lower) + ' to ' + str(bin_upper))




# We will need to create some new directories so that we can write
# each DNA sequence to the correct one based on the directory name
for bin_lower in list(range(100,1000,100)) :
    bin_upper = bin_lower + 99
    bin_directory_name = str(bin_lower) + '_' + str(bin_upper)
    os.mkdir(bin_directory_name)



# What directories do we have now? Don't want to use os.listdir() !
len(os.listdir())


# Could do this
for name in os.listdir(".") :
   if os.path.isdir(name) :
     print (name)




# Or use os.walk() to generate the names in a directory tree
mydirs = os.walk('.')
mydirs


mydirs[0]

# Access the contents with a for loop
sorted([x[0] for x in mydirs])


# Can we use the generator object again? No, would need to recreate it.
sorted([x[0] for x in mydirs])


# Let's run our programme now
# Re-create a new directory for each size bin, essentially deleting anything there previously
for bin_lower in list(range(100,1000,100)) : 
    bin_upper = bin_lower + 99 
    bin_directory_name = str(bin_lower) + '_' + str(bin_upper) 
    shutil.rmtree(bin_directory_name)
    os.mkdir(bin_directory_name) 

# create a variable to hold an arbitrary sequence identifier number
seq_number = 1


# Go through each file in the directory
for file_name in os.listdir('dna_files') :
# Check if the file name ends with .dna
  if file_name.endswith('.dna') :
    print('Reading sequences from ' + file_name)
# Open the file and process each line
    dna_file = open('dna_files/' + file_name)
    for line in dna_file :
# Calculate the sequence length
      dna = line.rstrip('\n')
      length = len(dna)
      print('\tsequence length is ' + str(length))
# Go through each size bin and check if the sequence belongs in it
      for bin_lower in list(range(100,1000,100)) :
        bin_upper = bin_lower + 99
        if length >= bin_lower and length <= bin_upper :
          print('\t\tbin is ' + str(bin_lower) + ' to ' + str(bin_upper))
          bin_directory_name = str(bin_lower) + '_' + str(bin_upper)
          output_path = bin_directory_name + '/' + str(seq_number) + '.dna'
          output = open(output_path, 'w')
          output.write(dna)
          output.close()
# Increment the sequence number
          seq_number = seq_number+1


# Finally, see what we been done
os.listdir('dna_files')


# Look in one of the size directories we made
sorted(os.listdir('100_199'))


# Open one of the files in a size directory we made
print(open('100_199/3.dna').read())

# Was it the right size for that directory?
len(open('100_199/3.dna').read())
