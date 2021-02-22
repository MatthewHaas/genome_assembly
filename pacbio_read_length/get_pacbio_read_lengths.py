import sys
import os
import matplotlib.pyplot as plt

infile = open(sys.argv[1], 'r')

data = infile.readlines()

read_lengths = list()

for line in data:
        line = line.split('_')
        length = int(line[3])
        read_lengths.append(length)

        # Find the sum/total of all the read lengths from PacBio sequencing
# Ultimate goal is to find the N50.
total = 0
for element in range(0, len(read_lengths)):
	total = total + read_lengths[element]
print("Sum of all PacBio read lengths:", total)

# Sort PacBio read lengths in descending order
# reverse = True puts list in descending order
sortedReadLengths = sorted(read_lengths, reverse = True)

# N50 is defined as the smallest length such that 50% of the assembly is contained in contigs or scaffolds equal to or larger than this number
# In this case, these are read lengths, not contigs or scaffolds.
runningTotal = 0
for element in range (0, len(sortedReadLengths)):
	runningTotal = runningTotal + sortedReadLengths[element]
	if runningTotal < (total/2):
		continue
	elif runningTotal >= (total/2):
		print("The N50 is: ", sortedReadLengths[element])
		break
	else:
		print("Something is wrong.")

plt.hist(read_lengths, color = 'blue', edgecolor = 'black', bins = int(500))
plt.xlabel('PacBio Read Length (bp)')
plt.ylabel('Frequency')
plt.savefig('pacbio_length_distr.png')
