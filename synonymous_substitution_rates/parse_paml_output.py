import sys
import os

infile = open(sys.argv[1], 'r')
#outfile = open(sys.argv[2], 'w')
# built-in outfile not being used because when launched from a shell script (with loop) the output is overwritten. instead, I opted to use print() so that output will go to the standard output of the shell script. I decided to name that file all_yn00_data.txt which you can see in the shell script that pairs with this python script.

pamlOutput = infile.readlines()

count = 0
for line in pamlOutput:
        line = line.rstrip('\r\n')
        count += 1
        if count == 80:
                print(line)
                break
#outfile.close()
