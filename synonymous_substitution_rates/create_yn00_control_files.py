import os
import sys
import re

infile = open(sys.argv[1], 'r') # this will always be the template file (yn00_template.ctl)
outfile = open(sys.argv[2], 'w') # this will always be in the form of orthogroup/yn00_orthogroup_control_file

orthogroupName = str(sys.argv[3]) # this is to get the orthogroup name from the command-line input for use in this script

controlFile = infile.readlines()

for line in controlFile:
        if 'alignment_file_in_phylip_format' in line:
                line = re.sub('alignment_file_in_phylip_format', str(orthogroupName + '_slim_backtranslated_with_space.phy'), line)
                outfile.write(line)
        elif 'yn00_output_file.txt' in line:
                line = re.sub('yn00_output_file.txt', str('yn00_' + orthogroupName + '_output_file.txt'), line)
                outfile.write(line)
        else:
                outfile.write(line)
outfile.close()
