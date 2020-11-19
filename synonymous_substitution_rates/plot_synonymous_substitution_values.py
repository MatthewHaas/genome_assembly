import sys
import os
import re
import numpy as np
from matplotlib import pyplot as plt


infile = open(sys.argv[1], 'r')
#outfile = open(sys.argv[2], 'w')

data = infile.readlines()

synonSubstValues = list()
for line in data:
        line = str(line)
        dS_value = re.findall(r'(0.[0-9]+)\)$', line)
        try:
                dS_value = float(dS_value[0])
        except IndexError:
                continue
        synonSubstValues.append(dS_value)
print(synonSubstValues)

binValues = np.arange(0, 0.6, 0.01).tolist()

plt.title('Distribution of synonymous substitution rates')
plt.xlabel('synonymous substitution rate')
plt.ylabel('frequency')
plt.hist(synonSubstValues, bins = binValues)
plt.savefig('synonymous_substitution_value_distribution.png')
