# README for alignment_scripts

These scripts were written to align orthologous shattering genes between NWR and _O. sativa_. Initially, the results were going to be included in the Supporting Information, but were ultimately cut from the final version of the manuscript.

## Workflow
1. [setup_orthogroup_directories.sh](alignment_scripts/setup_orthogroup_directories.sh) was used to set up the directory structure.
2. [run_filter_fasta.sh](run_filter_fasta.sh) was used to filter the orthogroups using the Python script [filter_fasta.py](filter_fasta.py) so that it would only contain sequences from NWR (_Z. palustris_), _Z. latifolia_, and _O. sativa_.
3. [run_clustal_omega.sh](run_clustal_omega.sh) was used to align the orthogroup sequences.
