# README for alignment_scripts

These scripts were written to align orthologous shattering genes between NWR and _O. sativa_. Initially, the results were going to be included in the Supporting Information, but were ultimately cut from the final version of the manuscript.

## Workflow
1. [setup_orthogroup_directories.sh](setup_orthogroup_directories.sh) was used to set up the directory structure.
2. [run_filter_fasta.sh](run_filter_fasta.sh) was used to filter the orthogroups using the Python script [filter_fasta.py](filter_fasta.py) so that it would only contain sequences from NWR (_Z. palustris_), _Z. latifolia_, and _O. sativa_.
3. [run_clustal_omega.sh](run_clustal_omega.sh) was used to align the orthogroup sequences.
4. [run_backtranslate.sh](run_backtranslate.sh) uses Tom Kono's [Backtranslate_Orthogroup_TK.py](Backtranslate_Orthogroup_TK.py) script to convert the protein sequence used in the orthogroup files back to DNA sequence (using the CDS FASTA files)
5. [convert_fasta_to_phlyip.sh](convert_fasta_to_phlyip.sh) uses the Python script [convert_fasta_to_phylip.py](convert_fasta_to_phylip.py) to convert the FASTA sequence that results from the above scripts to PHYLIP format for use in PAML.
