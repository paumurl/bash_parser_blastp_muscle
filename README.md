# bash_parser_blastp_muscle
A bash parser that given a  protein fasta/multifasta and a two column text file (first: species identifier, second: url link to the subject proteome), will perform a blastp analysis for the query against the proteomes, and will create a MUSCLE allignment and neighbor tree for each query protein.

Expected input files can be find in example_files/ for clarification.


### ejercicio_bash2.sh
This .sh script requires 5 arguments:
- query_sequences.fa: the protein multifasta
- ncbi_urls_file: the url links to the ncbi proteomes of each subject
- output_folder: the desired name for the folder in which the project will be stored
- blast-identity: sequence identity cut off (value 0-100), 70 recommended
- blast-coverage: sequence coverage cut off (value 0-100), 40 recommended

Upon script execution there will be a folder with the name output_folder, inside which you can find two more folders (and the log file): one folder for data (in which all the data used for the analysis can be found) and one for the results. Inside the results folder you will get the blastp results and a folder for each query protein in your query_sequences.fa, and each subfolder will contain the allignment, the neighbor joining tree and the fasta blastp hits.
The terminal will print a summary of the blastp hits: the total hits and the hits for each protein in the query_sequences.fa.

Expected files for output_folder/data are as following:
- the unzipped downloaded proteomes
- Identifiers.tsv: a tsv two column text file containing the species identifiers and their proteome directory
- proteome.fa: the proteome database of all subject proteomes
- query.fa: the given query_sequences.fa
- species_fasta_ID.tsv: a tsv two column text file containing the species identifiers and their fasta header


Expected files for output_folder/results are as following:
- a folder for each protein query: which will contain the alignment file (.aln), the tree file (.nx) and the fasta blastp hits for said protein.
- blast_result: the untouched blastp results
- blast_result_filtered: upon awk filter of the blast-identity and blast_coverage values
- blast_result final: the blastp_result_filtered after replacing the fasta header of the proteome hit with its species identifier instead (using species_fasta_ID.tsv)
- proteinquery: a text file with the names of each protein in query_sequences.fa
