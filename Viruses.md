# Alignment with MALT

MALT, an acronym for MEGAN alignment tool, is a sequence alignment and analysis tool designed for processing high-throughput sequencing data, especially in the context of metagenomics



MALT manual https://software-ab.cs.uni-tuebingen.de/download/malt/manual.pdf

Herbig, A., Maixner, F., Bos, K. I., Zink, A., Krause, J., & Huson, D. H. (2016). MALT: Fast alignment and analysis of metagenomic DNA sequence data applied to the Tyrolean Iceman. BioRxiv, 050559.

## 1.1 Convert BZ2 to GZ2

```bash for file in $(find /mnt/Cromosoma/mavila/falvarez/MALT/Viruses -type f -name "*bz2"); do name="${file%.bz2}" bzcat "$file" | gzip > "$name.gz" done ``` </pre>


## 2.2 Create the following script
```bash  

nano MALT.sh

#!/bin/bash
#
# Your job name
#$ -N My_Job
#
# Use current working directory
#$ -cwd
#
# Join stdout and stderr
#$ -j y
#
# Run job through bash shell
#$ -S /bin/bash
#
# Send an email after the job has finished
#$ -m e
#$ -M falvarez.alihuen@outlook.com

# Amount of memory requested
#$ -l vf=10G
#$ -pe openmp 12
#export omp_num_threads=12

#
# If modules are needed, source modules environment (Do not delete the next line):
. /etc/profile.d/modules.sh

module load anaconda3/2021.05
module load malt/0.5.2

# Define the necessary variables
index=/panfs/pan5/mavila/asanchezl/VIRALDATABASE/index_V_DB/
input_dir=your/path
outputdir=your/path

# Create output directories if they do not exist
mkdir -p "${outputdir}/alignments"
mkdir -p "${outputdir}/RMA_files"

# Iterate over each subdirectory in input_dir
for sub_dir in "$input_dir"/*/; do
    echo "Processing directory: $sub_dir"
    # Enter the subdirectory
    cd "$sub_dir" || continue
    #  Check for .unmapped.fastq files
    if ls *.unmapped.fastq.gz 1> /dev/null 2>&1; then
        # Execute the MALT command for each file .unmapped.fastq.gz
        for file in *.unmapped.fastq.gz; do
            malt-run -m BlastN -at SemiGlobal -i "${file}" -d "${index}" -o "${outputdir}/RMA_files" -a "${outputdir}/alignments" -oa "${outputdir}/alignments" --minPercentIdentity 85 -e 0.001
        done
    else
        echo "No .unmapped.fastq files found in $sub_dir"
    fi
    # Back to original directory
    cd - > /dev/null
done

```` ``` ````
#BlastN mode: This specifies the BlastN search mode, which is used to search for similar nucleotide sequences in a database.
#at SemiGlobal: Uses semi-global alignment, which allows a query sequence to be aligned with a database without the need for matching at the beginning or end.
# -minPercentIdentity 95: Sets the minimum identity threshold to 95%. This means that only matches with a similarity of 95% or higher will be considered.
# -index /your/path: Specifies the location of the MALT database index.
# inFile file.fastq.collapsed.unmapped.fastq.gz: Specifies the input file.


## 2.3 Download the 'rma6' output files from MALT and open them in Megan6 in your local machine.

MEGAN6 is an algorithm that classifies sequences taxonomically based on the lowest common ancestor of the taxa with which the sequence has the best alignments. This provides a more accurate and specific classification through multiple comparisons with reference databases.

Set the Megan6 options to change the LCA parameters:
-minimum percent identity = 80

If you find a virus, right-click on it and select 'Show alignment'.
-Files > Export > Each Read Separately.

Observe whether the reads are distributed throughout the genome, and whether they are repeated sequences.

## 2.4 tblastn

Paste each different read into the BLAST search  (https://blast.ncbi.nlm.nih.gov/Blast.cgi?PROGRAM=blastn&PAGE_TYPE=BlastSearch&LINK_LOC=blastsearch) in the 'Enter Query Sequence' box.

Pathogen is confirmed with low e-value and high identity percentage.


## 2.5 Check for repeated readings on next page

https://www.repeatmasker.org/cgi-bin/WEBRepeatMasker


















