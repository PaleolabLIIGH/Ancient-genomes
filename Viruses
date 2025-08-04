# Alignment with MALT
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
