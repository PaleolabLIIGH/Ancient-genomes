
Script for genotyping human aDNA using a processed BAM file and the MEGA SNP panel. The script generates a missing data report and a SNP count from the aDNA genotypes. It then integrates the ancient dataset with a human reference panel from the 1000 Genomes Project and conducts a Principal Component Analysis (PCA).
```
#!/bin/bash

# Job settings for SGE
#$ -N genotyping
#$ -cwd
#$ -j y
#$ -S /bin/bash
#$ -m e
#-l vf=10G

# Load required modules
module load gcc/5.1.0
module load htslib/1.2.1
module load plink/1.07
module load samtools/1.2
module load eigensoft/6.0.1.1

# Source modules environment if needed
. /etc/profile.d/modules.sh

# Enable extended globbing
shopt -s extglob

# File paths
BAM_FOLDER="/mnt/Cromosoma/mavila/xroca/data/HSJN_CON/trimmed" #CHANGE accordingly

REFERENCE="/panfs/pan5/mavila/vvilla/genome_ref/human/GRCh38.p13/GRCh38.p13_chr-only/GRCh38.p13.rCRS.fasta"
POS_FILE="/panfs/pan5/mavila/vvilla/Reference_panels/1000GP_MEGAsites/hg38/1000GP.Phase3.v5a.MEGA.sites.10_pops.filt_hg38_sites_flipped_NO-alt_geno.pos"
MAP_FILE="/panfs/pan5/mavila/vvilla/Reference_panels/1000GP_MEGAsites/hg38/1000GP.Phase3.v5a.MEGA.sites.10_pops.filt_hg38_sites_flipped_NO-alt_geno.rand.map"
PED_FILE="/panfs/pan5/mavila/vvilla/Reference_panels/1000GP_MEGAsites/hg38/1000GP.Phase3.v5a.MEGA.sites.10_pops.filt_hg38_sites_flipped_NO-alt_geno.rand.ped"
REFERENCE_PLINK="/panfs/pan5/mavila/vvilla/Reference_panels/1000GP_MEGAsites/hg38/1000GP.Phase3.v5a.MEGA.sites.10_pops.filt_hg38_sites_flipped_NO-alt_geno"
REFERENCE_NAME="1000GP.Phase3.v5a.MEGA.sites.10_pops.filt_hg38_sites_flipped_NO-alt_geno"

# Temporary folder for storing individual genotyping results
TEMP_DIR="./temp_genotypes"
mkdir -p "$TEMP_DIR"

# Iterate over all BAM files in the folder
for file in "$BAM_FOLDER"/*.bam; do
    # Get the base name of the BAM file
    base=$(basename "$file" .bam)
    # Get the base name of the individual
    base_short=$(basename "$file" _L001_R1-2_001.fastq.collapsed.bz2.realign.MD.trimmed.sorted.bam) #CHANGE accordingly

    # Index the BAM file (if not already indexed)
    if [ ! -f "$file.bai" ]; then
        samtools index "$file"
    fi

    # Generate mpileup and convert to TPED format
    samtools mpileup -f "$REFERENCE" -l "$POS_FILE" "$file" | \
    /mnt/Cromosoma/mavila/vvilla/scripts/script_PCA/get_random_allele_30.pl > "$TEMP_DIR/$base_short.Q30.rand.mpileup"

    /mnt/Cromosoma/mavila/vvilla/scripts/script_PCA/overlap_mpileup_tped.pl "$MAP_FILE" "$TEMP_DIR/$base_short.Q30.rand.mpileup" > "$TEMP_DIR/$base_short.tped"
    /mnt/Cromosoma/mavila/vvilla/scripts/script_PCA/overlap_mpileup_map.pl "$MAP_FILE" "$TEMP_DIR/$base_short.Q30.rand.mpileup" > "$TEMP_DIR/$base_short.map"

    echo "$base_short $base_short 0 0 3 1" > "$TEMP_DIR/$base_short.tfam"

    #Extract SNP IDs from the TPED file and store them in a file
    cut -f2 "$TEMP_DIR/$base_short.tped" >  "$TEMP_DIR/$base_short.tped.ids"

    #Use PLINK to extract the same SNPs from the ancient sample in the reference panel and create binary PLINK files (.bed, .bim, .fam).
    plink --bfile $REFERENCE_PLINK --extract "$TEMP_DIR/$base_short.tped.ids" --make-bed --out "$TEMP_DIR/$REFERENCE_NAME.$base_short.preFilt_sites" --allow-no-sex

    #Join the pre-filtered reference file with the TPED sample data, extract alleles that match between the reference panel and the ancient sample
    paste "$TEMP_DIR/$REFERENCE_NAME.$base_short.preFilt_sites.bim" "$TEMP_DIR/$base_short.tped" | \
    nice awk '{if ((($5 ~ /[A]/ || $6 ~ /[A]/) && $11 ~ /[A]/) || ( ($5 ~ /[C]/ || $6 ~ /[C]/) && $11 ~ /C/) || (($5 ~ /[G]/ || $6 ~ /[G]/) && $11 ~ /G/ ) || (($5 ~ /[T]/ || $6 ~ /[T]/) && $11 ~ /[T]/)) print $1,$2,$3,$4}' > "$TEMP_DIR/$base_short.filt_sites.map"

    #Get a list of overlapping SNP IDs in the reference panel and in the ancient sample that passed the allele consistency check.
    cut -f2 -d ' ' "$TEMP_DIR/$base_short.filt_sites.map" >  "$TEMP_DIR/$base_short.filt_sites.map.ids"

    #Filter the ancient sample to Keep Only the Consistent SNPs
    plink --tfile "$TEMP_DIR/$base_short" --extract "$TEMP_DIR/$base_short.filt_sites.map.ids" --recode --out "$TEMP_DIR/$base_short.filtered" --allow-no-sex

done

#Create a folder for intermediate files
mkdir "$TEMP_DIR/previous"
mv $TEMP_DIR/*.Q30.rand.mpileup $TEMP_DIR/*.tped $TEMP_DIR/*.tfam $TEMP_DIR/*ids $TEMP_DIR/*preFilt_sites* $TEMP_DIR/*filt_sites.map "$TEMP_DIR/previous"
cd $TEMP_DIR
mv !(*.filtered).map previous

#Create the merging file
ls *.ped | sed 's/.ped//' | while read f; do echo "$f.ped $f.map"; done > merging_list.txt

#Merge ancient individuals
plink --merge-list merging_list.txt --recode --out HSJN_CON_genotypes --allow-no-sex

#Get a missing report of the ancient dataset
plink --file HSJN_CON_genotypes --missing --out HSJN_CON_genotypes_missing_report

#Get a SNP count of the ancient dataset
awk 'NR > 1 {print $1, $5 - $4}' HSJN_CON_genotypes_missing_report.imiss > HSJN_CON_genotypes_snp_count.txt

mkdir filtered
mv *.filtered.* filtered/

#Sym link the reference
ln -s $PED_FILE .
ln -s $MAP_FILE .

#Create second merging file (better to do it separetly)
echo "HSJN_CON_genotypes.ped HSJN_CON_genotypes.map" >> merging_list2.txt
echo "$REFERENCE_NAME.rand.ped $REFERENCE_NAME.rand.map" >> merging_list2.txt

#Merge ancient dataset with reference
plink --merge-list merging_list2.txt --recode --out $REFERENCE_NAME.HSJN_CON_genotypes --allow-no-sex

#Change phenotype missigness value
sed -i 's/-9/1/g' $REFERENCE_NAME.HSJN_CON_genotypes.ped

#To project with lsqproject, we need to add a column to the merged PED file indicating whether the samples are ancient or reference.
awk '$5 = $5 FS "REFERENCE"' $REFERENCE_NAME.HSJN_CON_genotypes.ped > $REFERENCE_NAME.HSJN_CON_genotypes.pop.ped

#Sanity check
awk '{if($6 == "REFERENCE") print $1}' $REFERENCE_NAME.HSJN_CON_genotypes.pop.ped

#Create a list with the ancient samples IDs
awk '!array[$1]++' HSJN_CON_genotypes.ped | cut -f1-10 -d' ' | awk '{print $1}' > ancient_list.txt

# Loop through each line in ancient_list.txt
while IFS= read -r label; do
    # Replace 'REFERENCE' with 'ANCIENT'
    sed -i "/$label/s/REFERENCE/ANCIENT/g" $REFERENCE_NAME.HSJN_CON_genotypes.pop.ped
done < ancient_list.txt

#Sanity check
awk '{if($6 == "ANCIENT") print $1}' $REFERENCE_NAME.HSJN_CON_genotypes.pop.ped

#Create poplist file
echo "REFERENCE" > $REFERENCE_NAME.poplist.txt

#Delete the 7th Column
awk '!($7="")' $REFERENCE_NAME.HSJN_CON_genotypes.pop.ped  > $REFERENCE_NAME.HSJN_CON_genotypes.pop2.ped

#Run PCA
echo "genotypename:	$REFERENCE_NAME.HSJN_CON_genotypes.pop2.ped
snpname:	$REFERENCE_NAME.HSJN_CON_genotypes.map
indivname:	$REFERENCE_NAME.HSJN_CON_genotypes.pop2.ped
evecoutname:	$REFERENCE_NAME.HSJN_CON_genotypes.pca.evec.txt
evaloutname:	$REFERENCE_NAME.HSJN_CON_genotypes.pca.eval.txt
poplistname:	$REFERENCE_NAME.poplist.txt
lsqproject:	YES" > $REFERENCE_NAME.HSJN_CON_genotypes.pca.par.txt

smartpca -p $REFERENCE_NAME.HSJN_CON_genotypes.pca.par.txt > $REFERENCE_NAME.HSJN_CON_genotypes.lsqproject.out
```
