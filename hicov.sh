#!/bin/bash

# HiCov

# This script is a workhorse of hCoV project born out of collaboration between IMI MF, University of Ljubljana
# and Hiscox Lab, University of Liverpool.

# This procedure performs read trimming,  mapping on the provided reference genome, primer trimming,
# calculation of the mapping statistics, consensus calling and variant analysis.

# Author: Alen Suljic
# Inspiration: Martin Bosilj, Hannah Goldswain
# Maiden voyage: 10.07.2024

# Prerequisites:
######################################################################################################################################
# 1. hicov.sh script
# 2. Reference directory (reference.fasta, reference.gff3, reference.primer.bed)
# 3. Singulaity container hicov.sif
######################################################################################################################################

# Usage:
# 1. Copy the script to the working directory (this will be the output directory!). 
# 2. Reference directory must contain  "reference.fasta", "reference.gff3" and "reference.primer.bed" files!
# 3. The script is run inside dedicated singularity container "hicov.sif".
# 4. If required, adjust pipeline parameters listed below to your requirements.
# 5. Script takes reference name as positional argument 1 and input data (fastq.gz) directory as positional argument 2.

###################################################################################################
# Command example:
# bash hicov.sh hcov_229e /home/imi.mf/asuljic/experiments/seasonal_coronas/run_id201384/data/e229
###################################################################################################
# Parameters:
######################################
## thread count
thr=24
## minimum phread quality score
qqp=20
## minimum read length after trimming
lr=30
## minimum read depth for consensus
dc=10
######################################

# Routine:

# initiate log file
exec > >(tee -a "experiment.log") 2>&1

# Create sample list
ls "${2}" | cut -d "_" -f 1 | sort -V | uniq > samples

# Trim reads with fastp
echo "Trimming reads"

for i in $(cat samples); do

	sample="${i}"
	echo $i;
	mkdir -p trimmed qc

	fastp \
	--in1 "${2}"/${sample}_R1.fastq.gz \
	--in2 "${2}"/${sample}_R2.fastq.gz \
	--out1 trimmed/${sample}_trim_R1.fastq.gz \
	--out2 trimmed/${sample}_trim_R2.fastq.gz \
	--cut_front \
	--cut_tail \
	--cut_window_size 4 \
	--cut_mean_quality $qqp \
	--qualified_quality_phred $qqp \
	--length_required $lr \
	--correction \
	--detect_adapter_for_pe \
	--trim_poly_x \
	--trim_poly_g \
	--html qc/${sample}.fastp.html \
	--json qc/${sample}.fastp.json \
	--thread $thr \

done

# Index reference
echo "Indexing reference sequence"

bwa index reference/"${1}".fasta
samtools faidx reference/"${1}".fasta
samtools dict reference/"${1}".fasta > reference/"${1}".dict

# Alignment
echo "Mapping reads to reference"

for i in $(cat samples); do

	sample="${i}"
	echo $i;
	mkdir -p mappings

## Align PE reads and trim primer sequences
	bwa mem -t $thr -k 17 -R "@RG\tID:1\tSM:${sample}\tPL:ILLUMINA" reference/"${1}".fasta trimmed/${sample}_trim_R1.fastq.gz trimmed/${sample}_trim_R2.fastq.gz | \
		samtools view -uhS -F4 -@$thr - | \
		samtools sort -@$thr -n - | \
		ivar trim -b reference/"${1}".primer.bed -m 0 -q 0 -e | \
		samtools sort -@$thr -n - | \
		samtools fixmate -@$thr -m - mappings/${sample}_PE_fixmate.bam
	samtools sort -@$thr mappings/${sample}_PE_fixmate.bam | \
		samtools markdup -rsS -@$thr - mappings/${sample}_PE.bam
	samtools index mappings/${sample}_PE.bam

## Recalculate mismatches and deletions between aligned reads and reference sequence

	samtools calmd -Ar -u -@$thr mappings/${sample}_PE.bam reference/"${1}".fasta | \
        samtools sort -@$thr -o mappings/${sample}_final.bam 
	samtools index mappings/${sample}_final.bam 

	find mappings -type f -not -name "*final*" -delete

done

# Mapping statistics

echo "Calculating mapping statistics"

for i in $(cat samples); do

	sample="${i}"
	echo $i;
	mkdir -p stats

	echo -e "${sample}" > stats/${sample}_mapstats_name.log
	samtools flagstat mappings/${sample}_final.bam > stats/${sample}_allstats.log
	samtools coverage -q $qqp -Q $qqp mappings/${sample}_final.bam > stats/${sample}_coverage.log
	samtools depth -aa -H -q $qqp -Q $qqp mappings/${sample}_final.bam -o stats/${sample}.covdepth
	
	cat stats/${sample}_mapstats_name.log stats/${sample}_allstats.log stats/${sample}_coverage.log > stats/${sample}_stats.log
	rm  stats/${sample}_allstats.log stats/${sample}_mapstats_name.log stats/${sample}_coverage.log

done

## Generate mapping stats report

echo "Generating mapping statistics report"

mkdir -p results

echo -e "sample\trname\tstartpos\tendpos\tnumreads\tcovbases\tcoverage\tmeandepth\tmeanbaseq\tmeanmapq\tr1_nreads\tr2_nreads" > results/mapstats.tsv

for i in $(cat samples); do

  sample="${i}"

  first_line=$(head -n 1 stats/${sample}_stats.log)
  
  nineteenth_line=$(sed -n '19p' stats/${sample}_stats.log)

  r1_nreads=$(gunzip -c trimmed/${sample}_trim_R1.fastq.gz | awk '{s++}END{print s/4}')
  r2_nreads=$(gunzip -c trimmed/${sample}_trim_R2.fastq.gz | awk '{s++}END{print s/4}')
  
  echo -e "$first_line\t$nineteenth_line\t$r1_nreads\t$r2_nreads" >> results/mapstats.tsv
  
done

# Generate Consensus sequence
## threshold for consensus 0.5, depth 10
echo "Generating consensus sequence"

for i in $(cat samples); do

	sample="${i}"
	echo $i;
	mkdir -p consensus/ results

	samtools mpileup -aa -A -d 0 -E -f reference/"${1}".fasta -q $qqp -Q $qqp mappings/${sample}_final.bam | \
	ivar consensus -p consensus/${sample} -t 0.5 -q $qqp -m $dc

	seqtk seq consensus/${sample}.fa | \
		tr "?RYSWKMBDHVN.ryswkmbdhvn" "N" | \
		tr "-" "N" | \
		tr [:lower:] [:upper:] | grep -iv ">" - > consensus/tmp.seq

	echo ">"$i > consensus/${sample}.fa
	cat consensus/${sample}.fa consensus/tmp.seq  > consensus/${sample}.fasta
	find ./consensus -type f ! -name '*.fasta' -delete

done

cat consensus/*.fasta > results/consensus_sequences.fasta

# Variant calling 

echo "Variant analysis"

for i in $(cat samples); do

	sample="${i}"
	echo $i;
	mkdir -p variants

	samtools mpileup -aa -A -d 0 -E -f reference/"${1}".fasta -q 20 -Q 20 mappings/${sample}_final.bam | \
	ivar variants -p variants/${sample} -t 0.01 -q $qqp -m 10 -r reference/"${1}".fasta -g reference/"${1}".gff3

done

# Data transformtion
## Add sample name as first column

echo "Results consolidation and beautification"

for file in variants/*.tsv; do 

    filename_noext=$(basename -- "$file")
    filename_noext="${filename_noext%.*}"

    awk -v filename="$filename_noext" 'BEGIN{OFS="\t"} {$1=filename; print}' "$file" > "${file}.tmp"
    mv "${file}.tmp" "$file"

done


for file in stats/*.covdepth; do 

    filename_noext=$(basename -- "$file")
    filename_noext="${filename_noext%.*}"

    awk -v filename="$filename_noext" 'BEGIN{OFS="\t"} {$1=filename; print}' "$file" > "${file}.tmp"
    mv "${file}.tmp" "$file"

done

# Merge variant annotation and coverage files to single dataframe

awk 'FNR==1 && NR!=1{next;}{print}' variants/*.tsv > results/sleek_variants.tsv

awk 'FNR==1 && NR!=1{next;}{print}' stats/*.covdepth > results/coverage.tsv

tr '\t' ',' < results/coverage.tsv > results/coverage.csv

rm results/coverage.tsv

# Cleanup
mkdir -p logs
mv experiment.log logs/
mv samples logs/


echo "Analysis complete"
