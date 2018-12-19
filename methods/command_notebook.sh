#!/bin/bash

######## Create and activate conda environment for EMS analysis
conda create -n ems_env multiqc fastqc trim-galore samtools bedtools vcftools snpeff snpsift vcfkit vt bwa seqtk
source activate ems_env

###### Set a home base and go there.
bieserHome="/data/gpfs/assoc/inbre/richland/projects/projects_bieser/"
mkdir -p $bieserHome
cd $bieserHome

##### Download all 5 samples from basespace
# these bs download commands assume I did a `bs auth` when i installed `bs`

bs -V # should be ≥0.10.5 to be sure for this download syntax works.
# now we download each biosample
bs download biosample -n 1801-1 -o ./ &
bs download biosample -n 1801-2 -o ./ &
bs download biosample -n 1801-3 -o ./ &
bs download biosample -n 1801-4 -o ./ &
bs download biosample -n 1801-5 -o ./ &

###### Structuring the project part 1
mkdir fastq_inparts
mkdir fastq_joined
mkdir docs
mkdir docs/json_illumina
mv *.json docs/json_illumina/

######### Unzipping the lane fastqs for concatenation+rezipping
unpigz -vk 1801-*/*.gz
mv 1801-*/*.fastq fastq_joined/
mv 1801-* fastq_inparts/

ls fastq_joined/


########### Concatenation and re-zipping
cd fastq_joined
ls
for myset in `ls *L001_R1_001.fastq`
  do
  echo pattern based on $myset
  mypattern=${myset//L001_R1_001.fastq}
  echo pattern generalized to $mypattern
  # concatenate R1s
  myoutR1="${mypattern}"R1.fastq
  echo joining "${mypattern}"L00{1..4}_R1_001.fastq to $myoutR1
  cat "${mypattern}"L00{1..4}_R1_001.fastq > $myoutR1
	# concatenate R2s
  myoutR2="${mypattern}"R2.fastq
  echo joining "${mypattern}"L00{1..4}_R2_001.fastq to $myoutR2
  cat "${mypattern}"L00{1..4}_R2_001.fastq > $myoutR2
done

# checking that nothing broke
ls
wc -l *.fastq # should be divisible by 4 for all. the new file should equal the sum of its inputs. R1 should equal R2.
# OUTPUT OF wc -l here
# 19936576 A44_S2_L001_R1_001.fastq
# 19936576 A44_S2_L001_R2_001.fastq
# 18839864 A44_S2_L002_R1_001.fastq
# 18839864 A44_S2_L002_R2_001.fastq
# 20739280 A44_S2_L003_R1_001.fastq
# 20739280 A44_S2_L003_R2_001.fastq
# 19170280 A44_S2_L004_R1_001.fastq
# 19170280 A44_S2_L004_R2_001.fastq
# 78686000 A44_S2_R1.fastq
# 78686000 A44_S2_R2.fastq
# 26632400 Control_S5_L001_R1_001.fastq
# 26632400 Control_S5_L001_R2_001.fastq
# 26364960 Control_S5_L002_R1_001.fastq
# 26364960 Control_S5_L002_R2_001.fastq
# 27872388 Control_S5_L003_R1_001.fastq
# 27872388 Control_S5_L003_R2_001.fastq
# 26054596 Control_S5_L004_R1_001.fastq
# 26054596 Control_S5_L004_R2_001.fastq
# 106924344 Control_S5_R1.fastq
# 106924344 Control_S5_R2.fastq
# 33743500 cos2_S3_L001_R1_001.fastq
# 33743500 cos2_S3_L001_R2_001.fastq
# 33731732 cos2_S3_L002_R1_001.fastq
# 33731732 cos2_S3_L002_R2_001.fastq
# 35150356 cos2_S3_L003_R1_001.fastq
# 35150356 cos2_S3_L003_R2_001.fastq
# 33395140 cos2_S3_L004_R1_001.fastq
# 33395140 cos2_S3_L004_R2_001.fastq
# 136020728 cos2_S3_R1.fastq
# 136020728 cos2_S3_R2.fastq
# 36565920 H22_S4_L001_R1_001.fastq
# 36565920 H22_S4_L001_R2_001.fastq
# 36638576 H22_S4_L002_R1_001.fastq
# 36638576 H22_S4_L002_R2_001.fastq
# 38093588 H22_S4_L003_R1_001.fastq
# 38093588 H22_S4_L003_R2_001.fastq
# 36048784 H22_S4_L004_R1_001.fastq
# 36048784 H22_S4_L004_R2_001.fastq
# 147346868 H22_S4_R1.fastq
# 147346868 H22_S4_R2.fastq
# 35601920 L31_S1_L001_R1_001.fastq
# 35601920 L31_S1_L001_R2_001.fastq
# 35646164 L31_S1_L002_R1_001.fastq
# 35646164 L31_S1_L002_R2_001.fastq
# 37063132 L31_S1_L003_R1_001.fastq
# 37063132 L31_S1_L003_R2_001.fastq
# 35170452 L31_S1_L004_R1_001.fastq
# 35170452 L31_S1_L004_R2_001.fastq
# 143481668 L31_S1_R1.fastq
# 143481668 L31_S1_R2.fastq
# 2449838432 total

# delete the unzipped unjoined intermediates
rm *_001.fastq
ls
pigz -v -p 2 *.fastq
ls -lh

######## Building Drosophila reference
cd $bieserHome
mkdir -p ref_genome/dna
mkdir -p ref_genome/gff
mkdir -p ref_genome/gtf

cd ref_genome/dna
wget   ftp://ftp.flybase.net/genomes/Drosophila_melanogaster/dmel_r6.24_FB2018_05/fasta/dmel-all-chromosome-r6.24.fasta.gz
wget ftp://ftp.flybase.net/genomes/Drosophila_melanogaster/dmel_r6.24_FB2018_05/fasta/md5sum.txt

cd $bieserHome/ref_genome/gff
wget   ftp://ftp.flybase.net/genomes/Drosophila_melanogaster/dmel_r6.24_FB2018_05/gff/dmel-all-r6.24.gff.gz
wget   ftp://ftp.flybase.net/genomes/Drosophila_melanogaster/dmel_r6.24_FB2018_05/gff/dmel-all-no-analysis-r6.24.gff.gz
wget   ftp://ftp.flybase.net/genomes/Drosophila_melanogaster/dmel_r6.24_FB2018_05/gff/dmel-all-filtered-r6.24.gff.gz
wget   ftp://ftp.flybase.net/genomes/Drosophila_melanogaster/dmel_r6.24_FB2018_05/gff/problem_case_filter.pl
wget   ftp://ftp.flybase.net/genomes/Drosophila_melanogaster/dmel_r6.24_FB2018_05/gff/md5sum.txt

cd $bieserHome/ref_genome/gtf
wget   ftp://ftp.flybase.net/genomes/Drosophila_melanogaster/dmel_r6.24_FB2018_05/gtf/dmel-all-r6.24.gtf.gz
wget   ftp://ftp.flybase.net/genomes/Drosophila_melanogaster/dmel_r6.24_FB2018_05/gtf/md5sum.txt
# done downloading. now figure out differences in the GFF file choices?

####### What's in those GFFs?

####### QC
cd $bieserHome/fastq_joined
fastqc -t 10 *.gz -o ../qc_joined

cd $bieserHome/qc_joined
multiqc ./

cd $bieserHome/fastq_joined

# autodetect adapter from illumina, nextera, and small RNA
# did some of the samples get detected as small RNA adapter. if so,
# WHAT HOW WHAT?
for myfastqlist in `ls *R1.fastq.gz`; do
  echo pattern based on $myfastqlist
  mypattern=${myfastqlist//R1.fastq.gz}
  trim_galore -q 20 --length 36 --trim-n --fastqc --paired "${mypattern}"* &
done

# forcing illumina truseq adapter
# compare vs. the full-auto above
for myforcelist in `ls *R1.fastq.gz`; do
  echo pattern based on $myforcelist
  mypatternforce=${myforcelist//R1.fastq.gz}
  trim_galore -q 20 --length 36 --trim-n --illumina --fastqc --paired "${mypatternforce}"* &
done

# tests 3 and 4
# Test 3. removing polyG from the autos
for myfastqlist in `ls *R1_val_1.fq.gz`; do
  echo pattern based on $myfastqlist
  mypattern=${myfastqlist//R1_val_1.fq.gz}
  trim_galore -o ../fastq_trim_auto_G -a G{18} -q 20 --length 36 --fastqc --paired "${mypattern}"R1_val_1.fq.gz "${mypattern}"R2_val_2.fq.gz &
done

# Test 4. removing polyG from the forced
for myforcelist in `ls *R1_val_1.fq.gz`; do
  echo pattern based on $myforcelist
  mypatternforce=${myforcelist//R1_val_1.fq.gz}
  trim_galore -o ../fastq_trim_force_G -a G{18} -q 20 --length 36 --fastqc --paired "${mypatternforce}"R1_val_1.fq.gz "${mypatternforce}"R2_val_2.fq.gz &
done


A44_S2_R1_val_1_val_1.fq.gz  Control_S5_R1_val_1_val_1.fq.gz  cos2_S3_R1_val_1_val_1.fq.gz  H22_S4_R1_val_1_val_1.fq.gz  L31_S1_R1_val_1_val_1.fq.gz

srun -c 9 --mem 6000 --time 480 bwa mem -t 8 \
-o sam_bwa_mem/A44.sam \
ref_genome/Ensembl_BDGP6/Drosophila_melanogaster.BDGP6.dna.toplevel.fa.gz \
fastq_trim_auto_G/A44_S2_R1_val_1_val_1.fq.gz \
fastq_trim_auto_G/A44_S2_R2_val_2_val_2.fq.gz &> A44.bwa.log &

srun -c 9 --mem 6000 --time 480 bwa mem -t 8 \
-o sam_bwa_mem/Control.sam \
ref_genome/Ensembl_BDGP6/Drosophila_melanogaster.BDGP6.dna.toplevel.fa.gz \
fastq_trim_auto_G/Control_S5_R1_val_1_val_1.fq.gz \
fastq_trim_auto_G/Control_S5_R2_val_2_val_2.fq.gz &> Control.bwa.log &

srun -c 9 --mem 6000 --time 480 bwa mem -t 8 \
-o sam_bwa_mem/cos2.sam \
ref_genome/Ensembl_BDGP6/Drosophila_melanogaster.BDGP6.dna.toplevel.fa.gz \
fastq_trim_auto_G/cos2_S3_R1_val_1_val_1.fq.gz \
fastq_trim_auto_G/cos2_S3_R2_val_2_val_2.fq.gz &> cos2.bwa.log &

srun -c 9 --mem 6000 --time 480 bwa mem -t 8 \
-o sam_bwa_mem/H22.sam \
ref_genome/Ensembl_BDGP6/Drosophila_melanogaster.BDGP6.dna.toplevel.fa.gz \
fastq_trim_auto_G/H22_S4_R1_val_1_val_1.fq.gz \
fastq_trim_auto_G/H22_S4_R2_val_2_val_2.fq.gz &> H22.bwa.log &

srun -c 9 --mem 6000 --time 480 bwa mem -t 8 \
-o sam_bwa_mem/L31.sam \
ref_genome/Ensembl_BDGP6/Drosophila_melanogaster.BDGP6.dna.toplevel.fa.gz \
fastq_trim_auto_G/L31_S1_R1_val_1_val_1.fq.gz \
fastq_trim_auto_G/L31_S1_R2_val_2_val_2.fq.gz &> L31.bwa.log &

cd sam_bwa_mem
for i in A44 Control cos2 H22 L31; do
  srun -c 9 --mem 6000 --time 480 samtools view -@8 -b -o $i.bam $i.sam &
done

for i in A44 Control cos2 H22 L31; do
  srun -c 9 --mem 8000 --time 480 samtools sort -@8 -o $i.sort.bam $i.bam &
done

for i in A44 Control cos2 H22 L31; do
  srun -c 5 --mem 8000 --time 480 samtools rmdup -S $i.sort.bam $i.rmdup.sort.bam &
done

for i in A44 Control cos2 H22 L31; do
  srun -c 9 --mem 8000 --time 120 samtools index -@8 $i.rmdup.sort.bam &
done

## this for loop runs as a separate script that gets invoked like:
## srun -c 10 --mem 24000 --time 130 bash pileupscript.sh &
for i in Control A44 cos2 H22 L31; do
bcftools mpileup -Q 30 -C 50 -P Illumina \
-Ov -B -f ../ref_genome/Ensembl_BDGP6/Drosophila_melanogaster.BDGP6.dna.toplevel.fa \
${i}.rmdup.sort.bam | bcftools call -c -Ov -v -V 'indels' \
-o ${i}.calls.vcf &
done

## this for loop also runs as a separate script that gets invoked like:
## srun -c 10 --mem 24000 --time 130 bash callscript.sh &
## (and only after the pileup script finishes successfully)
for i in Control A44 cos2 H22 L31; do
bcftools mpileup -Q 30 -C 50 -P Illumina \
-Ov -B -f ../ref_genome/Ensembl_BDGP6/Drosophila_melanogaster.BDGP6.dna.toplevel.fa \
-o ${i}.piles.vcf \
${i}.rmdup.sort.bam &
done


# EMS does C to T transversions, 99% of the time. So, C->T or G->A changes.
