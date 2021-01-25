#Script 1
#Remove - from file names

#!/bin/bash
for file in *-*; do
  mv $file "${file//-/}"
done

#convert to unix before executing 

#Script 2
#FastQC Quality Check

#BSUB -L /bin/bash              # uses the bash login shell to initialize the job's execution environment.
#BSUB -J fastqc                 # job name
#BSUB -n 2                     # assigns 2 cores for execution
#BSUB -R "span[ptile=2]"       # assigns 20 cores per node
#BSUB -R "rusage[mem=2700]"     # reserves 2700MB memory per core
#BSUB -M 2700                   # sets to 2,700MB (~2.7GB) per process enforceable memory limit. Total job memory = (M * n)
#BSUB -W 24:00                  # sets to 24 hours the job's runtime wall-clock limit.
#BSUB -o stdout.%J              # directs the job's standard output to stdout.jobid
#BSUB -e stderr.%J              # directs the job's standard error to stderr.jobid
#BSUB -u venura.herath@tamu.edu       #Send all emails to email_address
#BSUB -B -N                  #Send email on job begin (-B) and end (-N)

module load FastQC/0.11.9-Java-11

########## INPUTS ##########
SAMPLES="PVX_2d_1 PVX_2d_2 PVX_2d_3 PVX_3d_1 PVX_3d_2 PVX_3d_3 VMOCK_2d_1 VMOCK_2d_2 VMOCK_2d_3 VMOCK_3d_1 VMOCK_3d_2 VMOCK_3d_3"


######## PARAMETERS ########
threads=2                       # make sure this is <= your BSUB -n value


################################### COMMANDS ###################################
# use -o <directory> to save results to <directory> instead of directory where reads are located
#   <directory> must already exist before using -o <directory> option
# --nogroup will calculate average at each base instead of bins after the first 50 bp
# fastqc runs one thread per file; using 20 threads for 2 files does not speed up the processing

for SAMPLE in $SAMPLES; do
fastqc -t $threads -o ./fastQC ${SAMPLE}_1.fq.gz ${SAMPLE}_2.fq.gz
done

#convert to unix before executing 


#Script 3
#HISAT2 Allignment and SAM to BAM convertion

#BSUB -L /bin/bash              # uses the bash login shell to initialize the job's execution environment.
#BSUB -J hisat2                 # job name
#BSUB -n 20                     # assigns 20 cores for execution
#BSUB -R "span[ptile=20]"       # assigns 20 cores per node
#BSUB -R "rusage[mem=2700]"     # reserves 2700MB memory per core
#BSUB -M 2700                   # sets to 2,700MB (~2.7GB) per process enforceable memory limit. Total job memory = (M * n)
#BSUB -W 24:00                  # sets to 24 hours the job's runtime wall-clock limit.
#BSUB -o stdout.%J              # directs the job's standard output to stdout.jobid
#BSUB -e stderr.%J              # directs the job's standard error to stderr.jobid

module load HISAT2/2.2.0-foss-2018b
module load Python/3.6.6-foss-2018b

#bash script for hisat2; align all .fq.gz files to indexed reference genome to generate .sam files

SAMPLES="PVX_2d_1 PVX_2d_2 PVX_2d_3 PVX_3d_1 PVX_3d_2 PVX_3d_3 VMOCK_2d_1 VMOCK_2d_2 VMOCK_2d_3 VMOCK_3d_1 VMOCK_3d_2 VMOCK_3d_3"
threads=20

for SAMPLE in $SAMPLES; do
    hisat2 -p $threads --dta --rna-strandness RF -x /scratch/datasets/genome_indexes/other_genomes/potato/hisat2/DM_1-3_516_R44_potato_genome_assembly.v6.1 -1 ${SAMPLE}_1.fq.gz -2 ${SAMPLE}_2.fq.gz -S ${SAMPLE}.sam
done

#Purge module to avoid conflicts
module purge

#bash script for samtools; convert .sam files to .bam files
module load SAMtools/1.9-intel-2018b

for SAMPLE in $SAMPLES; do
    samtools sort -@ $threads -o ${SAMPLE}.bam ${SAMPLE}.sam
done

#bash script for samtools; index our .bam files to obtain .bam.bai files using samtools

for SAMPLE in $SAMPLES; do
    samtools index ${SAMPLE}.bam ${SAMPLE}.bam.bai
done

#convert to unix before executing

#Script 4
#Stringtie to DESEQ2

#BSUB -L /bin/bash              # uses the bash login shell to initialize the job's execution environment.
#BSUB -J hisat2                 # job name
#BSUB -n 20                     # assigns 20 cores for execution
#BSUB -R "span[ptile=20]"       # assigns 20 cores per node
#BSUB -R "rusage[mem=2700]"     # reserves 2700MB memory per core
#BSUB -M 2700                   # sets to 2,700MB (~2.7GB) per process enforceable memory limit. Total job memory = (M * n)
#BSUB -W 24:00                  # sets to 24 hours the job's runtime wall-clock limit.
#BSUB -o stdout.%J              # directs the job's standard output to stdout.jobid
#BSUB -e stderr.%J              # directs the job's standard error to stderr.jobid
#BSUB -u venura.herath@tamu.edu       #Send all emails to email_address
#BSUB -B -N                  #Send email on job begin (-B) and end (-N)

#Loading New modules

module load StringTie/2.1.4-GCC-9.3.0
module load Python/2.7.18-GCCcore-9.3.0

#bash script for stringtie; assemble transcripts using stringtie

SAMPLES="PVX_2d_1 PVX_2d_2 PVX_2d_3 PVX_3d_1 PVX_3d_2 PVX_3d_3 VMOCK_2d_1 VMOCK_2d_2 VMOCK_2d_3 VMOCK_3d_1 VMOCK_3d_2 VMOCK_3d_3"
threads=20

for SAMPLE in $SAMPLES; do
    stringtie --rf -p $threads -G DM_1-3_516_R44_potato.v6.1.working_models.gff3 -o ${SAMPLE}.gtf -l ${SAMPLE} ${SAMPLE}.bam
done

#merge transcripts
stringtie --merge -p $threads  -G DM_1-3_516_R44_potato.v6.1.working_models.gff3 -o stringtie_merged.gtf mergelist.txt

#purge modules
module purge

#Loading GFFcompare 
module load GffCompare/0.10.6-GCCcore-7.3.0

#Using Gffcompare
gffcompare -r DM_1-3_516_R44_potato.v6.1.working_models.gff3 -G -o merged stringtie_merged.gtf

#purge modules
module purge

#bash script for stringtie; assemble transcripts using stringtie

SAMPLES="PVX_2d_1 PVX_2d_2 PVX_2d_3 PVX_3d_1 PVX_3d_2 PVX_3d_3 VMOCK_2d_1 VMOCK_2d_2 VMOCK_2d_3 VMOCK_3d_1 VMOCK_3d_2 VMOCK_3d_3"
threads=20


module load StringTie/2.1.4-GCC-9.3.0
module load Python/2.7.18-GCCcore-9.3.0

for SAMPLE in $SAMPLES; do
    mkdir ballgown/${SAMPLE}
    stringtie --rf -e -B -p $threads -G stringtie_merged.gtf -o ballgown/${SAMPLE}/${SAMPLE}.gtf ${SAMPLE}.bam
done

#converting data to DESEQ2 compatible format

prepDE.py -i ballgown
#end

