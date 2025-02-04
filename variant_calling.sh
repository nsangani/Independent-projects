#!/bin/bash
#SBATCH -J NGS_homework
#SBATCH -p general
#SBATCH -o filename_%j.txt
#SBATCH -e filename_%j.err
#SBATCH --mail-type=fail
#SBATCH --mail-user=nsangani@iu.edu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --mem=90G
#SBATCH --time=8:00:00

# the above commands will let the HPC know you're requesting for 1 node with 8 cores, memory 90 GB, and running time is 8 hours.

# the simplest variant calling pipeline with bwa-samtools-bcftools
# we only map to one chrosome X
module load bwa/0.7.12
module load samtools/1.9
module load bcftools/1.9

# set a workdirectory for your job
WORKDIR=/N/u/cpdong/Carbonate/NGS_hw

[ -d $WORKDIR ] && echo "Directory Exists" || mkdir -m 755 $WORKDIR
mkdir -p  $WORKDIR/refgenome/
mkdir -p  $WORKDIR/fastq/
mkdir -p  $WORKDIR/bam/
mkdir -p  $WORKDIR/vcf/

# fastq file download in GoogleDrive shared link. ~ 5 minutes
cd $WORKDIR/fastq/
# the following long syntax only for downloading google drive files using wget commands
wget --load-cookies /tmp/cookies.txt "https://docs.google.com/uc?export=download&confirm=$(wget --quiet --save-cookies /tmp/cookies.txt --keep-session-cookies --no-check-certificate 'https://docs.google.com/uc?export=download&id=12ncw3mMptP-5PmeZZo9wEmpOasbwxPK0' -O- | sed -rn 's/.*confirm=([0-9A-Za-z_]+).*/\1\n/p')&id=12ncw3mMptP-5PmeZZo9wEmpOasbwxPK0" -O S1_R1.fq.gz && rm -rf /tmp/cookies.txt
wget --load-cookies /tmp/cookies.txt "https://docs.google.com/uc?export=download&confirm=$(wget --quiet --save-cookies /tmp/cookies.txt --keep-session-cookies --no-check-certificate 'https://docs.google.com/uc?export=download&id=12g125t7LmeFRtGuIWFSL83ywwAd_jdAi' -O- | sed -rn 's/.*confirm=([0-9A-Za-z_]+).*/\1\n/p')&id=12g125t7LmeFRtGuIWFSL83ywwAd_jdAi" -O S1_R2.fq.gz && rm -rf /tmp/cookies.txt


# download reference genome file via Gencode and build alignment genome index
cd $WORKDIR/refgenome/
# about 10 minutes
# wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_32/GRCh38.primary_assembly.genome.fa.gz
wget ftp://ftp.ensembl.org/pub/release-98/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.X.fa.gz
gunzip Homo_sapiens.GRCh38.dna.chromosome.X.fa.gz
# about 1 hour or more
bwa index Homo_sapiens.GRCh38.dna.chromosome.X.fa



# alignment paired fastq files using built reference index
bwa mem -t 8 -M $WORKDIR/refgenome/Homo_sapiens.GRCh38.dna.chromosome.X.fa $WORKDIR/fastq/S1_R1.fq.gz $WORKDIR/fastq/S1_R2.fq.gz > $WORKDIR/bam/read.sam



# bam file precessing with samtools before calling variants
cd $WORKDIR/bam/
samtools view -bS read.sam >read.bam
samtools flagstat read.bam
samtools sort -@ 8 -o read.sorted.bam  read.bam
samtools view -h -F4  -q 5 read.sorted.bam |samtools view -bS |samtools rmdup -  read.filter.rmdup.bam
samtools index read.filter.rmdup.bam



# using bcftools calling viriants
samtools mpileup -ugf $WORKDIR/refgenome/Homo_sapiens.GRCh38.dna.chromosome.X.fa  $WORKDIR/bam/read.filter.rmdup.bam  |bcftools call -vmO z -o $WORKDIR/vcf/read.bcftools.vcf.gz

# end
