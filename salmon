#Starting with Chaoborus 
cd chaoborus-trimmed
#job script 
nano salmonC_trimmed.pbs
#!/bin/bash

#PBS -l walltime=120:00:00,select=1:ncpus=12:mem=20gb
#PBS -N trinity_chao_trim_CLEAN
#PBS -A #alloc-code
#PBS -m abe
#PBS -M #contact
#PBS -o 20220120-chao-output.txt
#PBS -e 20220120-chao-error.txt

################################################################################

module load CVMFS_CC
module load trinity/2.14.0
module load samtools
module load jellyfish
module load gcc/9.3.0  openmpi/4.0.3
module load salmon/1.7.0
module load python
module load scipy-stack
module load kallisto

PATH=/bin/salmon-latest_linux_x86_64/bin/:$PATH
echo "Job execution start: $(date)"

cd $TMPDIR

CPU=12
TXHOME=/lab_data/assemblies/transcriptome/chaoborus-trimmed/
TX=chaoborus.trimmed.Trinity.fasta

cp $TXHOME/$TX $TMPDIR
cp $TXHOME/${TX}.gene_trans_map $TMPDIR

$TRINITY_HOME/util/align_and_estimate_abundance.pl --transcripts $TX \
        --seqType fq --samples_file /lab_data/sequencing/RNAseq/chaoborus/20220321/chao-samples.txt \
        --est_method salmon --output_dir salmon.quant --thread_count $CPU --gene_trans_map ${TX}.gene_trans_map \

tar -czf $PBS_O_WORKDIR/results.tgz ./*



#Mochlonyx 
cd mochlonyx-trimmed
#job script 
nano salmonM_trimmed.pbs
#!/bin/bash

#PBS -l walltime=120:00:00,select=1:ncpus=12:mem=20gb
#PBS -N trinity_moch_trim_CLEAN
#PBS -A #alloc-code
#PBS -m abe
#PBS -M #contact
#PBS -o 2022xxxx-moch-output.txt
#PBS -e 2022xxxx-moch-error.txt

################################################################################

module load CVMFS_CC
module load trinity/2.14.0
module load samtools
module load jellyfish
module load gcc/9.3.0  openmpi/4.0.3
module load salmon/1.7.0
module load python
module load scipy-stack
module load kallisto

PATH=/bin/salmon-latest_linux_x86_64/bin/:$PATH
echo "Job execution start: $(date)"

cd $TMPDIR

CPU=12
TXHOME=/lab_data/assemblies/transcriptome/mochlonyx-trimmed/
TX=mochlonyx.trimmed.Trinity.fasta

cp $TXHOME/$TX $TMPDIR
cp $TXHOME/${TX}.gene_trans_map $TMPDIR

$TRINITY_HOME/util/align_and_estimate_abundance.pl --transcripts $TX \
        --seqType fq --samples_file /lab_data/sequencing/RNAseq/mochlonyx/2023xxxx/moch-samples.txt \
        --est_method salmon --output_dir salmon.quant --thread_count $CPU --gene_trans_map ${TX}.gene_trans_map \

tar -czf $PBS_O_WORKDIR/results.tgz ./*



#Eucorethra 
cd eucorethra-trimmed
#job script 
nano salmonE_trimmed.pbs
#!/bin/bash

#PBS -l walltime=120:00:00,select=1:ncpus=12:mem=20gb
#PBS -N trinity_euc_trim_CLEAN
#PBS -A #alloc-code
#PBS -m abe
#PBS -M #contact
#PBS -o 20220120-euc-output.txt
#PBS -e 20220120-euc-error.txt

################################################################################

module load CVMFS_CC
module load trinity/2.14.0
module load samtools
module load jellyfish
module load gcc/9.3.0  openmpi/4.0.3
module load salmon/1.7.0
module load python
module load scipy-stack
module load kallisto

PATH=/bin/salmon-latest_linux_x86_64/bin/:$PATH
echo "Job execution start: $(date)"

cd $TMPDIR

CPU=12
TXHOME=/lab_data/assemblies/transcriptome/eucorethra-trimmed/
TX=eucorethra.trimmed.Trinity.fasta

cp $TXHOME/$TX $TMPDIR
cp $TXHOME/${TX}.gene_trans_map $TMPDIR

$TRINITY_HOME/util/align_and_estimate_abundance.pl --transcripts $TX \
        --seqType fq --samples_file /lab_data/sequencing/RNAseq/eucorethra/20220817/euc-samples.txt \
        --est_method salmon --output_dir salmon.quant --thread_count $CPU --gene_trans_map ${TX}.gene_trans_map \

tar -czf $PBS_O_WORKDIR/results.tgz ./*

