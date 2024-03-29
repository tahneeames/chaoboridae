##Chaoborus assembly using Trinity RNAseq
#Create sequencing directory for Chaoborus samples 
cd /lab_data/sequencing/RNAseq/chaoborus/20220321

##create sample file for assembly from reads 
nano chao-samples.txt

#  Or,
#      --samples_file /lab_data/sequencing/RNAseq/chaoborus/20220321       tab-delimited text file indicating biological replicate relationships.
#                                  ex.



chao_as  chao_as_rep1 /lab_data/sequencing/RNAseq/chaoborus/20220321/NS.1829.002.NEBNext_dual_i7_A1---NEBNext_dual_i5_A1.as_1-17_1_R1.fastq.gz /lab_data/sequencing/RNAseq/chaoborus/20220321/NS.1829.002.NEBNext_dual_i7_A1---NEBNext_dual_i5_A1.as_1-17_1_R2.fastq.gz
chao_as  chao_as_rep2 /lab_data/sequencing/RNAseq/chaoborus/20220321/NS.1829.002.NEBNext_dual_i7_D1---NEBNext_dual_i5_D1.as_1-17_2_R1.fastq.gz /lab_data/sequencing/RNAseq/chaoborus/20220321/NS.1829.002.NEBNext_dual_i7_D1---NEBNext_dual_i5_D1.as_1-17_2_R2.fastq.gz  
chao_as  chao_as_rep3 /lab_data/sequencing/RNAseq/chaoborus/20220321/NS.1829.002.NEBNext_dual_i7_G1---NEBNext_dual_i5_G1.as_1-17_3_R1.fastq.gz /lab_data/sequencing/RNAseq/chaoborus/20220321/NS.1829.002.NEBNext_dual_i7_G1---NEBNext_dual_i5_G1.as_1-17_3_R2.fastq.gz
chao_mt  chao_mt_rep1 /lab_data/sequencing/RNAseq/chaoborus/20220321/NS.1829.002.NEBNext_dual_i7_B1---NEBNext_dual_i5_B1.mt_1-17_1_R1.fastq.gz /lab_data/sequencing/RNAseq/chaoborus/20220321/NS.1829.002.NEBNext_dual_i7_B1---NEBNext_dual_i5_B1.mt_1-17_1_R2.fastq.gz
chao_mt  chao_mt_rep2 /lab_data/sequencing/RNAseq/chaoborus/20220321/NS.1829.002.NEBNext_dual_i7_E1---NEBNext_dual_i5_E1.mt_1-17_2_R1.fastq.gz /lab_data/sequencing/RNAseq/chaoborus/20220321/NS.1829.002.NEBNext_dual_i7_E1---NEBNext_dual_i5_E1.mt_1-17_2_R2.fastq.gz
chao_mt  chao_mt_rep3 /lab_data/sequencing/RNAseq/chaoborus/20220321/NS.1829.002.NEBNext_dual_i7_H1---NEBNext_dual_i5_H1.mt_1-17_3_R1.fastq.gz /lab_data/sequencing/RNAseq/chaoborus/20220321/NS.1829.002.NEBNext_dual_i7_H1---NEBNext_dual_i5_H1.mt_1-17_3_R2.fastq.gz
chao_gi  chao_gi_rep1 /lab_data/sequencing/RNAseq/chaoborus/20220321/NS.1829.002.NEBNext_dual_i7_C1---NEBNext_dual_i5_C1.gi_1-17_1_R1.fastq.gz /lab_data/sequencing/RNAseq/chaoborus/20220321/NS.1829.002.NEBNext_dual_i7_C1---NEBNext_dual_i5_C1.gi_1-17_1_R2.fastq.gz
chao_gi  chao_gi_rep2 /lab_data/sequencing/RNAseq/chaoborus/20220321/NS.1829.002.NEBNext_dual_i7_F1---NEBNext_dual_i5_F1.gi_1-17_2_R1.fastq.gz /lab_data/sequencing/RNAseq/chaoborus/20220321/NS.1829.002.NEBNext_dual_i7_F1---NEBNext_dual_i5_F1.gi_1-17_2_R2.fastq.gz
chao_gi  chao_gi_rep3 /lab_data/sequencing/RNAseq/chaoborus/20220321/NS.1829.002.NEBNext_dual_i7_A2---NEBNext_dual_i5_A2.gi_1-17_3_R1.fastq.gz /lab_data/sequencing/RNAseq/chaoborus/20220321/NS.1829.002.NEBNext_dual_i7_A2---NEBNext_dual_i5_A2.gi_1-17_3_R2.fastq.gz


##Assemble with raw reads
#copy samples file to assembly directory 
cp chao-samples.txt /lab_data/assemblies/transcriptome/chaoborus-untrimmed

#move to assembly directory  
cd chaoborus-untrimmed

#job script  
nano trinityC_untrimmed.pbs
#!/bin/bash

#PBS -l walltime=120:00:00,select=1:ncpus=32:mem=187gb
#PBS -N trinity_chao_untrim_CLEAN
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

echo "Job execution start: $(date)"

cd $TMPDIR

${TRINITY_HOME}/Trinity --seqType fq        \
        --CPU 32 --max_memory 100G \
        --samples_file $PBS_O_WORKDIR/chao-samples.txt   \

rm $TMPDIR/trinity_out_dir/*.gz
rm -r $TMPDIR/trinity_out_dir/chrysalis
rm -r $TMPDIR/trinity_out_dir/insilico_read_normalization
rm -r $TMPDIR/trinity_out_dir/read_partitions
rm -r $TMPDIR/trinity_out_dir/__salmon_filt.chkpts

tar -czf $PBS_O_WORKDIR/results.tgz ./*



##Assemble with trimmomatic 
#copy samples file to assembly directory 
cp chao-samples.txt /lab_data/assemblies/transcriptome/chaoborus-trimmed

#move to assembly directory 
cd chaoborus-trimmed

#job script 
nano trinityC_trimmed.pbs
#!/bin/bash

#PBS -l walltime=120:00:00,select=1:ncpus=32:mem=375gb
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

echo "Job execution start: $(date)"

cd $TMPDIR


${TRINITY_HOME}/Trinity --seqType fq        \
        --CPU 32 --max_memory 80G --trimmomatic \
        --samples_file $PBS_O_WORKDIR/chao-samples.txt   \

cp trinity_out_dir/*Trinity.*.fasta $PBS_O_WORKDIR/

rm trinity_out_dir/*.fq trinity_out_dir/*.fq.gz trinity_out_dir/*.fa
tar -czf $PBS_O_WORKDIR/results.tgz ./*




##Mochlonyx assembly using Trinity RNAseq
#Create sequencing directory for Mochlonyx samples 
cd /lab_data/sequencing/RNAseq/mochlonyx/2023xxxx

##create sample file for assembly from reads 
nano moch-samples.txt

#  Or,
#      --samples_file /lab_data/sequencing/RNAseq/mochlonyx/2023xxxx       tab-delimited text file indicating biological replicate relationships.
#                                  ex.



chao_as  chao_as_rep1 /lab_data/sequencing/RNAseq/ /lab_data/sequencing/RNAseq/
chao_as  chao_as_rep2 /lab_data/sequencing/RNAseq/ /lab_data/sequencing/RNAseq/ 
chao_as  chao_as_rep3 /lab_data/sequencing/RNAseq/ /lab_data/sequencing/RNAseq/
chao_mt  chao_mt_rep1 /lab_data/sequencing/RNAseq/ /lab_data/sequencing/RNAseq/
chao_mt  chao_mt_rep2 /lab_data/sequencing/RNAseq/ /lab_data/sequencing/RNAseq/
chao_mt  chao_mt_rep3 /lab_data/sequencing/RNAseq/ /lab_data/sequencing/RNAseq/
chao_gi  chao_gi_rep1 /lab_data/sequencing/RNAseq/ /lab_data/sequencing/RNAseq/
chao_gi  chao_gi_rep2 /lab_data/sequencing/RNAseq/ /lab_data/sequencing/RNAseq/
chao_gi  chao_gi_rep3 /lab_data/sequencing/RNAseq/ /lab_data/sequencing/RNAseq/


##Assemble with raw reads
#copy samples file to assembly directory 
cp moch-samples.txt /lab_data/assemblies/transcriptome/mochlonyx-untrimmed

#move to assembly directory  
cd mochlonyx-untrimmed

#job script  
nano trinityM_untrimmed.pbs
#!/bin/bash

#PBS -l walltime=120:00:00,select=1:ncpus=32:mem=187gb
#PBS -N trinity_moch_untrim_CLEAN
#PBS -A #alloc-code
#PBS -m abe
#PBS -M #contact
#PBS -o 2023xxxx-moch-output.txt
#PBS -e 2023xxxx-moch-error.txt

################################################################################

module load CVMFS_CC
module load trinity/2.14.0
module load samtools
module load jellyfish
module load gcc/9.3.0  openmpi/4.0.3
module load salmon/1.7.0
module load python
module load scipy-stack

echo "Job execution start: $(date)"

cd $TMPDIR

${TRINITY_HOME}/Trinity --seqType fq        \
        --CPU 32 --max_memory 100G \
        --samples_file $PBS_O_WORKDIR/moch-samples.txt   \

rm $TMPDIR/trinity_out_dir/*.gz
rm -r $TMPDIR/trinity_out_dir/chrysalis
rm -r $TMPDIR/trinity_out_dir/insilico_read_normalization
rm -r $TMPDIR/trinity_out_dir/read_partitions
rm -r $TMPDIR/trinity_out_dir/__salmon_filt.chkpts

tar -czf $PBS_O_WORKDIR/results.tgz ./*



##Assemble with trimmomatic 
#copy samples file to assembly directory 
cp moch-samples.txt /lab_data/assemblies/transcriptome/mochlonyx-trimmed

#move to assembly directory 
cd mochlonyx-trimmed

#job script 
nano trinityM_trimmed.pbs
#!/bin/bash

#PBS -l walltime=120:00:00,select=1:ncpus=32:mem=375gb
#PBS -N trinity_moch_trim_CLEAN
#PBS -A #alloc-code
#PBS -m abe
#PBS -M #contact 
#PBS -o 2023xxxx-moch-output.txt
#PBS -e 2023xxxx-moch-error.txt

################################################################################

module load CVMFS_CC
module load trinity/2.14.0
module load samtools
module load jellyfish
module load gcc/9.3.0  openmpi/4.0.3
module load salmon/1.7.0
module load python
module load scipy-stack

echo "Job execution start: $(date)"

cd $TMPDIR


${TRINITY_HOME}/Trinity --seqType fq        \
        --CPU 32 --max_memory 80G --trimmomatic \
        --samples_file $PBS_O_WORKDIR/moch-samples.txt   \

cp trinity_out_dir/*Trinity.*.fasta $PBS_O_WORKDIR/

rm trinity_out_dir/*.fq trinity_out_dir/*.fq.gz trinity_out_dir/*.fa
tar -czf $PBS_O_WORKDIR/results.tgz ./*




##Eucorethra assembly using Trinity RNAseq
#Create sequencing directory for Eucorethra samples 
cd /lab_data/sequencing/RNAseq/eucorethra/20220817

##create sample file for assembly from reads 
nano euc-samples.txt
#  Or,
#      --samples_file /lab_data/sequencing/RNAseq/eucorethra/20220817        tab-delimited text file indicating biological replicate relationships.

#                                  ex.


euc_tr      euc_tr_rep1 /lab_data/sequencing/RNAseq/eucorethra/20220817/NS.1961.004.UDP0021_i7---UDP0021_i5.tr1_eu_R1.fastq.gz     /lab_data/sequencing/RNAseq/eucorethra/20220817/NS.1961.004.UDP0021_i7---UDP0021_i5.tr1_eu_R2.fastq.gz
euc_tr      euc_tr_rep2 /lab_data/sequencing/RNAseq/eucorethra/20220817/NS.1961.004.UDP0024_i7---UDP0024_i5.tr2_eu_R1.fastq.gz     /lab_data/sequencing/RNAseq/eucorethra/20220817/NS.1961.004.UDP0024_i7---UDP0024_i5.tr2_eu_R2.fastq.gz
euc_tr      euc_tr_rep3 /lab_data/sequencing/RNAseq/eucorethra/20220817/NS.1961.004.UDP0027_i7---UDP0027_i5.tr3_eu_R1.fastq.gz     /lab_data/sequencing/RNAseq/eucorethra/20220817/NS.1961.004.UDP0027_i7---UDP0027_i5.tr3_eu_R2.fastq.gz
euc_mt      euc_mt_rep1 /lab_data/sequencing/RNAseq/eucorethra/20220817/NS.1961.004.UDP0022_i7---UDP0022_i5.mt1_eu_R1.fastq.gz     /lab_data/sequencing/RNAseq/eucorethra/20220817/NS.1961.004.UDP0022_i7---UDP0022_i5.mt1_eu_R2.fastq.gz
euc_mt      euc_mt_rep2 /lab_data/sequencing/RNAseq/eucorethra/20220817/NS.1961.004.UDP0025_i7---UDP0025_i5.mt2_eu_R1.fastq.gz     /lab_data/sequencing/RNAseq/eucorethra/20220817/NS.1961.004.UDP0025_i7---UDP0025_i5.mt2_eu_R2.fastq.gz
euc_mt      euc_mt_rep3 /lab_data/sequencing/RNAseq/eucorethra/20220817/NS.1961.004.UDP0028_i7---UDP0028_i5.mt3_eu_R1.fastq.gz     /lab_data/sequencing/RNAseq/eucorethra/20220817/NS.1961.004.UDP0028_i7---UDP0028_i5.mt3_eu_R2.fastq.gz
euc_gi      euc_gi_rep1 /lab_data/sequencing/RNAseq/eucorethra/20220817/NS.1961.004.UDP0023_i7---UDP0023_i5.gi1_eu_R1.fastq.gz     /lab_data/sequencing/RNAseq/eucorethra/20220817/NS.1961.004.UDP0023_i7---UDP0023_i5.gi1_eu_R2.fastq.gz
euc_gi      euc_gi_rep2 /lab_data/sequencing/RNAseq/eucorethra/20220817/NS.1961.004.UDP0026_i7---UDP0026_i5.gi2_eu_R1.fastq.gz     /lab_data/sequencing/RNAseq/eucorethra/20220817/NS.1961.004.UDP0026_i7---UDP0026_i5.gi2_eu_R2.fastq.gz


##Assemble with raw reads
#copy samples file to assembly directory 
cp euc-samples.txt /lab_data/assemblies/transcriptome/eucorethra-untrimmed

#move to job directory  
cd eucorethra-untrimmed

#job script  
nano trinityEI_untrimmed.pbs
#!/bin/bash

#PBS -l walltime=120:00:00,select=1:ncpus=32:mem=187gb
#PBS -N trinity_eucor_notrim_CLEAN
#PBS -A #alloc-code
#PBS -m abe
#PBS -M #contact
#PBS -o 20220120-output.txt
#PBS -e 20220120-error.txt

################################################################################

module load CVMFS_CC
module load trinity/2.14.0
module load samtools
module load jellyfish
module load gcc/9.3.0  openmpi/4.0.3
module load salmon/1.7.0
module load python
module load scipy-stack


echo "Job execution start: $(date)"

cd $TMPDIR

cd $TMPDIR

${TRINITY_HOME}/Trinity --seqType fq        \
        --CPU 32 --max_memory 100G \
        --samples_file $PBS_O_WORKDIR/euc-samples.txt   \
        --SS_lib_type RF


cp trinity_out_dir/*Trinity.fasta $PBS_O_WORKDIR/

rm trinity_out_dir/*.fq trinity_out_dir/*.fq.gz

tar -czf $PBS_O_WORKDIR/results.tgz ./*



##Assemble with trimmomatic 
#copy samples file to assembly directory 
cp euc-samples.txt /lab_data/assemblies/transcriptome/eucorethra-trimmed

#move to assembly directory 
cd eucorethra-trimmed

#job script 
nano trinityEI_trimmed.pbs
#!/bin/bash

#PBS -l walltime=120:00:00,select=1:ncpus=32:mem=372gb
#PBS -N trinity_eucor_trim_CLEAN
#PBS -A #alloc-code
#PBS -m abe
#PBS -M #contact
#PBS -o 20220120-output.txt
#PBS -e 20220120-error.txt

################################################################################

module load CVMFS_CC
module load trinity/2.14.0
module load samtools
module load jellyfish
module load gcc/9.3.0  openmpi/4.0.3
module load salmon/1.7.0
module load python
module load scipy-stack

echo "Job execution start: $(date)"

cd $TMPDIR

${TRINITY_HOME}/Trinity --seqType fq        \
        --CPU 32 --max_memory 60G --trimmomatic \
        --samples_file $PBS_O_WORKDIR/euc-samples.txt   \
        --SS_lib_type RF

cp trinity_out_dir/*Trinity.*.fasta $PBS_O_WORKDIR/

rm trinity_out_dir/*.fq trinity_out_dir/*.fq.gz trinity_out_dir/*.fa

tar -czf $PBS_O_WORKDIR/results.tgz ./*
 
