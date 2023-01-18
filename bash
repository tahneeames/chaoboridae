
##this was all executed on a hpc cluster - some as batch jobs and some as interactive jobs using these paths - to run, make sure to change these paths 

### Chaoborus trivittatus assembly and annotation ###

##connect to cluster and set directory - using the project file to store sample files 
cd /arc/project/trivittatus

##transcriptome fasta files uploaded onto  cluster  
##these are the raw bp reads for each of the nine samples, which were then renamed by tissue, replicate, and read orientation
Ct_as1_left.fastq.gz
Ct_as1_right.fastq.gz
Ct_as2_left.fastq.gz
Ct_as2_right.fastq.gz
Ct_as3_left.fastq.gz
Ct_as3_right.fastq.gz
Ct_gi1_left.fastq.gz
Ct_gi1_right.fastq.gz
Ct_gi2_left.fastq.gz
Ct_gi2_right.fastq.gz
Ct_gi3_left.fastq.gz
Ct_gi3_right.fastq.gz
Ct_mt1_left.fastq.gz
Ct_mt1_right.fastq.gz
Ct_mt2_left.fastq.gz
Ct_mt2_right.fastq.gz
Ct_mt3_left.fastq.gz
Ct_mt3_right.fastq.gz


##create samplesFile for assembly 
nano ctrivTX1_samplesFile.txt

#  Or,
#      --samples_file /arc/project/trivittatus       tab-delimited text file indicating biological replicate relationships.
#                                  ex.



Ct_as  Ct_as_rep1 /arc/project/trivittatus/Ct_as1_left.fastq.gz /arc/project/trivittatus/Ct_as1_right.fastq.gz
Ct_as  Ct_as_rep2 /arc/project/trivittatus/Ct_as2_left.fastq.gz /arc/project/trivittatus/Ct_as2_right.fastq.gz    
Ct_as  Ct_as_rep3 /arc/project/trivittatus/Ct_as3_left.fastq.gz /arc/project/trivittatus/Ct_as3_right.fastq.gz
Ct_mt  Ct_mt_rep1 /arc/project/trivittatus/Ct_mt1_left.fastq.gz /arc/project/trivittatus/Ct_mt1_right.fastq.gz
Ct_mt  Ct_mt_rep2 /arc/project/trivittatus/Ct_mt2_left.fastq.gz /arc/project/trivittatus/Ct_mt2_right.fastq.gz
Ct_mt  Ct_mt_rep3 /arc/project/trivittatus/Ct_mt3_left.fastq.gz /arc/project/trivittatus/Ct_mt3_right.fastq.gz
Ct_gi  Ct_gi_rep1 /arc/project/trivittatus/Ct_gi1_left.fastq.gz /arc/project/trivittatus/Ct_gi1_right.fastq.gz
Ct_gi  Ct-gi_rep2 /arc/project/trivittatus/Ct_gi2_left.fastq.gz /arc/project/trivittatus/Ct_gi2_right.fastq.gz
Ct_gi  Ct_gi_rep3 /arc/project/trivittatus/Ct_gi3_left.fastq.gz /arc/project/trivittatus/Ct_gi3_right.fastq.gz


##make a chaoborus trivittatus directory in the scratch folder 
mkdir /scratch/trivittatusRaw
mkdir /scratch/trivittatusRaw/trinityAssembly
cd /scratch/trivittatusRaw/trinityAssembly

##run the trinity assembly in this directoty as a batch job - first running assembly with no trimming
##create a pbs file 
nano trinity_untrim.pbs

#!/bin/bash
 
#PBS -l walltime=120:00:00,select=1:ncpus=32:mem=100gb
#PBS -N trinity_ctrivTX1
#PBS -A [account name]
#PBS -m abe
#PBS -M [cluster name]
#PBS -o output.txt
#PBS -e error.txt
 
################################################################################
 

export TRINITY_HOME=/arc/project/bin/trinityrnaseq-v2.13.2/

#load modules 
module load Software_Collection/2019
module load samtools
module load boost
module load python
module load py-virtualenv

#set current cd as job directory 
cd $PBS_O_WORKDIR

source /scratch/trinity_env/bin/activate

#use shell script to run the assembly 
sh runUntrim.sh


##create execution file 
nano runUntrim.sh
#!/bin/bash -ve

#######################################################
##  Run Trinity to Generate Transcriptome Assemblies ##
#######################################################

export PATH=/arc/project/bin/bowtie2-2.4.5/:/arc/project/bin/jellyfish-2.3.0/bin/:/arc/project/bin/salmon-1.7.0_linux_x86_64/bin/:$PATH
export TRINITY_HOME=/arc/project/bin/trinityrnaseq-v2.13.2/

if [ -z ${TRINITY_HOME} ]; then
    echo "Must set env var TRINITY_HOME"
    exit 1
fi

#parameters for assembly 
${TRINITY_HOME}/Trinity --seqType fq        \
        --CPU 32 --max_memory 100G  \
        --samples_file /arc/project/trivittatus/ctrivTX1_samplesFile.txt  \
           

##### Done Running Trinity #####

if [ $* ]; then
    # check full-length reconstruction stats:

    ${TRINITY_HOME}/util/misc/illustrate_ref_comparison.pl __indiv_ex_sample_derived/refSeqs.fa trinity_out_dir.Trinity.fasta 90

    ./test_FL.sh --query  trinity_out_dir.Trinity.fasta --target __indiv_ex_sample_derived/refSeqs.fa --no_reuse

fi
           


##submit the job 
qsub trinity_test.pbs

##download the fasta files locally (trinity_out_dir.Trinity.fasta and trinity_out_dir.Trinity.fasta.gene_trans_map) - I used sftp 

##run assembly with trimmomatic to compare a trimmed and untrimmed assembly. To keep all files independent and trackable, it's best to make a new directory
mkdir /scratch/trivitatusTrim
mkdir /scratch/trivitatusTrim/trinityAssembly
cd /scratch/trivitatusTrim/trinityAssembly

nano trinity_trimmed.pbs
#!/bin/bash
 
#PBS -l walltime=120:00:00,select=4:ncpus=32:mem=100gb
#PBS -N trinity_ctrivTX2
#PBS -A [account name]
#PBS -m abe
#PBS -M [cluster name]
#PBS -o output.txt
#PBS -e error.txt
 
################################################################################

export TRINITY_HOME=/arc/project/bin/trinityrnaseq-v2.13.2/

module load Software_Collection/2019
module load samtools
module load boost
module load python
module load py-virtualenv
 
cd $PBS_O_WORKDIR

source /scratch/trinity_env/bin/activate

sh runTrim.sh

##write run file 
nano runTrim.sh
#!/bin/bash -ve

#######################################################
##  Run Trinity to Generate Transcriptome Assemblies ##
#######################################################

export PATH=/arc/project/bin/bowtie2-2.4.5/:/arc/project/bin/jellyfish-2.3.0/bin/:/arc/project/bin/salmon-1.7.0_linux_x86_64/bin/:$PATH
export TRINITY_HOME=/arc/project/bin/trinityrnaseq-v2.13.2/

if [ -z ${TRINITY_HOME} ]; then
    echo "Must set env var TRINITY_HOME"
    exit 1
fi


${TRINITY_HOME}/Trinity --seqType fq        \
        --CPU 32 --max_memory 100G --trimmomatic \
        --samples_file /arc/project/trivittatus/ctrivTX1_samplesFile.txt   \
           

##### Done Running Trinity #####

if [ $* ]; then
    # check full-length reconstruction stats:

    ${TRINITY_HOME}/util/misc/illustrate_ref_comparison.pl __indiv_ex_sample_derived/refSeqs.fa trinity_out_dir.Trinity.fasta 90

    ./test_FL.sh --query  trinity_out_dir.Trinity.fasta --target __indiv_ex_sample_derived/refSeqs.fa --no_reuse

fi

##submit this job 
qsub trinity_trimmed.pbs


##running fastqc to assess read quality (especially between trimmed and raw reads)
mkdir /scratch/trivittatus_readQuality
cd /scratch/trivittatus_readQuality
#!/bin/bash

#PBS -l walltime=72:00:00,select=1:ncpus=18:mem=100gb
#PBS -N fastQC
#PBS -A [account name]
#PBS -m abe
#PBS -M [cluster name]
#PBS -o output.txt
#PBS -e error.txt
 
################################################################################


module load CVMFS_CC
module load fastqc 

cd $PBS_O_WORKDIR

#execute fastqc 
fastqc /arc/project/trivittatus/*.gz -o -t 18 /scratch/trivittatus_readQuality/fastqc_report.txt


##assessing quality of reads 
##extracting contig length from trimmed assembly
##move to file directory
cd /scratch/trivittatusTrim/trinityAssembly
##make text file with contig lengths 
grep -oP '(?<=len=)(\s+)?\K([^ ]*)' trinity_out_dir.Trinity.fasta > contigLength_trimRead.txt

##get 5 prime and 3 prime partials, write into text files 
cd trinity_out_dir.Trinity.fasta.transdecoder_dir

grep -o "5prime" longest_orfs.cds > fivePrime_regionsTr.txt
grep -o "3prime" longest_orfs.cds > threePrime_regionsTr.txt
##all you really need are these numbers - don't need the actual txt files, but okay to have anyway
wc -l fivePrime_regionsTr.txt
wc -l threePrime_regionsTr.txt

##get trinity stats log
/arc/project/trinityrnaseq-v2.13.2/util/TrinityStats.pl trinity_out_dir.Trinity.fasta > trinityTrim_Stats.log

##do the same with the assembly without trimmomatic 
cd scratch/trivittatusRaw/trinityAssembly
grep -oP '(?<=len=)(\s+)?\K([^ ]*)' trinity_out_dir.Trinity.fasta > contigLength_Read.txt

cd trinity_out_dir.Trinity.fasta.transdecoder_dir

grep -o "5prime" longest_orfs.cds > fivePrime_regionsUn.txt
grep -o "3prime" longest_orfs.cds > threePrime_regionsUn.txt
wc -l fivePrime_regionsUn.txt
wc -l threePrime_regionsUn.txt

/arc/project/bin/trinityrnaseq-v2.13.2/util/TrinityStats.pl trinity_out_dir.Trinity.fasta > trinityUntrimmed_Stats.log

##download the contigLength txt files and trinity stats files locally
##plot in R (code available separately)

##move forward with the trimmed assembly for the annotation steps
#make a new directory for the annotation in trim parent directory
mkdir /scratch/trivittatusTrim/trinotateAssembly
cd /scratch/trivittatusTrim/trinotateAssembly

#download swissprot database 
wget  ftp://ftp.ncbi.nlm.nih.gov/blast/db/swissprot.tar.gz
gunnzip swissprot.tar.gz
tar -xvf swissprot.tar 

##make trinotate batch job file 
nano trivTrinotate.pbs 

#!/bin/bash

#PBS -l walltime=96:00:00,select=1:ncpus=40:mem=187gb
#PBS -N trinotate1
#PBS -A [account name]
#PBS -m abe
#PBS -M [cluster name]
#PBS -o output.txt
#PBS -e error.txt
 
################################################################################

cd $PBS_O_WORKDIR

module load CVMFS_CC; module load gcc; module load trinotate
module load hmmer; module load blast+; module load signalp; module load rnammer
module load transdecoder; module load tmhmm 



# HMMSCAN
hmmscan --cpu 40 --domtblout TrinotatePFAM.out /scratch/trinotate/Trinotate-Trinotate-v3.2.2/Pfam-A.hmm /scratch/trivittatusTrim/trinityAssembly/longest_orfs.pep > pfam.log


# SignalP
signalp -f short -n signalp.out /scratch/trivittatusTrim/trinityAssembly/longest_orfs.pep > signalp.out

# TMHMM
tmhmm --short < /scratch/trivittatusTrim/trinityAssembly/longest_orfs.pep > tmhmm.out


#RNAmmer
/scratch/trinotate/Trinotate-Trinotate-v3.2.2/util/rnammer_support/RnammerTranscriptome.pl --transcriptome /scratch/trivittatusTrim/trinityAssembly/trinity_out_dir.Trinity.fasta --path_to_rnammer /cvmfs/soft.computecanada.ca/easybuild/software/2020/avx2/Core/rnammer/1.2/rnammer

#blastx 
blastx -query /scratch/trivittatusTrim/trinityAssembly/trinity_out_dir.Trinity.fasta \
         -db /scratch/trivittatusTrim/trinotateAssembly/swissprot -num_threads 2 \
         -max_target_seqs 1 -outfmt 6 -evalue 1e-5 \
          > blastx.outfmt6

#blastp 
blastp -query /scratch/trivittatusTrim/trinityAssembly/longest_orfs.pep \
         -db /scratch/trivittatusTrim/trinotateAssembly/swissprot -num_threads 2 \
         -max_target_seqs 1 -outfmt 6 -evalue 1e-5 \
          > blastp.outfmt6



##write conf.txt for autotrinotae pipeline

##########################################################################
# Globals. Specify resource locations and other templated parameter values
# Use format {__token__} when using token values in command strings.
# Other templated parameters are defined by the parent script.
##########################################################################


[GLOBALS]

# progs
TRANSDECODER_DIR=/cvmfs/soft.computecanada.ca/easybuild/software/2020/avx2/Compiler/gcc9/transdecoder/5.5.0/bin/
BLASTX_PROG=/scratch/trivittatusTrim/trinotateAssembly/blastx
BLASTP_PROG=/scratch/trivittatusTrim/trinotateAssembly/blastp
SIGNALP_PROG=/scratch/trivittatusTrim/trinotateAssembly/signalp
TMHMM_PROG=/scratch/trivittatusTrim/trinotateAssembly/tmhmm
RNAMMER_TRANS_PROG=/cvmfs/soft.computecanada.ca/easybuild/software/2020/avx2/Compiler/gcc9/trinotate/3.2.2/util/rnammer_support/RnammerTranscriptome.pl
RNAMMER=/cvmfs/soft.computecanada.ca/easybuild/software/2020/avx2/Core/rnammer/1.2/rnammer
HMMSCAN_PROG=hmmscan

# dbs
SWISSPROT_PEP=/scratch/trivittatusTrim/trinotateAssembly/uniprot_sprot.pep
PFAM_DB=/scratch/trivittatusTrim/trinotateAssembly/Pfam-A.hmm


####################
#  BioIfx computes 
####################

[TRANSDECODER_LONGORF]
RANK=100
RUN=T
CMD={__TRANSDECODER_DIR__}/TransDecoder.LongOrfs -t {__TRANSCRIPTS_FASTA__}


[TRANSDECODER_PREDICT]
RANK=101
RUN=T
CMD={__TRANSDECODER_DIR__}/TransDecoder.Predict -t {__TRANSCRIPTS_FASTA__} --cpu {__CPU__}


[BLASTX_SPROT_TRANS]
RANK=200
RUN=T
CMD={__BLASTX_PROG__} -db {__SWISSPROT_PEP__} -query {__TRANSCRIPTS_FASTA__} -num_threads {__CPU__} -max_target_seqs 1 -outfmt 6 -evalue 1e-5 > swissprot.blastx.outfmt6


[BLASTX_SPROT_PEP]
RANK=300
RUN=T
CMD={__BLASTP_PROG__} -query {__TRANSCRIPTS_FASTA__}.transdecoder.pep -db {__SWISSPROT_PEP__} -num_threads {__CPU__}  -max_target_seqs 1 -outfmt 6 -evalue 1e-5 > swissprot.blastp.outfmt6


[PFAM]
RANK=400
RUN=T
CMD={__HMMSCAN_PROG__} --cpu {__CPU__} --domtblout TrinotatePFAM.out {__PFAM_DB__} {__TRANSCRIPTS_FASTA__}.transdecoder.pep  > pfam.log


[SIGNALP]
RANK=500
RUN=T
CMD={__SIGNALP_PROG__} -f short -n signalp.out {__TRANSCRIPTS_FASTA__}.transdecoder.pep > sigP.log


[TMHMM]
RANK=600
RUN=T
CMD={__TMHMM_PROG__} --short < {__TRANSCRIPTS_FASTA__}.transdecoder.pep > tmhmm.out

[RNAMMER]
RANK=700
RUN=T
CMD={__RNAMMER_TRANS_PROG__} --transcriptome {__TRANSCRIPTS_FASTA__} --path_to_rnammer {__RNAMMER__} 
# generates file: {__TRANSCRIPTS_FASTA__}.rnammer.gff


#############################
# Trinotate Database Loading
#############################

[TRINOTATE_INIT]
RANK=1100
RUN=T
CMD={__TRINOTATE_HOME__}/Trinotate {__TRINOTATE_SQLITE__} init --gene_trans_map {__GENE_TO_TRANS_MAP__} --transcript_fasta {__TRANSCRIPTS_FASTA__} --transdecoder_pep {__TRANSCRIPTS_FASTA__}.transdecoder.pep

[TRINOTATE_LOAD_SPROT_BLASTX]
RANK=1200
RUN=T
CMD={__TRINOTATE_HOME__}/Trinotate  {__TRINOTATE_SQLITE__} LOAD_swissprot_blastx swissprot.blastx.outfmt6 

[TRINOTATE_LOAD_SPROT_BLASTP]
RANK=1300
RUN=T
CMD={__TRINOTATE_HOME__}/Trinotate  {__TRINOTATE_SQLITE__} LOAD_swissprot_blastp swissprot.blastp.outfmt6 


[TRINOTATE_LOAD_PFAM]
RANK=1400
RUN=T
CMD={__TRINOTATE_HOME__}/Trinotate  {__TRINOTATE_SQLITE__} LOAD_pfam TrinotatePFAM.out

[TRINOTATE_LOAD_RNAMMER]
RANK=1500
RUN=T
CMD={__TRINOTATE_HOME__}/Trinotate  {__TRINOTATE_SQLITE__} LOAD_rnammer {__TRANSCRIPTS_FASTA__}.rnammer.gff

[TRINOTATE_LOAD_TMHMM]
RANK=1600
RUN=T
CMD={__TRINOTATE_HOME__}/Trinotate  {__TRINOTATE_SQLITE__} LOAD_tmhmm tmhmm.out

[TRINOTATE_LOAD_SIGNALP]
RANK=1700
RUN=T
CMD={__TRINOTATE_HOME__}/Trinotate  {__TRINOTATE_SQLITE__} LOAD_signalp signalp.out

[TRINOTATE_REPORT]
RANK=1800
RUN=T
CMD={__TRINOTATE_HOME__}/Trinotate  {__TRINOTATE_SQLITE__} report > Trinotate.xls


[EXTRACT_GO]
RANK=1900
RUN=T
CMD={__TRINOTATE_HOME__}/util/extract_GO_assignments_from_Trinotate_xls.pl  --Trinotate_xls Trinotate.xls -T -I > Trinotate.xls.gene_ontology

[NAME_UPDATES]
RANK=2000
RUN=T
CMD={__TRINOTATE_HOME__}/util/annotation_importer/import_transcript_names.pl {__TRINOTATE_SQLITE__} Trinotate.xls

 

##create boilerplate and run the compilation to an xls file 

cd /scratch/trivittatusTrim/trinotateAssembly

module load CVMFS_CC; module load gcc; module load trinotate;
module load transdecoder; module load tmhmm; 
module load blast+; module load signalp; module load rnammer; module load hmmer

/cvmfs/soft.computecanada.ca/easybuild/software/2020/avx2/Compiler/gcc9/trinotate/3.2.2/admin/Build_Trinotate_Boilerplate_SQLite_db.pl trinTrim

gunzip uniprot_sprot.dat.gz
gunzip Pfam-A.hmm.gz ; hmmpress Pfam-A.hmm
makeblastdb -in uniprot_sprot.pep -dbtype prot

##write pbs file for batch job
nano autoTrin.pbs

#!/bin/bash

#PBS -l walltime=96:00:00,select=1:ncpus=40:mem=187gb
#PBS -N trinotate_sql
#PBS -A [account name]
#PBS -m abe
#PBS -M [cluster name]
#PBS -o output.txt
#PBS -e error.txt
 
################################################################################

cd $PBS_O_WORKDIR


module load CVMFS_CC; module load gcc; module load trinotate;
module load transdecoder; module load tmhmm; 
module load blast+; module load signalp; module load rnammer; module load hmmer

/cvmfs/soft.computecanada.ca/easybuild/software/2020/avx2/Compiler/gcc9/trinotate/3.2.2/auto/autoTrinotate.pl --Trinotate_sqlite trinTrim.sqlite --transcripts /scratch/trivittatusTrim/trinityAssembly/trinity_out_dir.Trinity.fasta --gene_to_trans_map /scratch/trivittatusTrim/trinityAssembly/trinity_out_dir.Trinity.fasta.gene_trans_map --conf conf.txt --CPU 40


##run targeted blast searches 
## make a database to query 
mkdir /scratch/trivittatusTrim/dbQuery
cd  /scratch/trivittatusTrim/dbQuery
module load CVMFS_CC; module load gcc; module load blast+

#make database from trinotate output 
#prot database from longest orfs
makeblastdb -in /scratch/trivittatusTrim/trinityAssembly/longest_orfs.pep -title chaoDB_pro -dbtype prot -out chaoDB_pro 
#nucleotide database 
makeblastdb -in /scratch/trivittatusTrim/trinityAssembly/trinity_out_dir.Trinity.fasta -title chao_nt -dbtype nt -out chao_nt


##run targeted gene experssion on some expected genes from vectorbase 
module load CVMFS_CC; module load gcc; module load blast+

##add compiled drosophila gene protein sequences into a fasta file based on expected/known epithelial genes - dMelanogaster_fullProt.fasta file 

##blast - can do interactive on cluster 
##query assembled database for these expected genes from drosophila 
blastp -query /scratch/trivittatusTrim/dbQuery/dMelanogaster_prot.fasta\
         -db /scratch/trivittatusTrim/dbQuery/chaoDB_pro -out cTriv_aaTarget.outfmt6 \
         -evalue 1e-20 -max_target_seqs 1 -outfmt 6     

#additional stats from the targeted blast
module load samtools

samtools faidx dMelanogaster_prot.fasta

##download outputs locally, compile into spreadsheet (ctriv_targetedQuery.txt file on data dryad)
##use R to plot the summary results (see R code file availale on repository)

##differential expression analysis for the three tissue types

cd /scratch/trivittatusTrim/trinityAssembly

##abundance estimations for each sample
nano salmon.pbs
#!/bin/bash

#PBS -l walltime=96:00:00,select=1:ncpus=40:mem=187gb
#PBS -N abundEst
#PBS -A [account name]
#PBS -m abe
#PBS -M [cluster name]
#PBS -o outputAE.txt
#PBS -e errorAE.txt
 
################################################################################

cd $PBS_O_WORKDIR

module load CVMFS_CC; module load gcc; module load openmpi; module load samtools; module load salmon; module load trinity

## prep the reference and run the alignment/estimation (can use same sample file used in trinity assembly)
/cvmfs/soft.computecanada.ca/easybuild/software/2020/avx2/Core/trinity/2.14.0/trinityrnaseq-v2.14.0/util/align_and_estimate_abundance.pl  --transcripts /scratch/trivittatusTrim/trinityAssembly/trinity_out_dir.Trinity.fasta --seqType fq --samples_file /arc/project/trivittatus/ctrivTX1_samplesFile.txt --est_method salmon --trinity_mode --prep_reference --output_dir salmon_outdir


nano abunEst_matrix.pbs

#!/bin/bash

#PBS -l walltime=96:00:00,select=1:ncpus=40:mem=187gb
#PBS -N abundEst
#PBS -A [account name]
#PBS -m abe
#PBS -M [cluster name]
#PBS -o outputAE2.txt
#PBS -e errorAE2.txt
 
################################################################################

cd $PBS_O_WORKDIR

module load CVMFS_CC; module load gcc/9.3.0; module load openmpi/4.0.3; module load samtools; module load salmon/1.7.0; module load trinity; module load perl; module load r; module load cmake


# get a list of the salmon quant.sf files so we don't have to list them individually
find . -maxdepth 2 -name "quant.sf" | tee salmon.quant_files.txt


# generate the abundance matrix
/cvmfs/soft.computecanada.ca/easybuild/software/2020/avx2/Core/trinity/2.14.0/trinityrnaseq-v2.14.0/util/abundance_estimates_to_matrix.pl --est_method salmon --gene_trans_map /scratch/trivittatusTrim/trinityAssembly/trinity_out_dir.Trinity.fasta.gene_trans_map --quant_files salmon.quant_files.txt --name_sample_by_basedir


##download files locally to work with in R (trivittatus_geneCounts.matrix file)
## see R file in repository for assembly, annotation, and differential expression analysis 

### Eucorethra underwoodi assembly and annotation ###

#eucorethra samples directory 
mkdir /arc/project/eucorethra

#download and rename samples 
Eu_tr1_left.fastq.gz
Eu_tr1_right.fastq.gz
Eu_mt1_left.fastq.gz
Eu_mt1_right.fastq.gz
Eu_gi1_left.fastq.gz
Eu_gi1_right.fastq.gz
Eu_tr2_left.fastq.gz
Eu_tr2_right.fastq.gz
Eu_mt2_left.fastq.gz
Eu_mt2_right.fastq.gz
Eu_gi2_left.fastq.gz
Eu_gi2_right.fastq.gz
Eu_tr3_left.fastq.gz
Eu_tr3_right.fastq.gz
Eu_mt3_left.fastq.gz
Eu_mt3_right.fastq.gz
Eu_gi3_left.fastq.gz
Eu_gi3_right.fastq.gz

##create samplesFile
nano eucorTX_samplesFile.txt

#  Or,
#      --samples_file /arc/project/eucorethra      tab-delimited text file indicating biological replicate relationships.
#                                  ex.



Eu_tr  Eu_tr_rep1 /arc/project/eucorethra/Eu_tr1_left.fastq.gz /arc/project/eucorethra/Eu_tr1_right.fastq.gz
Eu_tr  Eu_tr_rep2 /arc/project/eucorethra/Eu_tr2_left.fastq.gz /arc/project/eucorethra/Eu_tr2_right.fastq.gz    
Eu_tr  Eu_tr_rep3 /arc/project/eucorethra/Eu_tr3_left.fastq.gz /arc/project/eucorethra/Eu_tr3_right.fastq.gz
Eu_mt  Eu_mt_rep1 /arc/project/eucorethra/Eu_mt1_left.fastq.gz /arc/project/eucorethra/Eu_mt1_right.fastq.gz
Eu_mt  Eu_mt_rep2 /arc/project/eucorethra/Eu_mt2_left.fastq.gz /arc/project/eucorethra/Eu_mt2_right.fastq.gz
Eu_mt  Eu_mt_rep3 /arc/project/eucorethra/Eu_mt3_left.fastq.gz /arc/project/eucorethra/Eu_mt3_right.fastq.gz
Eu_gi  Eu_gi_rep1 /arc/project/eucorethra/Eu_gi1_left.fastq.gz /arc/project/eucorethra/Eu_gi1_right.fastq.gz
Eu_gi  Eu_gi_rep2 /arc/project/eucorethra/Eu_gi2_left.fastq.gz /arc/project/eucorethra/Eu_gi2_right.fastq.gz
Eu_gi  Eu_gi_rep3 /arc/project/eucorethra/Eu_gi3_left.fastq.gz /arc/project/eucorethra/Eu_gi3_right.fastq.gz


##make a eucorethra scratch directory for raw (untrimmed) reads   
mkdir /scratch/eucorethraRaw
mkdir /scratch/eucorethraRaw/trinityAssembly
cd /scratch/eucorethraRaw/trinityAssembly

##run the trinity assembly in this directoty as a batch job - first running assembly untrimmed 
#create a pbs file 
nano trinityEI_raw.pbs

#!/bin/bash

#PBS -l walltime=120:00:00,select=1:ncpus=32:mem=100gb
#PBS -N trinity_eucTX1
#PBS -A [account name]
#PBS -m abe
#PBS -M [cluster name]
#PBS -o output.txt
#PBS -e error.txt

################################################################################


export TRINITY_HOME=/arc/project/bin/trinityrnaseq-v2.13.2/

module load Software_Collection/2019
module load gcc/5.4.0
module load samtools
module load boost/1.54.0
module load python
module load py-virtualenv


cd $PBS_O_WORKDIR

source /scratch/trinity_env/bin/activate

sh runRaw_ei.sh


##create execution file 
nano runRaw_ei.sh
#!/bin/bash -ve

#######################################################
##  Run Trinity to Generate Transcriptome Assemblies ##
#######################################################

export PATH=/arc/project/bin/bowtie2-2.4.5/:/arc/project/bin/jellyfish-2.3.0/bin/:/arc/project/bin/salmon-1.7.0_linux_x86_64/bin/:$PATH
export TRINITY_HOME=/arc/project/bin/trinityrnaseq-v2.13.2/

if [ -z ${TRINITY_HOME} ]; then
    echo "Must set env var TRINITY_HOME"
    exit 1
fi

#parameters for assembly 
${TRINITY_HOME}/Trinity --seqType fq        \
        --CPU 32 --max_memory 100G  \
        --samples_file /arc/project/eucorethra/eucorTX_samplesFile.txt   \
        --SS_lib_type RF

##### Done Running Trinity #####

if [ $* ]; then
    # check full-length reconstruction stats:

    ${TRINITY_HOME}/util/misc/illustrate_ref_comparison.pl __indiv_ex_sample_derived/refSeqs.fa eucorethraTrinity_raw.fasta 90

    ./test_FL.sh --query  eucorethraTrinity_raw.fasta --target __indiv_ex_sample_derived/refSeqs.fa --no_reuse

fi

           


##submit the job 
qsub trinityEI_raw.pbs

##make new directory for trimmed assembly
mkdir /scratch/eucorethraTrim
mkdir /scratch/eucorethraTrim/trinityAssembly
cd /scratch/eucorethraTrim/trinityAssembly


nano trinityEI_trimmed.pbs
#!/bin/bash
 
#PBS -l walltime=120:00:00,select=4:ncpus=32:mem=100gb
#PBS -N trinity_eucorTX2
#PBS -A [account name]
#PBS -m abe
#PBS -M [cluster name]
#PBS -o output.txt
#PBS -e error.txt
 
################################################################################

export TRINITY_HOME=/arc/project/bin/trinityrnaseq-v2.13.2/

module load samtools
module load boost/1.71.0
module load python
module load py-virtualenv

 
cd $PBS_O_WORKDIR

source /scratch/trinity_env/bin/activate

sh runTrim_ei.sh

##write run file 
nano runTrim_ei.sh
#!/bin/bash -ve

#######################################################
##  Run Trinity to Generate Transcriptome Assemblies ##
#######################################################

export PATH=/arc/project/bin/bowtie2-2.4.5/:/arc/project/bin/jellyfish-2.3.0/bin/:/arc/project/bin/salmon-1.7.0_linux_x86_64/bin/:$PATH
export TRINITY_HOME=/arc/project/bin/trinityrnaseq-v2.13.2

if [ -z ${TRINITY_HOME} ]; then
    echo "Must set env var TRINITY_HOME"
    exit 1
fi


${TRINITY_HOME}/Trinity --seqType fq        \
        --CPU 32 --max_memory 100G --trimmomatic \
        --samples_file /arc/project/eucorethra/eucorTX_samplesFile.txt   \
        --SS_lib_type RF \
           

##### Done Running Trinity #####

if [ $* ]; then
    # check full-length reconstruction stats:

    ${TRINITY_HOME}/util/misc/illustrate_ref_comparison.pl __indiv_ex_sample_derived/refSeqs.fa eucorethraTrinity_trim.fasta 90

    ./test_FL.sh --query  eucorethraTrinity_trim.fasta --target __indiv_ex_sample_derived/refSeqs.fa --no_reuse

fi

##submit this job 
qsub trinityEI_trimmed.pbs


##fastqc 
mkdir /scratch/eucorethraTrim/fastqcTrim 
cd /scratch/eucorethraTrim/fastqcTrim

nano eucorethraFQ.pbs

#!/bin/bash

#PBS -l walltime=72:00:00,select=1:ncpus=18:mem=100gb
#PBS -N fastQC
#PBS -A [account name]
#PBS -m abe
#PBS -M [cluster name]
#PBS -o output.txt
#PBS -e error.txt
 
################################################################################


module load CVMFS_CC
module load fastqc 

cd $PBS_O_WORKDIR

#execute fastqc 
fastqc /project/eucorethra/*.gz -o -t 18 /scratch/eucorethraTrim/fastqcTrim/fastqc_report.txt

#back into the assembly directory for assembly stats 
cd /scratch/eucorethraTrim/trinityAssembly

module load CVMFS_CC
module load gcc 
module load trinity

/cvmfs/soft.computecanada.ca/easybuild/software/2020/avx2/Core/trinity/2.14.0/trinityrnaseq-v2.14.0/util/TrinityStats.pl trinity_out_dir.Trinity.fasta > trinityTrim_Stats.log

##get more summary stats 
#contig lengths
grep -oP '(?<=len=)(\s+)?\K([^ ]*)' trinity_out_dir.Trinity.fasta > contigLength_Read.txt

#read end bias 
cd ../transDeco/trinity_out_dir.Trinity.fasta.transdecoder_dir

#five prime partials 
grep -o "5prime" longest_orfs.cds > fivePrime_regions.txt
#three prime partials
grep -o "3prime" longest_orfs.cds > threePrime_regions.txt
##all you really need are these numbers - don't need the actual txt files, but okay to have anyway
wc -l fivePrime_regions.txt
wc -l threePrime_regions.txt



#transdecoder step
#in the eucorethraTrim directory 
mkdir /scratch/eucorethraTrim/transDeco
cd /scratch/eucorethraTrim/transdeco 

module load CVMFS_CC
module load gcc
module load transedecoder

#identify ORFs
/cvmfs/soft.computecanada.ca/easybuild/software/2020/avx2/Compiler/gcc9/transdecoder/5.5.0/bin/TransDecoder.LongOrfs -t /scratch/eucorethraTrim/trinityAssembly/trinity_out_dir.Trinity.fasta

#predict coding ORFs
/cvmfs/soft.computecanada.ca/easybuild/software/2020/avx2/Compiler/gcc9/transdecoder/5.5.0/bin/TransDecoder.Predict -t /scratch/eucorethraTrim/trinityAssembly/trinity_out_dir.Trinity.fasta 

#stats log 
module load CVMFS_CC
module load trinity

/cvmfs/soft.computecanada.ca/easybuild/software/2020/avx2/Core/trinity/2.14.0/trinityrnaseq-v2.14.0/util/TrinityStats.pl trinity_out_dir.Trinity.fasta > trinityTrim_Stats.log


##trinotate pipeline 

#copy rhe signalp exe 
#in eucorethraTrim dir 
mkdir /scratch/eucorethraTrim/trnotate
cd /scratch/eucorethraTrim/trnotate
cp signalp /scratch/eucorethraTrim/trnotate/trnotate
cd trnotate 

nano trinotate.pbs
#!/bin/bash

#PBS -l walltime=96:00:00,select=1:ncpus=40:mem=187gb
#PBS -N trinotate
#PBS -A [account name]
#PBS -m abe
#PBS -M [cluster name]
#PBS -o output.txt
#PBS -e error.txt
 
################################################################################

cd $PBS_O_WORKDIR

module load CVMFS_CC; module load gcc; module load trinotate
module load hmmer; module load blast+; module load signalp; module load rnammer
module load transdecoder; module load tmhmm 



# HMMSCAN
hmmscan --cpu 40 --domtblout TrinotatePFAM.out /scratch/trinotate/Trinotate-Trinotate-v3.2.2/Pfam-A.hmm /scratch/eucorethraTrim/transDeco/trinity_out_dir.Trinity.fasta.transdecoder_dir/longest_orfs.pep > pfam.log


# SignalP
signalp -f short -n signalp.out /scratch/eucorethraTrim/transDeco/trinity_out_dir.Trinity.fasta.transdecoder_dir/longest_orfs.pep > signalp.log

# TMHMM
tmhmm --short < /scratch/eucorethraTrim/transDeco/trinity_out_dir.Trinity.fasta.transdecoder_dir/longest_orfs.pep > tmhmm.out


#RNAmmer
/scratch/trinotate/Trinotate-Trinotate-v3.2.2/util/rnammer_support/RnammerTranscriptome.pl --transcriptome /scratch/eucorethraTrim/trinityAssembly/trinity_out_dir.Trinity.fasta --path_to_rnammer /cvmfs/soft.computecanada.ca/easybuild/software/2020/avx2/Core/rnammer/1.2/rnammer

#blastp
blastp -query /scratch/transDeco/trinity_out_dir.Trinity.fasta.transdecoder_dir/longest_orfs.pep \
         -db /scratch/eucorethraTrim/trinotateTrim/swissprot -num_threads 2 \
         -max_target_seqs 1 -outfmt 6 -evalue 1e-5 \
          > blastp.outfmt6

#blastx
blastx -query /scratch/eucorethraTrim/transDeco/trinity_out_dir.Trinity.fasta.transdecoder_dir/longest_orfs.pep \
         -db /scratch/eucorethraTrim/trinotateTrim/swissprot -num_threads 2 \
         -max_target_seqs 1 -outfmt 6 -evalue 1e-5 \
          > blastx.outfmt6

signalp -f short -n signalp.out /scratch/eucorethraTrim/transDeco/trinity_out_dir.Trinity.fasta.transdecoder_dir/longest_orfs.pep > sigP.log


nano autoTrin.pbs
#!/bin/bash

#PBS -l walltime=96:00:00,select=1:ncpus=40:mem=187gb
#PBS -N autoTrin
#PBS -A [account name]
#PBS -m abe
#PBS -M [cluster name]
#PBS -o outputS.txt
#PBS -e errorS.txt

################################################################################

cd $PBS_O_WORKDIR

#load our modules from the compute Canada stack
module load CVMFS_CC; module load gcc/9.3.0; module load trinotate/3.2.2;
module load transdecoder; module load tmhmm;
module load blast+; module load signalp; module load rnammer; module load hmmer

# make the boilerplate sql database
/cvmfs/soft.computecanada.ca/easybuild/software/2020/avx2/Compiler/gcc9/trinotate/3.2.2/admin/Build_Trinotate_Boilerplate_SQLite_db.pl trimReads

# prepare the pfam and sprot databases
gunzip uniprot_sprot.dat.gz
gunzip Pfam-A.hmm.gz ; hmmpress Pfam-A.hmm
makeblastdb -in uniprot_sprot.pep -dbtype prot

# run the autoTrinotate pipeline: note! you must have updated the conf.txt file to match the sockeye paths etc.
/cvmfs/soft.computecanada.ca/easybuild/software/2020/avx2/Compiler/gcc9/trinotate/3.2.2/auto/autoTrinotate.pl --Trinotate_sqlite trimReads.sqlite --transcripts /scratch/eucorethraTrim/trinityAssembly/trinity_out_dir.Trinity.fasta --gene_to_trans_map /scratch/eucorethraTrim/trinityAssembly/trinity_out_dir.Trinity.fasta.gene_trans_map --conf conf.txt --CPU 40

 
##write conf.txt 
nano conf.txt

##########################################################################
# Globals. Specify resource locations and other templated parameter values
# Use format {__token__} when using token values in command strings.
# Other templated parameters are defined by the parent script.
##########################################################################


[GLOBALS]

# progs
TRANSDECODER_DIR=/cvmfs/soft.computecanada.ca/easybuild/software/2020/avx2/Compiler/gcc9/transdecoder/5.5.0/bin/
BLASTX_PROG=/cvmfs/soft.computecanada.ca/easybuild/software/2020/avx2/Compiler/gcc9/blast+/2.13.0/bin/blastx
BLASTP_PROG=/cvmfs/soft.computecanada.ca/easybuild/software/2020/avx2/Compiler/gcc9/blast+/2.13.0/bin/blastp
SIGNALP_PROG=/scratch/eucorethraTrim/trnotate/signalp
TMHMM_PROG=/cvmfs/soft.computecanada.ca/easybuild/software/2020/Core/tmhmm/2.0c/bin/tmhmm 
RNAMMER_TRANS_PROG=/cvmfs/soft.computecanada.ca/easybuild/software/2020/avx2/Compiler/gcc9/trinotate/3.2.2/util/rnammer_support/RnammerTranscriptome.pl
RNAMMER=/cvmfs/soft.computecanada.ca/easybuild/software/2020/avx2/Core/rnammer/1.2/rnammer
HMMSCAN_PROG=/cvmfs/soft.computecanada.ca/easybuild/software/2020/avx2/Core/hmmer/3.2.1/bin/hmmscan
TRINOTATE_HOME=/cvmfs/soft.computecanada.ca/easybuild/software/2020/avx2/Compiler/gcc9/trinotate/3.2.2

# dbs
SWISSPROT_PEP=/scratch/eucorethraTrim/trnotate/uniprot_sprot.pep
PFAM_DB=/scratch/eucorethraTrim/trnotate/Pfam-A.hmm


####################
#  BioIfx computes 
####################

[TRANSDECODER_LONGORF]
RANK=100
RUN=T
CMD={__TRANSDECODER_DIR__}/TransDecoder.LongOrfs -t {__TRANSCRIPTS_FASTA__}


[TRANSDECODER_PREDICT]
RANK=101
RUN=T
CMD={__TRANSDECODER_DIR__}/TransDecoder.Predict -t {__TRANSCRIPTS_FASTA__} --cpu {__CPU__}


[BLASTX_SPROT_TRANS]
RANK=200
RUN=T
CMD={__BLASTX_PROG__} -db {__SWISSPROT_PEP__} -query {__TRANSCRIPTS_FASTA__} -num_threads {__CPU__} -max_target_seqs 1 -outfmt 6 -evalue 1e-5 > swissprot.blastx.outfmt6


[BLASTX_SPROT_PEP]
RANK=300
RUN=T
CMD={__BLASTP_PROG__} -query {__TRANSCRIPTS_FASTA__}.transdecoder.pep -db {__SWISSPROT_PEP__} -num_threads {__CPU__}  -max_target_seqs 1 -outfmt 6 -evalue 1e-5 > swissprot.blastp.outfmt6


[PFAM]
RANK=400
RUN=T
CMD={__HMMSCAN_PROG__} --cpu {__CPU__} --domtblout TrinotatePFAM.out {__PFAM_DB__} {__TRANSCRIPTS_FASTA__}.transdecoder.pep  > pfam.log


[SIGNALP]
RANK=500
RUN=T
CMD={__SIGNALP_PROG__} -f short -n signalp.out {__TRANSCRIPTS_FASTA__}.transdecoder.pep > sigP.log


[TMHMM]
RANK=600
RUN=T
CMD={__TMHMM_PROG__} --short < {__TRANSCRIPTS_FASTA__}.transdecoder.pep > tmhmm.out

[RNAMMER]
RANK=700
RUN=T
CMD={__RNAMMER_TRANS_PROG__} --transcriptome {__TRANSCRIPTS_FASTA__} --path_to_rnammer {__RNAMMER__} 
# generates file: {__TRANSCRIPTS_FASTA__}.rnammer.gff


#############################
# Trinotate Database Loading
#############################

[TRINOTATE_INIT]
RANK=1100
RUN=T
CMD={__TRINOTATE_HOME__}/Trinotate {__TRINOTATE_SQLITE__} init --gene_trans_map {__GENE_TO_TRANS_MAP__} --transcript_fasta {__TRANSCRIPTS_FASTA__} --transdecoder_pep {__TRANSCRIPTS_FASTA__}.transdecoder.pep

[TRINOTATE_LOAD_SPROT_BLASTX]
RANK=1200
RUN=T
CMD={__TRINOTATE_HOME__}/Trinotate  {__TRINOTATE_SQLITE__} LOAD_swissprot_blastx swissprot.blastx.outfmt6 

[TRINOTATE_LOAD_SPROT_BLASTP]
RANK=1300
RUN=T
CMD={__TRINOTATE_HOME__}/Trinotate  {__TRINOTATE_SQLITE__} LOAD_swissprot_blastp swissprot.blastp.outfmt6 


[TRINOTATE_LOAD_PFAM]
RANK=1400
RUN=T
CMD={__TRINOTATE_HOME__}/Trinotate  {__TRINOTATE_SQLITE__} LOAD_pfam TrinotatePFAM.out

[TRINOTATE_LOAD_RNAMMER]
RANK=1500
RUN=T
CMD={__TRINOTATE_HOME__}/Trinotate  {__TRINOTATE_SQLITE__} LOAD_rnammer {__TRANSCRIPTS_FASTA__}.rnammer.gff

[TRINOTATE_LOAD_TMHMM]
RANK=1600
RUN=T
CMD={__TRINOTATE_HOME__}/Trinotate  {__TRINOTATE_SQLITE__} LOAD_tmhmm tmhmm.out

[TRINOTATE_LOAD_SIGNALP]
RANK=1700
RUN=T
CMD={__TRINOTATE_HOME__}/Trinotate  {__TRINOTATE_SQLITE__} LOAD_signalp signalp.out

[TRINOTATE_REPORT]
RANK=1800
RUN=T
CMD={__TRINOTATE_HOME__}/Trinotate  {__TRINOTATE_SQLITE__} report > Trinotate.xls


[EXTRACT_GO]
RANK=1900
RUN=T
CMD={__TRINOTATE_HOME__}/util/extract_GO_assignments_from_Trinotate_xls.pl  --Trinotate_xls Trinotate.xls -T -I > Trinotate.xls.gene_ontology

[NAME_UPDATES]
RANK=2000
RUN=T
CMD={__TRINOTATE_HOME__}/util/annotation_importer/import_transcript_names.pl {__TRINOTATE_SQLITE__} Trinotate.xls


##salmon isoform and gene counts 
mkdir /scratch/eucorethraTrim/difExp 
cd /scratch/eucorethraTrim/difExp

nano salmon.pbs 
#!/bin/bash

#PBS -l walltime=96:00:00,select=1:ncpus=40:mem=187gb
#PBS -N abundEst
#PBS -A [account name]
#PBS -m abe
#PBS -M [cluster name]
#PBS -o output.txt
#PBS -e error.txt
 
################################################################################

cd $PBS_O_WORKDIR

module load CVMFS_CC; module load gcc; module load openmpi; module load samtools; module load salmon; module load trinity

## prep the reference and run the alignment/estimation (can use same sample file used in trinity assembly)
/cvmfs/soft.computecanada.ca/easybuild/software/2020/avx2/Core/trinity/2.14.0/trinityrnaseq-v2.14.0/util/align_and_estimate_abundance.pl  --transcripts /scratch/eucorethraTrim/trinityAssembly/trinity_out_dir.Trinity.fasta --seqType fq --samples_file /project/eucorethra/eucorTX_samplesFile.txt --est_method salmon --trinity_mode --prep_reference --output_dir salmon_outdir


nano abunEst_matrix.pbs

#!/bin/bash

#PBS -l walltime=96:00:00,select=1:ncpus=40:mem=187gb
#PBS -N abundEst
#PBS -A [account name]
#PBS -m abe
#PBS -M [cluster name]
#PBS -o output2.txt
#PBS -e error2.txt
 
################################################################################

cd $PBS_O_WORKDIR

module load CVMFS_CC; module load gcc/9.3.0; module load openmpi/4.0.3; module load samtools; module load salmon/1.7.0; module load trinity; module load perl; module load r; module load cmake


# get a list of the salmon quant.sf files so we don't have to list them individually
find . -maxdepth 2 -name "quant.sf" | tee salmon.quant_files.txt


# generate the abundance matrix
/cvmfs/soft.computecanada.ca/easybuild/software/2020/avx2/Core/trinity/2.14.0/trinityrnaseq-v2.14.0/util/abundance_estimates_to_matrix.pl --est_method salmon --gene_trans_map /scratch/eucorethraTrim/trinityAssembly/trinity_out_dir.Trinity.fasta.gene_trans_map --quant_files salmon.quant_files.txt --name_sample_by_basedir

#targeted blast search
#make database 
mkdir /scratch/eucorethraTrim/blasts 
cd /scratch/eucorethraTrim/blasts 

module load CVMFS_CC 
module load gcc 
module load blast+

#nucleotide database
makeblastdb -in /scratch/eucorethraTrim/transDeco/trinity_out_dir.Trinity.fasta.transdecoder_dir/longest_orfs.cds -title eUnd_nlDB -dbtype nucl -out eUnd_nlDB 

#protein databade 
makeblastdb -in /scratch/eucorethraTrim/transDeco/trinity_out_dir.Trinity.fasta.transdecoder_dir/longest_orfs.pep -title eUnd_prDB -dbtype prot -out eUnd_prDB 

#the query files are the same ones we used for the chaoborus trivittatus assembly so point to those
##query assembled database for these expected genes from drosophila 
blastp -query /scratch/trinotateTrim/dMelanogaster_prot.fasta \
         -db /scratch/eucorethraTrim/blasts/eUnd_pro -out eucor_aaTarget.outfmt6 \
         -evalue 1e-20 -max_target_seqs 1 -outfmt 6


