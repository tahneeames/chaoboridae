
#based on the assembly statistics it seems like the trimmed assembly is the way to go - move forward with Trinotate on this assembly 
cd chaoborus-trimmed
#job script 
nano trinotate_Ctrimmed.pbs
#!/bin/bash

#PBS -l walltime=24:00:00,select=1:ncpus=12:mem=48gb
#PBS -N trinotate
#PBS -A #alloc-code
#PBS -m abe
#PBS -M #contact
#PBS -o 20230217-trinotate-output.txt
#PBS -e 20230217-trinotate-error.txt

################################################################################


###### load the compute canada stack, and then load the canu module
module load CVMFS_CC
module load gcc/9.3.0
module load trinotate/3.2.2
module load tmhmm
module load signalp
module load rnammer
module load blast+

### copy over working files
cd $PBS_O_WORKDIR

if [[ ! -f results.tgz ]]
then
        echo "results.tgz does not exist on your filesystem."
        echo "zip and send!"

        tar -czhf results.tgz *
fi

cp $PBS_O_WORKDIR/results.tgz $TMPDIR


BASE_tx=chaoborus.trimmed.Trinity.fasta
CPU=12

### change to local storage on compute node and extract working files
cd $TMPDIR
tar -xzf ./results.tgz
rm ./results.tgz

### print out tmpdir for reference, and let's go!
echo $TMPDIR


#### STEP 1 ##### transdecoder

if [[ ! -f ${BASE_tx}.transdecoder.pep ]]
then
        echo "<file> does not exist on your filesystem."
        echo "Time to do Transdecoder!"

        TransDecoder.LongOrfs -t ${BASE_tx}
        TransDecoder.Predict -t ${BASE_tx} --cpu ${CPU}

fi


#### STEP 2 #### PFAM

if [[ ! -f pfam.log ]]
then
        echo "<file> does not exist on your filesystem."
        echo "Time to do PFAM/HMMScan!"

        hmmscan --cpu ${CPU} --domtblout TrinotatePFAM.out \
                /lab_data/db/Pfam/Pfam-A.hmm \
                ${BASE_tx}.transdecoder.pep > pfam.log
fi



#### STEP 3 #### SignalP
if [[ ! -f sigP.log ]]
then
        echo "<file> does not exist on your filesystem."
        echo "Time to do signalP!"

        /bin/signalp -f short -n signalp.out \
                 ${BASE_tx}.transdecoder.pep > sigP.log

fi


#### STEP 4 #### tmhmm

if [[ ! -f tmhmm.out ]]
then
        echo "<file> does not exist on your filesystem."
        echo "Time to do tmhmm!"

        tmhmm --short  ${BASE_tx}.transdecoder.pep > tmhmm.out

fi


#### STEP 5 #### RNAmmer

if [[ ! -f eucorethra.trimmed.Trinity.fasta.rnammer.gffffffff ]]
then
        echo "<file> does not exist on your filesystem."
        echo "Time to do rnammer"

        /cvmfs/soft.computecanada.ca/easybuild/software/2020/avx2/Compiler/gcc9/trinotate/3.2.2/util/rnammer_support/RnammerTranscriptome.pl --transcriptome  ${BASE_tx} \
                --path_to_rnammer /arc/project/st-benmat01-1/bin/rnammer-1.2/rnammer
fi



###### if you use the local disk for computing, you *need to copy* your results files when you're done
tar -czf $PBS_O_WORKDIR/results.tgz ./*



#same with the Moxhlonyx assembly 
#move to assembly directory 
cd mochlonyx-trimmed
#job script 
nano trinotate_Mtrimmed.pbs
#!/bin/bash

#PBS -l walltime=24:00:00,select=1:ncpus=12:mem=48gb
#PBS -N trinotate
#PBS -A #alloc-code
#PBS -m abe
#PBS -M #contact
#PBS -o 2023xxxx-trinotate-output.txt
#PBS -e 2023xxxx-trinotate-error.txt

################################################################################


###### load the compute canada stack, and then load the canu module
module load CVMFS_CC
module load gcc/9.3.0
module load trinotate/3.2.2
module load tmhmm
module load signalp
module load rnammer
module load blast+

### copy over working files
cd $PBS_O_WORKDIR

if [[ ! -f results.tgz ]]
then
        echo "results.tgz does not exist on your filesystem."
        echo "zip and send!"

        tar -czhf results.tgz *
fi

cp $PBS_O_WORKDIR/results.tgz $TMPDIR


BASE_tx=mochlonyx.trimmed.Trinity.fasta
CPU=12

### change to local storage on compute node and extract working files
cd $TMPDIR
tar -xzf ./results.tgz
rm ./results.tgz

### print out tmpdir for reference, and let's go!
echo $TMPDIR


#### STEP 1 ##### transdecoder

if [[ ! -f ${BASE_tx}.transdecoder.pep ]]
then
        echo "<file> does not exist on your filesystem."
        echo "Time to do Transdecoder!"

        TransDecoder.LongOrfs -t ${BASE_tx}
        TransDecoder.Predict -t ${BASE_tx} --cpu ${CPU}

fi


#### STEP 2 #### PFAM

if [[ ! -f pfam.log ]]
then
        echo "<file> does not exist on your filesystem."
        echo "Time to do PFAM/HMMScan!"

        hmmscan --cpu ${CPU} --domtblout TrinotatePFAM.out \
                /lab_data/db/Pfam/Pfam-A.hmm \
                ${BASE_tx}.transdecoder.pep > pfam.log
fi



#### STEP 3 #### SignalP
if [[ ! -f sigP.log ]]
then
        echo "<file> does not exist on your filesystem."
        echo "Time to do signalP!"

        /bin/signalp -f short -n signalp.out \
                 ${BASE_tx}.transdecoder.pep > sigP.log

fi


#### STEP 4 #### tmhmm

if [[ ! -f tmhmm.out ]]
then
        echo "<file> does not exist on your filesystem."
        echo "Time to do tmhmm!"

        tmhmm --short  ${BASE_tx}.transdecoder.pep > tmhmm.out

fi


#### STEP 5 #### RNAmmer

if [[ ! -f eucorethra.trimmed.Trinity.fasta.rnammer.gffffffff ]]
then
        echo "<file> does not exist on your filesystem."
        echo "Time to do rnammer"

        /cvmfs/soft.computecanada.ca/easybuild/software/2020/avx2/Compiler/gcc9/trinotate/3.2.2/util/rnammer_support/RnammerTranscriptome.pl --transcriptome  ${BASE_tx} \
                --path_to_rnammer /arc/project/st-benmat01-1/bin/rnammer-1.2/rnammer
fi



###### if you use the local disk for computing, you *need to copy* your results files when you're done
tar -czf $PBS_O_WORKDIR/results.tgz ./*




#Same with Eucorethra - move to assembly directory 
cd eucorethra-trimmed
#job script 
nano trinotate_Eirimmed.pbs
#!/bin/bash

#PBS -l walltime=24:00:00,select=1:ncpus=12:mem=48gb
#PBS -N trinotate
#PBS -A #alloc-code
#PBS -m abe
#PBS -M #contact
#PBS -o 20230217-trinotate-output.txt
#PBS -e 20230217-trinotate-error.txt

################################################################################


###### load the compute canada stack, and then load the canu module
module load CVMFS_CC
module load gcc/9.3.0
module load trinotate/3.2.2
module load tmhmm
module load signalp
module load rnammer
module load blast+

### copy over working files
cd $PBS_O_WORKDIR

if [[ ! -f results.tgz ]]
then
        echo "results.tgz does not exist on your filesystem."
        echo "zip and send!"

        tar -czhf results.tgz *
fi

cp $PBS_O_WORKDIR/results.tgz $TMPDIR


BASE_tx=eucorethra.trimmed.Trinity.fasta
CPU=12

### change to local storage on compute node and extract working files
cd $TMPDIR
tar -xzf ./results.tgz
rm ./results.tgz

### print out tmpdir for reference, and let's go!
echo $TMPDIR


#### STEP 1 ##### transdecoder

if [[ ! -f ${BASE_tx}.transdecoder.pep ]]
then
        echo "<file> does not exist on your filesystem."
        echo "Time to do Transdecoder!"

        TransDecoder.LongOrfs -t ${BASE_tx}
        TransDecoder.Predict -t ${BASE_tx} --cpu ${CPU}

fi


#### STEP 2 #### PFAM

if [[ ! -f pfam.log ]]
then
        echo "<file> does not exist on your filesystem."
        echo "Time to do PFAM/HMMScan!"

        hmmscan --cpu ${CPU} --domtblout TrinotatePFAM.out \
                /lab_data/db/Pfam/Pfam-A.hmm \
                ${BASE_tx}.transdecoder.pep > pfam.log
fi



#### STEP 3 #### SignalP
if [[ ! -f sigP.log ]]
then
        echo "<file> does not exist on your filesystem."
        echo "Time to do signalP!"

        /bin/signalp -f short -n signalp.out \
                 ${BASE_tx}.transdecoder.pep > sigP.log

fi


#### STEP 4 #### tmhmm

if [[ ! -f tmhmm.out ]]
then
        echo "<file> does not exist on your filesystem."
        echo "Time to do tmhmm!"

        tmhmm --short  ${BASE_tx}.transdecoder.pep > tmhmm.out

fi


#### STEP 5 #### RNAmmer

if [[ ! -f eucorethra.trimmed.Trinity.fasta.rnammer.gffffffff ]]
then
        echo "<file> does not exist on your filesystem."
        echo "Time to do rnammer"

        /cvmfs/soft.computecanada.ca/easybuild/software/2020/avx2/Compiler/gcc9/trinotate/3.2.2/util/rnammer_support/RnammerTranscriptome.pl --transcriptome  ${BASE_tx} \
                --path_to_rnammer /arc/project/st-benmat01-1/bin/rnammer-1.2/rnammer
fi



###### if you use the local disk for computing, you *need to copy* your results files when you're done
tar -czf $PBS_O_WORKDIR/results.tgz ./*


