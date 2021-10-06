#!/bin/bash
#SBATCH --job-name=PACO       
#SBATCH --ntasks=1 
#SBATCH --cpus-per-task=12
#SBATCH --mem-per-cpu=10G
#SBATCH --time=240:00:00
#SBATCH --mail-user=alexander.piper@agriculture.vic.gov.au
#SBATCH --mail-type=ALL
#SBATCH --account=pathogens
#SBATCH --export=none

#Check if job is launched as array
if [ -z "$SLURM_ARRAY_TASK_COUNT" ]; then 
  echo SLURM_ARRAY_TASK_COUNT unset; 
  echo You must launch this job as an array
  echo see https://slurm.schedmd.com/job_array.html
  echo for info on how to run arrays
  exit 1
fi
Index=sequence_index.txt


# Make sure that sequence index file is there before we do anything
if [[ ! -f "${Index}" ]]; then
  echo "Error sequence index file ${Index} does not exist"
  exit 1
fi

#Get necessary files
JobName=$(sed -n ${SLURM_ARRAY_TASK_ID}p ${Index})
JobPath=$(dirname ${JobName})
Sample=$(basename ${JobName} .R)

# Double check that array index is valid
if [[ ! -f "${JobName}" ]]; then
  echo "Error array index doesnt match up with index file"
  echo "Array index is  ${SLURM_ARRAY_TASK_ID}"
  exit 1
fi

# Goto tmp to do our processing
echo $TMPDIR
cd $TMPDIR
tmp_dir=$(mktemp -d -t ci-XXXXXXXXXX)
cd $tmp_dir

#Copy data files, database and decompress all
cp ${JobName} .
cp /group/pathogens/Alexp/Metabarcoding/Psyllid_microbiome/paco_para/psyllid_beast_tree.nwk .
cp /group/pathogens/Alexp/Metabarcoding/Psyllid_microbiome/paco_para/ps2.rds .
cp /group/pathogens/Alexp/Metabarcoding/Psyllid_microbiome/paco_para/phytree.nwk .
cp /group/pathogens/Alexp/Metabarcoding/Psyllid_microbiome/paco_para/aligned_asv.fasta.gz .
cp /group/pathogens/Alexp/Metabarcoding/Psyllid_microbiome/paco_para/COI.csv .
pwd
ls

#Load modules
module purge
module load R/4.1.0-foss-2021a

###START###

echo ${Sample}

Rscript ${JobName} 

date

# put all output files back where we started
cp *.rds ${SLURM_SUBMIT_DIR}