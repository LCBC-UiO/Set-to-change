#!/bin/bash

#SBATCH -J FS732_base_S2C
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=8G
#SBATCH --time=24:00:00
#SBATCH --account=p23_lcbc
##SBATCH --partition=accel
#SBATCH --output logs_freesurfer/slurm-%j.txt

module purge
echo "LOADING SINGULARITY MODULE"
module load singularity/3.7.3
echo `which singularity`
#set -o errexit
unset PYTHONPATH

sub=${1}
fs_dir=${2}
tmp_dir=${3}
logname=${4}

mv logs_freesurfer/slurm-${SLURM_JOBID}.txt ${logname}

rm -r ${tmp_dir}
mkdir -p ${tmp_dir}
export TMPDIR=${tmp_dir}
env | grep TMPDIR

export FREESURFER_HOME=/cluster/projects/p23/tools/mri/freesurfer/freesurfer.7.3.2
source $FREESURFER_HOME/SetUpFreeSurfer.sh
export SUBJECTS_DIR=${fs_dir}/${sub}

rm -rf ${SUBJECTS_DIR}/fsaverage
ln -s $FREESURFER_HOME/subjects/fsaverage ${SUBJECTS_DIR}

echo ""
echo "####################################################"
echo "                  FS 7.3.2 - BASE                   "
echo "####################################################"

sessions=$(ls -d ${fs_dir}/${sub}/ses-??)
time_point_flag=""
for path_ses in ${sessions}; do
    ses=$(basename ${path_ses})
    time_point_flag+=" -tp ${ses}"
done

echo ${time_point_flag}
recon-all -base base ${time_point_flag} -all
