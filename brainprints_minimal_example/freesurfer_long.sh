#!/bin/bash

#SBATCH -J FS732_long_S2C
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=6G
#SBATCH --time=12:00:00
#SBATCH --account=p23
##SBATCH --partition=accel
#SBATCH --output logs_freesurfer/slurm-%j.txt

module purge
echo "LOADING SINGULARITY MODULE"
module load singularity/3.7.3
echo `which singularity`
#set -o errexit
unset PYTHONPATH

sub=${1}
ses=${2}
fs_dir=${3}
tmp_dir=${4}
logname=${5}

mv logs_freesurfer/slurm-${SLURM_JOBID}.txt ${logname}

rm -r ${tmp_dir}
mkdir -p ${tmp_dir}
export TMPDIR=${tmp_dir}
env | grep TMPDIR

export FREESURFER_HOME=/cluster/projects/p23/tools/mri/freesurfer/freesurfer.7.3.2
source $FREESURFER_HOME/SetUpFreeSurfer.sh
export SUBJECTS_DIR=${fs_dir}/${sub}

echo ""
echo "####################################################"
echo "                  FS 7.3.2 - LONG                   "
echo "####################################################"

# NOTE: important to use the raw T2 from the default stream as the long-stream automatically will 
# load and apply the .ltas coregistering this to default-orig.mgz + to base template. If not using 
# T2raw (ie. using something further down the pipeline, data would potentially be coregd twice)
path_T2w="${fs_dir}/${sub}/${ses}/mri/orig/T2raw.mgz"
if [ ! -f ${path_T2w} ]; then
    echo "WARNING: CANNOT FIND T2W-SPACE, RUNNING RECON-ALL WITHOUT T2pial"
    recon-all -long ${ses} base -hires -all
else
    recon-all -long ${ses} base -T2 ${path_T2w} -T2pial -hires -all
fi
