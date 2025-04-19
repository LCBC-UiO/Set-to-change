#!/bin/bash

#SBATCH -J FS732_def_S2C
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4G
#SBATCH --time=48:00:00
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
ses=${2}
bids_dir=${3}
fs_dir=${4}
tmp_dir=${5}
logname=${6}

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
echo "                    FS 7.3.2                        "
echo "####################################################"

path_T1w=$(ls ${bids_dir}/${sub}/${ses}/anat/*acq-mprage_T1w.nii.gz)
path_T2w=$(ls ${bids_dir}/${sub}/${ses}/anat/*acq-space_T2w.nii.gz)

# MANUAL FIX FOR PROBLEMATIC SUBJECT(S)
if [[ "${sub}" == "sub-1750721" && "${ses}" == "ses-03" ]]; then
    echo "WARNING: mprage for ${sub} --> ${ses} is slightly off in terms of voxel size. Conforming to 0.8mm isotropic"
    mri_convert \
        --conform_size 0.8 -rt cubic \
        -rl ${bids_dir}/${sub}/ses-01/anat/${sub}_ses-01_acq-mprage_T1w.nii.gz \
        ${path_T1w} \
        ${bids_dir}/${sub}/${ses}/extra_data/${sub}_${ses}_acq-mprageCONFORMED_T1w.nii.gz
            
    path_T1w=${bids_dir}/${sub}/${ses}/extra_data/${sub}_${ses}_acq-mprageCONFORMED_T1w.nii.gz
fi

run_recon_with_t2=true

if [ -z ${path_T2w} ]; then
    echo "WARNING: CANNOT FIND T2W-SPACE, RUNNING RECON-ALL WITHOUT T2pial"
    run_recon_with_t2=false
fi

if [ -z ${path_T1w} ]; then
    echo "ERROR: MISSING MPRAGE"
    exit 1
fi

if [ $(echo ${path_T1w} | wc -w) -gt 1 ]; then
    echo "ERROR: FOUND MORE THAN ONE MPRAGE"
    echo "PROOF: $(echo ${path_T1w})"
    exit 1
fi

if [ $(echo ${path_T2w} | wc -w) -gt 1 ]; then
    echo "ERROR: FOUND MORE THAN ONE T2-SPACE"
    echo "PROOF: $(echo ${path_T2w})"
    exit 1
fi

mkdir -p ${fs_dir}/${sub}

if [ ${run_recon_with_t2} = true ]; then
    recon-all \
        -all \
        -s ${fs_dir}/${sub}/${ses} \
        -i ${path_T1w} \
        -T2 ${path_T2w} \
        -T2pial \
        -hires
else
    recon-all \
        -all \
        -s ${fs_dir}/${sub}/${ses} \
        -i ${path_T1w} \
        -hires
fi