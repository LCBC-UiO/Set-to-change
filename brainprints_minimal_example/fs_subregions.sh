#!/bin/bash

#SBATCH -J S2C_fssubregions
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=8G
#SBATCH --time=12:00:00
#SBATCH --account=p23
##SBATCH --partition=accel
#SBATCH --output logs_fssubregions/slurm-%j.txt

module purge
module load MCR/R2019b.8
#set -o errexit
unset PYTHONPATH

sub=${1}
logname=${2}
tmp_dir=${3}

mv logs_fssubregions/slurm-${SLURM_JOBID}.txt ${logname}

project_dir="/ess/p23/cluster/data/bids_fmri/s2c_gridtask_prisma"
bids_dir="${project_dir}/bids"
fs_dir="${bids_dir}/derivatives/freesurfer_7.3.2"

rm -r ${tmp_dir}
mkdir -p ${tmp_dir}
export TMPDIR=${tmp_dir}
env | grep TMPDIR

# Loop over sessions
sessions=$(ls -d ${fs_dir}/${sub}/ses-??)
for path_ses in ${sessions}; do
    ses=$(basename ${path_ses})

    echo ""
    echo "SESSION: ${ses}"
    echo ""

    echo ""
    echo "###############################################################"
    echo "segment_subregions (FS 7.3.2), patched version (core.py) by MHS"
    echo "###############################################################"

    for pipeline in "long" "long_alt1" "long_alt2" "cross"; do
        
        if [ ${pipeline} == "long" ]; then
            sesdir=${ses}.long.base
            seg_argument="--long-base base"
        elif [ ${pipeline} == "long_alt1" ]; then         
            sesdir=${ses}.long.base_alt1
            seg_argument="--long-base base_alt1"     
        elif [ ${pipeline} == "long_alt2" ]; then         
            sesdir=${ses}.long.base_alt2
            seg_argument="--long-base base_alt2"                     
        elif [ ${pipeline} == "cross" ]; then
            sesdir=${ses}
            seg_argument="--cross ${ses}"
        fi
              
        thalseg_output=${fs_dir}/${sub}/${sesdir}/mri/ThalamicNuclei*.txt #wildcard as is differs between cross and long pipeline
        if ! ls ${thalseg_output} 1> /dev/null 2>&1; then

            echo "${sub} --> ${ses} --> ${pipeline} RUNNING SEGMENT_SUBREGIONS (THALAMUS)"

            singularity exec --cleanenv --contain \
            -B ${tmp_dir} \
            -B ${fs_dir} \
            -B freesurfer_license.txt:/usr/local/freesurfer/.license \
            -B core_patched_MHS.py:/usr/local/freesurfer/python/packages/freesurfer/subregions/core.py \
            --env TMPDIR=${tmp_dir} \
            --env SUBJECTS_DIR=${fs_dir}/${sub} \
                freesurfer_7.3.2.sif \
                    segment_subregions \
                        thalamus \
                        ${seg_argument} \
                        --threads 4
                                             
        else
            echo "${sub} --> ${ses} --> ${pipeline} SEGMENT_SUBREGIONS (THALAMUS) ALREADY RUN"
            echo -e "\t proof: $(ls ${thalseg_output})"
        fi
    done

    echo ""
    echo "####################################################"
    echo "mri_sclimbic_seg (FS 7.3.2), patched version by MHS"
    echo "####################################################"


    for pipeline in "long" "long_alt1" "long_alt2" "cross"; do
    
        if [ ${pipeline} == "long" ]; then
            sesdir=${ses}.long.base
        elif [ ${pipeline} == "long_alt1" ]; then
            sesdir=${ses}.long.base_alt1
        elif [ ${pipeline} == "long_alt2" ]; then
            sesdir=${ses}.long.base_alt2            
        elif [ ${pipeline} == "cross" ]; then
            sesdir=${ses}
        fi
        
        path_to_t1w=${fs_dir}/${sub}/${sesdir}/mri/nu.mgz
        path_to_conformed_t1w=${fs_dir}/${sub}/${sesdir}/mri/nu.1mm.mgz
        path_to_sclimbic_output=${fs_dir}/${sub}/${sesdir}/mri/sclimbic.mgz

        if [ ! -f "${path_to_sclimbic_output}" ]; then

            echo "${sub} --> ${ses} --> ${pipeline} RUNNING MRI_SCLIMBIC_SEG"

            singularity exec --cleanenv --contain \
            -B ${tmp_dir} \
            -B ${fs_dir} \
            -B ${project_dir}/scripts/freesurfer/mri_sclimbic_seg_patched_MHS:/usr/local/freesurfer/python/scripts/mri_sclimbic_seg \
            -B freesurfer_license.txt:/usr/local/freesurfer/.license \
            --env TMPDIR=${tmp_dir} \
            --env SUBJECTS_DIR=${fs_dir}/${sub} \
                freesurfer_7.3.2.sif \
                    /bin/bash -c \
                    "mri_convert \
                        --conform \
                        -rt cubic \
                        ${path_to_t1w} \
                        ${path_to_conformed_t1w} \
                    && \
                    mri_sclimbic_seg \
                        --i "${path_to_conformed_t1w}" \
                        --o "${path_to_sclimbic_output}" \
                        --write_volumes \
                        --write_qa_stats \
                        --threads 1"
        else
            echo "${sub} --> ${ses} --> ${pipeline} MRI_SCLIMBIC_SEG ALREADY RUN"
            echo -e "\t proof: $(ls ${path_to_sclimbic_output})"
        fi
        
    done


    echo ""
    echo "#################################"
    echo "hippocampal subfields (FS 7.3.2)"
    echo "#################################"

    for pipeline in "long" "long_alt1" "long_alt2" "cross"; do
    
        if [ ${pipeline} == "long" ]; then
            sesdir=${ses}.long.base
            lta_ses="--lta ${fs_dir}/${sub}/base/mri/transforms/${ses}_to_base.lta"
        elif [ ${pipeline} == "long_alt1" ]; then
            sesdir=${ses}.long.base_alt1
            lta_ses="--lta ${fs_dir}/${sub}/base_alt1/mri/transforms/${ses}_to_base_alt1.lta"
        elif [ ${pipeline} == "long_alt2" ]; then
            sesdir=${ses}.long.base_alt2
            lta_ses="--lta ${fs_dir}/${sub}/base_alt2/mri/transforms/${ses}_to_base_alt2.lta"                          
        elif [ ${pipeline} == "cross" ]; then
            sesdir=${ses}
            lta_ses="--regheader"
        fi
        
        path_T2=$(ls ${bids_dir}/${sub}/${ses}/extra_data/*acq-tse_T2w.nii.gz)
        if [ $(echo ${path_T2} | wc -w) -gt 1 ]; then
            echo "ERROR: FOUND MORE THAN ONE T2-TSE HIPPOSLAB"
            echo "proof: $(echo ${path_T2})"
            exit 1
        fi

        if [ -z ${path_T2} ]; then
            echo "WARNING: CANNOT FIND T2w-TSE hipposlab, RUNNING HCsubfield WITH T2-SPACE"
            path_T2=${fs_dir}/${sub}/${ses}/mri/orig/T2raw.mgz
            if [ ${pipeline} == "long" ]; then
                lta_ses="--lta ${fs_dir}/${sub}/${ses}.long.base/mri/transforms/T2raw.lta"
            elif [ ${pipeline} == "long_alt1" ]; then
                lta_ses="--lta ${fs_dir}/${sub}/${ses}.long.base_alt1/mri/transforms/T2raw.lta"
            elif [ ${pipeline} == "long_alt2" ]; then
                lta_ses="--lta ${fs_dir}/${sub}/${ses}.long.base_alt2/mri/transforms/T2raw.lta"                           
            fi
        fi

        if [ ! -f ${path_T2} ] ; then
            echo "ERROR: could not find hipposlab or T2-space for ${sub} ${ses}"
            echo "proof: "$(ls ${path_T2})""
            continue
        fi

        if [ ! -f "${fs_dir}/${sub}/${sesdir}/mri/rh.hippoSfVolumes-T2only.v22.txt" ]; then

            singularity exec --cleanenv --contain \
            -B ${bids_dir} \
            -B ${tmp_dir}:/tmp \
            -B freesurfer_license.txt:/usr/local/freesurfer/.license \
            -B /cluster/software/EL9/amd/zen/easybuild/software/MCR/R2019b.8/v97:/usr/local/freesurfer/MCRv97 \
            --env SUBJECTS_DIR=${fs_dir}/${sub} \
            --env FSLOUTPUTTYPE=NIFTI_GZ \
                freesurfer_7.3.2.sif \
                    /bin/bash -c \
                    "mri_vol2vol \
                        ${lta_ses} \
                        --mov ${path_T2} \
                        --targ ${path_T2} \
                        --o ${fs_dir}/${sub}/${sesdir}/mri/orig/T2hippo.mgz \
                        --no-resample \
                    && \
                    segmentHA_T2.sh \
                        ${sesdir} \
                        ${fs_dir}/${sub}/${sesdir}/mri/orig/T2hippo.mgz \
                        T2multispectral \
                        1 \
                    && \
                    segmentHA_T2.sh \
                        ${sesdir} \
                        ${fs_dir}/${sub}/${sesdir}/mri/orig/T2hippo.mgz \
                        T2only \
                        0"
        else
            echo "${sub} --> ${ses} --> ${pipeline} segmentHA_T2.sh ALREADY RUN"
            echo -e "\t proof: $(ls ${fs_dir}/${sub}/${sesdir}/mri/rh.hippoSfVolumes-T2only.v22.txt)"
        fi
    
    done    
done

echo ""
echo "##########################################"
echo "                  DONE                    "
echo "##########################################"

