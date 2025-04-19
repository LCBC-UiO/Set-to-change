#!/bin/bash

project_dir="/ess/p23/cluster/data/bids_fmri/s2c_gridtask_prisma"
bids_dir="${project_dir}/bids"
fs_dir="${bids_dir}/derivatives/freesurfer_7.3.2"

c=0
for bidssub_dir in $(ls -d ${bids_dir}/sub-???????); do
    sub=$(basename ${bidssub_dir})
    sessions=$(ls -d ${bids_dir}/${sub}/ses-??)

    for path_ses in ${sessions}; do
        ses=$(basename ${path_ses})
        tmp_dir=${project_dir}/scratch/freesurfer_7.3.2/${sub}_${ses}

        logname="logs_freesurfer/slurm.freesurfer-default.${sub}_${ses}.log"
        if [ -f ${logname} ]; then
            echo "${sub} --> ${ses}, default log exists"
            echo "proof: $(ls logs_freesurfer/slurm.freesurfer-default.${sub}_${ses}.log)"

            if [ ! -z "$(cat ${logname} | grep "exited with ERRORS")" ]; then
                echo "Log indicates that Freesurfer exited with ERRORS. Deleting and rerunning"
                rm ${logname}
                rm -r ${fs_dir}/${sub}/${ses}

            elif [ ! -z "$(cat ${logname} | grep "ERROR: You are trying to re-run an existing subject")" ]; then
                echo "Log indicates that Freesurfer didn't start because subject already existed. Deleting and rerunning"
                rm ${logname}
                rm -r ${fs_dir}/${sub}/${ses}

            elif [ ! -z "$(cat ${logname} | grep "finished without error")" ]; then
                echo -e "\t FS SUCCESSFUL"
                continue

            else
                echo -e "\t FS STILL RUNNING?"
                continue
            fi
        else
            rm -r ${fs_dir}/${sub}/${ses}
        fi

        echo "SUBMITTING ${sub} --> ${ses}, freesurfer_default.sh"
        sbatch freesurfer_default.sh ${sub} ${ses} ${bids_dir} ${fs_dir} ${tmp_dir} ${logname}
        echo ""
        c=$((c+1))

        # if [ ${c} -gt 63 ]; then
        #     echo ${c}
        #     exit 0
        # fi

    done

done

echo ${c}