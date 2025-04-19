#!/bin/bash

project_dir="/ess/p23/cluster/data/bids_fmri/s2c_gridtask_prisma"
fs_dir="${project_dir}/bids/derivatives/freesurfer_7.3.2"

rm logs_freesurfer/errors_base.txt
c=0
dryrun=false

for sub_dir in $(ls -d ${fs_dir}/sub-???????); do
    sub=$(basename ${sub_dir})
    logname="logs_freesurfer/slurm.freesurfer-base.${sub}.log"
    tmp_dir=${project_dir}/scratch/freesurfer_7.3.2/${sub}_base

    #First check that all cross-sectional steps ran successfully
    status_cross=good
    lognames_cross=$(ls logs_freesurfer/slurm.freesurfer-default.${sub}_ses-??.log)
    for logname_cross in ${lognames_cross}; do
        if [ -z "$(cat ${logname_cross} | grep "finished without error")" ]; then
            echo "${sub}"
            echo -e "\t FS CROSS (ONE OR SEVERAL SESSIONS) FINISHED WITH ERROR OR IS STILL RUNNING. FIX THIS FIRST \n"
            echo "FSCROSS ERROR: ${logname_cross}" >> logs_freesurfer/errors_base.txt
            status_cross=bad
        fi
    done
    if [ ${status_cross} = "bad" ]; then continue; fi

    #If base log already exist, check if it ran with errors. If success, check if new timepoints have arrived and need to be incorporated
    if [ -f ${logname} ]; then
        echo "${sub} --> base log exists"
        echo "proof: $(ls logs_freesurfer/slurm.freesurfer-base.${sub}.log)"

        if [ ! -z "$(cat ${logname} | grep "exited with ERRORS")" ]; then
            echo "Log indicates that Freesurfer exited with ERRORS. Deleting and rerunning"
            echo ${logname} >> logs_freesurfer/errors_base.txt

            if [ ${dryrun} = "false" ]; then
                rm ${logname}
                rm -r ${fs_dir}/${sub}/base
            fi

        elif [ ! -z "$(cat ${logname} | grep "ERROR: You are trying to re-run an existing subject")" ]; then
            echo "Log indicates that Freesurfer didn't start because subject already existed. Deleting and rerunning"
            echo "Re-run? ${logname}" >> logs_freesurfer/errors_base.txt

            if [ ${dryrun} = "false" ]; then
                rm ${logname}
                rm -r ${fs_dir}/${sub}/base
            fi

        elif [ ! -z "$(cat ${logname} | grep "finished without error")" ]; then

            if [ $(ls -d $fs_dir/$sub/ses-??  | wc -w) -gt $(cat $fs_dir/$sub/base/base-tps  | wc -w) ]; then
                echo "Higher number of cross-sectional sessions than base-tps"
                echo "New timepoint added for this participant"
                echo "Deleting base/long-folders and rerunning with new cross-sectional time point"
                echo "Also deleting samseg-folder+logs and  Deleting base/long-folders and rerunning with new cross-sectional time point"

                if [ ${dryrun} = "false" ]; then
                    rm ${logname}
                    rm logs_freesurfer/slurm.freesurfer-long.${sub}*.log
                    rm logs_samseg/*${sub}*
                    rm logs_fssubregions/*${sub}*
                    rm -r ${fs_dir}/${sub}/base
                    rm -r ${fs_dir}/${sub}/ses-??.long.base
                    rm -r ${fs_dir}/${sub}/samseg
                fi

            else
                echo -e "\t FS BASE SUCCESSFUL (and no new time points discovered)"
                continue
            fi

        else
            echo -e "\t FS BASE STILL RUNNING?"
            continue
        fi
    else
        echo -e "${sub} --> base NO LOG FOUND"
        if [ ${dryrun} = "false" ]; then
            rm -r ${fs_dir}/${sub}/base
        fi
    fi

    if [ ${dryrun} = "false" ]; then
        sbatch freesurfer_base.sh ${sub} ${fs_dir} ${tmp_dir} ${logname}
        echo "SUBMITTED ${sub} --> base, freesurfer_base.sh"
        echo ""
    fi

    c=$((c+1))

    # if [ ${c} -gt 6 ]; then
    #     echo ${c}
    #     exit 0
    # fi

done

echo ${c}