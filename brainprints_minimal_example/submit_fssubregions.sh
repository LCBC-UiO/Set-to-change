#!/bin/bash

project_dir="/ess/p23/cluster/data/bids_fmri/s2c_gridtask_prisma"
fs_dir="${project_dir}/bids/derivatives/freesurfer_7.3.2"

c=0
dryrun=false # false means failed stuff will be deleted and failed+new jobs will be started
rm logs_fssubregions/errors_fssubregions.txt

for sub_dir in $(ls -d ${fs_dir}/sub-???????); do
    sub=$(basename ${sub_dir})
    logname="logs_fssubregions/slurm.fssubregions.${sub}.log"
    tmp_dir=${project_dir}/scratch/fs_subregions/${sub}

    if [ -f ${logname} ]; then
        echo -e " \n ${sub}: fssubregions log exists"
        echo "proof: $(ls logs_fssubregions/slurm.fssubregions.${sub}.log)"
        if [ ! -z "$(cat ${logname} | tail | grep "DONE")" ]; then
            echo -e "SUCCESS (likely, at least looks like the script finished successfully) \n"
            continue
        elif [ ! -z "$(cat ${logname} | grep "ERROR")" ]; then
            if [ ${dryrun} = "false" ]; then
                #rm $(logname)            
                echo -e "\t ERROR detected. Skipping \n"
                echo "${sub}" >> logs_fssubregions/errors_fssubregions.txt
                continue
            fi
        else
            echo -e "\t STILL RUNNING? \n"      
            continue
        fi
    else
        echo "${sub} --> ${ses}: no fssubregions log found"
        echo "Removing some intermediate files that might have been produced if the job has failed"
        ls ${fs_dir}/${sub}/ses-??.long.base/scripts/IsRunningHP*.lh+rh
        if [ ${dryrun} = "false" ]; then
            rm ${fs_dir}/${sub}/ses-??.long.base/scripts/IsRunningHP*.lh+rh
            rm ${fs_dir}/${sub}/ses-??.long.base/scripts/*hippocampal-subfields*
            rm ${fs_dir}/${sub}/ses-??.long.base/mri/transforms/T1_to_T2*
        fi        
    fi    

    echo -e "\t SUBMITTING\n"

    if [ ${dryrun} = "false" ]; then
        sbatch fs_subregions.sh ${sub} ${logname} ${tmp_dir} 
        echo "SUBMITTED ${sub}, FS_SUBREGIONS: fs_subregion.sh"
        echo ""
    fi

    c=$((c+1))

    # if [ ${c} -gt 64 ]; then
    #     echo ${c}
    #     exit 0
    # fi

done
cat logs_fssubregions/errors_fssubregions.txt
echo ${c}
