#!/bin/bash

project_dir="/ess/p23/cluster/data/bids_fmri/s2c_gridtask_prisma"
fs_dir="${project_dir}/bids/derivatives/freesurfer_7.3.2"

c=0
dryrun=false # false means failed stuff will be deleted and failed+new jobs will be started
rm logs_samseg/errors_samseg.txt
echo ""
for sub_dir in $(ls -d ${fs_dir}/sub-???????); do
    sub=$(basename ${sub_dir})
    logname="logs_samseg/slurm.samseg.${sub}.log"
    tmp_dir=${project_dir}/scratch/samseg/${sub}

    if [ -f ${logname} ]; then
        echo "${sub}: samseg log exists"
        echo "proof: $(ls logs_samseg/slurm.samseg.${sub}.log)"
        if [ ! -z "$(cat ${logname} | grep "complete: ")" ]; then
            echo "SUCCESS (likely, at least looks like the script finished successfully)"
            echo ""
            continue
        elif [ ! -z "$(cat ${logname} | grep "ERROR")" ]; then
            if [ ${dryrun} = "false" ]; then
                rm $(logname)
            fi
            echo -e "\n\t ERROR detected. Check error before rerunning \n"
            echo ${logname} >> logs_samseg/errors_samseg.txt
            continue       
        else
            echo -e "\n\t STILL RUNNING? \n"
            continue
        fi
    else
        echo "${sub}; no samseg log found"
    fi

    echo -e "\t SUBMITTING\n"

    if [ ${dryrun} = "false" ]; then
        sbatch samseg.sh ${sub} ${fs_dir} ${tmp_dir} ${logname} 
        echo -e "\t SUBMITTED ${sub}, SAMSEG: samseg.sh"
        echo ""
    fi

    c=$((c+1))
#    if [ ${c} -gt 28 ]; then
#        echo ${c}
#        exit 0
#    fi

done

echo ${c}
