#!/bin/bash

project_dir="/ess/p23/cluster/data/bids_fmri/s2c_gridtask_prisma"
fs_dir="${project_dir}/bids/derivatives/freesurfer_7.3.2"

c=0
dryrun=false # false means failed stuff will be deleted and failed+new jobs will be started
rm logs_freesurfer/errors_long.txt

for sub_dir in $(ls -d ${fs_dir}/sub-???????); do
    sub=$(basename ${sub_dir})
    if [ -d ${sub_dir}/base ]; then
        logname_base="logs_freesurfer/slurm.freesurfer-base.${sub}.log"

        sessions=$(ls -d ${fs_dir}/${sub}/ses-??)
        for path_ses in ${sessions}; do
            ses=$(basename ${path_ses})

            tmp_dir=${project_dir}/scratch/freesurfer_7.3.2/${sub}_${ses}_long
            logname="logs_freesurfer/slurm.freesurfer-long.${sub}_${ses}.log"

            if [ -f ${logname} ]; then
                echo "${sub} --> ${ses}: long log exists"
                #echo "proof: $(ls logs_freesurfer/slurm.freesurfer-long.${sub}_${ses}.log)"

                if [ ! -z "$(cat ${logname} | grep "exited with ERRORS")" ]; then
                    echo "Log indicates that Freesurfer exited with ERRORS. Deleting and rerunning"
                    echo ${logname} >> logs_freesurfer/errors_long.txt

                    if [ ${dryrun} = "false" ]; then
                        rm ${logname}
                        rm -r ${fs_dir}/${sub}/${ses}.long.base
                    fi

                elif [ ! -z "$(cat ${logname} | grep "ERROR: You are trying to re-run an existing subject")" ]; then
                    echo "Log indicates that Freesurfer didn't start because subject already existed. Deleting and rerunning"
                    echo ${logname} >> logs_freesurfer/errors_long.txt

                    if [ ${dryrun} = "false" ]; then
                        rm ${logname}
                        rm -r ${fs_dir}/${sub}/${ses}.long.base
                    fi

                elif [ ! -z "$(cat ${logname} | grep "finished without error")" ]; then
                    echo -e "FS LONG SUCCESSFUL \n"
                    continue

                else
                    echo -e "\t FS LONG STILL RUNNING? \n"
                    continue
                fi
            else
                echo "${sub} --> ${ses}: NO LONG LOG FOUND"
                if [ -z "$(cat ${logname_base} | grep "finished without error")" ]; then
                    echo -e "\t FS BASE FINISHED WITH ERROR OR IS STILL RUNNING. FIX THIS FIRST \n"
                    echo "FSBASE ERROR: ${logname}" >> logs_freesurfer/errors_long.txt

                    continue
                else
                    echo "Removing ${sub}/${ses}.long.base if it exists"
                    if [ ${dryrun} = "false" ]; then
                        rm -r ${fs_dir}/${sub}/${ses}.long.base
                    fi
                fi
            fi

            if [ ${dryrun} = "false" ]; then
                sbatch freesurfer_long.sh ${sub} ${ses} ${fs_dir} ${tmp_dir} ${logname}
                echo "SUBMITTED ${sub} --> ${ses}, LONG: freesurfer_long.sh"
                echo ""
            fi

            c=$((c+1))

        done
    fi
    # if [ ${c} -gt 64 ]; then
    #     echo ${c}
    #     exit 0
    # fi
done
echo ${c}