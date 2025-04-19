#!/bin/bash

#SBATCH -J samseg_S2C
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=7G
#SBATCH --time=24:00:00
#SBATCH --account=p23
##SBATCH --partition=accel
##SBATCH --gres=gpu:0
#SBATCH --output logs_samseg/slurm-%j.txt

module purge
#echo "LOADING SINGULARITY MODULE"
#module load singularity/3.7.3
#echo `which singularity`
#set -o errexit
unset PYTHONPATH

sub=${1}
fs_dir=${2}
tmp_dir=${3}
logname=${4}

mv logs_samseg/slurm-${SLURM_JOBID}.txt ${logname}

rm -r ${tmp_dir}
mkdir -p ${tmp_dir}
export TMPDIR=${tmp_dir}
env | grep TMPDIR

export FREESURFER_HOME=/ess/p23/cluster/tools/mri/freesurfer/freesurfer.7.3.2
source $FREESURFER_HOME/SetUpFreeSurfer.sh
export SUBJECTS_DIR=${fs_dir}/${sub}

echo ""
echo "####################################################"
echo "                  SAMSEG - LONG                    "
echo "####################################################"

sessions=$(ls -d ${fs_dir}/${sub}/ses-??.long.base)

if [ $(echo $sessions | wc -w) -gt 1 ]; then
    echo "MULTIPLE VISITS DETECTED"            
    echo "CHECKING IF SAMSEG-LONG HAS ALREADY BEEN RUN FOR ALL VISITS"
    counter=1
    all_tps_done=true
    for session in ${sessions}; do
        tp_var=$(printf "tp%03d" "${counter}")
        if [ ! -f "${fs_dir}/${sub}/samseg-long/${tp_var}/sbtiv.stats" ]; then
            all_tps_done=false
            break
        else
            echo $(ls -l "${fs_dir}/${sub}/samseg-long/${tp_var}/sbtiv.stats")
        fi
        ((counter++))
    done

    if [ "${all_tps_done}" = false ]; then

        echo "SAMSEG-LONG NOT RUN FOR ALL VISITS. RUNNING NOW"
                        
        useT2=true
        for path_ses in ${sessions}; do
            if [ ! -f ${path_ses}/mri/T2.prenorm.mgz ]; then
                echo "WARNING: at least one session missing T2.prenorm.mgz. Running SAMSEG_LONG with mprages only"
                useT2=false
                break
            fi
        done

        time_point_flag=""
        for path_ses in ${sessions}; do
            if [ ${useT2} = "true" ]; then
                time_point_flag+=" --timepoint ${path_ses}/mri/orig.mgz ${path_ses}/mri/T2.prenorm.mgz"
            else
                time_point_flag+=" --timepoint ${path_ses}/mri/orig.mgz"
            fi
        done

        echo ${time_point_flag}

        if [ ${useT2} = "true" ]; then
            run_samseg_long ${time_point_flag} \
                --pallidum-separate \
                --output ${fs_dir}/${sub}/samseg-long \
                --threads 4
        else
            run_samseg_long ${time_point_flag} \
                --output ${fs_dir}/${sub}/samseg-long \
                --threads 4
        fi
        
    else

        echo "SAMSEG-LONG ALREADY RUN FOR ALL VISITS. MOVING ON"    
    fi    
fi


echo ""
echo "####################################################"
echo "                  SAMSEG - LONG-ALT2                "
echo "####################################################"

sessions=$(ls -d ${fs_dir}/${sub}/ses-??.long.base_alt2)

if [ $(echo $sessions | wc -w) -gt 1 ]; then
    echo "MULTIPLE VISITS DETECTED"        
    echo "CHECKING IF SAMSEG-LONG_ALT2 HAS ALREADY BEEN RUN FOR ALL VISITS"
    counter=1
    all_tps_done=true
    for session in ${sessions}; do
        tp_var=$(printf "tp%03d" "${counter}")
        if [ ! -f "${fs_dir}/${sub}/samseg-long_alt2/${tp_var}/sbtiv.stats" ]; then
            all_tps_done=false
            break
        else
            echo $(ls -l "${fs_dir}/${sub}/samseg-long_alt2/${tp_var}/sbtiv.stats")
        fi
        ((counter++))
    done

    if [ "${all_tps_done}" = false ]; then

        echo "SAMSEG-LONG_ALT2 NOT RUN FOR ALL VISITS. RUNNING NOW"
                        
        useT2=true
        for path_ses in ${sessions}; do
            if [ ! -f ${path_ses}/mri/T2.prenorm.mgz ]; then
                echo "WARNING: at least one session missing T2.prenorm.mgz. Running SAMSEG_LONG_ALT2 with mprages only"
                useT2=false
                break
            fi
        done

        time_point_flag=""
        for path_ses in ${sessions}; do
            if [ ${useT2} = "true" ]; then
                time_point_flag+=" --timepoint ${path_ses}/mri/orig.mgz ${path_ses}/mri/T2.prenorm.mgz"
            else
                time_point_flag+=" --timepoint ${path_ses}/mri/orig.mgz"
            fi
        done

        echo ${time_point_flag}

        if [ ${useT2} = "true" ]; then
            run_samseg_long ${time_point_flag} \
                --pallidum-separate \
                --output ${fs_dir}/${sub}/samseg-long_alt2 \
                --threads 4
        else
            run_samseg_long ${time_point_flag} \
                --output ${fs_dir}/${sub}/samseg-long_alt2 \
                --threads 4
        fi
        
    else

        echo "SAMSEG-LONG_ALT2 ALREADY RUN FOR ALL VISITS. MOVING ON"    
    fi    
else
    echo "ONLY ONE VISIT DETECTED - FAKING LONG_ALT2 PIPELINE FOR THIS ONE"        
    useT2=true
    path_ses=${sessions}
    if [ ! -f ${path_ses}/mri/T2.prenorm.mgz ]; then
        echo "WARNING: missing T2.prenorm.mgz. Running SAMSEG_LONG_ALT2 with mprage only"
        useT2=false
    fi
    
    if [ ! -f "${fs_dir}/${sub}/samseg-long_alt2/tp001/sbtiv.stats" ]; then    
        echo "SAMSEG-LONG_ALT2 NOT RUN. RUNNING NOW"    
        if [ ${useT2} = "true" ]; then
            run_samseg \
            --input ${path_ses}/mri/orig.mgz ${path_ses}/mri/T2.prenorm.mgz \
            --pallidum-separate \
            --output ${fs_dir}/${sub}/samseg-long_alt2/tp001 \
            --threads 4
        else
            run_samseg \
            --input ${path_ses}/mri/orig.mgz \
            --output ${fs_dir}/${sub}/samseg-long_alt2/tp001 \
            --threads 4
        fi
    else
        echo "SAMSEG-LONG_ALT2 ALREADY RUN. MOVING ON"        
    fi
fi



echo ""
echo "####################################################"
echo "                  SAMSEG - CROSS                    "
echo "####################################################"

sessions=$(ls -d ${fs_dir}/${sub}/ses-??)
echo "CHECKING IF SAMSEG-CROSS HAS ALREADY BEEN RUN FOR ALL VISITS"
all_tps_done=true

for session in ${sessions}; do
    tp_var=$(basename ${session})
    if [ ! -f "${fs_dir}/${sub}/samseg-cross/${tp_var}/sbtiv.stats" ]; then
        all_tps_done=false
        break
    else
        echo $(ls -l "${fs_dir}/${sub}/samseg-cross/${tp_var}/sbtiv.stats")
    fi
done

if [ "${all_tps_done}" = false ]; then

    echo "SAMSEG-CROSS NOT RUN FOR ALL VISITS. RUNNING NOW"

    useT2=true
    for path_ses in ${sessions}; do
        if [ ! -f ${path_ses}/mri/T2.prenorm.mgz ]; then
            echo "WARNING: missing T2.prenorm.mgz in at least one session. Running SAMSEG with mprage only"
            useT2=false
            break
        fi
    done

    for path_ses in ${sessions}; do
        ses=$(basename ${path_ses})
        if [ ${useT2} = "true" ]; then
            run_samseg \
            --input ${path_ses}/mri/orig.mgz ${path_ses}/mri/T2.prenorm.mgz \
            --pallidum-separate \
            --output ${fs_dir}/${sub}/samseg-cross/${ses} \
            --threads 4
        else
            run_samseg \
            --input ${path_ses}/mri/orig.mgz \
            --output ${fs_dir}/${sub}/samseg-cross/${ses} \
            --threads 4
        fi
    done
        
else
    echo "SAMSEG-CROSS ALREADY RUN FOR ALL VISITS. MOVING ON"    
fi    
