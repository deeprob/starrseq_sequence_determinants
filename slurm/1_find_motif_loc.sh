#!/bin/bash
#SBATCH --account=girirajan
#SBATCH --partition=girirajan
#SBATCH --job-name=starr_motifloc
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --time=400:0:0
#SBATCH --mem-per-cpu=10G
#SBATCH --chdir /data5/deepro/starrseq/main_library/6_sequence_determinants/src
#SBATCH -o /data5/deepro/starrseq/main_library/6_sequence_determinants/slurm/logs/out_find_motif.log
#SBATCH -e /data5/deepro/starrseq/main_library/6_sequence_determinants/slurm/logs/err_find_motif.log
#SBATCH --nodelist laila

# >>> conda initialize >>>
# !! Contents within this block are managed by 'conda init' !!
__conda_setup="$('/data5/deepro/miniconda3/bin/conda' 'shell.bash' 'hook' 2> /dev/null)"
if [ $? -eq 0 ]; then
    eval "$__conda_setup"
else
    if [ -f "/data5/deepro/miniconda3/etc/profile.d/conda.sh" ]; then
        . "/data5/deepro/miniconda3/etc/profile.d/conda.sh"
    else
        export PATH="/data5/deepro/miniconda3/bin:$PATH"
    fi
fi
unset __conda_setup
# <<< conda initialize <<<

conda activate starrseq_mea

echo `date` starting job on $HOSTNAME

python /data5/deepro/starrseq/main_library/6_sequence_determinants/src/1_find_motif_loc.py

echo `date` ending job
