#!/bin/bash
#SBATCH --time=04:00:00
#SBATCH --nodes=1
#SBATCH --account=rrg-mchakrav-ab
#SBATCH --ntasks=40

# this script is to find a working seed for AROMA-based ICA in difficult datasts
# by default, it will test 100 different seeds and return in working_seed.txt a list of potential seeds to be used. 

module load gcc/8.3.0; module load openblas/0.3.7; module load fsl/6.0.4

cd /scratch/m/mchakrav/gjoanes/multirat/test/scratch/ica_seed

ls -d sub* | while read line
do
cd $line
touch working_seed.txt

echo $line 
for i in {1..100}
do
seed=$((1 + RANDOM))
nii=$(ls *input.nii.gz)

melodic -i ${nii} -o ${seed} -m mask.nii.gz -d 10 --Ostats --nobet --mmthresh=0.5 --report --seed=${seed} 

grep -q "No convergence" $seed/log.txt || echo $seed >> working_seed.txt
done

cd ..
done


# comands to find cp problematic input images to aroms into a test folder. 
#line='0100209_ses-1'
#cp ../confound/aromas/sub-${line}/confound_correction_wf_main_post_wf/confound_correction_wf/_split_name_sub-${line}_run-1_bold/regress/*aroma_input.nii.gz sub-${line}/
#cp ../confound/aromas/sub-${line}/confound_correction_wf_main_post_wf/confound_correction_wf/_split_name_sub-${line}_run-1_bold/regress/aroma_out/mask.nii.gz sub-${line}/

