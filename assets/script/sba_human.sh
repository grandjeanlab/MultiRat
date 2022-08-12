!#bin/bash

# declare relevant env variables
FSLDIR='/opt/fsl/6.0.0'
template='/opt/fsl/6.0.0/data/standard/MNI152_T1_2mm.nii.gz'
template_mask='/opt/fsl/6.0.0/data/standard/MNI152_T1_2mm_brain_mask.nii.gz'
root_dir='/project/4180000.19/human'
roi_folder=$root_dir'/ROI'

#location of the S1200 HCP dataset
data_dir='/project_freenas/3022017.01/S1200/'


mkdir -p $root_dir'/script'
mkdir -p $root_dir'/output'

mkdir -p $roi_folder

#make the ROIs based on the 2mm MNI template
cd $roi_folder

sphere=4

fslmaths ${template} -mul 0 -add 1 -roi 43 1 87 1 37 1 0 1 tmp
fslmaths tmp -kernel sphere ${sphere} -fmean -bin ACA
rm tmp.nii.gz

fslmaths ${template} -mul 0 -add 1 -roi 21 1 58 1 61 1 0 1 tmp
fslmaths tmp -kernel sphere ${sphere} -fmean -bin S1_r
rm tmp.nii.gz

fslmaths ${template} -mul 0 -add 1 -roi 71 1 58 1 61 1 0 1 tmp
fslmaths tmp -kernel sphere ${sphere} -fmean -bin S1_l
rm tmp.nii.gz


cd $root_dir


counter=0


#loop through the processed scans, for each scans, loop through the region-of-interests and estimate seed-based analysis
ls $data_dir | while read subject
do
echo "now doing scan "$subject
echo "counters "$counter

for scan_name in rfMRI_REST1_7T_PA #rfMRI_REST2_7T_AP rfMRI_REST3_7T_PA rfMRI_REST4_7T_AP #rfMRI_REST2_LR rfMRI_REST1_RL rfMRI_REST2_RL 
do
file_path=$data_dir'/'$subject'/MNINonLinear/Results/'$scan_name'/'$scan_name'_hp2000_clean.nii.gz'


if [ -f "$file_path" ]; then
counter=$((counter+1))

qscript='SBA_job_'$counter'.sh'


# create a directory on scratch
echo "mkdir -p /scratch/data/joagra/"$subject > $root_dir/script/$qscript
echo "cd /scratch/data/joagra/"$subject >> $root_dir/script/$qscript
# load afni
echo "module load afni" >> $root_dir/script/$qscript
#apply blur and bandpass. Afni is the easiest. 
echo "3dTproject -automask  -blur 2.5 -passband 0.01 0.1 -input "$file_path" -prefix tproject.nii.gz" >> $root_dir/script/$qscript
# resample to 2mm because that is the space in which the seeds are
echo "3dresample -input tproject.nii.gz -prefix resample.nii.gz -master /project/4180000.19/human/ROI/S1_r.nii.gz" >> $root_dir/script/$qscript
#extract the timeseries for each ROI
echo "fslmeants -i resample.nii.gz -m /project/4180000.19/human/ROI/S1_r.nii.gz -o S1_r" >> $root_dir/script/$qscript
echo "fslmeants -i resample.nii.gz -m /project/4180000.19/human/ROI/S1_l.nii.gz -o S1_l" >> $root_dir/script/$qscript
echo "fslmeants -i resample.nii.gz -m /project/4180000.19/human/ROI/ACA.nii.gz -o ACA" >> $root_dir/script/$qscript
#perform correlation analysis and output to relevant directory
echo "1ddot -terse  S1_r  S1_l ACA > "$root_dir"/output/"$subject >> $root_dir/script/$qscript
#clean scratch. 
echo "cd .." >> $root_dir/script/$qscript
echo "rm -rf "$subject >> $root_dir/script/$qscript

#make script executable
chmod +x $root_dir/script/$qscript
#because we have access to a cluster system, we submit jobs. 
qsub -l 'procs=1,mem=16gb,walltime=01:00:00' $root_dir/script/$qscript
#if a cluster is not available, comment the line above and use the following instead. Might be a tad slow to run. 
#bash ${PWD}/script/$qscript
  
 
fi
done
done




##the estimation of specificity is done in R, because it was faster to do it this way
#setwd('/project/4180000.19/human/output/')

#file_list<-dir('.')

#sp<-c()

#for(i in file_list){
#  df<-read.table(i)
#  if(df[1,2]>=0.1 & df[1,3]<0.1){sp<-c(sp,1)}else{sp<-c(sp,0)}
#}

