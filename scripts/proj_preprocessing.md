MultiRAT analysis code
================
Joanes Grandjean

![rat art](../assets/img/rat_art.png) \# Foreword The preprocessing is
performed in four steps.

1.  Functional preprocessing  
2.  [Qality control](scripts/proj_qa.md)  
3.  Confound removal
4.  Analysis

Currently this is achieved all through a single script submitted as a
sbatch job on a SLURM-running super computer (niagara.computecanada.ca),
and submitted with the following command `sbatch
--export=bidsID='02004',TR='1' script/job_mr_stim.sh`

Importantly, in this script, RABIES is run in a singularity. This was
built using `singularity build rabies.0.2.0-dev.simg
docker://gabdesgreg/rabies:0.2.0-dev`.

See the [qality control](scripts/proj_qa.md). See the analysis for
stimulus evoked. (pending) See the analysis for resting-state. (pending)

### Here is the script for stimulus-evoked preprocessing and confound regression

``` bash
#!/bin/bash
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --account=rrg-mchakrav-ab
#SBATCH --ntasks=80


# load modules necessary to run FSL on niagara
module load gcc/8.3.0
module load openblas/0.3.7
module load fsl/6.0.4

# provide path to templates and masks. 
template=/template/SIGMA_Wistar_Rat_Brain_TemplatesAndAtlases_Version1.1/SIGMA_Rat_Anatomical_Imaging/SIGMA_Rat_Anatomical_InVivo_Template/SIGMA_InVivo_Brain_Template.nii
mask=/template/SIGMA_Wistar_Rat_Brain_TemplatesAndAtlases_Version1.1/SIGMA_Rat_Anatomical_Imaging/SIGMA_Rat_Anatomical_InVivo_Template/SIGMA_InVivo_Brain_Mask.nii
wm=/template/SIGMA_Wistar_Rat_Brain_TemplatesAndAtlases_Version1.1/SIGMA_Rat_Anatomical_Imaging/SIGMA_Rat_Anatomical_InVivo_Template/SIGMA_InVivo_WM_bin.nii.gz
csf=/template/SIGMA_Wistar_Rat_Brain_TemplatesAndAtlases_Version1.1/SIGMA_Rat_Anatomical_Imaging/SIGMA_Rat_Anatomical_InVivo_Template/SIGMA_InVivo_CSF_bin.nii.gz
atlas=/template/SIGMA_Wistar_Rat_Brain_TemplatesAndAtlases_Version1.1/SIGMA_Rat_Brain_Atlases/SIGMA_Anatomical_Atlas/SIGMA_Anatomical_Brain_Atlas.nii
roi=/template/roi/


# define main folder inputs/outputs
tar_file=$PROJECT/multirat/ds/${bidsID}.tar.gz
bids_folder=$SCRATCH/multirat/ds/${bidsID}
preprocess_folder=$SCRATCH/multirat/preprocess/${bidsID}


# clean and make directories
rm -rf ${bids_folder}
rm -rf ${preprocess_folder}
mkdir -p ${bids_folder}
mkdir -p ${preprocess_folder}

# make output directories
mkdir -p $PROJECT/multirat/qa/epi2anat
mkdir -p $PROJECT/multirat/qa/anat2template
mkdir -p $PROJECT/multirat/qa/template2std
mkdir -p $PROJECT/multirat/qa/temporal
mkdir -p $PROJECT/multirat/snr/atlas
mkdir -p $PROJECT/multirat/snr/tsnr/pre
mkdir -p $PROJECT/multirat/motion/confound
mkdir -p $PROJECT/multirat/motion/FD
mkdir -p $PROJECT/multirat/log
mkdir -p $PROJECT/multirat/tar/confound
mkdir -p $PROJECT/multirat/tar/analysis
mkdir -p $PROJECT/multirat/stim

#untar the tarball containing the specified dataset
tar -xvzf ${tar_file} -C ${bids_folder} 

singularity run -B ${PROJECT}/template:/template -B ${bids_folder}:/rabies_in:ro \
-B ${preprocess_folder}:/rabies_out \
/project/m/mchakrav/gjoanes/rabies.0.2.0-dev.simg -p MultiProc preprocess /rabies_in /rabies_out \
--template_reg_script multiRAT \
--coreg_script light_SyN \
--no_STC \
--anat_template ${template} \
--brain_mask ${mask} \
--WM_mask ${wm} \
--CSF_mask ${csf} \
--vascular_mask ${csf} \
--labels ${atlas} \
--TR ${TR}s \
--commonspace_resampling 0.2x0.2x0.2


# copy preprocessing outputs of interest to designated folders
cp ${preprocess_folder}/rabies_preprocess.log $PROJECT/multirat/log/${bidsID}.log
cp ${preprocess_folder}/QC_report/EPI2Anat/* $PROJECT/multirat/qa/epi2anat/
cp ${preprocess_folder}/QC_report/Anat2Template/* $PROJECT/multirat/qa/anat2template/
cp ${preprocess_folder}/QC_report/temporal_diagnosis/* $PROJECT/multirat/qa/temporal/
cp ${preprocess_folder}/QC_report/Template2Commonspace/_registration.png $PROJECT/multirat/qa/template2std/${bidsID}.png

cp ${preprocess_folder}/bold_datasink/bold_labels/*/*/*.nii.gz $PROJECT/multirat/snr/atlas/
cp ${preprocess_folder}/confounds_datasink/confounds_csv/*/*/*.csv $PROJECT/multirat/motion/confound/
cp ${preprocess_folder}/confounds_datasink/FD_csv/*/*/*.csv $PROJECT/multirat/motion/FD/


# extract tSNR from atlas. 
ls ${preprocess_folder}/bold_datasink/commonspace_bold_labels | while read label
do

label_path=$(find ${preprocess_folder}/bold_datasink/commonspace_bold_labels/${label}/*/*.nii.gz)
tSNR_path=$(find ${preprocess_folder}/bold_datasink/tSNR_filename/${label}/*/*.nii.gz)
mask_path=$(find ${preprocess_folder}/bold_datasink/commonspace_bold_mask/${label}/*/*.nii.gz)

fslmeants -i ${tSNR_path} \
-m ${mask_path} \
--label=${label_path} \
-o $PROJECT/multirat/snr/tsnr/pre/${label}.txt

done


denoise=smooth

current_confound=$SCRATCH/multirat/confound/${denoise}/${bidsID}


rm -rf $current_confound
mkdir -p $current_confound

singularity run -B ${PROJECT}/template:/template -B ${bids_folder}:/rabies_in \
-B ${preprocess_folder}:/rabies_out -B ${current_confound}:/confound_out \
/project/m/mchakrav/gjoanes/rabies.0.2.0-dev.simg -p MultiProc confound_regression /rabies_out /confound_out \
--TR ${TR}s --commonspace_bold \
--smoothing_filter 0.5 



tar -czvf $PROJECT/multirat/tar/confound/${bidsID}_${denoise}.tar.gz ${current_confound}

cp ${current_confound}/confound_regression_datasink/cleaned_timeseries/*/*/*.nii.gz $PROJECT/multirat/stim/
```