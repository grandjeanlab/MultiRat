MultiRAT asset preparation
================
Joanes Grandjean

![rat art](../assets/img/rat_art.png)

1.  [Functional preprocessing](#functional-preprocessing)  
2.  [Confound removal and Analysis](#confound-removal-and-analysis)   
3.  [Confound remove (minimal)](#confound-remove-minimal)

Currently this is achieved all through a single script submitted as a
sbatch job on a SLURM-running super computer (niagara.computecanada.ca). The exact call my differ on your environement. 

Importantly, in this script, RABIES is run in a singularity. This was
built using `singularity build rabies-0.3.5.simg docker://gabdesgreg/rabies:0.3.5`.

Each session was processed separately to increase the odds of success. Individualized scripts were generated using [this script](proj_submitjobs.ipynb). Failed confound_correction runs were detected using [this script](proj_resubmitjobs.ipnb)

For latest updates on RABIES, consult the [Github](https://github.com/CoBrALab/RABIES)

In some instances, running RABIES in parallel using the `-p MultiProc` is known to cause some crashes. If so, try running in linear (default) mode. 

Please report issues [here](https://github.com/CoBrALab/RABIES/issues)

## Preregistration
From the preregistration: 
"All scans will be preprocessed using RABIES (https://github.com/CoBrALab/RABIES), a BIDS-based software based on the fMRIprep pipeline. Data will be co-registered into Waxholm space (https://scalablebrainatlas.incf.org/rat/PLCJB14). Denoising will be performed using motion regression concurently with either: global signal regression, white matter+CSF signal regression, or  ICA-AROMA automatic method (https://github.com/maartenmennes/ICA-AROMA) adapted for the rodent and RABIES (https://github.com/Gab-D-G/conf_reg_pkg). Temporal filtering will be applied at 0.01-0.1 Hz for all scans (3dbandpass). To account for anesthesia requiring larged bandpass filter, we will also try a 0.01-0.25 Hz filter. Smoothing will be applied at 0.5mm for all scans (3dBlurInMask)."   


## Functional preprocessing  
The script belows initiate directories and run the preprocessing. This is common to both resting-state and stimulus-evoked datasets. 

 --commonspace_resampling 0.3x0.3x0.3 --anat_inho_cor_method N4_reg --bold_inho_cor_method N4_reg --anat_autobox --coreg_script Rigid 

The following options are justified:   
`--fast_commonspace`, each session is treated individually, there is no need for the generation of a common space.
`--anat_inho_cor_method N4_reg`, In my trials, I found N4_reg generally more reliable. Notably, for failed runs, I replaced N4_reg with no_reg which allowed completion of additional scans. 
`--anat_autobox`, I find autobox helps greatly. In some instances, either this wasn't enough or autobox cropped too much. In which case, I cropped manually using [this script](script/resize_anat_func.sh). 
`--coreg_script Rigid`, I find that rigid EPI to anatomical registration is more reliable using rigid registration. The downside is that this registration does not correct for (minor) deformations. Higher quality datasets may want to use Syn registrations instead. 
`--commonspace_resampling 0.3x0.3x0.3`, This is exclusively to try to save disk space. This is selected on the basis of the resolutions found among most datasets, and my need to limit disk usage.  

Quality metrics were carefully examined for each dataset before proceeding to the confound regression and analysis steps. See [QA steps](proj_qa.ipynb). 

``` bash
#!/bin/bash
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --account=rrg-mchakrav-ab
#SBATCH --ntasks=40


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


TR=1.0

# define input/output folders for subject sub-0101000 ses-1
orig_folder=/scratch/m/mchakrav/gjoanes/multirat/test/bids/sub-0101000/ses-1
template_folder=/scratch/m/mchakrav/gjoanes/multirat/test/template
bids_folder=/scratch/m/mchakrav/gjoanes/multirat/test/scratch/bids/sub-0101000_ses-1
preprocess_folder=/scratch/m/mchakrav/gjoanes/multirat/test/scratch/preprocess/sub-0101000_ses-1

# clean and make directories
rm -rf $bids_folder; rm -rf $preprocess_folder; mkdir -p $bids_folder/sub-0101000; mkdir -p $preprocess_folder

# move scans into individualized BIDS directories 
cp -r /scratch/m/mchakrav/gjoanes/multirat/test/bids/sub-0101000/ses-1 $bids_folder/sub-0101000/

# the RABIES call with optional arguments
singularity run -B ${template_folder}:/template -B ${bids_folder}:/rabies_in:ro -B ${preprocess_folder}:/rabies_out \
/project/m/mchakrav/gjoanes/rabies-0.3.5.simg -p MultiProc \
preprocess /rabies_in /rabies_out  \
--anat_template ${template} \
--brain_mask ${mask} \
--WM_mask ${wm} \
--CSF_mask ${csf} \
--vascular_mask ${csf} \
--labels ${atlas} \
--TR ${TR}s \
--fast_commonspace \
--commonspace_resampling 0.3x0.3x0.3 \
--anat_inho_cor_method N4_reg \
--bold_inho_cor_method N4_reg \
--anat_autobox \
--coreg_script Rigid 


# copy preprocessing outputs of interest to designated folders
cp ${preprocess_folder}/rabies_preprocess.log /scratch/m/mchakrav/gjoanes/multirat/test/export/log/preprocess/sub-0101000_ses-1.log
cp ${preprocess_folder}/preprocess_QC_report/EPI2Anat/* /scratch/m/mchakrav/gjoanes/multirat/test/export/qa/epi2anat/
cp ${preprocess_folder}/preprocess_QC_report/temporal_features/* /scratch/m/mchakrav/gjoanes/multirat/test/export/qa/temporal/
cp ${preprocess_folder}/preprocess_QC_report/commonspace_reg_wf.Native2Atlas/*png /scratch/m/mchakrav/gjoanes/multirat/test/export/qa/anat2template/
cp ${preprocess_folder}/bold_datasink/commonspace_labels/*/*/*.nii.gz /scratch/m/mchakrav/gjoanes/multirat/test/export/snr/atlas/
cp ${preprocess_folder}/confounds_datasink/confounds_csv/*/*/*.csv /scratch/m/mchakrav/gjoanes/multirat/test/export/motion/confound/
cp ${preprocess_folder}/confounds_datasink/FD_csv/*/*/*.csv  /scratch/m/mchakrav/gjoanes/multirat/test/export/motion/FD/


# extract tSNR from atlas. 
ls ${preprocess_folder}/bold_datasink/commonspace_labels | while read label
do 
label_path=$(find ${preprocess_folder}/bold_datasink/commonspace_labels/${label}/*/*.nii.gz)
tSNR_path=$(find ${preprocess_folder}/bold_datasink/tSNR_map_preprocess/${label}/*/*.nii.gz)
mask_path=$(find ${preprocess_folder}/bold_datasink/commonspace_mask/${label}/*/*.nii.gz) 
fslmeants -i ${tSNR_path} -m ${mask_path} --label=${label_path} -o /scratch/m/mchakrav/gjoanes/multirat/test/export/snr/tsnr/${label}.txt
done

```

## Confound removal  and Analysis   
The script belows execute confound regression using 5 different methods, then runs seed-based analysis. 
The confound regression methods and seeds are defined in the pregistration.

#### Confound regression list
1. aroma + 0.01 - 0.1 Hz bandpass filter (aromas)
2. aroma + 0.01 - 0.20 Hz bandpass filter (aromal)
3. aroma + 0.01 - 0.1 Hz bandpass filter + scrubbing  (aromasr)
4. Global signal regression + 0.01 - 0.1 Hz bandpass filter  (GSRs)
5. White matter + CSF + 0.01 - 0.1 Hz bandpass filter (WMCSFs)

#### Seed list
1. S1 barrel field (S1bf)
2. Caudate putamen (CPu)
3. Anterior cingulate area (ACA)
Other lists were excluded amidst concerns that incomplete field-of-view coverage would make the analysis script crash. This turned out later to be not the case (but was instead memory issue with multiprocess plugin).  


``` bash

denoise_list=("aromas" "aromal" "GSRs" "WMCSFs" "aromasr")
denoise_parameters=("--lowpass 0.1 --run_aroma --aroma_dim 10" "--lowpass 0.20 --run_aroma --aroma_dim 10" "--lowpass 0.1 --conf_list global_signal mot_6" "--lowpass 0.1 --conf_list WM_signal CSF_signal mot_6" "--lowpass 0.1 --run_aroma --aroma_dim 10 --apply_scrubbing")

for denoise in ${denoise_list[*]}; do

id=$(echo ${denoise_list[@]/$denoise//} | cut -d/ -f1 | wc -w | tr -d ' ')
param=${denoise_parameters[$id]}



current_confound=$SCRATCH/multirat/confound/${denoise}/${bidsID}
current_analysis=$SCRATCH/multirat/analysis/${denoise}/${bidsID}

mkdir -p $current_confound
mkdir -p $current_analysis

singularity run -B ${PROJECT}/template:/template -B ${bids_folder}:/rabies_in \
-B ${preprocess_folder}:/rabies_out -B ${current_confound}:/confound_out \
/project/m/mchakrav/gjoanes/rabies-0.3.5.simg confound_correction /rabies_out /confound_out \
--read_datasink \
--TR ${TR}s --commonspace_bold \
--highpass 0.01 \
--smoothing_filter 0.5 ${param}

# save important output
cp ${current_confound}/rabies_confound_correction.log /scratch/m/mchakrav/gjoanes/multirat/test/export/log/confound/sub-0101000_ses-1_${denoise}.log

# if using aromas denoising, we want to save the aroma step for validation, as well as the fully cleaned images. 
if [ "$denoise" == "aromas" ]
then 
cp -r ${current_confound}/confound_correction_wf_datasink/cleaned_timeseries/*/*.nii.gz /scratch/m/mchakrav/gjoanes/multirat/test/export/aromas_scan

ls ${current_confound}/confound_correction_wf_datasink/aroma_out | while read label
do 
mkdir -p /scratch/m/mchakrav/gjoanes/multirat/test/export/aromas_qa/${label} 
cp -r ${current_confound}/confound_correction_wf_datasink/aroma_out/${label}/aroma_out/melodic.ica/report/* /scratch/m/mchakrav/gjoanes/multirat/test/export/aromas_qa/${label}
cp -r ${current_confound}/confound_correction_wf_datasink/aroma_out/${label}/aroma_out/classified_motion_ICs.txt /scratch/m/mchakrav/gjoanes/multirat/test/export/aromas_qa/${label}/classification.txt 
done
fi



singularity run -B ${PROJECT}/template/roi:/roi_folder -B ${PROJECT}/template:/template -B ${bids_folder}:/rabies_in \
-B ${preprocess_folder}:/rabies_out -B ${current_confound}:/confound_out -B ${current_analysis}:/analysis_out \
/project/m/mchakrav/gjoanes/rabies-0.3.5.simg analysis /confound_out /analysis_out \
--FC_matrix --seed_list /roi_folder/ACA_l.nii.gz /roi_folder/CPu_l.nii.gz /roi_folder/MOp_l.nii.gz /roi_folder/S1bf_l.nii.gz

# save important output
cp ${current_analysis}/rabies_analysis.log /scratch/m/mchakrav/gjoanes/multirat/test/export/log/analysis/sub-0101000_ses-1_${denoise}.log
seed_dir=/scratch/m/mchakrav/gjoanes/multirat/test/export/seed/${denoise}; mkdir -p ${seed_dir}; cp -r ${current_analysis}/analysis_wf_datasink/seed_correlation_maps/*/*/* $seed_dir
NA_dir=/scratch/m/mchakrav/gjoanes/multirat/test/export/NA/${denoise}; mkdir -p ${NA_dir}; cp -r ${current_analysis}/analysis_wf_datasink/matrix_data_file/*/* $NA_dir

done
```

## Confound remove (minimal) 
The script below is a minimal confound regression script for stimulus-evoked processing track. 

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


# define folder outputs
tar_file=$PROJECT/multirat/ds/${bidsID}.tar.gz
bids_folder=$SCRATCH/multirat/ds/${bidsID}
preprocess_folder=$SCRATCH/multirat/preprocess/${bidsID}


# clean and make directories
rm -rf ${bids_folder}
rm -rf ${preprocess_folder}
mkdir -p ${bids_folder}
mkdir -p ${preprocess_folder}


mkdir -p $SCRATCH/multirat/tar/confound
mkdir -p $SCRATCH/multirat/stim



denoise=stim

current_confound=$SCRATCH/multirat/confound/${denoise}/${bidsID}


rm -rf $current_confound
mkdir -p $current_confound

singularity run -B ${PROJECT}/template:/template -B ${bids_folder}:/rabies_in \
-B ${preprocess_folder}:/rabies_out -B ${current_confound}:/confound_out \
/project/m/###/###/rabies-0.2.1.simg -p MultiProc confound_regression /rabies_out /confound_out \
--TR ${TR}s --commonspace_bold \
--smoothing_filter 0.5 



tar -czvf $SCRATCH/multirat/tar/confound/${bidsID}_${denoise}.tar.gz ${current_confound}

cp ${current_confound}/confound_regression_datasink/cleaned_timeseries/*/*/*.nii.gz $SCRATCH/multirat/stim/

```