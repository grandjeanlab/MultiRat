MultiRAT asset preparation
================
Joanes Grandjean

![rat art](../assets/img/rat_art.png)

1.  [Functional preprocessing](#functional-preprocessing)  
2.  [Confound removal and Analysis](#confound-removal-and-analysis)   
3.  [Confound remove (minimal)](#confound-remove-minimal)

Currently this is achieved all through a single script submitted as a
sbatch job on a SLURM-running super computer (niagara.computecanada.ca),
and submitted with the following command `sbatch
--export=bidsID='02004',TR='1' script/job_mr_stim.sh`

Importantly, in this script, RABIES is run in a singularity. This was
built using `singularity build rabies-0.2.1.simg docker://gabdesgreg/rabies:0.2.1`.

For latest updates on RABIES, consult the [Github](https://github.com/CoBrALab/RABIES)

In some instances, running RABIES in parralel using the `-p MultiProc` is known to cause some crashes. If so, try running in linear (default) mode. 

Please report issues [here](https://github.com/CoBrALab/RABIES/issues)

## Preregistration
From the preregistration: 
"All scans will be preprocessed using RABIES (https://github.com/CoBrALab/RABIES), a BIDS-based software based on the fMRIprep pipeline. Data will be co-registered into Waxholm space (https://scalablebrainatlas.incf.org/rat/PLCJB14). Denoising will be performed using motion regression concurently with either: global signal regression, white matter+CSF signal regression, or  ICA-AROMA automatic method (https://github.com/maartenmennes/ICA-AROMA) adapted for the rodent and RABIES (https://github.com/Gab-D-G/conf_reg_pkg). Temporal filtering will be applied at 0.01-0.1 Hz for all scans (3dbandpass). To account for anesthesia requiring larged bandpass filter, we will also try a 0.01-0.25 Hz filter. Smoothing will be applied at 0.5mm for all scans (3dBlurInMask)."   


## Functional preprocessing  
The script belows initiate directories and run the preprocessing. This is common to both resting-state and stimulus-evoked datasets. 

The following options are justified:   
`--template_reg_script multiRAT`, Coregistrations are criticial. I've optimized the anatomical to template registration script in  a dedicated script for this study. In my hands, it is more robust and fast than Syn or light_SyN registrations. 
`--coreg_script light_SyN`, This is the hardest trade-off. Indeed, most of the excluded scans are done so due to failed EPI to anatomical registration. Using purely rigid registrations would ensure minimal exclusions, but I find that light non-linear registration helps enhance the fit between the (slightly) distorted EPI and their anatomical images. 
`--no_STC`, I find it impossible to know the slice acquisition order. Most users would know approximately, especially for older datasets. Slice timing correction has been shown to be necessary for longer TR, but only have minimal impact on low TR. Since the maximum TR is in the ballpark of 2 s, I've decided to omit this step to exclude the risk of inducing artifacts. 
`--commonspace_resampling 0.3x0.3x0.3`, This is exclusively to try to save disk space. This is selected on the basis of the resolutions found among most datasets, and my need to limit disk usage.  

Quality metrics are outputed and carefully examined for each dataset before proceeding to the condound regression and analysis steps. 

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


# define folder outputs
tar_file=$PROJECT/multirat/ds/${bidsID}.tar.gz
bids_folder=$SCRATCH/multirat/ds/${bidsID}
preprocess_folder=$SCRATCH/multirat/preprocess/${bidsID}


# clean and make directories
rm -rf ${bids_folder}
rm -rf ${preprocess_folder}
mkdir -p ${bids_folder}
mkdir -p ${preprocess_folder}

# make output directories
mkdir -p $SCRATCH/multirat/qa/epi2anat
mkdir -p $SCRATCH/multirat/qa/anat2template
mkdir -p $SCRATCH/multirat/qa/template2std
mkdir -p $SCRATCH/multirat/qa/temporal
mkdir -p $SCRATCH/multirat/snr/atlas
mkdir -p $SCRATCH/multirat/snr/tsnr/pre
mkdir -p $SCRATCH/multirat/motion/confound
mkdir -p $SCRATCH/multirat/motion/FD
mkdir -p $SCRATCH/multirat/log
mkdir -p $SCRATCH/multirat/tar/confound
mkdir -p $SCRATCH/multirat/tar/analysis
mkdir -p $SCRATCH/multirat/stim



#untar the tarball containing the specified dataset
tar -xvzf ${tar_file} -C ${bids_folder} 

singularity run -B ${PROJECT}/template:/template -B ${bids_folder}:/rabies_in:ro \
-B ${preprocess_folder}:/rabies_out \
/project/m/###/###/rabies-0.2.1.simg -p MultiProc preprocess /rabies_in /rabies_out \
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
--commonspace_resampling 0.3x0.3x0.3


# copy preprocessing outputs of interest to designated folders
cp ${preprocess_folder}/rabies_preprocess.log $SCRATCH/multirat/log/${bidsID}.log
cp ${preprocess_folder}/QC_report/EPI2Anat/* $SCRATCH/multirat/qa/epi2anat/
cp ${preprocess_folder}/QC_report/Anat2Template/* $SCRATCH/multirat/qa/anat2template/
cp ${preprocess_folder}/QC_report/temporal_diagnosis/* $SCRATCH/multirat/qa/temporal/
cp ${preprocess_folder}/QC_report/Template2Commonspace/_registration.png $SCRATCH/multirat/qa/template2std/${bidsID}.png

cp ${preprocess_folder}/bold_datasink/bold_labels/*/*/*.nii.gz $SCRATCH/multirat/snr/atlas/
cp ${preprocess_folder}/confounds_datasink/confounds_csv/*/*/*.csv $SCRATCH/multirat/motion/confound/
cp ${preprocess_folder}/confounds_datasink/FD_csv/*/*/*.csv $SCRATCH/multirat/motion/FD/


# extract tSNR from atlas. 
ls ${preprocess_folder}/bold_datasink/commonspace_bold_labels | while read label
do

label_path=$(find ${preprocess_folder}/bold_datasink/commonspace_bold_labels/${label}/*/*.nii.gz)
tSNR_path=$(find ${preprocess_folder}/bold_datasink/tSNR_filename/${label}/*/*.nii.gz)
mask_path=$(find ${preprocess_folder}/bold_datasink/commonspace_bold_mask/${label}/*/*.nii.gz)

fslmeants -i ${tSNR_path} \
-m ${mask_path} \
--label=${label_path} \
-o $SCRATCH/multirat/snr/tsnr/pre/${label}.txt

done

```

## Confound removal  and Analysis   
The script belows execute confound regression using 5 different methods, then runs seed-based analysis. 
The confound regression methods and seeds are defined in the pregistration.

#### Confound regression list
1. aroma + 0.01 - 0.1 Hz bandpass filter (aromas)
2. aroma + 0.01 - 0.25 Hz bandpass filter (aromal)
3. aroma + 0.01 - 0.1 Hz bandpass filter + scrubbing  (aromasr)
4. Global signal regression + 0.01 - 0.1 Hz bandpass filter  (GSRs)
5. White matter + CSF + 0.01 - 0.1 Hz bandpass filter (WMCSFs)

#### Seed list
1. S1 barrel field (S1bf)
2. Caudate putamen (CPu)
3. Anterior cingulate area (ACA)
Other lists were excluded amidst concerns that incomplete field-of-view coverage would make the analysis script crash. This turned out later to be not the case (but was instead memory issue with multiprocess plugin).  


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


# define folder outputs
tar_file=$PROJECT/multirat/ds/${bidsID}.tar.gz
bids_folder=$SCRATCH/multirat/ds/${bidsID}
preprocess_folder=$SCRATCH/multirat/preprocess/${bidsID}


denoise_list=("aromas" "aromal" "GSRs" "WMCSFs" "aromasr")
denoise_parameters=("--lowpass 0.1 --run_aroma --aroma_dim 10" "--lowpass 0.25 --run_aroma --aroma_dim 10" "--lowpass 0.1 --conf_list global_signal mot_6" "--lowpass 0.1 --conf_list WM_signal CSF_signal mot_6" "--lowpass 0.1 --run_aroma --aroma_dim 10 --apply_scrubbing")

for denoise in ${denoise_list[*]}; do

id=$(echo ${denoise_list[@]/$denoise//} | cut -d/ -f1 | wc -w | tr -d ' ')
param=${denoise_parameters[$id]}



current_confound=$SCRATCH/multirat/confound/${denoise}/${bidsID}
current_analysis=$SCRATCH/multirat/analysis/${denoise}/${bidsID}

mkdir -p $current_confound
mkdir -p $current_analysis

singularity run -B ${PROJECT}/template:/template -B ${bids_folder}:/rabies_in \
-B ${preprocess_folder}:/rabies_out -B ${current_confound}:/confound_out \
/project/m/###/###/rabies-0.2.1.simg -p MultiProc confound_regression /rabies_out /confound_out \
--TR ${TR}s --commonspace_bold \
--highpass 0.01 \
--smoothing_filter 0.5 ${param}


tar -czvf $SCRATCH/multirat/tar/confound/${bidsID}_${denoise}.tar.gz ${current_confound}


echo '' > ${PROJECT}/roi_list_burner.txt 
ls ${PROJECT}/template/roi/*_l.nii.gz | while read roi_file
do

roi_output='/roi_folder/'$(basename $roi_file)
echo -n $roi_output" " >> ${PROJECT}/roi_list_burner.txt 
done 
roi_file_final=$(cat ${PROJECT}/roi_list_burner.txt)

rm ${PROJECT}/roi_list_burner.txt 

singularity run -B ${PROJECT}/template/roi:/roi_folder -B ${PROJECT}/template:/template -B ${bids_folder}:/rabies_in \
-B ${preprocess_folder}:/rabies_out -B ${current_confound}:/confound_out -B ${current_analysis}:/analysis_out \
/project/m/###/###/rabies-0.2.1.simg -p MultiProc analysis /confound_out /analysis_out \
--seed_list ${roi_file_final} --TR ${TR}s 

tar -czvf $SCRATCH/multirat/tar/analysis/${bidsID}_${denoise}.tar.gz ${current_analysis}

seed_dir=$SCRATCH/multirat/seed/${denoise}
seed_txt_dir=$SCRATCH/multirat/seed_txt/${denoise}
#matrix_dir=$SCRATCH/multirat/matrix/${denoise}


mkdir -p ${seed_dir}
mkdir -p ${seed_txt_dir}
#mkdir -p ${matrix_dir}



cp -r ${current_analysis}/analysis_datasink/seed_correlation_maps/*/*/*/* $seed_dir
#cp -r ${current_analysis}/analysis_datasink/matrix_data_file/*/*/* $matrix_dir


#for every right hemisphere roi, extract the FC parameter from its correponding left hemisphere seed map
ls ${PROJECT}/template/roi/*_r.nii.gz | while read roi_file
do
roi_r=$(remove_ext $(basename $roi_file))
roi_l=$(echo "$roi_r" | tr r l)
ls ${current_analysis}/analysis_datasink/seed_correlation_maps/*/*/_seed_name_${roi_l}/* | while read seed_map
do
seed_name=$(remove_ext $(basename $seed_map))
seed_output=${seed_txt_dir}/${seed_name}'_'${roi_r}'.txt'
flirt -in ${roi_file} -ref ${seed_map} -out tmp.nii.gz -applyxfm -interp nearestneighbour
fslmeants -i ${seed_map} -m tmp.nii.gz -o ${seed_output}
rm -rf tmp.nii.gz
done
done

#for every S1bf seed map, extract the FC parameter within the ACA
ls ${current_analysis}/analysis_datasink/seed_correlation_maps/*/*/_seed_name_S1bf_l/* | while read seed_map
do
roi_file=${PROJECT}/template/roi/ACA_l.nii.gz
roi_r='ACA_l'
seed_name=$(remove_ext $(basename $seed_map))
seed_output=${seed_txt_dir}/${seed_name}'_'${roi_r}'.txt'
flirt -in ${roi_file} -ref ${seed_map} -out tmp.nii.gz -applyxfm -interp nearestneighbour
fslmeants -i ${seed_map} -m tmp.nii.gz -o ${seed_output}
rm -rf tmp.nii.gz
done

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