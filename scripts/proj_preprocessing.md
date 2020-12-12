MultiRAT analysis code
================
Joanes Grandjean

![rat art](../assets/img/rat_art.png) \# Foreword The preprocessing is
performed in three broad steps.

1.  Anatomical image preparation  
2.  Dataset-specific anatomical template  
3.  Functional preprocessing, confound removal, and analysis

This separation allows separating the time-consuming functional
preprocessing from the anatomical preprocessingm, and the implementation
of QA and checks at the different levels of the preprocessing.

## Anatomical image preparation

This step consists of bias-field correction in ANTs using the N4
function, denoising, and rigid affine registration to the template.

for every anatomical scan I run the following:  
`N4BiasFieldCorrection -d 3 -i $anat.id -o $n4.id`  
`DenoiseImage -d 3 -i $n4.id -o $denoise.id` `antsIntroduction.sh -d 3
-t RA -i $denoise.id -r $template -o $reg.id,`

The chunk below outputs a `master.sh` file allowing to run all
anatomical preprocessing in parallel using `qsub`, or sequentially if
`use.qsub <- FALSE`.

``` r
library(stringr)

# Load bash environment variable, including asset location needed to make our RABIES call. 
readRenviron("../bash_env.sh")
analysis_folder <- Sys.getenv("analysis_folder")
template <- Sys.getenv("template")
template_mask <- Sys.getenv("template_mask")
template_WM <- Sys.getenv("template_WM")
template_CSF <- Sys.getenv("template_CSF")
atlas <- Sys.getenv("atlas")
ROI <- Sys.getenv("ROI")

study <- read.csv('../assets/table/meta_data.tsv',sep='\t') # Load the meta data (currently not available on public repository.)


ds.select <- c(1001,1002,1003) # Select which dataset to preprocess
use.singularity <- TRUE #select if it is run within a singularity environment (not currently implemented)
use.qsub <- TRUE #Select if you want to submit calls using qsub job submission command (for HPC). 

# Prepare output directories 
dir.create(path = file.path(analysis_folder,'preprocess','anatomical','tmp'), recursive = TRUE, showWarnings = FALSE)
dir.create(path = file.path(analysis_folder,'preprocess','anatomical','indiv_reg'), recursive = TRUE, showWarnings = FALSE)
dir.create(path = file.path(analysis_folder,'preprocess','anatomical','indiv_script'), recursive = TRUE, showWarnings = FALSE)

# Environment to be exported to PATH within the singularity. 
singularity.env <- 'export FSLDIR=/opt/fsl/6.0.1; export PATH=$PATH:$FSLDIR/bin; export ANTSPATH=/home/rabies/ants-v2.3.1/bin; export RABIES_VERSION=0.2.0-dev; export RABIES=/home/rabies/RABIES-0.2.0-dev; export PYTHONPATH="${PYTHONPATH}:$RABIES"; export PATH=$PATH:$RABIES/bin; export PATH=$PATH:$RABIES/rabies/shell_scripts; export PATH=$RABIES/twolevel_ants_dbm:$PATH;'

# The detail for a singularity call. I need to call FSL commands from oustide singularity, hence the /opt/fsl argument
singularity.call <- paste('/opt/singularity/3.5.2/bin/singularity exec -B /opt/fsl -B ',analysis_folder,' /opt/rabies/0.2.0/rabies-0.2.0-dev.simg bash -c \'',sep='')

# Sets the duration, number of processor and memory for HPC job submission. 
qsub.short <- '| qsub -l \'procs=1,mem=8gb,walltime=1:00:00\''
qsub.long <- '| qsub -l \'procs=1,mem=12gb,walltime=72:00:00\''



study.sub<-study[study$rat.ds %in% ds.select, ]

for(i in 1:dim(study.sub[1])){

sub.id <- study.sub$rat.sub[i]
ses.id <- study.sub$rat.ses[i]

output.name <- paste(sub.id, ses.id,sep='_')
output.script <- file.path(analysis_folder,'preprocess','anatomical','indiv_script',paste('run_anat-',Sys.Date(),'-',output.name,'.sh',sep=''))


anat.id <- file.path(analysis_folder,'bids',paste('sub-0',sub.id,sep=''),paste('ses-',ses.id,sep=''),'anat',paste('sub-0',sub.id,'_','ses-',ses.id,'_T2w.nii.gz',sep=''))
n4.id <- file.path(analysis_folder,'preprocess','anatomical','tmp',paste(output.name,'_N4.nii.gz',sep=''))
denoise.id <- file.path(analysis_folder,'preprocess','anatomical','tmp',paste(output.name,'_dn.nii.gz',sep=''))
reg.id <- file.path(analysis_folder,'preprocess','anatomical','tmp',output.name)
deformed.id <- file.path(analysis_folder,'preprocess','anatomical','tmp',paste(output.name,'deformed.nii.gz',sep=''))
rm.id <- file.path(analysis_folder,'preprocess','anatomical','tmp',paste(output.name,'*',sep=''))

anat.preprocess.call <- paste('analysis_folder=',analysis_folder,
                              '; N4BiasFieldCorrection -d 3 -i ',anat.id,' -o ',n4.id,
                              '; DenoiseImage -d 3 -i ',n4.id,' -o ', denoise.id, 
                              '; antsIntroduction.sh -d 3 -t RA -i ',denoise.id, ' -r ', template, ' -o ',reg.id,
                              '; mv ',deformed.id,' ',file.path(analysis_folder,'preprocess','anatomical','indiv_reg'),
                              '; rm ', rm.id
                              , sep='')


 
sink(output.script) 
cat(paste('cd ', file.path(analysis_folder,'preprocess','anatomical','tmp'),sep=''))
cat("\n")
if(use.singularity){cat(c(singularity.call,singularity.env))}
cat(anat.preprocess.call)
if(use.singularity){cat("\'")}
sink()
} 




output.script.master<-file.path(analysis_folder,'preprocess','anatomical','indiv_script',paste('run_anat-',Sys.Date(),'-master.sh',sep=''))

sink(output.script.master) 


for(i in 1:dim(study.sub[1])){

sub.id <- study.sub$rat.sub[i]
ses.id <- study.sub$rat.ses[i]

output.name <- paste(sub.id, ses.id,sep='_')
output.script <- file.path(analysis_folder,'preprocess','anatomical','indiv_script',paste('run_anat-',Sys.Date(),'-',output.name,'.sh',sep=''))

if(use.qsub){cat(paste('echo \"',output.script,'\" ',sep=''),qsub.short)}else{cat(output.script)}
cat("\n")

}
sink()
```

## Dataset-specific anatomical template

The second step consists of creating a study anatomical template for
each dataset and register it to the SIGMA template to use its assets.

at its core the scripts runs the following:  
<br> 1. An study template based on the preprocessed anatomicals from the
chunk above.  
`antsMultivariateTemplateConstruction.sh -d 3 -m 30x20x10 -t SY -s MI
-c 0 -i 3 -j 4 -n 0 -o ./00 $anat.list`  
<br> 2. A linear registration to the SIGMA template.  
`antsRegistration --dimensionality 3 --float 0 -a 0 -v 1 --output
reg/anat2std --interpolation Linear --winsorize-image-intensities
[0.005,0.995] --use-histogram-matching 0 --initial-moving-transform
[std.nii.gz,00template0.nii.gz,1] --transform Rigid[0.1] --metric
MI[std.nii.gz,00template0.nii.gz,1,32,Regular,0.25] --convergence
[10000x5000x2000x1000,1e-6,10] --shrink-factors 8x4x2x1
--smoothing-sigmas 3x2x1x0vox --transform Affine[0.1] --metric
MI[std.nii.gz,00template0.nii.gz,1,32,Regular,0.25] --convergence
[10000x5000x2000x1000,1e-6,10] --shrink-factors 8x4x2x1
--smoothing-sigmas 3x2x1x0vox -x [std_mask.nii.gz]`  
`antsApplyTransforms -i 00template0.nii.gz -r std.nii.gz -t
reg/anat2std0GenericAffine.mat -o 00template0_lin.nii.gz` <br> 3. A
non-linear registration to the SIGMA template.  
`antsRegistration --dimensionality 3 --float 0 -a 0 -v 1 --output
reg/anat2std --interpolation Linear --winsorize-image-intensities
[0.005,0.995] --use-histogram-matching 0 --transform SyN[0.1,3,0]
--metric CC[std.nii.gz,00template0_lin.nii.gz,1,4] --convergence
[500x100x50,1e-6,10] --shrink-factors 4x2x1 --smoothing-sigmas 2x1x0vox
-x [std_mask.nii.gz]`  
`antsApplyTransforms -i 00template0.nii.gz -r std.nii.gz -t
reg/anat2std0GenericAffine.mat -t reg/anat2std0Warp.nii.gz
-o 00template0_nlin.nii.gz`  
<br> 4. A diagnostic image indicative of the registration accuracy.  
`slicer 00template0_nlin.nii.gz std.nii.gz -s 2 -x 0.35 sla.png -x 0.45
slb.png -x 0.55 slc.png -x 0.65 sld.png -y 0.35 sle.png -y 0.45 slf.png
-y 0.55 slg.png -y 0.65 slh.png -z 0.35 sli.png -z 0.45 slj.png -z 0.55
slk.png -z 0.65 sll.png ; pngappend sla.png + slb.png + slc.png +
sld.png + sle.png + slf.png + slg.png + slh.png + sli.png + slj.png +
slk.png + sll.png highres2standard1.png ; slicer
$template 00template0_nlin.nii.gz -s 2 -x 0.35 sla.png -x 0.45 slb.png
-x 0.55 slc.png -x 0.65 sld.png -y 0.35 sle.png -y 0.45 slf.png -y 0.55
slg.png -y 0.65 slh.png -z 0.35 sli.png -z 0.45 slj.png -z 0.55 slk.png
-z 0.65 sll.png ; pngappend sla.png + slb.png + slc.png + sld.png +
sle.png + slf.png + slg.png + slh.png + sli.png + slj.png + slk.png +
sll.png highres2standard2.png ; pngappend highres2standard1.png -
highres2standard2.png anat2standard_nlin.png; rm -f sl?.png
highres2standard2.png; rm highres2standard1.png`

``` r
library(stringr)

# Load bash environment variable, including asset location needed to make our RABIES call. 
readRenviron("../bash_env.sh")
analysis_folder <- Sys.getenv("analysis_folder")
template <- Sys.getenv("template")
template_mask <- Sys.getenv("template_mask")
template_WM <- Sys.getenv("template_WM")
template_CSF <- Sys.getenv("template_CSF")
atlas <- Sys.getenv("atlas")
ROI <- Sys.getenv("ROI")

study <- read.csv('../assets/table/meta_data.tsv',sep='\t') # Load the meta data (currently not available on public repository.)


ds.select <- c(1001,1002,1003) # Select which dataset to preprocess
use.singularity <- TRUE #select if it is run within a singularity environment (not currently implemented)
use.qsub <- TRUE #Select if you want to submit calls using qsub job submission command (for HPC). 

# Prepare output directories 
dir.create(path = file.path(analysis_folder,'preprocess','anatomical','tmp'), recursive = TRUE, showWarnings = FALSE)
dir.create(path = file.path(analysis_folder,'preprocess','anatomical','group_reg'), recursive = TRUE, showWarnings = FALSE)
dir.create(path = file.path(analysis_folder,'preprocess','anatomical','group_script'), recursive = TRUE, showWarnings = FALSE)
dir.create(path = file.path(analysis_folder,'preprocess','qa','group_reg'), recursive = TRUE, showWarnings = FALSE)

# Environment to be exported to PATH within the singularity. 
singularity.env <- 'export FSLDIR=/opt/fsl/6.0.1; export PATH=$PATH:$FSLDIR/bin; export ANTSPATH=/home/rabies/ants-v2.3.1/bin; export RABIES_VERSION=0.2.0-dev; export RABIES=/home/rabies/RABIES-0.2.0-dev; export PYTHONPATH="${PYTHONPATH}:$RABIES"; export PATH=$PATH:$RABIES/bin; export PATH=$PATH:$RABIES/rabies/shell_scripts; export PATH=$RABIES/twolevel_ants_dbm:$PATH;'

# The detail for a singularity call. I need to call FSL commands from oustide singularity, hence the /opt/fsl argument
singularity.call <- paste('/opt/singularity/3.5.2/bin/singularity exec -B /opt/fsl -B ',analysis_folder,' /opt/rabies/0.2.0/rabies-0.2.0-dev.simg bash -c \'',sep='')

# Sets the duration, number of processor and memory for HPC job submission. 
qsub.short <- '| qsub -l \'procs=1,mem=8gb,walltime=1:00:00\''
qsub.long <- '| qsub -l \'procs=1,mem=12gb,walltime=72:00:00\''


ds<-1001
for(ds in ds.select){
study.sub <- study[study$rat.ds %in% ds, ]
  

output.name <- ds
output.script <- file.path(analysis_folder,'preprocess','anatomical','group_script',paste('run_group_anat-',Sys.Date(),'-',output.name,'.sh',sep=''))


anat.list <- paste(file.path(analysis_folder,'preprocess','anatomical','indiv_reg/'),study.sub$rat.sub,'_',study.sub$rat.ses,'deformed.nii.gz',sep='', collapse = ' ')
tmp.dir <- file.path(analysis_folder,'preprocess','anatomical','tmp',ds)
group.reg.dir <- file.path(analysis_folder,'preprocess','anatomical','group_reg',ds)

anat.preprocess.call <- paste('analysis_folder=',analysis_folder,
                              '; mkdir -p ', tmp.dir, '; mkdir -p ', group.reg.dir,'; cd ',tmp.dir,
                              '; cp ',template,' std.nii.gz',
                              '; cp ',template_mask,' std_mask.nii.gz',
                              '; antsMultivariateTemplateConstruction.sh -d 3 -m 30x20x10 -t SY -s MI -c 0 -i 3 -j 4 -n 0 -o ./00 ',  anat.list,
                              '; mkdir -p reg', 
                              '; antsRegistration --dimensionality 3 --float 0 -a 0 -v 1 --output reg/anat2std --interpolation Linear --winsorize-image-intensities [0.005,0.995] --use-histogram-matching 0 --initial-moving-transform [std.nii.gz,00template0.nii.gz,1] --transform Rigid[0.1] --metric MI[std.nii.gz,00template0.nii.gz,1,32,Regular,0.25] --convergence [10000x5000x2000x1000,1e-6,10] --shrink-factors 8x4x2x1 --smoothing-sigmas 3x2x1x0vox --transform Affine[0.1] --metric MI[std.nii.gz,00template0.nii.gz,1,32,Regular,0.25] --convergence [10000x5000x2000x1000,1e-6,10] --shrink-factors 8x4x2x1 --smoothing-sigmas 3x2x1x0vox -x [std_mask.nii.gz]',
                              '; antsApplyTransforms -i 00template0.nii.gz -r std.nii.gz -t reg/anat2std0GenericAffine.mat  -o 00template0_lin.nii.gz',
                              '; slicer 00template0_lin.nii.gz std.nii.gz -s 2 -x 0.35 sla.png -x 0.45 slb.png -x 0.55 slc.png -x 0.65 sld.png -y 0.35 sle.png -y 0.45 slf.png -y 0.55 slg.png -y 0.65 slh.png -z 0.35 sli.png -z 0.45 slj.png -z 0.55 slk.png -z 0.65 sll.png ; pngappend sla.png + slb.png + slc.png + sld.png + sle.png + slf.png + slg.png + slh.png + sli.png + slj.png + slk.png + sll.png highres2standard1.png ; slicer std.nii.gz 00template0_lin.nii.gz -s 2 -x 0.35 sla.png -x 0.45 slb.png -x 0.55 slc.png -x 0.65 sld.png -y 0.35 sle.png -y 0.45 slf.png -y 0.55 slg.png -y 0.65 slh.png -z 0.35 sli.png -z 0.45 slj.png -z 0.55 slk.png -z 0.65 sll.png ; pngappend sla.png + slb.png + slc.png + sld.png + sle.png + slf.png + slg.png + slh.png + sli.png + slj.png + slk.png + sll.png highres2standard2.png ; pngappend highres2standard1.png - highres2standard2.png anat2standard_lin.png; rm -f sl?.png highres2standard2.png; rm highres2standard1.png',
                              '; antsRegistration --dimensionality 3 --float 0 -a 0 -v 1 --output reg/anat2std --interpolation Linear --winsorize-image-intensities [0.005,0.995] --use-histogram-matching 0 --transform SyN[0.1,3,0] --metric CC[std.nii.gz,00template0_lin.nii.gz,1,4] --convergence [500x100x50,1e-6,10] --shrink-factors 4x2x1 --smoothing-sigmas 2x1x0vox -x [std_mask.nii.gz]',
                              '; antsApplyTransforms -i 00template0.nii.gz -r std.nii.gz -t reg/anat2std0GenericAffine.mat -t reg/anat2std0Warp.nii.gz -o 00template0_nlin.nii.gz',
                              '; slicer 00template0_nlin.nii.gz std.nii.gz -s 2 -x 0.35 sla.png -x 0.45 slb.png -x 0.55 slc.png -x 0.65 sld.png -y 0.35 sle.png -y 0.45 slf.png -y 0.55 slg.png -y 0.65 slh.png -z 0.35 sli.png -z 0.45 slj.png -z 0.55 slk.png -z 0.65 sll.png ; pngappend sla.png + slb.png + slc.png + sld.png + sle.png + slf.png + slg.png + slh.png + sli.png + slj.png + slk.png + sll.png highres2standard1.png ; slicer $template 00template0_nlin.nii.gz -s 2 -x 0.35 sla.png -x 0.45 slb.png -x 0.55 slc.png -x 0.65 sld.png -y 0.35 sle.png -y 0.45 slf.png -y 0.55 slg.png -y 0.65 slh.png -z 0.35 sli.png -z 0.45 slj.png -z 0.55 slk.png -z 0.65 sll.png ; pngappend sla.png + slb.png + slc.png + sld.png + sle.png + slf.png + slg.png + slh.png + sli.png + slj.png + slk.png + sll.png highres2standard2.png ; pngappend highres2standard1.png - highres2standard2.png anat2standard_nlin.png; rm -f sl?.png highres2standard2.png; rm highres2standard1.png',
                              ';mv 00template0_nlin.nii.gz ', group.reg.dir,'template.nii.gz',
                              ';mv anat2standard_nlin.png ', file.path(analysis_folder,'preprocess','qa','group_reg'),'/',ds,'.png',
                              '#;rm -r ',tmp.dir
                              , sep='')


 
sink(output.script) 
cat(paste('cd ', file.path(analysis_folder,'preprocess','anatomical','tmp'),sep=''))
cat("\n")
if(use.singularity){cat(c(singularity.call,singularity.env))}
cat(anat.preprocess.call)
if(use.singularity){cat("\'")}
sink()
} 




output.script.master<-file.path(analysis_folder,'preprocess','anatomical','group_script',paste('run_group_anat-',Sys.Date(),'-master.sh',sep=''))

sink(output.script.master) 


for(ds in ds.select){
study.sub <- study[study$rat.ds %in% ds, ]
  

output.name <- ds
output.script <- file.path(analysis_folder,'preprocess','anatomical','group_script',paste('run_group_anat-',Sys.Date(),'-',output.name,'.sh',sep=''))

if(use.qsub){cat(paste('echo \"',output.script,'\" ',sep=''),qsub.long)}else{cat(output.script)}
cat("\n")

}
sink()
```

# RABIES argument justification

RABIES is run in 3 steps (preprocessing, confound\_regression,
analysis). Here I justify my selection of arguments. Below is a full
`preprocessing` command with the augments I used.

`rabies --plugin MultiProc preprocess --no_STC --anat_template $template
--brain_mask $template_mask --WM_mask $template_WM --CSF_mask
$template_CSF --vascular_mask $template_CSF --labels $atlas --autoreg
--commonspace_resampling 0.35x0.35x0.35
--anatomical_resampling 0.25x0.25x0.25 $analysis_folder/bids
$analysis_folder/preprocess`

`--plugin MultiProc` : I run this software on a high-performance cluster
(HPC). Due to software installation limitation, I’ve built RABIES into a
singularity. While using the singularity mode, I cannot use the
traditional job submission (PBS) mode. Users interesting in running the
software on their environment will need to adapt this section.

‘–no\_STC’ : Unfortunately, it is not possible to ensure the slice
acquisition sequence for every dataset reliably. For short-TR sequences
(TR \< 2s), slice timing correction is through to only improve
marginally.

`-anat_template $template --brain_mask $template_mask --WM_mask
$template_WM --CSF_mask $template_CSF --vascular_mask $template_CSF
--labels $atlas` : The templates and masks arguments assume the assets
have been prepared as described in [2. Asset
preparation](../proj_asset.md). Notably, there was not vascular map
available with the template. In [Grandjean et
al. 2020](https://www.sciencedirect.com/science/article/pii/S1053811919308699),
I generated a vascular mask from the data by identifying overlapping ICA
components. It is unclear before preprocessing if I will find such
components in rat. Also, a vascular mask is a prerequisite for running
the `confound_regression` command.

`--autoreg` : as of the version used, this is the recommended
registration argument.

`--commonspace_resampling 0.4x0.4x0.4
--anatomical_resampling 0.25x0.25x0.25` : Due to space / computer time
constrains, I have to run the preprocessing with relative low
resolution. Indeed, I am limited to 72h walltime per script, and some
datasets would not complete running when using high resolutions.

`$analysis_folder/bids $analysis_folder/preprocess` : The standard input
/ output arguments.

# Script preparation.

Because of the variable TR and other parameter, individual dataset are
preprocessed separately. To accommodate that and generate the
preprocessing script, I’ve prepared a R script to generate the
preprocessing code.

``` r
library(stringr)

# Load bash environment variable, including asset location needed to make our RABIES call. 
readRenviron("../bash_env.sh")
analysis_folder <- Sys.getenv("analysis_folder")
template <- Sys.getenv("template")
template_mask <- Sys.getenv("template_mask")
template_WM <- Sys.getenv("template_WM")
template_CSF <- Sys.getenv("template_CSF")
atlas <- Sys.getenv("atlas")
ROI <- Sys.getenv("ROI")


output.name <- '1004-07'
study <- read.csv('../assets/table/meta_data.tsv',sep='\t') # Load the meta data (currently not available on public repository.)
output.script <- file.path(analysis_folder,'script',paste('run_rabies-',Sys.Date(),'-',output.name,'.sh',sep=''))


ds.select <- c(1004,1006,1007) # Select which dataset to preprocess
is.rs <- TRUE # Select if resting-state (TRUE) or stimulus-evoked (FALSE). Important because this will impact confound regression and analysis
use.singularity <- TRUE #select if it is run within a singularity environment (not currently implemented)
use.qsub <- TRUE #Select if you want to submit calls using qsub job submission command (for HPC). 

# resting state confound
confound.list <- c('aroma_s','aroma_l','WMCSF','GSR')
confound.arguments <- c('--lowpass 0.1 --run_aroma', '--lowpass 0.25 --run_aroma','--lowpass 0.1 --conf_list WM_signal CSF_signal mot_6','--lowpass 0.1 --conf_list global_signal mot_6')

seed.list <- c('S1bf_l', 'ACA_l', 'RSP_l','AI_l','MOp_l','CPu_l','dHC_l','AMG_l','TH_l')

# Prepare output directories 
dir.create(path = file.path(analysis_folder,'script'), recursive = TRUE, showWarnings = FALSE)
dir.create(path = file.path(analysis_folder,'tmp'), recursive = TRUE, showWarnings = FALSE)
dir.create(path = file.path(analysis_folder,'preprocess'), recursive = TRUE, showWarnings = FALSE)
dir.create(path = file.path(analysis_folder,'confound'), recursive = TRUE, showWarnings = FALSE)
dir.create(path = file.path(analysis_folder,'analysis'), recursive = TRUE, showWarnings = FALSE)
dir.create(path = file.path(analysis_folder,'qa'), recursive = TRUE, showWarnings = FALSE)
dir.create(path = file.path(analysis_folder,'snr/atlas'), recursive = TRUE, showWarnings = FALSE)
dir.create(path = file.path(analysis_folder,'snr/epi'), recursive = TRUE, showWarnings = FALSE)
dir.create(path = file.path(analysis_folder,'log'), recursive = TRUE, showWarnings = FALSE)

# Environment to be exported to PATH within the singularity. 
singularity.env <- 'export FSLDIR=/opt/fsl/6.0.1; export PATH=$PATH:$FSLDIR/bin; export ANTSPATH=/home/rabies/ants-v2.3.1/bin; export RABIES_VERSION=0.2.0-dev; export RABIES=/home/rabies/RABIES-0.2.0-dev; export PYTHONPATH="${PYTHONPATH}:$RABIES"; export PATH=$PATH:$RABIES/bin; export PATH=$PATH:$RABIES/rabies/shell_scripts; export PATH=$RABIES/twolevel_ants_dbm:$PATH;'

# The detail for a singularity call. I need to call FSL commands from oustide singularity, hence the /opt/fsl argument
singularity.call <- paste('/opt/singularity/3.5.2/bin/singularity exec -B /opt/fsl -B ',analysis_folder,' /opt/rabies/0.2.0/rabies-0.2.0-dev.simg bash -c \'',sep='')

# Sets the duration, number of processor and memory for HPC job submission. 
qsub.short <- ' |  qsub -l \'procs=1,mem=8gb,walltime=1:00:00\''
qsub.long <- ' |  qsub -l \'procs=4,mem=92gb,walltime=72:00:00\''

sink(output.script)

for( ds in ds.select){
  
#ds<-1
study.sub<-study[study$rat.ds == ds, ]
TR <- study.sub$func.TR[1]

mkdir.tmp.call <- paste('mkdir -p ',file.path(analysis_folder,'tmp',ds),sep='')
mv2tmp.call <- paste('cp -r ',file.path(analysis_folder,'bids',paste('sub-0',ds,'*',sep='')),' ', file.path(analysis_folder,'tmp',ds),sep='')
mkdir.qa.call<- paste('mkdir -p ',file.path(analysis_folder,'qa',ds),sep='')
mv2qa.call <- paste('mv ',file.path(analysis_folder,'preprocess',ds,'QC_report/*'),' ', file.path(analysis_folder,'qa',ds),sep='')

mv2atlas.call <- paste('mv ',file.path(analysis_folder,'preprocess',ds,'bold_datasink/bold_labels/*'),' ', file.path(analysis_folder,'snr/atlas'),sep='')
mv2epi.call <- paste('mv ',file.path(analysis_folder,'preprocess',ds,'bold_datasink/corrected_bold/*'),' ', file.path(analysis_folder,'snr/epi'),sep='')
mv2log.call <- paste('mv ',file.path(analysis_folder,'preprocess',ds,'rabies_preprocess.log'),' ', file.path(analysis_folder,'log',paste(ds,'.log',sep='')),sep='')
rm.tmp.call <- paste('rm -r ',file.path(analysis_folder,'tmp',ds),sep='')
rm.preprocess.call <- paste('rm -r ',file.path(analysis_folder,'preprocess',ds),sep='')


cat(paste('# now processing DS ',ds,sep=''))
cat("\n")

cat(mkdir.tmp.call)
if(use.qsub){cat(qsub.short)}
cat("\n")
cat(mkdir.qa.call)
if(use.qsub){cat(qsub.short)}
cat("\n")

cat(mv2tmp.call)
if(use.qsub){cat(qsub.short)}
cat("\n")
cat("\n")

# RABIES preprocessing call
if(use.singularity){
rabies.preprocess.call <- paste('analysis_folder=',analysis_folder,'; rabies --plugin MultiProc --scale_min_memory 10 preprocess --no_STC --anat_template ',template,' --brain_mask ', template_mask, ' --WM_mask ', template_WM, ' --CSF_mask ', template_CSF, ' --vascular_mask ',template_CSF, ' --labels ', atlas, ' --autoreg --commonspace_resampling 0.4x0.4x0.4 --anatomical_resampling 0.25x0.25x0.25 --TR ', TR, 's ', file.path(analysis_folder,'tmp',ds), ' ', file.path(analysis_folder,'preprocess',ds),sep='')
}else{
rabies.preprocess.call <- paste('analysis_folder=',analysis_folder,'; rabies preprocess --no_STC --anat_template ',template,' --brain_mask ', template_mask, ' --WM_mask ', template_WM, ' --CSF_mask ', template_CSF, ' --vascular_mask ',template_CSF, ' --labels ', atlas, ' --autoreg --commonspace_resampling 0.4x0.4x0.4 --anatomical_resampling 0.25x0.25x0.25 --TR ', TR, 's ', file.path(analysis_folder,'tmp',ds), ' ',file.path(analysis_folder,'preprocess',ds),sep='')}

if(use.singularity){cat(c(singularity.call,singularity.env))}
cat(rabies.preprocess.call)
if(use.singularity){cat("\'")}
if(use.qsub){cat(qsub.long)}
cat("\n")
cat("\n")

# Analysis path if dataset is resting-state
if(is.rs){
  for(confound in confound.list){
    
    # Call RABIES for confound regression. Loop across the different options specified
    rabies.confound.call <- paste('analysis_folder=',analysis_folder,'; rabies confound_regression ', file.path(analysis_folder,'preprocess',ds),' ',file.path(analysis_folder,'confound',ds,confound), ' --commonspace_bold --highpass 0.01 ', confound.arguments[which(confound == confound.list)], ' --smoothing_filter 0.5 --diagnosis_output  --TR ', TR, 's',sep='')

    if(use.singularity){cat(c(singularity.call,singularity.env))}
    cat(rabies.confound.call)
    if(use.singularity){cat("\'")}
    if(use.qsub){cat(qsub.long)}
    cat("\n")
    cat("\n")
    
    for(seed in seed.list){
      
      # RABIES analysis call
      rabies.analysis.call <- paste('analysis_folder=',analysis_folder,'; rabies analysis ', file.path(analysis_folder,'confound',ds, confound),' ',file.path(analysis_folder,'analysis',ds, confound), ' --seed_list ', ROI,seed,'.nii.gz',sep='')
      
      if(use.singularity){cat(c(singularity.call,singularity.env))}
      cat(rabies.analysis.call)
      if(use.singularity){cat("\'")}
      if(use.qsub){cat(qsub.long)}
      cat("\n")
      cat("\n")
    }
  }
}

# confound regression for stimulus evoked only high-pass filter. 
if(!is.rs){
  confound<-'stim'
  rabies.confound.call <- paste('analysis_folder=',analysis_folder,'; rabies confound_regression ', file.path(analysis_folder,'preprocess',ds),' ',file.path(analysis_folder,'confound',ds,confound), ' --commonspace_bold --highpass 0.01  --run_aroma --smoothing_filter 0.5 --diagnosis_output  --TR ', TR, 's',sep='')

    if(use.singularity){cat(c(singularity.call,singularity.env))}
    cat(rabies.confound.call)
    if(use.singularity){cat("\'")}
    if(use.qsub){cat(qsub.long)}
    cat("\n")
    cat("\n")
}

cat(mv2qa.call)
if(use.qsub){cat(qsub.short)}
cat("\n")

cat(mv2atlas.call)
if(use.qsub){cat(qsub.short)}
cat("\n")

cat(mv2epi.call)
if(use.qsub){cat(qsub.short)}
cat("\n")

cat(mv2log.call)
if(use.qsub){cat(qsub.short)}
cat("\n")

cat(rm.tmp.call)
if(use.qsub){cat(qsub.short)}
cat("\n")
cat("\n")

cat(rm.preprocess.call)
if(use.qsub){cat(qsub.short)}
cat("\n")
cat("\n")



}

sink()
```
