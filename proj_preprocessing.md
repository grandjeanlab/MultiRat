MultiRAT analysis code
================
Joanes Grandjean

![rat art](assets/img/rat_art.png)

# Foreword

This is a R markdown file which contains all the code for reproducing my
analysis. The code is meant to be followed step-wise. The raw fMRI
dataset will not be publicly available before the project preprint
publication on BioArxiv. The raw fMRI dataset can be made available
prior to publication upon request and review from the authors.

If re-using some of the scripts, please follow citations guidelines for
the software used. I’ve provided the links to the software wherever
possible. See also the [license](LICENSE.md) for this software.

The code is executed in `bash` (fMRI preprocessing) and `R` (analysis
and plots).

See the [environement](#Environement) section for the required software
libraries. The `R` libraries are organized using
[renv](https://rstudio.github.io/renv/).

To reproduce the code contained within this software, please follow
these steps  
1\. get the required dependencies for [bash](#Bash_environement) and
[R](#R_environement)  
2\. Update the variables `init_folder` and `analysis_folder` for both
bash and R environments

# Environement

The following was written in  
\- [R Version: 3.5.1](https://cran.r-project.org/)  
\- [rstudio Version: 1.3.959](https://rstudio.com/)

## Bash environement

  - [RABIES Version: 0.1.3](https://github.com/CoBrALab/RABIES) build
    from Dockerfile and converted into a Singulariy

Additional software include  
\- FSL Version: 6.0.1  
\- [Bruker2NIfTI
Version: 1.0.20180303](https://github.com/neurolabusc/Bru2Nii)  
\- curl  
\- unzip  
\- rm

Rstudio does not transfer variables between `bash` chunks. Hence, each
chunks needs to reload the environment. To achieve this seamlessly, I
use a [bash\_env.sh](bash_env.sh) file to re-initialize the environment
within each chunk. The content of it should be adapted by the user for
re-use.

``` bash

# Update this section to indicate where the code is kept (init_folder) and where the analysis is performed/stored (analysis_folder). 
init_folder="/home/traaffneu/joagra/code/MultiRat"
analysis_folder="/project/4180000.19/multiRat"


# no need to update the lines below. 
echo 'init_folder='$init_folder > bash_env.sh
echo 'analysis_folder='$analysis_folder >> bash_env.sh

echo 'template=$analysis_folder"/template/WHS_SD_rat_T2star_v1.01.nii.gz"'  >> bash_env.sh
echo 'template_mask=$analysis_folder"/template/WHS_SD_v2_brainmask_bin.nii.gz"' >> bash_env.sh
echo 'template_WM=$analysis_folder"/template/WHS_SD_v2_WM.nii.gz"' >> bash_env.sh
echo 'template_GM=$analysis_folder"/template/WHS_SD_v2_GM.nii.gz"' >> bash_env.sh
echo 'template_CSF=$analysis_folder"/template/WHS_SD_v2_CSF.nii.gz"' >> bash_env.sh
echo 'atlas=$analysis_folder"/template/WHS_SD_rat_atlas_v3.nii.gz"'  >> bash_env.sh
```

# Asset preparation

## Dowload and prepare the template

For this project, I will use the WHS rat template available
[here](https://www.nitrc.org/projects/whs-sd-atlas). The reference for
the template is: Papp EA, Leergaard TB, Calabrese E, Johnson GA, Bjaalie
JG (2014) “Waxholm Space atlas of the Sprague Dawley rat brain”
NeuroImage 97:374-386.
[doi 10.1016/j.neuroimage.2014.04.001](https://doi.org/10.1016/j.neuroimage.2014.04.001)

In this chunk, I download the template and atlas, and generate separate
white (WM), gray (GM), and cerebrospinal fluid (CSF) binary maps. An
additional xml file to make the atlas compatible with FSL was downloaded
from this [board](https://www.nitrc.org/forum/message.php?msg_id=29057),
and made available in the
[assets](assets/atlas/WHS_SD_rat_atlas_v3-FSL.xml)

``` bash
source bash_env.sh

mkdir -p $analysis_folder'/template'

curl https://www.nitrc.org/frs/download.php/9441/WHS_SD_rat_atlas_v2_pack.zip --output $analysis_folder'/template/WHS_SD.zip'
curl https://www.nitrc.org/frs/download.php/9746/WHS_SD_v2_white_gray_mask_clipped.nii.gz --output $analysis_folder'/template/WHS_SD_v2_white_gray_mask_clipped.nii.gz'
curl https://www.nitrc.org/frs/download.php/9748/WHS_SD_v2_brainmask_bin.nii.gz --output $analysis_folder'/template/WHS_SD_v2_brainmask_bin.nii.gz'
curl https://www.nitrc.org/frs/download.php/11404/WHS_SD_rat_atlas_v3.label --output $analysis_folder'/template/WHS_SD_rat_atlas_v3.label'
curl https://www.nitrc.org/frs/download.php/11403/WHS_SD_rat_atlas_v3.nii.gz --output $analysis_folder'/template/WHS_SD_rat_atlas_v3.nii.gz'


unzip -d $analysis_folder'/template/' $analysis_folder'/template/WHS_SD.zip' 
rm $analysis_folder'/template/WHS_SD.zip' 

fslmaths $analysis_folder'/template/WHS_SD_v2_white_gray_mask_clipped.nii.gz' -thr 1 -uthr 1 -bin $analysis_folder'/template/WHS_SD_v2_WM.nii.gz'
fslmaths $analysis_folder'/template/WHS_SD_v2_white_gray_mask_clipped.nii.gz' -thr 2 -uthr 2 -bin $analysis_folder'/template/WHS_SD_v2_GM.nii.gz'
fslmaths $analysis_folder'/template/WHS_SD_v2_white_gray_mask_clipped.nii.gz' -thr 3 -uthr 3 -bin $analysis_folder'/template/WHS_SD_v2_CSF.nii.gz'
```

## Dataset preparation

Datasets included in this study were accepted in any format (bruker,
dicom, nifti, minc). The first step consists of arranging all datasets
within the same convention. I opted for true voxel size and
**A**nterior-**P**osterior axis defined as the rostro-caudal axis. Some
datasets were provided with x10 inflated voxeld and the
**S**uperior-**I**nferior axis defined as the rostro-caudal axis
instead, e.g.:

![raw structrual image](assets/img/orient_pre.png)

These had to be corrected, and organized into
[BIDS](https://bids.neuroimaging.io/) format manually. To do so, I wrote
scripts using a combination of the following FSL and AFNI commands,
`fslinfo`, `fslmerge`, `fslorient`, `fslchpixdim`, `fslswapdim`, and
`3dresample`.

Two scripts used to convert datasets are provided as examples. [Convert
raw Bruker data](assets/script/convert_bruker.sh) and [convert nifti
data](assets/script/convert_nifti.sh). Raw Bruker data were converted
using the [Bruker2NIfTI](https://github.com/neurolabusc/Bru2Nii)
v1.0.20180303 package, written by Matthew Brett, Andrew Janke, Mikaël
Naveau, Chris Rorden. Please note that this software is no longer
supported. New users are invited to try
[BrkRaw](https://github.com/BrkRaw/bruker) instead.

Below is an example of a corrected structural image. Note how the
**S**uperior, **I**nferior, **A**nterior, **P**osterior axis labels are
indicated in `fsleyes`.

![corrected structrual image](assets/img/orient_post.png)

## Dataset preparation limiation

Unfortunately, I cannot ensure the **L**eft / **R**ight axis are
represented correctly across all datasets. While this is less of an
issue for resting-state fMRI, this is a caveat in the stimulus-evoked
fMRI arm of this study, and should be acknowledged as a limitation.
Similarly, I cannot ensure the slicing acquisition order, hence,
preprocessing is performed without slice timing correction.

# Testing RABIES

This is an example of a RABIES preprocessing call for dataset `ds01002`.
Because this is called within a Singularity, the call might differ from
platform to platform.

``` bash
source /home/traaffneu/joagra/code/MultiRat/bash_env.sh

export FSLDIR=/opt/fsl/6.0.1
export PATH=$PATH:$FSLDIR/bin


export ANTSPATH=/home/rabies/ants-v2.3.1/bin 

export RABIES_VERSION=0.1.3-dev
export RABIES=/home/rabies/RABIES-0.1.3-dev
export PYTHONPATH="${PYTHONPATH}:$RABIES"
export PATH=$PATH:$RABIES/bin
export PATH=$PATH:$RABIES/rabies/shell_scripts
export PATH=$RABIES/twolevel_ants_dbm:$PATH


rabies --plugin MultiProc preprocess --no_STC --anat_template $template --brain_mask $template_mask --WM_mask $template_WM --CSF_mask $template_CSF --labels $atlas -r light_SyN --template_reg_script light_SyN --commonspace_resampling 0.35x0.35x0.35 --anatomical_resampling 0.25x0.25x0.25 --autoreg $analysis_folder/data_small/ds01002 $analysis_folder/preprocess/ds01002
```

Practically, within my environment, I need to call the script above as
such.

``` bash
/opt/singularity/3.5.2/bin/singularity exec -B /opt/fsl -B /project/4180000.19/multiRat/ /opt/rabies/0.1.3/rabies-0.1.3-dev.simg bash -c /project/4180000.19/multiRat/script/run/rabies_2.sh 
```

### Quality control

RABIES outputs several QA/QC images that are directly relevant to assess
image registration (Func to Anat, Anat to template, template to
commonspace), and motion parameters.

Below are several examples of Func (top row) to Anat (bottom row) for
several datasets being preprocessed. While the majority of currently
processed scans passed quality control on the basis of Func to Anat
registration, some did not. Either RABIES will need optimization to
improve the registration generalization across datasets, or scans will
need to be excluded. Importantly, the study preregistration did not make
contingency in case some scans must be excluded.

\#\#\#\#Passed QC

![func2anat](assets/QC/sub-0100100_ses-1_run-1_bold_EPI2Anat.png)
![func2anat](assets/QC/sub-0100300_ses-1_run-1_bold_EPI2Anat.png)
![func2anat](assets/QC/sub-0100401_ses-1_run-1_bold_EPI2Anat.png)
![func2anat](assets/QC/sub-0100800_ses-1_run-1_bold_EPI2Anat.png)
![func2anat](assets/QC/sub-0100900_ses-1_run-1_bold_EPI2Anat.png)

\#\#\#\#Failed QC

![func2anat](assets/QC/sub-0100202_ses-1_run-1_bold_EPI2Anat.png)
![func2anat](assets/QC/sub-0100500_ses-1_run-1_bold_EPI2Anat.png)