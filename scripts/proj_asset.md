MultiRAT asset preparation
================
Joanes Grandjean

![rat art](../assets/img/rat_art.png)

## Dowload and prepare the template

For this project, I will use the SIGMA rat template available
[here](https://www.nitrc.org/projects/sigma_template) The reference for
the template is: Barrière, D.A., Magalhães, R., Novais, A. et al. The
SIGMA rat brain templates and atlases for multimodal MRI data analysis
and visualization. Nat Commun 10, 5699 (2019).
<https://doi.org/10.1038/s41467-019-13575-7>

Originally, I intended to use the [WHS
atlas](doi%2010.1016/j.neuroimage.2014.04.001) but differences in
contrast with the majority of the dataset, ex vivo template with empty
ventricles, and image artifacts rendered image registration complicated.
This is a deviation from the [preregistration](https://osf.io/emq4b).


From the preregistration: 
"All scans will be converted to NIFTI with original voxel size. Axis labels will be swapped so that NIFTI SI / AP / LR labels correspond to the right orientation. Scans will be organized using the Brain Imaging Data Structure (BIDS) format (https://bids.neuroimaging.io/).

Meta-data (see measured variables) will be kept in a master table in a tab-separated format, with corresponding JSON file according to BIDS format."   


```python
# init variables
init_folder='/home/traaffneu/joagra/code/MultiRat'
analysis_folder='/project/4180000.19/multiRat'

```


```python
# Create a template directory where the template will be stored. 
! mkdir -p $analysis_folder'/template'

# the following requires NITRC login credentials
! wget https://www.nitrc.org/frs/download.php/11708/SIGMA_Wistar_Rat_Brain_TemplatesAndAtlases_Version1.1.zip -O $analysis_folder'/template/tmp.zip'
! unzip -d $analysis_folder'/template/' $analysis_folder'/template/tmp.zip' 
! rm $analysis_folder'/template/tmp.zip' 
```

    --2021-04-05 00:03:35--  https://www.nitrc.org/frs/download.php/11708/SIGMA_Wistar_Rat_Brain_TemplatesAndAtlases_Version1.1.zip
    Resolving www.nitrc.org (www.nitrc.org)... 52.3.190.103
    Connecting to www.nitrc.org (www.nitrc.org)|52.3.190.103|:443... connected.
    HTTP request sent, awaiting response... 302 Found
    Location: /account/login.php?return_to=%2Ffrs%2Fdownload.php%2F11708%2FSIGMA_Wistar_Rat_Brain_TemplatesAndAtlases_Version1.1.zip&feedback=The+tool%2Fresource+administrator+has+requested+that+you+log+in+to+download+this+file. [following]
    --2021-04-05 00:03:35--  https://www.nitrc.org/account/login.php?return_to=%2Ffrs%2Fdownload.php%2F11708%2FSIGMA_Wistar_Rat_Brain_TemplatesAndAtlases_Version1.1.zip&feedback=The+tool%2Fresource+administrator+has+requested+that+you+log+in+to+download+this+file.
    Reusing existing connection to www.nitrc.org:443.
    HTTP request sent, awaiting response... 200 OK
    Length: unspecified [text/html]
    Saving to: ‘/project/4180000.19/multiRat/template/tmp.zip’
    
        [ <=>                                   ] 39.796      --.-K/s   in 0,09s   
    
    2021-04-05 00:03:35 (441 KB/s) - ‘/project/4180000.19/multiRat/template/tmp.zip’ saved [39796]
    
    unzip:  cannot find or open /project/4180000.19/multiRat/template/SIGMA_Wistar_Rat_Brain_TemplatesAndAtlases_Version1.1.zip, /project/4180000.19/multiRat/template/SIGMA_Wistar_Rat_Brain_TemplatesAndAtlases_Version1.1.zip.zip or /project/4180000.19/multiRat/template/SIGMA_Wistar_Rat_Brain_TemplatesAndAtlases_Version1.1.zip.ZIP.


## ROI generation

Here, I create regions-of-interest for seed-based analysis. The
preregistration specifies the seeds should be 0.9 mm3 and placed on both
hemispheres. It specifies the following ROIs: S1 barrel field area,
Cingulate area, Retrosplenial area, Insula area, motor area,
caudate-putamen, dorsal hippocampus, amygdala, thalamus.

However, because some datasets are provided with restricted FOV in the AP axis, I had to restrict the analysis to: 
S1 barrel field area, Cingulate area, Motor area, and Caudate-putamen



```python
from os import path, makedirs

template = path.join(path.sep, analysis_folder, "template/SIGMA_Wistar_Rat_Brain_TemplatesAndAtlases_Version1.1/SIGMA_Rat_Anatomical_Imaging/SIGMA_Rat_Anatomical_InVivo_Template/SIGMA_InVivo_Brain_Template.nii")
template_mask = path.join(analysis_folder,"template/SIGMA_Wistar_Rat_Brain_TemplatesAndAtlases_Version1.1/SIGMA_Rat_Anatomical_Imaging/SIGMA_Rat_Anatomical_InVivo_Template/SIGMA_InVivo_Brain_Mask.nii")
ROI = path.join(analysis_folder,"template/roi2/")

```


```python
# JG 11.12.2020, update ROI to new template specs. 
makedirs(ROI)

roi_name=path.join(ROI,'S1bf_l')
! fslmaths $template -mul 0 -add 1 -roi 33 9 237 9 134 9 0 1 $roi_name -odt int
roi_name=path.join(ROI,'S1bf_r')
! fslmaths $template -mul 0 -add 1 -roi 147 9 237 9 134 9 0 1 $roi_name -odt int

roi_name=path.join(ROI,'ACA_l')
! fslmaths $template -mul 0 -add 1 -roi 92 9 256 9 143 9 0 1 $roi_name -odt int
roi_name=path.join(ROI,'RSP_l')
! fslmaths $template -mul 0 -add 1 -roi 89 9 208 9 157 9 0 1 $roi_name -odt int

roi_name=path.join(ROI,'AI_l')
! fslmaths $template -mul 0 -add 1 -roi 45 9 268 9 99 9 0 1 $roi_name -odt int
roi_name=path.join(ROI,'AI_r')
! fslmaths $template -mul 0 -add 1 -roi 139 9 268 9 99 9 0 1 $roi_name -odt int

roi_name=path.join(ROI,'MOp_l')
! fslmaths $template -mul 0 -add 1 -roi 65 9 268 9 143 9 0 1 $roi_name -odt int
roi_name=path.join(ROI,'MOp_r')
! fslmaths $template -mul 0 -add 1 -roi 126 9 268 9 143 9 0 1 $roi_name -odt int

roi_name=path.join(ROI,'CPu_l')
! fslmaths $template -mul 0 -add 1 -roi 62 9 250 9 112 9 0 1 $roi_name -odt int
roi_name=path.join(ROI,'CPu_r')
! fslmaths $template -mul 0 -add 1 -roi 129 9 250 9 112 9 0 1 $roi_name -odt int

roi_name=path.join(ROI,'dHC_l')
! fslmaths $template -mul 0 -add 1 -roi 63 9 209 9 142 9 0 1 $roi_name -odt int
roi_name=path.join(ROI,'dHC_r')
! fslmaths $template -mul 0 -add 1 -roi 122 9 209 9 142 9 0 1 $roi_name -odt int

roi_name=path.join(ROI,'AMG_l')
! fslmaths $template -mul 0 -add 1 -roi 42 9 215 9 80 9 0 1 $roi_name -odt int
roi_name=path.join(ROI,'AMG_r')
! fslmaths $template -mul 0 -add 1 -roi 145 9 215 9 80 9 0 1 $roi_name -odt int

roi_name=path.join(ROI,'TH_l')
! fslmaths $template -mul 0 -add 1 -roi 63 9 215 9 113 9 0 1 $roi_name -odt int
roi_name=path.join(ROI,'TH_r')
! fslmaths $template -mul 0 -add 1 -roi 122 9 215 9 113 9 0 1 $roi_name -odt int
```

## Dataset preparation

Datasets included in this study were accepted in any format (bruker,
dicom, nifti, minc). The first step consists of arranging all datasets
within the same convention. I opted for true voxel size and
**A**nterior-**P**osterior axis defined as the rostro-caudal axis. Some
datasets were provided with x10 inflated voxels and the
**S**uperior-**I**nferior axis defined as the rostro-caudal axis
instead, e.g.:

![raw structrual image](../assets/img/orient_pre.png)

These had to be corrected, and organized into
[BIDS](https://bids.neuroimaging.io/) format manually. To do so, I wrote
scripts using a combination of the following FSL and AFNI commands,
`fslinfo`, `fslmerge`, `fslorient`, `fslchpixdim`, `fslswapdim`, and
`3dresample`.

Two scripts used to convert datasets are provided as examples. [Convert
raw Bruker data](../assets/script/convert_bruker.sh) and [convert nifti
data](../assets/script/convert_nifti.sh). Raw Bruker data were converted
using the [Bruker2NIfTI](https://github.com/neurolabusc/Bru2Nii)
v1.0.20180303 package, written by Matthew Brett, Andrew Janke, Mikaël
Naveau, Chris Rorden. Please note that this software is no longer
supported. New users are invited to try
[BrkRaw](https://github.com/BrkRaw/bruker) instead.

Below is an example of a corrected structural image. Note how the
**S**uperior, **I**nferior, **A**nterior, **P**osterior axis labels are
indicated in `fsleyes`.

![corrected structrual image](../assets/img/orient_post.png)

## Dataset preparation limiation

Unfortunately, I cannot ensure the **L**eft / **R**ight axis are
represented correctly across all datasets. While this is less of an
issue for resting-state fMRI, this is a caveat in the stimulus-evoked
fMRI arm of this study, and should be acknowledged as a limitation.
Similarly, I cannot ensure the slicing acquisition order, hence,
preprocessing is performed without slice timing correction.
