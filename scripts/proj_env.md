MultiRAT environement preparation
================
Joanes Grandjean

![rat art](../assets/img/rat_art.png)

# Foreword

This and the follow are R markdown files which contains all the code for
reproducing my analysis and detail the process. The code is meant to be
followed step-wise. The raw fMRI dataset will not be publicly available
before the project preprint publication on BioRxiv. The raw fMRI dataset
can be made available prior to publication upon request and review from
the authors.

If re-using some of the scripts, please follow citations guidelines for
the software used. I’ve provided the links to the software wherever
possible. See also the [license](../LICENSE.md) for this software.

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

  - [RABIES Version: 0.2.0-dev](https://github.com/CoBrALab/RABIES)
    build from and converted into a Singularity

Additional software include  
\- FSL Version: 6.0.1  
\- [Bruker2NIfTI
Version: 1.0.20180303](https://github.com/neurolabusc/Bru2Nii)  
\- bash binaries (curl, unzip, rm, mv, cp)

Rstudio does not transfer variables between `bash` chunks. Hence, each
chunks needs to reload the environment. To achieve this seamlessly, I
use a [bash\_env.sh](../bash_env.sh) file to re-initialize the
environment within each chunk. The content of it should be adapted by
the user for re-use.

``` bash

# Update this section to indicate where the code is kept (init_folder) and where the analysis is performed/stored (analysis_folder). 
init_folder="/home/traaffneu/joagra/code/MultiRat"
analysis_folder="/project/4180000.19/multiRat"


# no need to update the lines below. 
echo 'init_folder='$init_folder > bash_env.sh
echo 'analysis_folder='$analysis_folder >> bash_env.sh

echo 'template=$analysis_folder"/template/SIGMA_Wistar_Rat_Brain_TemplatesAndAtlases_Version1.1/SIGMA_Rat_Anatomical_Imaging/SIGMA_Rat_Anatomical_InVivo_Template/SIGMA_InVivo_Brain_Template.nii"'  >> bash_env.sh
echo 'template_mask=$analysis_folder"/template/SIGMA_Wistar_Rat_Brain_TemplatesAndAtlases_Version1.1/SIGMA_Rat_Anatomical_Imaging/SIGMA_Rat_Anatomical_InVivo_Template/SIGMA_InVivo_Brain_Mask.nii"' >> bash_env.sh
echo 'template_WM=$analysis_folder"/template/SIGMA_Wistar_Rat_Brain_TemplatesAndAtlases_Version1.1/SIGMA_Rat_Anatomical_Imaging/SIGMA_Rat_Anatomical_InVivo_Template/SIGMA_InVivo_WM.nii"' >> bash_env.sh
echo 'template_GM=$analysis_folder"/template/SIGMA_Wistar_Rat_Brain_TemplatesAndAtlases_Version1.1/SIGMA_Rat_Anatomical_Imaging/SIGMA_Rat_Anatomical_InVivo_Template/SIGMA_InVivo_GM.nii"' >> bash_env.sh
echo 'template_CSF=$analysis_folder"/template/SIGMA_Wistar_Rat_Brain_TemplatesAndAtlases_Version1.1/SIGMA_Rat_Anatomical_Imaging/SIGMA_Rat_Anatomical_InVivo_Template/SIGMA_InVivo_CSF.nii"' >> bash_env.sh
echo 'atlas=$analysis_folder"/template/SIGMA_Wistar_Rat_Brain_TemplatesAndAtlases_Version1.1/SIGMA_Rat_Brain_Atlases/SIGMA_Anatomical_Atlas/SIGMA_Anatomical_Brain_Atlas.nii"'  >> bash_env.sh
echo 'ROI=$analysis_folder"/template/roi/"'  >> bash_env.sh
```
