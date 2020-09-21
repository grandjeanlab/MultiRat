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

## R environement

  - ggplot2 Version: ggplot2\_3.3.2  
  - dplyr Version: dplyr\_1.0.1  
  - here Version: here\_0.1

<!-- end list -->

``` r
# Update this section to indicate where the code is kept (init_folder) and where the analysis is performed/stored (analysis_folder). 
init_folder<-"/home/traaffneu/joagra/code/MultiRat"  
analysis_folder<-"/project/4180000.19/multiRat"  


#load the R libraries 
library('ggplot2')  
library('dplyr')  
library('here')  
```