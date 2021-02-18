MultiRAT Quality Control
================
Joanes Grandjean

![rat art](../assets/img/rat_art.png)

### Quality control

RABIES outputs several QA/QC images that are directly relevant to assess
image registration (func to anat, anat to template, template to
standard), and motion parameters.

Below are several examples of QA/QC registration plots that have failed
the test. I have re-ran preprocessing for datasets with at least 2
failed QA/QC scans, or datasets that entirely failed QA/QC, as
registrations steps are non-deterministic. Failed registrations happen
despite consequential time was devoted to optimize the procedure. These
are often due to poor image quality or strong image artifacts.
Importantly, the study preregistration did not make contingency in case
some scans must be excluded.

``` r
df <- read.csv2('../assets/table/meta_data.tsv',sep='\t')
df.sub <- df[df$exclude == 'yes',]
summary(df.sub$exclude.reason)
```

    ##                     anat2template      empty files         epi2anat 
    ##                0                2                1                8 
    ## inconsistant FOV     template2std 
    ##                3               20

``` r
summary(as.factor(df.sub$rat.ds))
```

    ## 2001 2003 2004 2005 2006 2007 2008 
    ##    2    2   10   10    5    3    2

Here is a summary of the excluded datasets, grouped by exclusion reasons
and datasets. Currently including only stimulus-evoked datasets.

Below are the detailed failed QA/QC tests

#### Failed template2std registration

![func2anat](../assets/QC/template2std/02004.png)
![func2anat](../assets/QC/template2std/02005.png)

#### Failed anat2tempalte registration

![func2anat](../assets/QC/anat2tempalte/sub-0200307_ses-1_T2w_registration.png)

![func2anat](../assets/QC/anat2tempalte/sub-0200309_ses-1_T2w_registration.png)

#### Failed epi2anat registration

![func2anat](../assets/QC/epi2anat/sub-0200103_ses-1_run-1_bold_registration.png)

![func2anat](../assets/QC/epi2anat/sub-0200106_ses-1_run-1_bold_registration.png)

![func2anat](../assets/QC/epi2anat/sub-0200602_ses-1_run-1_bold_registration.png)

![func2anat](../assets/QC/epi2anat/sub-0200604_ses-1_run-1_bold_registration.png)

![func2anat](../assets/QC/epi2anat/sub-0200606_ses-1_run-1_bold_registration.png)

![func2anat](../assets/QC/epi2anat/sub-0200608_ses-1_run-1_bold_registration.png)

![func2anat](../assets/QC/epi2anat/sub-0200609_ses-1_run-1_bold_registration.png)

![func2anat](../assets/QC/epi2anat/sub-0200806_ses-1_run-1_bold_registration.png)
