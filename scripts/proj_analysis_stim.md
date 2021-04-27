MultiRAT analysis code
================
Joanes Grandjean

![rat art](../assets/img/rat_art.png)

# Foreword

This is a R markdown file which contains all the code for reproducing my
analysis. The code is meant to be followed step-wise. The raw fMRI
dataset will not be publicly available before the project preprint
publication on BioArxiv. The raw fMRI dataset can be made available
prior to publication upon request and review from the authors.

If re-using some of the scripts, please follow citations guidelines for
the software used. I’ve provided the links to the software wherever
possible. See also the [license](../LICENSE.md) for this software.



```python
# init variables
init_folder='/home/traaffneu/joagra/code/MultiRat'
analysis_folder='/project/4180000.19/multiRat'
```


```python
import os
import glob
import pandas as pd
import numpy as np

df = pd.read_csv('../assets/table/meta_data_20210411_snr.tsv', sep='\t')
```


```python
stim_list = glob.glob((os.path.join(analysis_folder, 'scratch', 'stim'))+'/*')
```


```python
from nilearn.input_data import NiftiMasker
import numpy as np
nifti_mask='/project/4180000.19/multiRat/template/roi/S1bf_r.nii.gz'
img='/project/4180000.19/multiRat/scratch/seed/aromas/sub-0100101_ses-1_run-1_bold_RAS_combined_aroma_cleaned_S1bf_l_corr_map.nii.gz'

NiftiMasker(nifti_mask).fit_transform(img).mean()

```




    0.49358446711758747




```python
from nilearn.image import resample_to_img
resample_to_img(stat_img, template)
```


```python
ts = ImageMeants()
ts.in_file = img
ts.mask = '/project/4180000.19/multiRat/template/roi/S1bf_r_RS.nii.gz'
ts.out_file = '/project/4180000.19/multiRat/template/roi/S1bf_'
```




    array([[1., 0., 0.],
           [0., 1., 0.],
           [0., 0., 1.]])


