![rat art](assets/img/rat_art.png)

# A collaborative rat functional MRI multi-center study.
A Rat fMRI multi-center study
*updated 2021_06_22*

### Executive summary
In this international collaborative project, we seek to gather the rodent imaging community toward performing a rat fMRI multi-center comparison, within the same template as Grandjean et al. NIMG 2020. We want to examine functional connectivity (FC) parameter distribution at the population level within key networks (somatosensory and default-mode network) of the rat brain, as well as establish connectivity sensitivity and specificity in the collected datasets. To do so, we will gather rat BOLD fMRI datasets from individual labs (n=10, any protocol).
In a desire to push further, we seek to obtain additional datasets pertaining to within-laboratory **test-retest** (i.e. scans that were collected in the same control animals within a period of 6 months or less), as well as datasets in **sensory-stimulation** fMRI (limb or whisker stimulation). Finally, we would like to organize a standardized data collection arm of this study for laboratories interested. This would include the **de-novo acquisition** of a dataset (n=10) using predefined, mutually-agreed, and standardized parameters consistent across all participating laboratories.
The end-goal is to make this as an available resource to researchers and to publish an extensive description of the collective dataset in a peer-reviewed journal.  

### Goal
1. To provide a large rat fMRI collective dataset  
2. To describe reference population parameters distributions (e.g. FC, motion, SNR) 
3. To provide evidence-based recommendations for rat fMRI acquisition 
4. To foster collaborations and discussion within the community  
5. Demonstrate test-retest reliability  
6. Investigate the human-factor variability during data acquisition using a standardized protocol.   
7. Organize an international seminar on rodent data standardization at the Donders Institute (NWO Scientific Meetings and Consultations 8. Domain Science) 


### Deliverables
* OSF.io preregistered study
* Publicly available collective dataset on openneuro.org consisting of individual lab dataset
* Publicly available code to replicate the study on github.com
* Bioarxiv preprint
* Journal publication
* Symposium on rodent fMRI standards held at the Donders Institute (if funding is available and COVID permit)

### Analysis
[1. Environement preparation](scripts/proj_env.md)  
[2. Asset preparation](scripts/proj_asset.md)  
[3. Dataset description](scripts/proj_dataset.md)   
[4. Preprocessing code](scripts/proj_preprocessing.md)   
[5. Qality control](scripts/proj_qa.md)   
[6. Analysis tSNR and motion](scripts/proj_analysis_snr.md)  
[7. Analysis seed-based analysis](scripts/proj_analysis_sba.md)    
8. Analysis ica analsysis    
[9. Analysis stimulus evoked](scripts/proj_analysis_stim.md)    

### Links
[License and permissions](LICENSE.md)  
[Collaborative model and project details](scripts/proj_detail.md)  
[Preregistration DOI: 10.17605/OSF.IO/EMQ4B](https://osf.io/emq4b)  
[Lab webpage](https://grandjeanlab.github.io/)  
[Twitter](https://twitter.com/grandjeanlab)  

### Usefull tooboxes
[RABIES](https://github.com/CoBrALab/RABIES), rodent fMRI preprocessing and analysis   
[BkrRaw](https://github.com/BrkRaw/bruker), convert bruker data to [BIDS](https://bids.neuroimaging.io/) format  
[SIGMA template](https://www.nature.com/articles/s41467-019-13575-7)   
[SAMRI](https://github.com/IBT-FMI/SAMRI), another rodent fMRI preprocessing pipeline   
[nirodent](https://github.com/nipreps/nirodents), a toolbox for rodent MRI processing   

### Deviations from preregistration
11.12.2020 - Use SIGMA template instead of WHS  
11.04.2021 - Changed analysis to Python   
11.04.2021 - Analysis using tSNR instead of SNR, because former is readily output in RABIES.  
11.04.2021 - reduced number volumes -> 1200 for ds 1001 (too long preprocessing time)   
11.04.2021 - cropped FOV for ds 1029, 1030, 1036 (improve registrations)   
12.04.2021 - Reduced number of seeds to S1bf, MOp, CPu, ACA because not all dataset had converage along A-P axis, and this seemed to cause RABIES crashes. (spoiler, it wasn't the reason) 
06.06.2021 - cropped FOV for ds 1023, 1038, 1039


