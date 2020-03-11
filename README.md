# A collaborative rat functional MRI multi-center study.
A Rat fMRI multi-center study
*updated 2020_03_10*

### 1.1 Executive summary
In this international collaborative project, we seek to gather the rodent imaging community toward performing a rat fMRI multi-center comparison, within the same template as Grandjean et al. NIMG 2020. We want to examine functional connectivity (FC) parameter distribution at the population level within key networks (somatosensory and default-mode network) of the rat brain, as well as establish connectivity sensitivity and specificity in the collected datasets. To do so, we will gather rat BOLD fMRI datasets from individual labs (n=10, any protocol).
In a desire to push further, we seek to obtain additional datasets pertaining to within-laboratory **test-retest** (i.e. scans that were collected in the same control animals within a period of 6 months or less), as well as datasets in **sensory-stimulation** fMRI (limb or whisker stimulation). Finally, we would like to organize a standardized data collection arm of this study for laboratories interested. This would include the **de-novo acquisition** of a dataset (n=10) using predefined, mutually-agreed, and standardized parameters consistent across all participating laboratories.
The end-goal is to make this as an available resource to researchers and to publish an extensive description of the collective dataset in a peer-reviewed journal.  

### 1.2 Goal
1. To provide a large rat fMR collective dataset
2. To describe reference population parameters distributions (e.g. FC, motion, SNR) 
3. To provide evidence-based recommendations for rat fMRI acquisition 
4. To foster collaborations and discussion within the community
5. Demonstrate test-retest reliability
6. Investigate the human-factor variability during data acquisition using a standardized protocol. 
7. Organize an international seminar on rodent data standardization at the Donders Institute (NWO Scientific Meetings and Consultations 8. Domain Science) 

### 1.3 Collaborative model
Every laboratory is invited to participate with 1 or multiple rat fMRI dataset (see **1.4 Data specification**). Invitation to collaborate is made through email, announced at conferences (ESMI, ISMRM), announced on [Twitter](https://twitter.com/grandjeanlab). See the list of currently listed laboratories (**section 2**). 

The collaborative model is as follows: 
- Every contributing laboratory can nominate any reasonable number of junior and senior staff. Junior staff who have been involved in data collection. Senior staff who have been involved in project coordination and funding acquisition. 
- Joanes Grandjean ([JG](https://grandjeanlab.github.io/pages/contact.html)) will put the collective dataset together and perform the primary analysis. 
- Andy Hess (AH) will oversee and coordinate the project. 
- Gabriel Desrosiers-Grégoire (GD) and M. Mallar Chakravarty (MMC) will provide the preprocessing software (RABIES)
- Every collaborator is invited to further contribute to the preregistration proposal, the analysis, and manuscript preparation and will be recognized in their contributions accordingly. 

In the ensuing journal publication, author list will be as such: 
JG, GD, junior with an extra contribution, [junior collaborator in alphabetical order], [senior collaborator in alphabetical order], senior with an extra contribution, MMC, AH. 
 
In addition to the main study, collaborators are invited to propose spin-off studies relying on the dataset, which will be provided by JG. The spin-off studies should respect an embargo period on publication, as agreed together with JG+AH. E.g. a comparison of connectivity between mice and rats using graph theory (AH). 

### 1.4 Data specification
Individual datasets each consist of an n=10 resting-state fMRI and corresponding anatomical scan obtained in rats, any strain, any gender, any age, any weight, any anesthesia or awake condition, acquired at any field strength, any vendor, any coil, with gradient-echo or spin-echo fMRI sequence. The anatomical scans are preferably T2 or T2* weighted. In addition to providing the scans in a nifti-convertible format, the dataset also consists of relevant meta-data (see **section 3.3.2 Measured variables**). Specific data pre-processing requirements (fieldmap correction or top-up need to be discussed separately with JG). Similar to [Grandjean et al 2020](https://www.sciencedirect.com/science/article/pii/S1053811919308699), the individual dataset will be anonymized with respect to the laboratory where it was produced.

In addition, we would like to examine test-retest, that is datasets (n=10) where the same rats are imaged within 6 months with the same imaging protocol. We also seek to gather stimulus-evoked fMRI dataset (n=10) from sensory stimulation, to compare localization and effect amplitude. 

Finally, we seek laboratories that would like to participate in a standardized acquisition in n=10 rat to estimate the between-lab variability when other parameters are kept constant. The parameters need to be agreed upon within participating laboratories. 

Material transfer will be performed under tacit e-mail agreement, or if the data provider host institution requires it, a Material Transfer Agreement will be signed between the provider host institution and the Radboud University Medical Center, Nijmegen, The Netherlands (JG’s host institution). 

### 1.5 Deliverables
* OSF.io preregistered study
* Publicly available collective dataset on openneuro.org consisting of individual lab dataset
* Publicly available code to replicate the study on github.com
* Bioarxiv preprint
* Journal publication
* Symposium on rodent fMRI standards held at the Donders Institute (if funding is available)

### 1.6 Timeline
- December 2019: AH contacted potential collaborators
- March 2020: Initial study proposal is submitted to collaborators, effort to recruit additional labs
- April-September 2020: Data collection/transfer from collaborators
- June 2020: Study preregistration on OSF.io
- July-December 2020: Analysis as per preregistration specification
- January-April 2021: Manuscript preparation. 

### 1.7 Communication
Important milestones: Group e-mail (bcc). 
Preregistration preparation: A google doc template. Every collaborator is invited to comment/edit. Please [contact me](https://grandjeanlab.github.io/pages/contact.html) to get access.
Analysis update: This Github repository github.com/grandjealab/MultiRat. Every collaborator is invited to comment/raise an issue. 
Manuscript preparation: Google doc. 

### 1.8 Data storage
Short-term data storage on google drive/dropbox to transfer data. 
Mid-term data storage at Donders Institute high-performance computer in BIDS format
Long-term data storage at the openneuro.org repository in BIDS format, __publicly available when the first preprint is submitted__.

### 1.9 Data processing
Data will be preprocessed using [RABIES](https://github.com/CoBrALab/RABIES), a BIDS-based software based on the fMRIprep pipeline. Data will be co-registered into [Waxholm space](https://scalablebrainatlas.incf.org/rat/PLCJB14). Seed-based analysis, group-ICA with automatic component estimation, and dual-regression analysis will be performed using FSL 6.0.1. Statistical analysis will take place in R (ROI) or Randomize (voxel-wise). For a detailed procedure, see **section 4 Detailed analysis**.

### 1.10 Retraction
At any time and without justification, collaborators can decide to be removed from the study. Their dataset will be deleted and their results not used in the final results. This needs to be communicated to JG or AH. 

## 2 List of contacted/participating laboratory (as of 10/03/2020)
Group | Institute | Country
--- | --- | ---
Shih, Yen-Yu Ian | Chapel Hill | USA
Eike Budinger, Jürgen Goldschmidt | LIN, Magdeburg | Germany
Cornelius Faber | Münster | Germany
Valerio Zerbi | Zürich | Switzerland
Kai-Hsiang Chuang | Brisbane | Australia
Jason Lerch |  Oxford | UK
Alessandro Gozzi | Italian Institute of Technology | Italy
Diana Cash | King’s College London | UK
Caitlin Fowler, Jamie Near, M. Mallar Chakravarty | McGill | Canada
Marc Dhenain | cea |  France
Shella Keilholz | Georgia Tech | USA
Dijkhuizen, R.M. |  Utrecht University | Netherlands
Marleen Verhoye, vanDer Linden Lab | U of Antwerp | Belgium


## 3 Study preregistration draft 
as per OSF.io template
### 3.1 Design Plan
#### 3.1.1 Blinding
None. Analysis tools are automatic and require minimal user inputs. The user inputs will be as follows: visual inspection of QA output, manual drawing of noise mask for signal-to-noise analysis. 
#### 3.1.2 Study design
Cross-sectional, cross individual-dataset comparisons. In the test-retest arm of the study, a longitudinal design will be applied (Timepoint 1 vs Timepoint 2) 
#### 3.1.3 Randomization
None.
### 3.2 Sampling Plan
#### 3.2.1 Explanation of existing data
Rat fMRI dataset consisting of n=10 fMRI and matching anatomical scan. We will accept rat data from any strain, gender, age, weight. Acquired on any field strength, any coil with GE-EPI or SE-EPI, any parameter.
#### 3.2.2 Data collection procedures
Invitation to collaborate and to make a dataset available is made through email, announced at conferences (ESMI, ISMRM), announced on Twitter (@grandjeanlab). In the de novo collection arm of the study, participating centers will agree on a standardized protocol prior to data acquisition. 
### 3.3 Variables
#### 3.3.1 Manipulated variables
Lab-specific acquisition parameters. 
#### 3.3.2 Measured variables
**Rat**: Strain, Gender, Age, Weight, Vendor(?)
**Anesthesia**: maintenance anesthesia, dose, time post-induction, average breathing rate, average heart rate. 
**Acquisition system**: Vendor, Field strength, Coil setup(T/R, T+R), Coil design (Cryo, RT-phased array, RT-single loop).
**Acquisition sequence**: sequence type (GE-EPI, SE-EPI, other) TR, TE, FA, in-plane resolution, slice thickness, slice gap, number of volumes, bandwidth(?), FOV sat(?).
**During analysis**: Average cortical signal-to-noise ratio, average cortical temporal cortical signal-to-noise ratio, average framewise displacement, functional connectivity (z-stat) between S1-S1 and Cg-Rsp seeds, functional connectivity (z-stat) within S1 and Cg extracted from dual-regression analysis following group-independent component analysis with melodic automatic component estimation. 
   
### 3.4 Analysis Plan
#### 3.4.1 Statistical models
##### 3.4.1.1 Sensitivity: parameters associated with connectivity strength.
Functional connectivity (FC) between selected region of interest (S1 - S1, Cg - Rsp), or within ICA component (S1, Cg)  will be modeled into a linear model in statistical software R. Functional connectivity parameter will be modeled as a function of strain (factorial), gender (factorial), weight (ordered factorial), cortical signal-to-noise ratio (continuous), temporal cortical signal-to-noise (continuous, if not correlated to signal-to-noise ratio), mean framewise displacement. Significance will be assessed using an analysis of variance. P-value threshold will be set at p<0.05 without additional correction.  
##### 3.4.1.2 Specificity: parameter associated with connectivity specificity. 
Assuming S1 and Cg belong to distinct anti- (or minimally-) correlated networks, FC specificity will be determined the 4 quadrant system in Grandjean et al 2020. A 𝛸2 test will be used to determine which factors (field strength, coil design, anesthesia, strain, gender) have a skewed distribution of specific FC. 
##### 3.4.1.3 Seed-based analysis voxel-wise analysis across the collective dataset.
A one-sample t-test will be performed across the collective dataset for seeds in the S1 barrel field area, Cingulate area, Retrosplenial area, Insula area, motor area, dorsal hippocampus, caudate putamen, amygdala, thalamus. Non-parametric statistical test (Randomize) will be used to estimate one-sample t-test, with TFCE cluster correction. Because of the high anticipated degrees of freedom (n>200), we will use a p-value threshold of 0.0001. The map will be indicated as a thresholded z-statistics map overlaid on the template. 
##### 3.4.1.4 Seed-based analysis voxel-wise analysis across the individual datasets
The seed-based analysis from the individual datasets (each consisting of n=10) will be examined with a more lenient parametric one-sample t-test (fsl_glm), without cluster correction and p-value threshold 0.05. This is to ensure that no FC is rejected (low false negative), but at the expense of a higher false-positive rate. The analysis across individual datasets will be summarized in an overlap map denoting the percentage of datasets reaching significance for each voxel.  

## 4 Detailed analysis
### 4.1 Dataset preprocessing
All scans will be converted to NIFTI with original voxel size. Axis labels will be swapped so that NIFTI SI / AP / LR labels correspond to the right orientation. Scans will be organized using the [Brain Imaging Data Structure (BIDS) format](https://bids.neuroimaging.io/).

Meta-data (see **3.3.2 measured variables**) will be kept in a master table in a tab-separated format, with corresponding JSON file according to BIDS format.   

All scans will be preprocessed using [RABIES](https://github.com/CoBrALab/RABIES), a BIDS-based software based on the fMRIprep pipeline. Data will be co-registered into [Waxholm space](https://scalablebrainatlas.incf.org/rat/PLCJB14). Denoising will be performed using an adapted [ICA-AROMA automatic method](https://github.com/maartenmennes/ICA-AROMA) adapted for [the rodent and RABIES](https://github.com/Gab-D-G/conf_reg_pkg). Temporal filtering will be applied at 0.01-0.1 Hz for all scans (3dbandpass). Bandpass filter will be applied at 0.5mm for all scans (3dblurinmask).

RABIES outputs diagnostic files for registration accuracy and motion correction. The diagnostic output will be examined for each scan to ensure that all datasets are properly preprocessed. Diagnostic output samples will be presented in the manuscript as supplementary information. 

ICA-AROMA classification accuracy will be verified in 2 scans/dataset. Acceptable classification accuracy will be set at 95% sensitivity (true noise component detected as such), and 5% false positive (signal component detected as noise). Classification accuracy will be reported as a table as supplementary information. Examples of components classified as noise and signal will be represented in the supplementary information. 

**No other denoising strategies will be applied (e.g. global signal regression).** 

Average framewise displacement will be estimated from motion correction parameters. 
Average cortical signal-to-noise ratio and temporal signal-to-noise ratio will be determined using an ROI across the whole cortex, back-projected from the template, and hand-drawn ROIs in the corner of EPIs (making sure no ghost or artifacts are captured). 

Average framewise displacement and signal-to-noise ratios will be added as columns to the meta-data table, and plotted as a function of datasets in the main figure, and added in the dataset description supplementary table. 


### 4.2 rsfMRI dataset comparison. 
Seed-based analysis for the following seeds (S1 barrel field area, Cingulate area, Retrosplenial area, Insula area, motor area, dorsal hippocampus, caudate putamen, amygdala, thalamus. ) placed on the left-hemisphere will be performed using FSL 6.0.1 (*fslmeants, fsl_glm*), and output as z-statistic maps. Seeds will be defined as 0.9 mm3 spheres. See the statistics section for group comparisons (**3.4.1.3**). FC extracted from contralateral seeds (or from Rsp relative to Cg) will be plotted per dataset. A sensitivity/specificity analysis will be performed using seeds in the S1 barrel field (both hemisphere) and in the cingulate and retrosplenial cortex. The Retrosplenial seed will be specific to the cingulate (part of the same network), whereas the seed in the S1 will be unspecific to the cingulate cortex (part of a different network). We hypothesize that FC specificity is characterized by high FC between specific seeds and low FC between unspecific ROI. The FC relative to the cingulate extracted from the Retrosplenial and S1 seeds will be plotted on the x and y-axis. Quadrants will be applied to determine the percentage of scans with specific/unspecific or no FC (determined as low FC in both seeds).  

Independent component analysis (ICA) will be performed in melodic using automatic dimension estimation. If this is impossible to compute due to a too large sample, we will revert to d=30 components. Selected group ICA will be represented in a main figure, with the full list in the supplement. To estimate individual-level FC within these components, a [dual-regression approach](https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/DualRegression) will be applied (*fslmeants, fsl_glm*). FC parameters from selected components will be extracted using the seed defined above and plotted as a function of the individual datasets. A sub-analysis per strain, e.g. Winstar, (and gender) will be performed if sufficient individual datasets (n>4) have been collected. 

### 4.3  rsfMRI test-retest.
If a sufficient dataset (n=4 dataset) has been collected, we will examine the reliability of the parameters from the analysis in 4.2 within the same animals. FC parameters will be plotted during the test and retest phase. We will determine the variance between scan 1 and scan 2 across all animals. 

### 4.4  stimulus-evoked fMRI comparison
If sufficient dataset (n=4 dataset per stimulated area) is collected, we will examine the spatial specificity and response amplitude distribution of stimulus-evoked fMRI experiments. Briefly, beta parameter estimates for the stimulation will be obtained using a GLM approach (*fsl_glm*). One-sample t-test group spatial maps (*Randomize*, TFCE correction, p-value threshold = 0.05) will be compared for spatial delineation. Beta parameter estimate within both the ipsi- and contralateral cortical region being stimulated will be extracted (*fslmeants*) and plotted per dataset. The time-locked BOLD course will be extracted from the same ROI and plotted as the average per dataset. 

### 4.5  A de-novo acquisition with standardized parameters
To assess the repeatability of a standardized protocol and the human factor inducing variability, datasets will be acquired with a mutually agreed fMRI protocol. This will only be performed if there is at least a n=4 dataset. 
The tentative parameter list is as follows: 
**Rat**: Wistar, male, 300g
**Anesthesia**: Free-breathing, induction with isoflurane 5% in ¼ O2 to medical Air (2ml/min) for 4 min. Transfer to cradle (isoflurane 3%), positioning using ear bar and face mask. Reduction to isoflurane 1.5% in ¼ O2 to medical air (2ml/min), and insertion into bore. fMRI is to take place 30min following induction. (precise timing to be determined). 
**Acquisition**: Any field strength, any coil, GE-EPI sequence, TR = 1000ms, TE = 17ms, FA = 50, FOV = 25 x 25 mm, matrix size 100 x 100, number of slice 20, slice thickness 0.5mm, slice gap 0.1mm, repetitions = 1000, shimming with MAPSHIM following B0 maps, with square voxel into cerebrum.  
The same parameters as in **4.2** will be extracted and compared across dataset. 


