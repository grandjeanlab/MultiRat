MultiRAT details
================
Joanes Grandjean

![rat art](../assets/img/rat_art.png)

### Collaborative model
Every laboratory is invited to participate with 1 or multiple rat fMRI dataset (see **Data specification**). Invitation to collaborate is made through email, announced at conferences (ESMI, ISMRM), announced on [Twitter](https://twitter.com/grandjeanlab). See the list of currently listed laboratories (**section 2**). 

The collaborative model is as follows: 
- Every contributing laboratory can nominate any reasonable number of junior and senior staff. Junior staff who have been involved in data collection. Senior staff who have been involved in project coordination and funding acquisition. 
- Joanes Grandjean ([JG](https://grandjeanlab.github.io/pages/contact.html)) will put the collective dataset together and perform the primary analysis. 
- Andy Hess (AH) will oversee and coordinate the project. 
- Gabriel Desrosiers-Gré?goire (GD) and M. Mallar Chakravarty (MMC) will provide the preprocessing software (RABIES)  
- Every collaborator is invited to further contribute to the preregistration proposal, the analysis, and manuscript preparation and will be recognized in their contributions accordingly. 

In the ensuing journal publication, author list will be as such: 
JG, GD, junior with an extra contribution, [junior collaborator in alphabetical order], [senior collaborator in alphabetical order], senior with an extra contribution, MMC, AH. 
 
In addition to the main study, collaborators are invited to propose spin-off studies relying on the dataset, which will be provided by JG. The spin-off studies should respect an embargo period on publication, as agreed together with JG+AH. E.g. a comparison of connectivity between mice and rats using graph theory (AH). 

### Data specification
Individual datasets each consist of an n=10 resting-state fMRI and corresponding anatomical scan obtained in rats, any strain, any gender, any age, any weight, any anesthesia or awake condition, acquired at any field strength, any vendor, any coil, with gradient-echo or spin-echo fMRI sequence. The anatomical scans are preferably T2 or T2* weighted. In addition to providing the scans in a nifti-convertible format, the dataset also consists of relevant meta-data. Specific data preprocessing requirements (fieldmap correction or top-up need to be discussed separately with JG). Similar to [Grandjean et al 2020](https://www.sciencedirect.com/science/article/pii/S1053811919308699), the individual dataset will be anonymized with respect to the laboratory where it was produced.

In addition, we would like to examine test-retest, that is datasets (n=10) where the same rats are imaged within 6 months with the same imaging protocol. We also seek to gather stimulus-evoked fMRI dataset (n=10) from sensory stimulation, to compare localization and effect amplitude. 

Finally, we seek laboratories that would like to participate in a standardized acquisition in n=10 rat to estimate the between-lab variability when other parameters are kept constant. The parameters need to be agreed upon within participating laboratories. 

Material transfer will be performed under tacit e-mail agreement, or if the data provider host institution requires it, a Material Transfer Agreement will be signed between the provider host institution and the Radboud University Medical Center, Nijmegen, The Netherlands (JG??s host institution). 

### Tentative timeline
- December 2019: AH contacted potential collaborators  
- March 2020: Initial study proposal is submitted to collaborators, effort to recruit additional labs  
- April-September 2020: Data collection/transfer from collaborators  
- June 2020: Study preregistration on OSF.io  
- July-December 2020: Analysis as per preregistration specification  
- January-April 2021: Manuscript preparation. 

### Communication
Important milestones: Group e-mail (bcc). 
Preregistration preparation: A google doc template. Every collaborator is invited to comment/edit. Please [contact me](https://grandjeanlab.github.io/pages/contact.html) to get access.
Analysis update: This Github repository github.com/grandjealab/MultiRat. Every collaborator is invited to comment/raise an issue. 
Manuscript preparation: Google doc. 

### Data storage
Short-term data storage on google drive/dropbox to transfer data. 
Mid-term data storage at Donders Institute high-performance computer in BIDS format
Long-term data storage at the openneuro.org repository in BIDS format, __publicly available when the first preprint is submitted__.

### Data processing
Data will be preprocessed using [RABIES](https://github.com/CoBrALab/RABIES), a BIDS-based software based on the fMRIprep pipeline. Data will be co-registered into [Waxholm space](https://scalablebrainatlas.incf.org/rat/PLCJB14). Seed-based analysis, group-ICA with automatic component estimation, and dual-regression analysis will be performed using FSL 6.0.1. Statistical analysis will take place in R (ROI) or Randomize (voxel-wise).

### Retraction
At any time and without justification, collaborators can decide to be removed from the study. Their dataset will be deleted and their results not used in the final results. This needs to be communicated to JG or AH. 

