Title: Accurate DNA Methylation Predictor for C9orf72 Repeat Expansion Alleles in the Pathogenic Range 

Authors: Naren Ramesh, Alexandria Evans, Kevin Wojta, Zhongan Yang, Marco P Boks, René S. Kahn, Sterre C.M. de Boer, Sven J. van der Lee, Yolande A.L. Pijnenburg, Lianne M. Reus, Roel A. Ophoff

Overview

This repository contains code and data for analysis of C9orf72 DNA methylation (DNAm) profiles and pathogenic expansion status. The main focus is on:
	1.	Differential Methylation Analysis (DMA) – using CpG methylation data.
	2.	LASSO-based prediction modeling – using CpG methylation values to predict pathogenic repeat expansion status.

The repository contains only de-identified or non-patient-specific data.
Note: Raw methylation array data and clinical source files required for DMA are not included due to patient privacy concerns.

1. Differential Methylation Analysis.R
	•	Purpose: Performs differential methylation analysis to identify CpGs associated with C9orf72 repeat expansion status.
	  •	This script cannot be run with only the provided files because the original methylation array .Rdata objects and patient metadata are excluded.
	•	Requirements: Required R packages are listed at the top of the script (e.g., minfi, limma, DMRcate, ggplot2).
	•	Data: Requires normalized methylation array data (not provided in this repository).
	•	Output: Statistical results and plots highlighting differentially methylated CpGs.

3. LASSO Regression By Array Type.R
	•	Purpose: Performs LASSO regression to predict pathogenic C9orf72 repeat expansion status across EPICv2, EPICv1, Methyl450K, and Methyl27K compatible CpGs.
	•	Requirements: Required R packages are listed at the top of the script.
	•	Data: Requires initial training cohort CpG data (C9orf72 CpGs.csv) and annotations of included CpGs (C9orf72 CpG Annotations.csv).
	•	Output: LASSO models and accuracy of predictions across array technology.

4. C9orf72 CpGs.csv
	•	Purpose: Dataset for LASSO-based prediction modeling.
	•	Format:
  	•	Each row = one patient/sample.
  	•	exp column: Binary variable (1 = pathogenic expansion present, 0 = no expansion).
  	•	Remaining columns = M-values for CpG sites included in the model.
   
5. C9orf72 CpG Annotations.csv
	•	Purpose: Metadata describing the CpGs used in the LASSO prediction model.
	•	Fields include:
  	•	CpG identifier (e.g., cg#######)
  	•	Genomic location (chromosome, position)
  	•	Associated gene region
  	•	Illumina methylation array platforms on which the CpG is available.

Note: Quality Control and Filtering for raw methylation data was done using methodoology outlined here: 
https://www.bioconductor.org/packages/release/workflows/vignettes/methylationArrayAnalysis/inst/doc/methylationArrayAnalysis.html

Privacy Notice:
All included datasets have been stripped of identifiable patient information.
Source methylation array files and clinical data are excluded in compliance with patient privacy requirements.

⸻

Citation:
If you use this repository or its methods in your research, please cite appropriately. This study is available on preprint here:
https://www.biorxiv.org/content/10.1101/2025.03.20.643775v1

Correspondence:
Any questions or concerns can be forwarded to naren.ramesh@gmail.com
