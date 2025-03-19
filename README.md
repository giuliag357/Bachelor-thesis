# Bachelor-thesis

This thesis explores using machine learning for tumor classification (prostate and breast cancer) from multi-omics data (mRNA-seq, DNA methylation, and CNA). While Fatima and Rueda’s iSOM-GSN method (DOI: 10.1093/bioinformatics/btaa500) used 2D transformations and CNNs, this work shows that a simpler random forest approach, combined with thorough data cleaning, feature selection (MutSigCV), and data augmentation (SMOTE), can achieve similarly impressive accuracy (95–99%).
The data used to run this project are from TCGA Prostate Adenocarcinoma (https://www.cbioportal.org/study/summary?id=prad_tcga) and TCGA Breast Invasive Carcinoma (https://www.cbioportal.org/study/summary?id=brca_tcga_pub2015).
