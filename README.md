# Final Year Project

This research focuses on Alzheimer's disease, caused by progressive neuron dysfunction. Using the GSE118553 dataset from NCBI GEO DATA, which includes brain tissue data from 27 control, 33 AsymAD, and 52 Alzheimer subjects, the study analyzed affected (entorhinal, temporal, frontal cortex) and spared (cerebellum) tissues. Employing Weighted Gene Co-expression Network Analysis (WGCNA), hub genes were identified through pairwise gene correlation. Machine learning algorithms validated the results, ensuring accuracy. The study uncovered novel genes involved in Alzheimer's pathogenesis, offering prospects for biomarkers and therapeutic targets, and showcasing machine learning's role in bioinformatics for treating neurodegenerative disorders.

## Objective
-Neurodegenerative diseases involve progressive neuron degeneration in the brain.
-Symptoms include memory loss, cognitive decline, and behavioral abnormalities.
-Alzheimer's and Parkinson's are the most common neurodegenerative disorders.
-Increasing prevalence due to aging population, posing significant healthcare challenges.
-Crucial to understand mechanisms like protein misfolding and oxidative stress.
-Loss of neuronal integrity and function leads to neurodegeneration.
-Molecular biology, bioinformatics, and computational biology can unravel disease complexities.
-Integration of large-scale data to identify genetic networks and signaling pathways.
-Key genes identification is promising for early diagnosis and disease-modifying -therapies.
-Computational approaches can enhance quality of life for affected individuals.

## Research Problem Statement
-Neurodegenerative diseases, including Alzheimer's, have a progressive and devastating impact on individuals and society.
-Understanding the molecular drivers is crucial for developing effective treatments.
-The study uses a computational approach to analyze genomics data and identify key genes related to Alzheimer's.
-The GSE118553 dataset includes brain samples from control subjects, individuals with early signs of Alzheimer's, and diagnosed patients.
-The research aims to extract insights into the molecular basis of Alzheimer's by examining affected brain regions and tissues.
-Goal: Identify genes involved in Alzheimer's pathogenesis using Weighted Gene Co-expression Network Analysis (WGCNA).
-Validate findings with machine learning algorithms for reliability.
-Aim to contribute to early diagnosis and development of therapeutic targets for Alzheimer's disease.

**Link To Dataset**
https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE118553

## Methodology

Data Preprocessing-R:
Preprocess raw datasets using R programming language.
Normalize data using quantile normalization to eliminate technical variability.
Transform and align column names across data frames for accurate merging.
Merge datasets to consolidate gene information and map to specific identifiers.
Join brain part information with normalized data to visualize gene interaction using Cytoscape.

Data Preprocessing-Python:
Import necessary libraries including sys, numpy, and pandas.
Parse command line arguments for input and output file paths.
Read matrix and family files, extract gene symbol and ID information, and remove duplicate probes.
Output trait and expression data to CSV files.

Weighted Gene Coexpression Analysis (WGCNA):
Load required R libraries for WGCNA and configure environment.
Prepare clean expression data for network construction, removing non-numeric data and transposing the data.
Detect and inspect outliers using dendrogram.
Visualize clusters and trait data using hierarchical clustering and color-coded heatmap.
Create network, select soft threshold power, and construct adjacency and Topological Overlay Matrix (TOM).
Group genes into modules using hierarchical clustering and dynamic tree cut.
Merge closely related modules and associate modules with phenotype traits.
Analyze module membership and gene significance, and generate CSV files for further analysis.

Machine Learning:
We utilized three algorithms: Random Forest, Gradient Boost, and XGBoost, all following a similar workflow for data preparation, model training, and prediction.
Importing Libraries: Import pandas, Scikit Learn, and classifiers like GradientBoostingClassifier, xgboost, and RandomForestClassifier.
Data Preparation: Separate features (X), encode target variable, and split data into training and testing sets.
Data Imputation: Impute missing values using the median.
Model Training: Initialize classifiers and train on imputed data.
Model Evaluation: Make predictions, test accuracy, and obtain feature importance.
Results: Feature importance to identify important hub genes.

**Requirements**
The requirement of this project involves high-end machines with GPU due to the requirements of processing large amounts of data. Along with the machines some other requirements are mentioned below. 

-Software:
R-studio version 1.4.1106
Visual Studio Code
Python version 3.9.x
Jupyter Notebook/Google Colaboratory

-Hardware:
RAM 64 GB
NVIDIA GeForce RTX 3090
Driver version:	31.0.15.3713
DirectX version:	12 (FL 12.1)
GPU Memory	55.9 GB

## Usage
-Environment Setup:
Ensure you have a high-end machine with a GPU, preferably NVIDIA GeForce RTX 3090, along with the required drivers installed (version 31.0.15.3713).
Install R-studio version 1.4.1106, Visual Studio Code, Python version 3.9.x, and either Jupyter Notebook or Google Colaboratory.
Ensure your machine has at least 64 GB of RAM.

-Dataset Preparation:
Access the dataset from NCBI GEO Database GSE118553.

-Data Preprocessing (R):
Install the following R packages: Affy, Tidyverse, Openxlsx, Dplyr, PreprocessCore.

-Data Preprocessing (Python):
Install the following Python libraries: Numpy, Pandas, Sys.

-WGCNA Setup:
Install the following R packages: BiocManager, GO.db, WGCNA, RColorBrewer, dynamicTreeCut.

-Machine Learning Setup:
Install the following Python libraries: Pandas, Scikit-learn, XGBoost, RandomForestClassifier.

## References
[1] Y. Lagisetty et al., “Identification of risk genes for Alzheimer’s disease by gene embedding,” Cell Genomics, vol. 2, no. 9, p. 100162, Sep. 2022, doi: 10.1016/j.xgen.2022.100162.

[2] J. J. Palop, J. Chin, and L. Mucke, “A network dysfunction perspective on neurodegenerative diseases,” Nature, vol. 443, no. 7113, pp. 768–773, Oct. 2006, doi: 10.1038/nature05289.

[3] G. G. Kovacs, “Concepts and classification of neurodegenerative diseases,” in Handbook of Clinical Neurology, 2018, pp. 301–307. doi: 10.1016/b978-0-12-802395-2.00021-3.

[4] M. A. Myszczynska et al., “Applications of machine learning to diagnosis and treatment of neurodegenerative diseases,” Nature Reviews Neurology, vol. 16, no. 8, pp. 440–456, Jul. 2020, doi: 10.1038/s41582-020-0377-8.

[5] W. Mandemakers, V. A. Morais, and B. De Strooper, “A cell biological perspective on mitochondrial dysfunction in Parkinson disease and other neurodegenerative diseases,” Journal of Cell Science, vol. 120, no. 10, pp. 1707–1716, May 2007, doi: 10.1242/jcs.03443.

[6] S. B. Prusiner, “Neurodegenerative diseases and prions,” The New England Journal of Medicine, vol. 344, no. 20, pp. 1516–1526, May 2001, doi: 10.1056/nejm200105173442006.

[7] D. P. Purohit, D. P. Perl, V. Haroutunian, P. Powchik, M. Davidson, and K. L. Davis, “Alzheimer disease and related neurodegenerative diseases in elderly patients with schizophrenia,” Archives of General Psychiatry, vol. 55, no. 3, p. 205, Mar. 1998, doi: 10.1001/archpsyc.55.3.205.

[8] J.-C. Lambert et al., “Meta-analysis of 74,046 individuals identifies 11 new susceptibility loci for Alzheimer’s disease,” Nature Genetics, vol. 45, no. 12, pp. 1452–1458, Oct. 2013, doi: 10.1038/ng.2802.

[9] G. Jun et al., “A novel Alzheimer disease locus located near the gene encoding tau protein,” Molecular Psychiatry, vol. 21, no. 1, pp. 108–117, Mar. 2015, doi: 10.1038/mp.2015.23.

[10] C. M. Karch, C. Cruchaga, and A. M. Goate, “Alzheimer’s Disease Genetics: From the bench to the clinic,” Neuron, vol. 83, no. 1, pp. 11–26, Jul. 2014, doi: 10.1016/j.neuron.2014.05.041.

[11] J. Hardy and G. A. Higgins, “Alzheimer’s Disease: The Amyloid Cascade hypothesis,” Science, vol. 256, no. 5054, pp. 184–185, Apr. 1992, doi: 10.1126/science.1566067.

[12] V. S. Marde et al., “Neurodegenerative disorders associated with genes of mitochondria,” Future Journal of Pharmaceutical Sciences, vol. 7, no. 1, Mar. 2021, doi: 10.1186/s43094-021-00215-5.

[13] V. Escott-Price et al., “Common polygenic variation enhances risk prediction for Alzheimer’s disease,” Brain, vol. 138, no. 12, pp. 3673–3684, Oct. 2015, doi: 10.1093/brain/awv268.

[14] Guan et al. "Weighted gene coexpression network analysis and machine learning reveal oncogenome associated microbiome plays an important role in tumor immunity and prognosis in pan-cancer." Journal of Translational Medicine 21.537 (2023): 1-21.

[15] Hu, Y., Yu, X., Zhou, Y., Yin, X., Hu, J., Lu, X., and Hu, C. (2020). Identification of novel genes related to Alzheimer's disease pathology using genome-wide co-expression network analysis and supervised machine learning. Frontiers in Aging Neuroscience, 12, 605961. doi: 10.3389/fnagi.2020.605961

[16] S. Balne and A. Elumalai, "Machine learning and deep learning algorithms used for the diagnosis of Alzheimer's" in Materials Today: Proceedings, vol. 47, no. 15, pp. 5151-5156, 2021, doi: 10.1016/j.matpr.2021.05.499.

[17] P. Kishore, Ch. Usha Kumari, M.N.V.S.S. Kumar, and T. Pavani, "Detection and analysis of Alzheimer’s disease using various machine learning algorithms," in Materials Today: Proceedings, vol. 45, no. 2, pp. 1502-1508, 2021, doi: 10.1016/j.matpr.2020.07.645.

[18] A. Amrutesh, G. B. C G, A. A, A. R. KP and G. S, "Alzheimer's Disease Prediction using Machine Learning and Transfer Learning Models," 2022 6th International Conference on Computation System and Information Technology for Sustainable Solutions (CSITSS), Bangalore, India, 2022, pp. 1-6, doi: 10.1109/CSITSS57437.2022.10026365.

[19] Diogo, V.S., Ferreira, H.A., Prata, D. et al. Early diagnosis of Alzheimer’s disease using machine learning: a multi-diagnostic, generalizable approach. Alz Res Therapy 14, 107 (2022). https://doi.org/10.1186/s13195-022-01047-y

[20] S. Harika, T. Yamini, T. Nagasaikamesh, S. H. Basha, S. Santosh Kumar, and Mrs. S. Sri DurgaKameswari, "Alzheimer's Disease Detection Using Different Machine Learning Algorithms," *International Journal of Research in Advent Technology (IJRASET)*, DOI: 10.22214/ijraset.2022.46937.

[21] H. Sun, J. Yang, X. Li, Y. Lyu, Z. Xu, H. He, X. Tong, T. Ji, S. Ding, C. Zhou, P. Han, J. Zheng, "Identification of feature genes and pathways for Alzheimer's disease via WGCNA and LASSO regression," Frontiers in Computational Neuroscience, vol. 16, p. 1001546, Sep. 2022, doi: 10.3389/fncom.2022.1001546.

[22] Langfelder, P., Horvath, S. WGCNA: an R package for weighted correlation network analysis. BMC Bioinformatics 9, 559 (2008), https://doi.org/10.1186/1471-2105-9-559

[23] Hong Yue, Bo Yang, Fang Yang, Xiao-Li Hu, Fan-Bin Kong, "Co-expression network-based analysis of hippocampal expression data associated with Alzheimer's disease using a novel algorithm," *Experimental and Therapeutic Medicine*, vol. 11, no. 5, pp. 1707-1715, 2016.

[24]  B. M. Tijms, A. M. Wink, W. de Haan, W. M. van der Flier, C. J. Stam, P. Scheltens, and F. Barkhof, "Alzheimer's disease: connecting findings from graph theoretical studies of brain networks," Neurobiology of Aging, vol. 34, no. 8, pp. 2023-2036, Aug. 2013. https://doi.org/10.1016/j.neurobiolaging.2013.02.020

[25] Brier, M. R., Thomas, J. B., Fagan, A. M., Hassenstab, J., Holtzman, D. M., Benzinger, T. L., Morris, J. C., and Ances, B. M. (2014). Functional connectivity and graph theory in preclinical Alzheimer's disease. *Neurobiology of Aging*, 35(4), 757-768. [https://doi.org/10.1016/j.neurobiolaging.2013.10.081](https://doi.org/10.1016/j.neurobiolaging.2013.10.081)

[26]Q. Jing, H. Zhang, X. Sun, Y. Xu, S. Cao, Y. Fang, X. Zhao, and C. Li, "A Comprehensive Analysis 
Identified Hub Genes and Associated Drugs in Alzheimer’s Disease," vol. 2021, Article ID 8893553, 
2021, https://doi.org/10.1155/2021/8893553.

[27]C. Zou, L. Su, M. Pan, L. Chen, H. Li, C. Zou, J. Xie, X. Huang, and M. Lu, "Exploration of novel 
biomarkers in Alzheimer’s disease based on four diagnostic models," Front. Aging Neurosci., vol. 15, p. 
1079433, Feb. 16, 2023. [Online]. Available: https://doi.org/10.3389/fnagi.2023.1079433

[28]Zhang, Y., Kiryu, H. Identification of oxidative stress-related genes differentially expressed in 
Alzheimer’s disease and construction of a hub gene-based diagnostic model. Sci Rep 13, 6817 (2023). 
https://doi.org/10.1038/s41598-023-34021-1

[29]A. Sharma and P. Dey, "A machine learning approach to unmask novel gene signatures and 
prediction of Alzheimer's disease within different brain regions," Genomics, vol. 113, no. 4, pp. 
1778-1789, 2021, ISSN 0888-7543, https://doi.org/10.1016/j.ygeno.2021.04.028.

[30]Y. Chen, Z. Li, X. Ge, H. Lv, and Z. Geng, "Identification of novel hub genes for Alzheimer’s disease associated with the hippocampus using WGCNA and differential gene analysis," Frontiers in Neuroscience, vol. 18, no. 1359631, Mar. 2024.

[31]H. Sun, J. Yang, X. Li, Y. Lyu, Z. Xu, H. He, X. Tong, T. Ji, S. Ding, C. Zhou, and P. Han, "Identification of feature genes and pathways for Alzheimer's disease via WGCNA and LASSO regression," Frontiers in Computational Neuroscience, vol. 16, p. 1001546, 2022, https://doi.org/10.3389/fncom.2022.1001546

[32]J. W. Liang, Z. Y. Fang, Y. Huang, Z. Y. Liuyang, X. L. Zhang, J. L. Wang, H. Wei, J. Z. Wang, X. C. Wang, J. Zeng, and R. Liu, "Application of weighted gene co-expression network analysis to explore the key genes in Alzheimer’s disease," Journal of Alzheimer's Disease, vol. 65, no. 4, pp. 1353-1364, 2018, DOI: 10.3233/JAD-180400.

[33]H. Wang, X. Han, and S. Gao, "Identification of potential biomarkers for pathogenesis of Alzheimer’s disease," Hereditas, vol. 158, pp. 1-8, 2021.

[34]J. Li, Y. Zhang, T. Lu, R. Liang, Z. Wu, M. Liu, L. Qin, H. Chen, X. Yan, S. Deng, and J. Zheng, "Identification of diagnostic genes for both Alzheimer’s disease and Metabolic syndrome by the machine learning algorithm," Frontiers in Immunology, vol. 13, p. 1037318, 2022, https://doi.org/10.3389/fimmu.2022.1037318

[35]J. Ren, B. Zhang, D. Wei, and Z. Zhang, "Identification of methylated gene biomarkers in patients with Alzheimer’s disease based on machine learning," BioMed Research International, 2020, https://doi.org/10.1155/2020/8348147

[36]B. Monk, A. Rajkovic, S. Petrus, A. Rajkovic, T. Gaasterland, and R. Malinow, "A machine learning method to identify genetic variants potentially associated with Alzheimer’s disease," Frontiers in Genetics, vol. 12, p. 647436, 2021, https://doi.org/10.3389/fgene.2021.647436

[37]H. Chen, Y. He, J. Ji, and Y. Shi, "A machine learning method for identifying critical interactions between gene pairs in Alzheimer's disease prediction," Frontiers in Neurology, vol. 10, p. 490162, 2019, https://doi.org/10.3389/fneur.2019.01162

[38]H. Alamro, M. A. Thafar, S. Albaradei, T. Gojobori, M. Essack, and X. Gao, "Exploiting machine learning models to identify novel Alzheimer’s disease biomarkers and potential targets," Scientific Reports, vol. 13, no. 1, p. 4979, 2023.

[39]Q. Jing, H. Zhang, X. Sun, Y. Xu, S. Cao, Y. Fang, X. Zhao, and C. Li, "A comprehensive analysis identified hub genes and associated drugs in Alzheimer’s disease," BioMed Research International, 2021, https://doi.org/10.1155/2021/8893553

[40]A. Sharma and P. Dey, "A machine learning approach to unmask novel gene signatures and prediction of Alzheimer's disease within different brain regions," Genomics, vol. 113, no. 4, pp. 1778-1789, 2021, https://doi.org/10.1016/j.ygeno.2021.04.028

