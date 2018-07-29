Progress Report
================
Abdullah Farouk
2018-03-17

What has changed based on the final proposal?
=============================================

### Did your dataset change? If so, why?

The dataset that we will be using has not changed. We are still using the datasets GSE57767 and GSE44667, and a private preterm birth placental DNA methylation dataset from the Robinson Lab for validation.

### Have you decided to do a different analysis than what was mentioned in your proposal? If so, Why?

In the initial proposal of our methodology, we intended to use beta regression analysis to identify differentially methylated sites between experimental groups. The objective was to compare the results achieved with beta regression to the results obtained from standard linear regression. We have identified that there are no R packages readily available to conduct beta regression on methylation array data. The most relevant R package (“BiSeq”) for beta regression on methylation data is designed to handle bisulfite sequencing data. Although it may be possible to implement this package for our data, it is not without extensive understanding of its functions. Therefore, we have decided not to use beta regression, and instead to only use standard linear regression. After normalization, our average M values appear normally distributed, and thus satisfy one of the major assumptions for using standard linear regression, supporting this choice.

Additionally, in the proposal we planned to use linear grouping algorithms (LGA) as a novel method of clustering our data. Having learned more about this procedure, we have realized that the complexity of such a method is not compatible with the scope of this project due to our limited timeline.

The final change involves examining batch effects present in the data. Reference to batch in the respective GEO dataset publications gave us the impression that addressing the effect of batch would be possible. However, after exploring the metadata very carefully, we have recognized that proper batch notation was not provided in either GEO dataset. Over the multiple generations of Illumina methylation arrays, it has become clear to researchers that batch effects exist between microarray chips, as well as between rows of the microarray chips. Since our datasets are from 2013 and 2014, it is probable that the authors did not include batch information, as the potential confounding effect of chip and row had not been widely realized.

### Are there any changes in task assignments of group members?

We have made no changes in the group member task assignments.

What is the progress of the analyses?
=====================================

### Since your initial proposal, you should have decided more concretely on what methods to use for each step of your analyses, and employed some of those methods. Briefly and concisely explain your methodology and progress for the aims you have investigated so far. Which parts were modified and which parts remained the same?

**We decided to combine the first two questions of this section in the progress report rubric to provide a view of the full scope of our project. The following are steps that have already been completed.**

1.  **Pre-processing** - (2 GEO datasets pre-processed, individually) During preprocessing, we filtered out probes based on the following criteria:
    -   Probes with poor detection P values ( &gt; 0.05) from the raw signal intensity data sets.
    -   Probes targeting the sex chromosomes.
    -   Probes targeting the 65 control SNP sites used by [Illumina 2010](https://www.illumina.com/documents/products/technotes/technote_cpg_loci_identification.pdf) for sample identification.
    -   Probes in silico predicted to be nonspecific and capable of cross-hybridizing to more than one genomic location according to previous studies, or in close proximity to polymorphic loci [Price et al. 2013](https://epigeneticsandchromatin.biomedcentral.com/articles/10.1186/1756-8935-6-4).

    After filtering the dataset, we calculated beta values from the raw methylated/unmethylated bead signal intensities using “M/(U+M)”, where M and U refer to the methylated and unmethylated signal intensities, respectively. We calculated M values from our beta values, using log(beta/1-beta). These formulas were sourced from [Du et al. 2010](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-11-587). Standard probe metadata provided by Illumina (probe-specific) were attached to the filtered beta and M value probesets.

    The two pre-processed GEO datasets were merged, on the basis of probes common to both after probe filtering. The merged metadata can also be found in the “Data” folder of our repository, which includes the factors “Disease Class” and “Gestational Age” (the point in time of the pregnancy in which the placenta was sampled).

2.  **Normalization** - After considering a few normalization methods, we decided that quantile normalization was appropriate for our data. Merged datasets were normalized, for both our beta and M value datasets.

3.  **Data exploration** - We performed multiple checks to ensure that the sample numbers and probes numbers were correct and missing data was not present. Densities of beta and M values were plotted before and after normalization.

4.  **Principal Component Analysis (PCA)** - PCA was conducted on methylation probes to identify the variability explained by factors available in the metadata. The first 10 principal components capture 37% of the variability in our data. The first 20 PC’s capture 48% of the variance. Gestational age (term vs. preterm) and disease status (preeclampsia vs. healthy) appear to be associated with the first two principal components. We had been hoping to explore batch effects with PCA, however it appears that batch is not a factor available in either of our GEO datasets; a major limitation of our datasets.

5.  **Agglomerative Hierarchical Clustering (AHC)** - We first filtered for probes found within CpG islands (CGI) and then aggregated CpGs to the CGI level (by taking the mean M value), prior to clustering. The clustering method we applied was Ward’s method and our distance metric was “euclidean”. We also verified the clustering using multiscale bootstrap resampling (although this was computationally intensive).

6.  **Standard linear regression on all CGIs** - Standard linear regression was performed (using the “limma” package) to identify the top differentially methylated CGIs between normal and preeclampsia birth cases. We identified 87 CGIs that are differently methylated, given a p-value cutoff of 0.01. AHC was performed again on those 87 CGIs.

**The following steps are still in progress. The methodologies we are planning to implement are listed below.**

1.  **Supervised learning methods to predict preterm birth (KNN or SVM)** - For machine learning, more exploratory work needs to be done. We plan to test out methods such as K-means, KNN, SVM and random forest. We will select the method that provided us with the maximum area under the curve (AUC) for the corresponding receiver operating characteristic curve (ROC curve). Ultimately, we plan to use supervised learning methods to predict the disease status of samples based on methylation data. We will accomplish this by filtering the CGI data based on the results of Limma. This should improve classification accuracy.

2.  **Analysis of differentially methylated hits** - Gene annotation analysis will be performed for our top differentially methylated hits based on CGI (hypo/hyper methylated). We plan to find the location of the differentially methylated CGIs previously identified. We will focus on the CGIs that are located within promoters and gene bodies. We will then use various databases (e.g.: Gene Ontology Consortium, OMIM) and literature search results to determine the function of gene and potentially deduce the association between the differential methylation and the pathology observed. Further, we will rely on literature search to identify a few potential candidate genes that we find to be differentially methylated that could be interesting for further research. The criteria we plan to use in identify such targets include novelty, biological relevance, and if the genes are targetable.

### What R packages or other tools are you using for your analyses? You do not need to provide your scripts in your report.

In addition to other standard R packages, specific packages required by our project are: limma, GEOquery, tidyverse, pvclust, e1071, randomForest, pheatmap, and cluster.

### Provide the links to any markdown reports within your repo to refer to the relevant analysis. Provide references.

[Pre-processing and normalization script](https://github.com/STAT540-UBC/Repo_team_GenX/blob/master/Code/pre_processing_data.Rmd)

[Initial analysis](https://github.com/STAT540-UBC/Repo_team_GenX/blob/master/Code/initial_analysis.md)

Results
=======

What are your primary results? Were you able to answer your hypothesis? Did you have any positive results? If no, postulate a discussion as to why that may be. Provide plots and/or tables to present your results. - List some challenges that you encountered? How will you address them?
--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

We identified 87 differentially methylated CpG islands between normal and preeclampsia cases (p-value cutoff = 0.01). Additionally, our PCA results suggest that gestational age and preeclampsia status (case vs. control) are related to some of the variability we are observing in our data. Our hypothesis was that differential methylation would be detectable in comparing preeclampsia cases to healthy pregnancies. So far, our results indicate several methylation differences between the two groups. We have not yet analyzed the genes associated with differentially methylated CGIs, or the position of those CGIs.

Please see the above linked reports for our preliminary results and figures.

We have encountered three major challenges thus far.

1.  We are unable to account for batch effects present in our data. Unfortunately, despite reference to batch in the text of the GEO data source publications, no metadata related to batch was actually available for either of our GEO datasets. As such, we cannot currently address batch effects related to microarray chip, chip row, or experimental batch. We may need to investigate potential batch-effect correction techniques to address this concern.

2.  As discussed above regarding changes in methodology, we have decided to drop beta regression as an analytical method, as the only known R package available for doing this type of regression is designed for DNA bisulfite sequencing data. Other package options (like “betareg”) do not necessarily account for multiple testing or poor variance estimates.

3.  Finally, we had intended to use the Linear Grouping Algorithm (LGA) on the methylation data, but have found it will be difficult to implement considering the scope of this project, and have thus focused our efforts in other areas.

References
==========

Du P, Zhang X, Huang CC, Jafari N, Kibbe WA, Hou L, Lin SM. 2010. **Comparison of Beta-value and M-value methods for quantifying methylation levels by microarray analysis**. *BMC Bioinformatics*. 11:587.

Illumina. 2010. **Technical Note: Epigenetics**. Illumina, Inc (US). Pub. No. 270-2007-006. Available from: <https://www.illumina.com/documents/products/technotes/technote_cpg_loci_identification.pdf>

Price EM, Cotton AM, Lam LL, Farre P, Emberly E, Brown CJ, Robinson WP, Kobor MS. 2013. **Additional annotation enhances potential for biologically-relevant analysis of the Illumina Infinium HumanMethylation450 BeadChip array**. *Epigenetics & Chromatin*. 6:4.
