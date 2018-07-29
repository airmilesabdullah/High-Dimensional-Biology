Technical\_Report
================
Nikolas Krstic
April 4, 2018

Differential methylation of preeclamptic placental tissue using the Illumina HumanMethylation450K array
=======================================================================================================

Introduction
------------

Preeclampsia (PE) is a hypertensive disorder that affects 2-8% of all pregnancies \[1\]. The symptoms range from maternal hypertension and proteinuria to more severe complications when the disease progresses to eclampsia, such as renal/liver dysfunction, and seizure \[2\]. Multiple studies have demonstrated that PE is associated with a particular placental DNA methylation pattern \[3, 4\]. It is important to conduct further analysis to confirm these findings, and to potentially train a clinical tool capable of predicting PE disease status by placnetla methylation signature. Thus, we performed a systematic analysis with 2 publicly available datasets (GSE57767 \[5\] and GSE44667 \[6\]) containing the Illumina 450K methylation data from placental tissue in a total of 53 PE cases and 32 controls. The purposes of our study was to 1) use methylation pattern to predict patients’ PE status; and 2) identify targets and pathways associated with PE for further study.

Methods
-------

First, we gathered methylation data from two publicly available datasets from the Gene Expression Omnibus (GEO) database (GEO asscession numbers GSE44667 \[5\] and GSE57767 \[6\]). We first handled the raw signal intensities in order to preprocess the two datasets identically. For further details of this step, please refer to [pre\_processing\_data.Rmd](https://github.com/STAT540-UBC/Repo_team_GenX/blob/master/Code/pre_processing_data.Rmd). From this pre-processing step, we generated the normalized beta and M values for all 85 samples.

Afterward, we performed our exploratory analysis on the data,in addition to modelling with Limma. Ultimately, we decided to aggregate the CpG site probe data (approx. 480 000 probes) to CpG-Islands (CGI). Several analytical approaches we applied include PCA (on both the CpG sites and CGIs), agglomerative hierarchical clustering on the CGIs, modelling with limma, and a feature selection step to prepare for our supervised learning approaches. The details for each of these steps can be found in [exploratory\_and\_limma\_analyses.md](https://github.com/STAT540-UBC/Repo_team_GenX/blob/master/Code/exploratory_and_limma_analyses.md). To assess whether the M values of the CGIs are normally distributed, we also performed two tests to check for normality: (1) Anderson-Darling and (2) Shapiro-Wilks tests. These tests can be found within [m\_value\_normality\_check.md](https://github.com/STAT540-UBC/Repo_team_GenX/blob/master/Code/m_value_normality_check.md). Lastly, after obtaining the top 87 differentially methylated CGIs from our above limma analyses, we performed a gene enrichment analysis to identify gene functions associated with these CGIs. This analysis can be found within [GO\_analysis.md](https://github.com/STAT540-UBC/Repo_team_GenX/blob/master/Code/GO_analysis.md)

For the supervised learning component of our analysis, our objective was to create a model that could predict the disease class of a sample ("normal" or "preeclampsia"). To do so, we split our data into a training set and a heldout test set (80% and 20% of the samples). This was done previously within [exploratory\_and\_limma\_analyses.md](https://github.com/STAT540-UBC/Repo_team_GenX/blob/master/Code/exploratory_and_limma_analyses.md) to conduct feature selection (so that the test set would not influence our choices). We compared several different approaches using nested cross-validation and then performed parameter tuning using repeated cross-validation, all on the **training set** only. The hold-out test set finally came into play when we predicted the disease classes, after having trained our model using the full training set as a final step. The details of these methods can be found within [Machine\_Learning\_Methods.md](https://github.com/STAT540-UBC/Repo_team_GenX/blob/master/Code/Machine_Learning_Methods.md). At one point in this analysis, we intended to use a separate dataset (GEO accession number GSE73375) as our test set, and use our original combined dataset as the training set. However, after implementing the exact same procedure above, the class probabilities of the test set were all the same. The creation of this test dataset can be found within [test\_data\_set\_creation.md](https://github.com/STAT540-UBC/Repo_team_GenX/blob/master/Code/test_data_set_creation.md).

Results
-------

We used Limma to identify the top differentially methylated CGIs. To do so we fit each CGI with a linear regression model as follows:

*M* *V**a**l**u**e*<sub>CGI</sub> = *D**i**s**e**a**s**e* *C**l**a**s**s*<sub>CGI</sub>  +  *G**e**s**t**a**t**i**o**n**a**l* *A**g**e*<sub>CGI</sub>

To correct for multiple hypothesis testing we set an FDR of 0.01. We found the coefficients for Disease class for 87 CGIs to be statistically different from 0 at this threshold.

Since traditional linear regression models assume normally distributed errors, we decided to test this assumption by checking the normality of our response variable; the CGI M values. We used the Anderson-Darling and Shapiro Wilks tests to do so by performing the test on the rows of our M value object (which correspond to CGIs). After using a bonferroni correction we found insufficient evidence to reject the hypothesis that the rows DO NOT follow a normal distribution for approximately 96% of them (CGIs).

We then performed hierarchical clustering on the samples using the top 87 differentially methylated (DM) CGIs (using the ward.D2 clustering method and euclidean distance metric). We do so to determine if we need to include additional sub-classes of preeclempsia or if retaining two would be sufficient. For Disease Class, the clusters seem to be relatively well-defined. This was verified by multiscale bootstrap resampling, which resulted in two stable clusters. Thus, we decided to retain two classes of Disease Class.

Next we performed hierarchical clustering on the samples using the top 87 differentially methylated CGIs (using the ward.D2 clustering method and euclidean distance metric). We did so to determine if we needed to include additional sub-classes of preeclempsia. For Disease Class, the clusters seem to be relatively well-defined. This was verified by multiscale bootstrap resampling, which resulted in two stable clusters. Thus, we decided to work with only two classes of Disease Class.

After that we sought to test if CGIs could predict the disease class of a placental cell sample. The first step was performing feature selection (ie pick a subset of 19315 CGIs in our dataset). We initially split our data (randomly) into a training and test set. We then carried out 5 fold cross validation to pick CGIs with the greatest explanatory power. This was done by using limma, with Disease class as our response and the 19315 CGIs as our covariates. Among each of the splits, the top 50 CGIs were chosen. This led to a selection of 162 unique CGIs as our "informative" features which were then used in our supervised learning methods.

The tuned SVM model (supervised learning method with the lowest misclassification error rate of the three tested) performed very well on our holdout set (17 samples the model was never trained on) with an AUC of 0.971

Additionally we employed reference free deconvolution to account for batch and cell heterogeneity (characteristic of placental tissue). We did so by using the first three PCs as covariates within limma to attempt to account for these effects. We obtained 6 differentially methylated (DM) CGIs (with an FDR of 0.01). Five of these CGIs are also found within the original 87 CGIs obtained without using the top 3 PCs as covariates. We also observed an association between Gestational Age and PC2 and a slight association between Disease Class and PC1 (based on plots we made). This could imply that some of the original 87 DM CGIs may have been false discoveries.

Our final analytical step in this project was the gene ontology analysis of the 87 differentially methylated CpG islands. It identified two processes as involving the most differentially methylated CGIs from our list. They are "biological regulation", which could relate to improper regulation of many processes related to the pathology of preeclampsia (PE) and "biosynthetic process"

Discussion
----------

In our initial linear modelling analysis (using limma), with disease status (preeclampsia/normotensive) and gestational age (term/preterm) as covariates, we identified 87 differentially methylated CpG islands at a statistical threshold of FDR &lt; 0.01.

One of the limitations of our investigation was our inability to use take batch effects from the two datasets into account, such as the effects that would arise from samples being run on different Illumina arrays, or at different times experimentally.Reference to batch in the respective GEO dataset publications gave us the impression that addressing the effect of batch would be possible. However, after exploring the metadata very carefully, it became clear that proper batch notation was not provided in either GEO dataset. Over the multiple generations of Illumina methylation arrays, researchers have recognized that batch effects exist between microarray chips, and even between sample rows of the microarray chips. Since our datasets are from 2013 and 2014, it is possible that the authors did not include batch information as the potential confounding effect of chip and row had not been widely realized.

Due to our inability to directly address batch effects, we applied a method of reference-free deconvolution by including the top 3 principal components as covariates in our linear model. When these were accounted for, only 6 sites were identified to be differentially methylated; reassuringly, 5/6 were also present in our list of 87 differentially methylated CGIs generated prior to reference-free deconvolution. These 5 CpG islands were associated with the following genes: CYP19A1, KRT19, GPR37L1, RNF217-AS1, and ZNF814. In previous studies, CGIs in genes from the protein families of GPR, RNF and ZNF have been identified as being differentially methylated in preeclamptic placentae, an interesting finding that lends some confidence to our data \[4\]. The disadvantage associated with this step is that we risk overcorrectio,n as PC2 does have a moderate association with the PE status.

Finally, we had planned to apply our trained supervised learning predictive model to a test dataset, GSE73375. With supervised machine learning, after performing nested cross-validation (CV) by 3 methods (KNN, SVM, and random forest) to classify disease status in our 68-sample training set, SVM was found to have the lowest error rate. The tuned SVM model was tested and found to perform well on our 17-sample test data with an AUC of 0.971. However, upon trying to apply this approach to a new GEO test dataset, we realized that it would be essential to pre-process this dataset in the same way that we had pre-processed the other two - this would include normalizing all 3 datasets together. Given the scope of our project we did not come up with a suitable solution to this problem, however have discussed the roadblock with both Emma and Paul and have devised appropriate ways of dealing with this that may be attempted in the future. At the present time, the generalizability of our linear model has yet to be tested.

Conclusion
----------

Overall, this project consisted of a differential methylation analysis of the placental DNA of preeclamptic and normotensive women. We made some interesting findings, and our results indicate that DNA methylation is indeed altered in this pathology. Future investigations into these datasets and other ways of dealing with batch effects (such as subtraction of potentially relevant principal components prior to linear modelling) may be conducted to continue our analyses. Additionally, it would be biologically interesting to both test our trained predictive model, as well as to more thoroughly investigate the genic locations of differentially methylated CGIs to potentially extend the findings of our study to possible inferences regarding effects on gene function and transcription alterations associated with preeclampsia. Ultimately, a clinical tool with the ability to predict preeclampsia based on DNA methylation analysis from a placental sample collected by chorionic villi sampling would be a very valuable application of our research, one which our initial results suggest might be possible.

References
----------

1.  Duley L. 2009. The Global impact of pre-eclampsia and eclampsia. Semin Perinatol. 33:130–137.

2.  \[ACOG\] The American College of Obstetricians and Gynecologists’ Task Force on Hypertension in Pregnancy. 2013. Hypertension in Pregnancy. Obstet Gynecol. 122:1122-31.

3.  Yeung KR, Chiu CL, Pidsley R, Makris A, Hennessy A, Lind JM. 2016. DNA methylation profiles in preeclampsia and healthy control placentas. Am J Physiol Heart Circ Physiol. 310:H1295-303.

4.  Anderson CM, Ralph JL, Wright ML, Linggi B, Ohm JE. 2014. DNA methylation as a biomarker for preeclampsia. Biol Res Nurs. 16:409-20.

5.  Anton L, Brown AG, Bartolomei MS, Elovitz MA. 2014. Differential methylation of genes associated with cell adhesion in preeclamptic placentas. PLoS One. 9:e100148.

6.  Blair JD, Yuen RK, Lim BK, McFadden DE, von Dadelszon P, Robinson WP. 2013. Widespread DNA hypomethylation at gene enhancer regions in placentas associated with early-onset pre-eclampsia. Mol Hum Reprod. 19:697-708.
