Proposal of statistical analyses to identify differentially methylated DNA sites associated with pre-term births
================
Nikolas Krstic, Amy Inkster, Abdullah Farouk, Yue Yang Shen
February 15, 2018

Motivation and background work
------------------------------

In 2013, 7.8% of all babies born in Canada were premature (Statistics Canada 2016). Investigating methylation differences between term (delivery &gt; 37 weeks) and preterm placentas may give us insight into the complex pathology of preterm birth. Studying DNA methylation in the placenta has proven to be a valid approach. Previous studies have shown associations between preterm births and both hyper- and hypo- methylated loci in placental DNA (Yeung et al. 2016). In this project, we will analyze placental methylation data collected from women affected by preeclampsia, one of the most common causes of preterm birth. Our hypothesis is that there will be differential methylation between placental tissue of preeclamptic births and normal pregnancies.

Division of labour
------------------

<table style="width:100%;">
<colgroup>
<col width="18%" />
<col width="33%" />
<col width="48%" />
</colgroup>
<thead>
<tr class="header">
<th align="left">Group Member</th>
<th align="left">Degree/Affiliations</th>
<th align="center">Projected Contributions/Tasks</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td align="left">Yue Yang Shen <br/><br/> BMLSc. Medical Laboratory Science</td>
<td align="left">Graduate student in Medical Genetics <br/><br/> Lab affiliation: Dr. Xiaoyan Jiang, at British Columbia Cancer Research Center (BCCRC).</td>
<td align="center">I am involved in the overall planning of the study. I will contribute to the statistical analysis as much as I can, and I will also perform literature search to bridge explain our findings. I hope to elucidate the biological implication of our data analysis and I will participate in the write-up of the final project.</td>
</tr>
<tr class="even">
<td align="left">Amy Inkster <br/><br/> BSc. Chemistry (Minor in Biology)</td>
<td align="left">Graduate student in Medical Genetics <br/><br/> Lab affiliation: Dr. Wendy Robinson, BCCHRI.</td>
<td align="center">Research biological context and contribute to the group’s understanding of the biology, as well as to interpretation of the data. Contribute to both the statistical analysis with curiosity and limited background in methylation data analysis, as well as contribute to the shared planning/organizing and writing-up of the final project.</td>
</tr>
<tr class="odd">
<td align="left">Abdullah Farouk <br/><br/> BA Economics and Finance</td>
<td align="left">MSc. graduate student in Statistics <br/><br/> Affiliations: Currently working on detecting outliers in financial data with Professors Natalia Nolde and Harry Joe.</td>
<td align="center">Hope to analyze the data with Nikolas to find statistically meaningful relationships in the data set. In addition to testing our hypothesis I plan to aid the visualization of our data set and assist in developing novel techniques of analyzing our data.</td>
</tr>
<tr class="even">
<td align="left">Nikolas Krstic <br/><br/> BSc. Statistics (Minor in Biology)</td>
<td align="left">MSc. graduate student in Statistics. <br/><br/> Affiliations: Currently a part-time analyst working at the British Columbia Centre for Disease Control (BCCDC).</td>
<td align="center">Investigate and apply several different statistical approaches to analyze the data (with Abdullah). I will also help bridge the gap between the statistical analysis and the biological interpretation of our results. Therefore, I will heavily participate in the methods, results and discussion of our final report.</td>
</tr>
</tbody>
</table>

Dataset
-------

We will combine two DNA methylation datasets publicly available on GEO. These datasets were generated using the Illumina HumanMethylation450 BeadChip, and both contain fluorescence intensity information corresponding to the CpG methylation status at 485,512 sites. These data sets are available in their raw forms (no pre-processing of the methylation signal intensity).

<table style="width:96%;">
<colgroup>
<col width="6%" />
<col width="69%" />
<col width="19%" />
</colgroup>
<thead>
<tr class="header">
<th align="left">Accession Number</th>
<th align="left">Data Description</th>
<th align="center">Data Source</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td align="left">GSE57767</td>
<td align="left">20 placental samples from normal pregnancies delivered at term, 19 samples delivered at term with preeclampsia (≥37 weeks) and 12 preterm preeclamptic samples. (51 columns, 485,512 rows)</td>
<td align="center">(Anton et al. 2014)</td>
</tr>
<tr class="even">
<td align="left">GSE44667</td>
<td align="left">20 chorionic villi samples of early-onset preeclampsia, and 20 gestational-age-matched controls. (40 columns, 485,512 rows)</td>
<td align="center">(Blair et al. 2013)</td>
</tr>
</tbody>
</table>

A third dataset that we may use to supplement our work is a preterm birth dataset from the Robinson lab at the BCCHRI. This is a dataset containing the methylation data from 44 patients (22 preterm cases affected by chorioamnionitis, 22 preterm controls). The chorioamnionitis cases may be excluded since our focus is on preeclampsia. This is an Illumina MethylationEPIC dataset, in an IDAT file format containing the raw fluorescence intensity data at 866, 895 CpG sites across the genome (44 columns, 866,895 rows). 90% of the sites from the Illumina 450K array are common to the Illumina EPIC probeset. To compare the 450K data to the EPIC data, we would select for the probes common to both datasets. This is made possible by the annotated probe identifiers in the Illumina array design. Other demographic factors included within the datasets are sex, ethnicity of the mother, gestational age, etc.

Aims and methodology
--------------------

With this project, our aim is to determine how the global DNA methylation landscape is different in the placenta of preeclamptic preterm birth cases versus normal pregnancies. Also, we would like to investigate whether subgroups exist among preterm birth cases. Lastly, by classifying the genic region that differentially methylated sites occupy, their effect on relevant biological pathways can be hypothesized. It is commonly accepted that promoter-CpG methylation is correlated with decreased gene expression, while gene-body methylation is often correlated with increased gene expression (Yang et al. 2014).

Since multiple datasets will be pooled together, identical preprocessing will be conducted on the raw data forms to ensure our analysis is consistent. This preprocessing will produce “beta” values; the standard metric that represents the proportion of methylation at interrogated sites. The exact method of preprocessing we intend to implement will be further investigated, but will involve such steps as normalization and filtering probes containing polymorphic sites or predicted to cross-hybridize to more than one genomic location (Wilhelm-Benartzi 2013). Subsequently, we will conduct principal component analysis (PCA) to identify potential batch effects within our dataset.

Clustering techniques will also help us establish subgroups among the samples. Alternatively, these techniques can be used to find clusters of co-methylated sites or possibly assign sites to genes or gene regions. Appropriate techniques include hierarchical agglomerative clustering and Recursive-Partitioning-Mixture Model (RPMM), although there are many others we plan to explore (Houseman et al. 2008). A novel way we plan to do so is to make use of linear grouping algorithms (Harrington 2010). They form linear clusters of the data and could help us identify biological functions or pathways that the genes containing methylated loci cluster into.

For the analysis to address our objectives, multiple approaches will be applied to simultaneously make comparisons and verify our conclusions. First, we will model the data using both Beta regression and standard linear regression methods (Siegmund 2011). This is because our response (DNA Methylation) is continuous and bounded between 0 and 1. We will then use random forests to predict if a particular genetic subgrouping belongs to a preterm or normal term placental sample, using methylation sites as covariates (Xi Chen 2012). This should help verify our results from clustering.

References
----------

1.  Statistics Canada. 2016. Preterm live births in Canada, 2000 to 2013. \[Catalogue number 82-625-X\]. Retrieved February 14, 2018 from Statistics Canada: <http://www.statcan.gc.ca/pub/82-625-x/2016001/article/14675-eng.htm>

2.  Anton L, Brown AG, Bartolomei MS, Elovitz MA. Differential methylation of genes associated with cell adhesion in preeclamptic placentas. PLoS One 2014;9(6):e100148

3.  Blair JD, Yuen RK, Lim BK, McFadden DE et al. Widespread DNA hypomethylation at gene enhancer regions in placentas associated with early-onset pre-eclampsia. Mol Hum Reprod 2013 Oct;19(10):697-708.

4.  Yeung KR, Chui CL, Ridlsey R, Makris A, Hennessy A, Lind JM. 2016. DNA methylation profiles in preeclampsia and healthy control placentas. Am J Physiol Heart Circ Physiol. 310:H1295-H1303.

5.  Yang X, Han H, De Carvalho DD, Lay FD, Jones PA, Liang G. 2014. Gene body methylation can alter gene expression and is a therapeutic target in cancer. Cancer Cell. 26(4):577-590.

6.  Wilhelm-Benartzi, C., Koestler, D., Karagas, M., Flanagan, J., Christensen, B., Kelsey, K., et al. (2013). Review of processing and analysis methods for DNA methylation array data. British Journal of Cancer, 109(6), 1394-1402.

7.  Houseman, E., Christensen, B., Yeh, R., Marsit, C., Karagas, M., Wrensch, M., et al. (2008). Model-based clustering of DNA methylation array data: A recursive-partitioning algorithm for high-dimensional data arising as a mixture of beta distributions. Bmc Bioinformatics, 9(1), 365-365.

8.  Harrington, J. and Salibian-Barrera, M. (2010). Finding approximate solutions to combinatorial problems with very large data sets using BIRCH.Computational Statistics and Data Analysis, 54, 655-667.

9.  Siegmund, K. D. (2011). Statistical approaches for the analysis of DNA methylation microarray data. Human Genetics, 129(6), 585-595.

10. Chen X, Ishwaran H. Random forests for genomic data analysis. Genomics. 2012 Jun 30;99(6):323-9.
