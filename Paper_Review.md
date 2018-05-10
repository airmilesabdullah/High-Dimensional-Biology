Paper Review
================
Abdullah Farouk
2018-02-19

[Gene Expression by Mouse Inner Ear Hair Cells during Development](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4405555/)
-------------------------------------------------------------------------------------------------------------------------

### Goals of the Study:

The authors try to identify the genes involved in inner ear mouse development. This is done to help study genes that cause deafness. They do so by creating profiles of genes expressed in hair cells (HC) relative to surrounding cells (SC); other cells in the sensory epithelium, from different tissue types across different stages of cell development. They find hereditary deafness genes highly expressed in hair cells (like Grxcr2 whose mutations cause recessive hearing loss). They conclude on the suggestion that the differential expression database they have compiled could help understand the causes of complex inner ear disorders.

### Data Collection and Analysis

RNA-Seq was used to generate HC and SC mRNAs. The mRNA expression data was stored in the Shared Harvard Inner Ear Laboratory Database (SHIELD). It primarily contains information on cell type (HCs and SCs), tissue source (cochlea and utricle), developmental stage (E16, P0, P4, and P7), Fold change values (FC) and FDR values.

Data on processes carried out like QPCR, immunocytochemistry, in situ hybridization and a second mouse strain at P0 were collected but not made available.

The experiment was carried out on a mouse transgenic strains. It expressed green fluorescence protein (GFP). A second strain from another mouse that expressed tdTomato was also used to validate results.

Cells from cochleae and utricles were separated using enzymatic dissociation followed by FACS. Transcriptomes of the purified HCs and SCs were then studied by High Throughput Sequencing (HTS). Each of 16 samples was sequenced to a depth of 40-100 million reads.

They compute a statistic known as Fold change (FC). It’s a ratio of mean gene expression levels pooled across 8 samples (GFP+/GFP-). This is used to identify gene enrichments in HCs and SCs. However the thresholds chosen seem to be arbitrarily selected with little explanation on what those numbers mean.

The authors used Principal Component Analysis (PCA) to understand the primary sources of variation in genes expressed. PCA also validated the reproducibility of results across samples.

PC1 is seen to be strongly associated with cell type (Fig 2a shows all the SC samples are to the left and all the HC samples to the right). PC2 on the other hand is loosely related to developmental stage (points corresponding to stages P4 and P7 are all at the top). PCA analysis was conducted on 24 samples from both mouse strains. 1/3 of the samples came from the second mouse. An interesting comparison would have been to remove data from the second mouse and see if the same factors explained most of the variation.

Hierarchical clustering and heat maps were used to visually identify new HC and SC genes. They first computed log transformations of the normalized expression values. These values were then standardized for each gene and averaged, for each gene, across samples. To assign genes to gene clusters the authors first defined the distance between two genes as correlation of the expression of the 3′-tags of two genes across samples. They then computed the distance between two gene clusters as the difference in the average expression values of the genes in each cluster (centroid linkage).

They found genes critical for mechanotransduction expressed in HCs but not SCs (eg:Cdh23 and Tmc1). Further, these genes were expressed early in the utricle but not in the cochlea as found in previous studies.

To identify deafness genes (DGs) they ranked the 18,199 genes expressed by enrichment level in HCs. They discovered that new DGs known to be HC specific (Usher syndrome genes) were enriched in HCs.

### Conclusion

This paper focuses on genes in HCs that are responsible for deafness. They identify new DGs and those mentioned in earlier literature. There is mention of multiple comparison tests run on normalized expression levels but there is no information of what they are. Further there is no mention of the results of these tests which makes it difficult to interpret the significance of the FDR values reported. Another limitation is that genes identified as highly enriched in HCs are enriched only relative to SCs.
