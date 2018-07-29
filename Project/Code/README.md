README
================

This folder is our codebase. It contains the different scripts used during different stages in our analysis. Please find a brief description of each script and a corresponding link to them in the table below:

Table of Contents
-----------------

<table style="width:100%;">
<colgroup>
<col width="18%" />
<col width="33%" />
<col width="48%" />
</colgroup>
<thead>
<tr class="header">
<th align="left">Script Name</th>
<th align="left">Description</th>
<th align="center">Link</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td align="left">pre_processing_data.Rmd</td>
<td align="left">This script details how we cleaned, computed beta and M values and merged the two GEO datasets used in this study.</td>
<td align="center"><a href="pre_processing_data.md">Click here</a></td>
</tr>
<tr class="even">
<td align="left">exploratory_and_limma_analyses.md</td>
<td align="left">This script builds on the pre processing script. It contains the code used to aggregate CpG probes to CpG Islands. <br/><br/> It also details the exploratory data analysis procedures undertaken to familiarize ourselves with the data at hand and the initial steps of our analysis (e.g.identifying the top differentially methylated cites )</td>
<td align="center"><a href="exploratory_and_limma_analyses.md">Click here</a></td>
</tr>
<tr class="odd">
<td align="left">Machine_Learning_Methods.Rmd</td>
<td align="left">This script describes the different supervised learning models tested to classify placental cells by disease status using CpG Islands as covariates. <br/><br/> It also highlights the cross validation techniques used to train, tune and choose the best performing model.</td>
<td align="center"><a href="Machine_Learning_Methods.Rmd">Click here</a></td>
</tr>
<tr class="even">
<td align="left">GO_analysis.md</td>
<td align="left">This script details the association of CGI sites to genes to then protiens which were ultimately used to identify a particular biological function.</td>
<td align="center"><a href="GO_analysis.md">Click here</a></td>
</tr>
<tr class="odd">
<td align="left">m_value_normality_check.md</td>
<td align="left">This script verifies the normality of the m values used as the response variable in the linear regression model fit to identiy the top differentially methylated CpG Islands.</td>
<td align="center"><a href="m_value_normality_check.md">Click here</a></td>
</tr>
<tr class="even">
<td align="left">test_data_set_creation.md</td>
<td align="left">This script details the procedures used to clean and normalize a publicly available methylation dataset on preeclempsia (GSE73375). <br/><br/> The test set was to be used to asses the generalizability of our chosen classification model</td>
<td align="center"><a href="test_data_set_creation.md">Click here</a></td>
</tr>
</tbody>
</table>
