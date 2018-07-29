README
================

This folder contains the data used in our analysis, as well as data objects created and saved after key analysis steps. Please find in the table below a brief description of each object stored in this folder, and link corresponding to each:

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
<td align="left">Beta_CGI_data.RDS</td>
<td align="left">This object contains the unnormalized beta values for all CGIs in the 85 samples used in our analyses.</td>
<td align="center"><a href="Beta_CGI_data.RDS">Click here</a></td>
</tr>
<tr class="even">
<td align="left">CGI_SL_Set.rdata</td>
<td align="left">This object contains the names of the CpG island sites selected by CV limma.</td>
<td align="center"><a href="CGI_SL_Set.rdata">Click here</a></td>
</tr>
<tr class="odd">
<td align="left">GO.xml</td>
<td align="left">This file contains the latest GO terms, and was used in the gene ontology analysis.</td>
<td align="center"><a href="GO.xml">Click here</a></td>
</tr>
<tr class="even">
<td align="left">M_CGI_data.RDS</td>
<td align="left">This object contains the normalized M values corresponding to all 19,394 CGI sites analyzed in our 85 samples.</td>
<td align="center"><a href="M_CGI_data.RDS">Click here</a></td>
</tr>
<tr class="odd">
<td align="left">Test_Set_Metadata.RDS</td>
<td align="left">This object contains the metadata for the 36 samples from GSE73375 that we were intending to use to test our trained machine learning approach (disease class with levels &quot;normal&quot; and &quot;preeclampsia&quot;, and gestational age with levels &quot;preterm&quot; and &quot;term&quot;).</td>
<td align="center"><a href="Test_Set_Metadata.RDS">Click here</a></td>
</tr>
<tr class="even">
<td align="left">m_val_test_set.rds</td>
<td align="left">This object contains the M values corresponding to each of the 36 samples from GSE73375, the dataset we intended to test our SVM approach with.</td>
<td align="center"><a href="m_val_test_set.RDS">Click here</a></td>
</tr>
<tr class="odd">
<td align="left">meta_data.csv</td>
<td align="left">This object contains the metadata for the 85 samples used in our analyses, and contains a disease class column with 2 levels: &quot;normal&quot; and &quot;preeclampsia&quot;, and a gestational age column with 2 levels: &quot;term&quot; and &quot;preterm&quot;, as well as sample identifiers.</td>
<td align="center"><a href="meta_data.csv">Click here</a></td>
</tr>
<tr class="even">
<td align="left">norm_m_values_test_set.rds</td>
<td align="left">This object contains the normalized M values for the 36 samples from GSE73375.</td>
<td align="center"><a href="norm_m_values_test_set.rds">Click here</a></td>
</tr>
<tr class="odd">
<td align="left">unnorm_m_val_test_set.rds</td>
<td align="left">This object contains the unnormalized M values for the 36 samples from GSE73376.</td>
<td align="center"><a href="unnorm_m_val_test_set.rds">Click here</a></td>
</tr>
</tbody>
</table>
