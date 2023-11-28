## Import XICRA counts

# Example to import XICRA data

In the following file `import_example.Rmd` find an example in R Markdown format to retrieve and load isomiR counts from XICRA pipeline. Additionally, find the html version in `import_example.html`.

XICRA creates the count of the reads at the isomiR level. In order to store information, XICRA generates unique entries using the parent miRNA, the variant and UID separated by character `&`.

This information is stored within the report/miRNA/ folder under the XICRA project generated named as: 
- `miRNA_expression-miraligner.csv`: that contains the count matrix for each isomiRs as rows per sample (columns)
- `miRNA_expression-miraligner_seq.csv` that contain the UID and the corresponding DNA sequence
- `miRNA_expression-miraligner_dup.csv` that contain duplicated entries with the same UID.

## Get data

Basically, we will need to load the package `XICRA.stats`
```
To get counts at the isomiR level, you can use the function within named `XICRA.stats::prepare_data` that basically uses `tidyr` to split strings by `&`

As an example:

```{r}

## XICRA data
isomir_data <- XICRA.stats::prepare_data(XICRA_csv_file, 'XICRA', checkNames = TRUE)
```

You might be interested in getting counts at the miRNA level. In order to do the analysis at the miRNA level, we need to aggregate isomiR counts by parent miRNA. We can use the `XICRA.stats::prepare_counts_by_miRNA` that basically uses `dplyr::group_by` and `dplyr::summarize_if` to aggregate accordingly.

As an example:

```{r}
## miRNA and variant data
miRNA_data <- XICRA.stats::prepare_counts_by_miRNA(isomir_data)
```

Finally, you can also get counts at the variant level using the isomiR counts.

```{r}
miRNA_variant_data <- XICRA.stats::prepare_counts_by_variant(isomir_data)
head(miRNA_variant_data)
dim(miRNA_variant_data)
```
