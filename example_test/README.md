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

```{r}
library(XICRA.stats)
```

```{r}
#############
## Get data
#############
XICRA_csv_file <- "./example_miRNA_files/miRNA_expression-miraligner.csv"
```
## Get counts at the isomiR level

You can use the function within named `XICRA.stats::prepare_data` that basically uses `tidyr` to split strings by `&`

```{r}
XICRA.stats::prepare_data
```

As an example:

```{r}

## XICRA data
isomir_data <- XICRA.stats::prepare_data(XICRA_csv_file, 'XICRA', checkNames = TRUE)
head(isomir_data)
names(isomir_data)
dim(isomir_data)
## set rownames
row.names(isomir_data) <- isomir_data$ID
```

## Get counts at the miRNA level

In order to do the analysis at the miRNA level, we need to aggregate counts by parent miRNA. We can use the `XICRA.stats::prepare_counts_by_miRNA` that basically uses `dplyr::group_by` and `dplyr::summarize_if` to aggregate accordingly.

```{r}
XICRA.stats::prepare_counts_by_miRNA
```
As an example:

```{r}
## miRNA and variant data
miRNA_data <- XICRA.stats::prepare_counts_by_miRNA(isomir_data)
head(miRNA_data)
dim(miRNA_data)
## set rownames
rownames(miRNA_data) <- miRNA_data$parent 
```

## Get counts at the variant level

Finally, we can also aggregate counts by variant, using the `XICRA.stats::prepare_counts_by_variant` that again groups and summarises using variant type.


```{r}
XICRA.stats::prepare_counts_by_variant
```


```{r}
miRNA_variant_data <- XICRA.stats::prepare_counts_by_variant(isomir_data)
head(miRNA_variant_data)
dim(miRNA_variant_data)
```

## Clean and prepare the final data

Before finishing we can discard columns that are not counts and will interfere with the analysis

```{r}
miRNA_data$parent <- NULL
isomir_data$parent <- NULL
isomir_data$variant <- NULL
isomir_data$UID <- NULL
isomir_data$ID <- NULL
```

We can create a list with the dataframes generated and removed them to increase readability of the code.
```{r}
my_data <- list("isomir_data" = isomir_data, 
                "miRNA_data"=miRNA_data, 
                "miRNA_variant_data" = miRNA_variant_data)


```
