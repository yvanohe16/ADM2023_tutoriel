R Notebook
================

``` r
download.file("https://github.com/ANF-MetaBioDiv/course-material/archive/refs/heads/main.zip", destfile = "main.zip")
```

``` r
refdb_folder <- here::here("data", "refdb")
refdb_folder
```

    ## [1] "/Users/yvanohemartinier/ADM_tutoriel_dada2/data/refdb"

``` r
if (!dir.exists(refdb_folder)) dir.create(refdb_folder, recursive = TRUE)
```

``` r
getOption("timeout")
```

    ## [1] 60

``` r
# so we change timeout to be 20 minutes
options(timeout = 1200)

# we save in variable the path to the refdb
# in the working space
silva_train_set <- file.path(refdb_folder,
                             "silva_nr99_v138.1_train_set.fa.gz")

silva_species_assignment <- file.path(refdb_folder,
                                      "silva_species_assignment_v138.1.fa.gz")

# then we download the files if they don't already exist

if (!file.exists(silva_train_set)) {
  download.file(
    "https://zenodo.org/record/4587955/files/silva_nr99_v138.1_train_set.fa.gz",
    silva_train_set,
    quiet = TRUE
  )
}

if (!file.exists(silva_species_assignment)) {
  download.file(
    "https://zenodo.org/record/4587955/files/silva_species_assignment_v138.1.fa.gz",
    silva_species_assignment,
    quiet = TRUE
  )
}
```

``` r
# devtools::load_all()
```

``` r
path_to_fastqs <- here::here("data", "raw")
```

``` r
fnFs <- sort(list.files(path_to_fastqs,
                        pattern = "_R1.fastq.gz",
                        full.names = TRUE))
```

``` r
fnRs <- sort(list.files(path_to_fastqs,
                        pattern = "_R2.fastq.gz",
                        full.names = TRUE))
```

## 2.2 Extract sample names

``` r
sample_names <- basename(fnFs) |>
  strsplit(split = "_") |>
  sapply(head, 1)
```

``` r
basename(fnFs) |>
  head()
```

    ## [1] "S11B_R1.fastq.gz" "S1B_R1.fastq.gz"  "S2B_R1.fastq.gz"  "S2S_R1.fastq.gz" 
    ## [5] "S3B_R1.fastq.gz"  "S3S_R1.fastq.gz"

``` r
basename(fnFs) |>
  strsplit(split = "_") |>
  head()
```

    ## [[1]]
    ## [1] "S11B"        "R1.fastq.gz"
    ## 
    ## [[2]]
    ## [1] "S1B"         "R1.fastq.gz"
    ## 
    ## [[3]]
    ## [1] "S2B"         "R1.fastq.gz"
    ## 
    ## [[4]]
    ## [1] "S2S"         "R1.fastq.gz"
    ## 
    ## [[5]]
    ## [1] "S3B"         "R1.fastq.gz"
    ## 
    ## [[6]]
    ## [1] "S3S"         "R1.fastq.gz"

``` r
basename(fnFs) |>
  strsplit(split = "_") |>
  sapply(head, 1) |>
  head()
```

    ## [1] "S11B" "S1B"  "S2B"  "S2S"  "S3B"  "S3S"

``` r
gsub("^.+/|_.+$", "", fnFs) |> head()
```

    ## [1] "S11B" "S1B"  "S2B"  "S2S"  "S3B"  "S3S"

## 3.Sequence quality check

\#\`\`\`{r} \# create a directory for the outputs quality_folder \<-
here::here(“outputs”, “dada2”, “quality_plots”)

if (!dir.exists(quality_folder)) { dir.create(quality_folder, recursive
= TRUE) }

qualityprofile(fnFs, fnRs, file.path(quality_folder,
“quality_plots.pdf”)) \`\`\`
