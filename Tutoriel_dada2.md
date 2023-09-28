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