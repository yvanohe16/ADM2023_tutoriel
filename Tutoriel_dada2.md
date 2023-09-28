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
devtools::load_all(path="/Users/yvanohemartinier/ADM_tutoriel_dada2/course-material-main/R")
```

    ## ℹ Loading ANF_metaB

``` r
path_to_fastqs <- here::here("data", "raw")
```

``` r
fnFs <- sort(list.files(path_to_fastqs,
                        pattern = "_R1.fastq.gz",
                        full.names = TRUE))
```

\#avec cette commande va chercher tous les noms de fichier
“\_R1.fastaq.gz” Mais on ne les veut pas tous

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

\#utilisation d’un plit, sample_names \<- basename(fnFs) \|\> Dans
basename on va rajouter une conditions, la suivante strsplit(split =
“*”) peremt de decouper les noms de base qui sont séparé avec le ”*”

``` r
basename(fnFs) |>
  head()
```

    ## [1] "S11B_R1.fastq.gz" "S1B_R1.fastq.gz"  "S2B_R1.fastq.gz"  "S2S_R1.fastq.gz" 
    ## [5] "S3B_R1.fastq.gz"  "S3S_R1.fastq.gz"

# avant découpage

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

# après découpage

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

``` r
# create a directory for the outputs
quality_folder <- here::here("outputs",
                             "dada2",
                             "quality_plots")

if (!dir.exists(quality_folder)) {
  dir.create(quality_folder, recursive = TRUE)
}

qualityprofile(fnFs,
               fnRs,
               file.path(quality_folder, "quality_plots.pdf"))
```

    ## quartz_off_screen 
    ##                 2

``` r
path_to_trimmed_reads <- here::here(
  "outputs",
  "dada2",
  "trimmed"
)

if (!dir.exists(path_to_trimmed_reads)) dir.create(path_to_trimmed_reads, recursive = TRUE)
```

``` r
primer_fwd  <- "CCTACGGGNBGCASCAG"
primer_rev  <- "GACTACNVGGGTATCTAAT"
```

``` r
Biostrings::readDNAStringSet(
  fnFs[1],
  format = "fastq",
  nrec = 10
)
```

    ## DNAStringSet object of length 10:
    ##      width seq                                              names               
    ##  [1]   293 CCTACGGGGGGCAGCAGTAGGGA...ACATCGGCTTAACCGATGAAGT M01522:260:000000...
    ##  [2]   293 CCTACGGGTGGCACCAGTAGGGA...CGGGGCTTAACCTCGGAACTGC M01522:260:000000...
    ##  [3]   292 CCTACGGGGCGCAGCAGGCGCGA...GGGACCGGGAGAGGTGTGAGGT M01522:260:000000...
    ##  [4]   293 CCTACGGGGTGCAGCAGTAGGGA...TCAAAACTCCCAGTCTAGAGTT M01522:260:000000...
    ##  [5]   291 CCTACGGGTGGCAGCAGTGGGGA...GCAGTGGAAACTGTTGGGCTTG M01522:260:000000...
    ##  [6]   293 CCTACGGGATGCAGCAGGCGCGA...GGGACCGGGAGAGGTGTGGGGG M01522:260:000000...
    ##  [7]   292 CCTACGGGATGCAGCAGTGGGGA...TTTAATCCTGATGAGCTAGAAA M01522:260:000000...
    ##  [8]   293 CCTACGGGGCGCAGCAGTAGGGA...TTAAAACTTTTGTTCTGGAATT M01522:260:000000...
    ##  [9]   292 CCTACGGGTTGCAGCAGTGGGGA...ATTAAAACTTTTCAGCTAGAGT M01522:260:000000...
    ## [10]   293 CCTACGGGAGGCAGCAGTGGGGA...CCCGGGCTCAACCTGGGAACGG M01522:260:000000...

``` r
Biostrings::readDNAStringSet(
  fnRs[1],
  format = "fastq",
  nrec = 10
)
```

    ## DNAStringSet object of length 10:
    ##      width seq                                              names               
    ##  [1]   301 GACTACCAGGGTATCTAATCCTG...GGCTGCTGGCACGAAGTTCGCC M01522:260:000000...
    ##  [2]   301 GACTACCGGGGTATCTAATCCTG...GGCTGCTGGCACGGAGTTAGCC M01522:260:000000...
    ##  [3]   300 AATCCGGTTCGTGCCCCTAGGCT...TCTTTCCCAGCCCTTATTCCAA M01522:260:000000...
    ##  [4]   301 GACTACCGGGGTATCTAATCCTG...GGCTGCTGGCACGGAGTTAGCC M01522:260:000000...
    ##  [5]   301 GACTACCGGGGTATCTAATCCCT...GGCTGCTGGCCCGGAATTAGCC M01522:260:000000...
    ##  [6]   301 GGTATCTAATCCGGTTCGTGCCC...CACCGTCCTTACCCCCCCCTTT M01522:260:000000...
    ##  [7]   301 GGTATCTAATCTTGTTTGCTCCC...CCCGACGTTAGCCGGGGCTTCT M01522:260:000000...
    ##  [8]   301 GACTACGAGGGTATCTAATCCCG...GGCTGCTGGCACGGAATTAGCC M01522:260:000000...
    ##  [9]   301 GGTATCTAATCCTCTTCGCTACC...CACGAAGTTAGCCGGACCTTCT M01522:260:000000...
    ## [10]   301 GACTACGGGGGTATCTAATCCTG...GGCTGCCGGCACGGGGTTAGCC M01522:260:000000...
