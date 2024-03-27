02-locus-HWE-evaluation
================
Diana Baetscher
2023-03-30

Using the complete dataset to evaluate HWE across loci and populations
in the baseline.

Updates: use the NAs explicit genos file. Don’t filter loci based on
missing data at this point. Just look at reference pops for evaluating
HWE.

Evaluating the loci that we’re using for self-assignment and bycatch
assignment. These were originally derived from lcWGS data and the
microhaplotype results suggest that the real variation present is
potentially quite different. Initially, I should probably just focus on
the baseline genotypes rather than including everything.

To evaluate the markers, I need to subset the reference baseline
genotypes for just the minimal amount of missing data.

``` r
library(tidyverse)
```

    ## ── Attaching core tidyverse packages ──────────────────────── tidyverse 2.0.0 ──
    ## ✔ dplyr     1.1.3     ✔ readr     2.1.4
    ## ✔ forcats   1.0.0     ✔ stringr   1.5.0
    ## ✔ ggplot2   3.4.4     ✔ tibble    3.2.1
    ## ✔ lubridate 1.9.3     ✔ tidyr     1.3.0
    ## ✔ purrr     1.0.2     
    ## ── Conflicts ────────────────────────────────────────── tidyverse_conflicts() ──
    ## ✖ dplyr::filter() masks stats::filter()
    ## ✖ dplyr::lag()    masks stats::lag()
    ## ℹ Use the conflicted package (<http://conflicted.r-lib.org/>) to force all conflicts to become errors

``` r
library(adegenet)
```

    ## Loading required package: ade4
    ## 
    ##    /// adegenet 2.1.10 is loaded ////////////
    ## 
    ##    > overview: '?adegenet'
    ##    > tutorials/doc/questions: 'adegenetWeb()' 
    ##    > bug reports/feature requests: adegenetIssues()

``` r
library(radiator)
library(DescTools)
```

``` r
# read in rds file with genotypes
genos_long <- read_rds("../data/processed/called_genos_na_explicit.rds")

# just select the reference samples based on metadata
samplesheets <- read_rds("../data/processed/metadata_bycatch_and_reference_20230413.rds")

reference_samples <- samplesheets %>%
  filter(Sample_Plate %in% c("BFAL001", "BFAL002") &
           !is.na(Location),
         !`Loc-Abbr` %in% c("Oahu", "Kauai, Kilauea")) %>%
  arrange(sampleID) %>%
  select(sampleID, gtseq_run, id, Location, `Loc-Abbr`)

# combine the genotypes with the reference
ref_genos <- reference_samples %>%
  left_join(., genos_long, by = c("gtseq_run", "id")) 
```

Deal explicitly with duplicate samples:

## Some initial filters

### Take highest read-depth call for multiply-genotyped DNA_IDs

I’m not sure if there are any of these, but best to leave it in here…

Now, here is a harder operation: if an individual is multiply-genotyped,
take the genotype with the highest total read depth.

``` r
# slow-ish function to get the total read depth column
tdepth <- function(a, d) {
  if(any(is.na(a))) {
    return(NA)
  }
  if(a[1]==a[2]) {
    return(d[1])
  } else {
    return(d[1] + d[2])
  }
  
}
# this takes the highest read-depth instance of each duplicately-genotyped individual.
geno_one_each <- ref_genos %>%
  group_by(sampleID, locus, gtseq_run) %>%
  mutate(total_depth = tdepth(allele, depth)) %>%
  ungroup() %>%
  arrange(sampleID, locus, total_depth, gtseq_run, depth) %>%
  group_by(sampleID, locus) %>%
  mutate(rank = 1:n()) %>% # this ranks all alleles for a given individual by depth
  ungroup() %>%
  filter(rank <= 2) # this takes the top ranked alleles and should remove duplicates
  

# how many samples now?
geno_one_each %>%
  filter(!is.na(sampleID)) %>%
  group_by(sampleID) %>% 
  select(sampleID, gtseq_run, id) %>%
  unique() %>%
  tally() 
```

    ## # A tibble: 148 × 2
    ##    sampleID     n
    ##    <chr>    <int>
    ##  1 12-0155      1
    ##  2 12-0270      1
    ##  3 12-0271      1
    ##  4 12-0371      1
    ##  5 12-0372      1
    ##  6 12-0381      1
    ##  7 13-0991      1
    ##  8 13-1187      1
    ##  9 13-1198      1
    ## 10 14-0454      1
    ## # ℹ 138 more rows

Take a look at missing data overall:

``` r
geno_one_each %>%
  ggplot(aes(x = reorder(locus, depth), y = reorder(sampleID, depth), fill = log10(depth))) +
  geom_tile()
```

![](02-locus-HWE-evaluation_files/figure-gfm/unnamed-chunk-2-1.png)<!-- -->

Deal with missing data:

## Missing data in loci

How many loci and how many alleles?

``` r
# alleles
geno_one_each %>%
  filter(!is.na(allele)) %>%
  select(locus, allele) %>%
  unique()
```

    ## # A tibble: 458 × 2
    ##    locus                allele           
    ##    <chr>                <chr>            
    ##  1 scaffold_0_4543438   GGG              
    ##  2 scaffold_0_4543438   GGA              
    ##  3 scaffold_1050_36334  T                
    ##  4 scaffold_1053_177405 A                
    ##  5 scaffold_1078_116737 GC               
    ##  6 scaffold_10_1333563  G                
    ##  7 scaffold_10_1333563  A                
    ##  8 scaffold_115_2200290 C                
    ##  9 scaffold_1164_119437 TATTATACTGGTATTAG
    ## 10 scaffold_116_1932131 CG               
    ## # ℹ 448 more rows

``` r
# loci
geno_one_each %>%
  filter(!is.na(allele)) %>%
  select(locus, allele) %>%
  unique() %>%
  group_by(locus) %>%
  tally() %>%
  arrange(desc(n))
```

    ## # A tibble: 189 × 2
    ##    locus                    n
    ##    <chr>                <int>
    ##  1 scaffold_1226_115783     6
    ##  2 scaffold_1589_18141      5
    ##  3 scaffold_316_437138      5
    ##  4 scaffold_33_2811656      5
    ##  5 scaffold_48_3128896      5
    ##  6 scaffold_75_1184494      5
    ##  7 scaffold_0_4543438       4
    ##  8 scaffold_1164_119437     4
    ##  9 scaffold_127_1105814     4
    ## 10 scaffold_15_2923659      4
    ## # ℹ 179 more rows

In the reference, 458 alleles across 189 loci with 1-6 alleles per
locus.

Some individuals are missing all their data, and some loci are missing
all their data. Beginning with loci:

``` r
# missing data across loci
locs_to_toss <- geno_one_each %>%
  group_by(locus) %>%
  mutate(missingness = ifelse(is.na(allele), 1, 0)) %>%
  summarise(sum(missingness)) %>% # 189 loci x2 = total
  filter(`sum(missingness)`>151) %>% # more than 50% missing data (given 151 individuals)
  select(locus) # drop those loci for now and see how the assignment goes

# just the keepers
genos_locs_filtered <- geno_one_each %>%
  anti_join(., locs_to_toss)
```

    ## Joining with `by = join_by(locus)`

Remove these 4 loci because of too much missing data in the reference
samples.

## Missing data in individuals

Total number of loci = 185

Total number of samples = 151

``` r
inds_to_toss <- genos_locs_filtered %>%
  group_by(sampleID) %>%
  mutate(missingness = ifelse(is.na(allele), 1, 0)) %>%
  summarise(sum(missingness)) %>%
  arrange(desc(`sum(missingness)`)) %>%
  #filter(`sum(missingness)` > 47) # remove samples with >25% missing data
 filter(`sum(missingness)` > 19) # remove samples with >10% missing data

# just the keepers
genos_locs_ind_filtered <- genos_locs_filtered %>%
  anti_join(., inds_to_toss)
```

    ## Joining with `by = join_by(sampleID)`

23 samples removed with \>10% missing data.

Now look at missing data:

``` r
genos_locs_ind_filtered %>%
  ggplot(aes(x = reorder(locus, -depth), y = reorder(sampleID, -depth), fill = log10(depth))) +
  geom_tile()
```

![](02-locus-HWE-evaluation_files/figure-gfm/unnamed-chunk-6-1.png)<!-- -->

Much better. Now move forward from here.

## Prep data for analysis

``` r
# first make integers of the alleles
alle_idxs <- genos_locs_ind_filtered %>% 
  filter(!is.na(sampleID)) %>%
  dplyr::select(sampleID, locus, gene_copy, allele) %>%
  group_by(locus) %>%
  mutate(alleidx = as.integer(factor(allele, levels = unique(allele)))) %>%
  ungroup() %>%
  arrange(sampleID, locus, alleidx) # rubias can handle NA's, so no need to change them to 0's

# reformat
reference <- alle_idxs %>%
  inner_join(., reference_samples) %>%
  select(-allele, -gtseq_run, -id) %>%
  select(`Loc-Abbr`, sampleID, everything()) %>%
  rename(collection = `Loc-Abbr`, indiv = sampleID)
```

    ## Joining with `by = join_by(sampleID)`

I need to finagle the dataset into a genind object.

``` r
# read in data
reference_genos <- reference

# number of missing alleles (max = 185*2)
inds_to_toss_missing_data <- reference_genos %>%
  group_by(indiv) %>%
  filter(is.na(alleidx)) %>%
  group_by(indiv, collection) %>%
  tally() %>%
  arrange(desc(n)) %>%
  filter(n > 9)
  
# if I removed all the indivs with > 10 missing alleles, that's only 5 indivs.
# what would the distribution of indivs per population be? 
# still ok, I'm pretty sure.
reference_genos %>%
  anti_join(., inds_to_toss_missing_data) %>%
  group_by(collection, indiv) %>%
  tally() %>%
  tally()
```

    ## Joining with `by = join_by(collection, indiv)`

    ## # A tibble: 6 × 2
    ##   collection           n
    ##   <chr>            <int>
    ## 1 Japan, Torishima    40
    ## 2 Kauai, Lehua         2
    ## 3 NWHI, FFS           34
    ## 4 NWHI, Kure          22
    ## 5 NWHI, Laysan         8
    ## 6 NWHI, Midway         6

I feel fine about that. Let’s see if that makes any of the HWE
calculations easier because less missing data?

``` r
# SKIP THIS FOR NOW.
# ref_genos_low_missing_data <- reference_genos %>%
#   anti_join(., inds_to_toss_missing_data)

ref_genos_low_missing_data <- reference_genos
```

Make the df match the requirements for tidy_genomic_data

``` r
long_df <- ref_genos_low_missing_data  %>%
  select(-gene_copy) %>%
  select(collection, everything()) %>%
  rename(INDIVIDUALS = indiv, STRATA = collection, MARKERS = locus, GT = alleidx)
```

Genotypes should be coded with 3 integers for each alleles. 6 integers
in total for the genotypes. e.g. 001002 or 111333 (for heterozygote
individual). 6 integers WITH separator: e.g. 001/002 or 111/333 (for
heterozygote individual). The separator can be any of these: “/”, “:”,
“\_“,”-“,”.”, and will be removed.

``` r
# create 3 digit integers from the genotypes
long_df$GT3 <- Format(long_df$GT, ldigits = 3, digits = 0)

# fix NAs
long_df0s <- long_df %>%
  mutate(GT3 = ifelse(is.na(GT3), "000", GT3)) # I don't love that this creates potential artifacts!
```

Now combine the GT3 column per indiv/marker:

``` r
# make the genos characters and then try pasting them as strings
long_df0s$GT3 <- as.character(long_df0s$GT3)

long_df3digit <- long_df0s %>%
  group_by(INDIVIDUALS, MARKERS) %>% 
  arrange(GT3, .by_group = TRUE) %>% 
  summarise(GENOTYPE = toString(GT3))
```

    ## `summarise()` has grouped output by 'INDIVIDUALS'. You can override using the
    ## `.groups` argument.

``` r
# paste strings together
long_df3digit$GENOTYPE <- gsub(", ","",long_df3digit$GENOTYPE)


# add back on species identity as strata
df_for_conversion <- long_df0s %>% 
  select(-GT, -GT3) %>%
  left_join(., long_df3digit) %>%
  unique() %>%
  rename(GT = GENOTYPE) %>%
  mutate(GT = ifelse(GT == "000000", NA, GT)) %>%
  filter(!str_detect(STRATA, "Kauai"))
```

    ## Joining with `by = join_by(INDIVIDUALS, MARKERS)`

``` r
df_for_conversion$STRATA <- as.factor(df_for_conversion$STRATA)
```

Double check - how many samples per population?

``` r
df_for_conversion %>%
  select(INDIVIDUALS, STRATA) %>%
  unique() %>%
  group_by(STRATA) %>%
  tally()
```

    ## # A tibble: 5 × 2
    ##   STRATA               n
    ##   <fct>            <int>
    ## 1 Japan, Torishima    43
    ## 2 NWHI, FFS           34
    ## 3 NWHI, Kure          24
    ## 4 NWHI, Laysan         9
    ## 5 NWHI, Midway        13

Let’s remove Lehua and Kilauea at this point because of so few samples.

``` r
convert_df_wo_2pops <- df_for_conversion 
```

``` r
# use the radiator package for this conversion
genind_df <- write_genind(convert_df_wo_2pops)
```

Basic population genetic evaluation of the markers, following, in part,
this: <https://bookdown.org/hhwagner1/LandGenCourse_book/WE_3.html>

``` r
sum_output <- summary(genind_df)

names(sum_output)
```

    ## [1] "n"         "n.by.pop"  "loc.n.all" "pop.n.all" "NA.perc"   "Hobs"     
    ## [7] "Hexp"

``` r
expected_hz <- rownames_to_column(as.data.frame(sum_output$Hexp), var = "locus") %>%
  rename(Hexp = `sum_output$Hexp`)
observed_hz <- rownames_to_column(as.data.frame(sum_output$Hobs), var = "locus") %>%
  rename(Hobs = `sum_output$Hobs`)

hz_df <- expected_hz %>%
  left_join(., observed_hz)
```

    ## Joining with `by = join_by(locus)`

Expected heterozygosity (here: Hexp) is the heterozygosity expected in a
population under HWE, and observed heterozygosity (here: Hobs) is the
observed number of heterozygotes at a locus divided by the total number
of genotyped individuals. Here are the global values (pooled across all
populations):

FIS \[= (HS - HI)/HS\] HI based on observed heterozygosities in
populations HS based on expected heterozygosities in populations

``` r
hz_df %>%
  group_by(locus) %>%
  mutate(Fis = (Hexp - Hobs)/ Hexp) %>%
  ggplot(aes(x = locus, y = Fis))+
  geom_point()
```

    ## Warning: Removed 11 rows containing missing values (`geom_point()`).

![](02-locus-HWE-evaluation_files/figure-gfm/calculate-Fis-1.png)<!-- -->
Ooof. There are some big outliers. Need to revisit the expectations
given multiple pops and small pop sizes and/or re-run with the full
dataset?

``` r
# expected heterozygosity per population
adegenet::Hs(genind2genpop(genind_df))
```

    ## 
    ##  Converting data from a genind to a genpop object... 
    ## 
    ## ...done.

    ## Japan,_Torishima        NWHI,_FFS       NWHI,_Kure     NWHI,_Laysan 
    ##        0.3155225        0.3617960        0.3531330        0.3728518 
    ##     NWHI,_Midway 
    ##        0.3436854

``` r
Hobs <- t(sapply(seppop(genind_df), function(ls) summary(ls)$Hobs))
  Hexp <- t(sapply(seppop(genind_df), function(ls) summary(ls)$Hexp))
  {cat("Expected heterozygosity (Hexp):", "\n")
  round(Hexp, 2)
  cat("\n", "Observed heterozygosity (Hobs):", "\n")
  round(Hobs, 2)}
```

    ## Expected heterozygosity (Hexp): 
    ## 
    ##  Observed heterozygosity (Hobs):

    ##                  scaffold_0_4543438 scaffold_102_694013 scaffold_1050_36334
    ## Japan,_Torishima               0.28                0.19                0.14
    ## NWHI,_FFS                      0.53                0.53                0.24
    ## NWHI,_Kure                     0.58                0.46                0.46
    ## NWHI,_Laysan                   0.78                0.44                0.44
    ## NWHI,_Midway                   0.46                0.31                0.31
    ##                  scaffold_1053_177405 scaffold_1078_116737 scaffold_10_1333563
    ## Japan,_Torishima                 0.37                 0.49                0.49
    ## NWHI,_FFS                        0.32                 0.56                0.32
    ## NWHI,_Kure                       0.50                 0.46                0.50
    ## NWHI,_Laysan                     0.33                 0.56                0.44
    ## NWHI,_Midway                     0.54                 0.46                0.46
    ##                  scaffold_115_2200290 scaffold_1164_119437 scaffold_116_1932131
    ## Japan,_Torishima                    0                 0.02                 0.33
    ## NWHI,_FFS                           0                 0.03                 0.24
    ## NWHI,_Kure                          0                 0.12                 0.50
    ## NWHI,_Laysan                        0                 0.22                 0.22
    ## NWHI,_Midway                        0                 0.15                 0.54
    ##                  scaffold_116_386790 scaffold_11_3057945 scaffold_121_1829963
    ## Japan,_Torishima                0.05                   0                    0
    ## NWHI,_FFS                       0.24                   0                    0
    ## NWHI,_Kure                      0.12                   0                    0
    ## NWHI,_Laysan                    0.11                   0                    0
    ## NWHI,_Midway                    0.15                   0                    0
    ##                  scaffold_121_565678 scaffold_1226_115783 scaffold_123_1115977
    ## Japan,_Torishima                0.47                 0.09                 0.72
    ## NWHI,_FFS                       0.47                 0.00                 0.74
    ## NWHI,_Kure                      0.54                 0.04                 0.71
    ## NWHI,_Laysan                    0.33                 0.11                 0.67
    ## NWHI,_Midway                    0.31                 0.85                 0.46
    ##                  scaffold_1273_1937 scaffold_127_1105814 scaffold_127_901882
    ## Japan,_Torishima               0.09                 0.51                   0
    ## NWHI,_FFS                      0.06                 0.50                   0
    ## NWHI,_Kure                     0.00                 0.58                   0
    ## NWHI,_Laysan                   0.00                 0.33                   0
    ## NWHI,_Midway                   0.00                 0.54                   0
    ##                  scaffold_12_5002902 scaffold_1327_103421 scaffold_140_259881
    ## Japan,_Torishima                0.47                 0.09                0.30
    ## NWHI,_FFS                       0.21                 0.32                0.38
    ## NWHI,_Kure                      0.46                 0.29                0.33
    ## NWHI,_Laysan                    0.33                 0.44                0.22
    ## NWHI,_Midway                    0.31                 0.23                0.38
    ##                  scaffold_143_1756997 scaffold_146_16700 scaffold_148_727111
    ## Japan,_Torishima                 0.00               0.35                0.51
    ## NWHI,_FFS                        0.09               0.44                0.47
    ## NWHI,_Kure                       0.08               0.50                0.50
    ## NWHI,_Laysan                     0.22               0.44                0.22
    ## NWHI,_Midway                     0.00               0.46                0.46
    ##                  scaffold_14_2621198 scaffold_14_5075430 scaffold_155_1707144
    ## Japan,_Torishima                0.21                0.40                 0.40
    ## NWHI,_FFS                       0.26                0.47                 0.26
    ## NWHI,_Kure                      0.42                0.54                 0.21
    ## NWHI,_Laysan                    0.44                0.44                 0.22
    ## NWHI,_Midway                    0.46                0.38                 0.31
    ##                  scaffold_155_468225 scaffold_157_815403 scaffold_1589_18141
    ## Japan,_Torishima                0.28                0.28                0.53
    ## NWHI,_FFS                       0.29                0.24                0.44
    ## NWHI,_Kure                      0.25                0.38                0.38
    ## NWHI,_Laysan                    0.22                0.33                0.56
    ## NWHI,_Midway                    0.31                0.23                0.46
    ##                  scaffold_15_2923659 scaffold_15_5063205 scaffold_164_1134575
    ## Japan,_Torishima                0.30                0.23                 0.14
    ## NWHI,_FFS                       0.21                0.26                 0.06
    ## NWHI,_Kure                      0.62                0.42                 0.21
    ## NWHI,_Laysan                    0.56                0.56                 0.22
    ## NWHI,_Midway                    0.23                0.54                 0.46
    ##                  scaffold_166_410622 scaffold_16_1050955 scaffold_16_1254862
    ## Japan,_Torishima                0.02                0.19                0.37
    ## NWHI,_FFS                       0.12                0.26                0.32
    ## NWHI,_Kure                      0.12                0.29                0.38
    ## NWHI,_Laysan                    0.33                0.00                0.44
    ## NWHI,_Midway                    0.08                0.31                0.23
    ##                  scaffold_16_31673 scaffold_16_32172 scaffold_177_673090
    ## Japan,_Torishima              0.44              0.44                0.12
    ## NWHI,_FFS                     0.47              0.35                0.44
    ## NWHI,_Kure                    0.42              0.42                0.33
    ## NWHI,_Laysan                  0.44              0.56                0.22
    ## NWHI,_Midway                  0.62              0.54                0.46
    ##                  scaffold_184_724429 scaffold_184_734991 scaffold_185_1133045
    ## Japan,_Torishima                0.53                0.72                 0.35
    ## NWHI,_FFS                       0.56                0.59                 0.21
    ## NWHI,_Kure                      0.67                0.67                 0.46
    ## NWHI,_Laysan                    0.44                0.44                 0.33
    ## NWHI,_Midway                    0.62                0.62                 0.46
    ##                  scaffold_185_1154507 scaffold_190_1605668 scaffold_199_875998
    ## Japan,_Torishima                 0.35                 0.14                0.58
    ## NWHI,_FFS                        0.26                 0.21                0.44
    ## NWHI,_Kure                       0.46                 0.33                0.38
    ## NWHI,_Laysan                     0.33                 0.33                0.67
    ## NWHI,_Midway                     0.62                 0.15                0.54
    ##                  scaffold_1_5556176 scaffold_204_1432955 scaffold_204_239685
    ## Japan,_Torishima               0.30                 0.19                0.23
    ## NWHI,_FFS                      0.53                 0.26                0.44
    ## NWHI,_Kure                     0.42                 0.12                0.33
    ## NWHI,_Laysan                   0.56                 0.22                0.22
    ## NWHI,_Midway                   0.31                 0.23                0.54
    ##                  scaffold_209_721065 scaffold_20_1133858 scaffold_210_1478805
    ## Japan,_Torishima                0.28                0.47                 0.09
    ## NWHI,_FFS                       0.47                0.38                 0.41
    ## NWHI,_Kure                      0.33                0.38                 0.25
    ## NWHI,_Laysan                    0.22                0.11                 0.11
    ## NWHI,_Midway                    0.23                0.31                 0.23
    ##                  scaffold_214_606303 scaffold_223_94277 scaffold_224_319624
    ## Japan,_Torishima                0.07                  0                0.19
    ## NWHI,_FFS                       0.26                  0                0.38
    ## NWHI,_Kure                      0.21                  0                0.21
    ## NWHI,_Laysan                    0.33                  0                0.44
    ## NWHI,_Midway                    0.46                  0                0.15
    ##                  scaffold_227_1219626 scaffold_229_770334 scaffold_234_1146621
    ## Japan,_Torishima                 0.21                0.70                 0.00
    ## NWHI,_FFS                        0.12                0.62                 0.18
    ## NWHI,_Kure                       0.08                0.54                 0.12
    ## NWHI,_Laysan                     0.00                0.56                 0.11
    ## NWHI,_Midway                     0.00                0.62                 0.23
    ##                  scaffold_237_341989 scaffold_238_668548 scaffold_245_674441
    ## Japan,_Torishima                0.19                0.58                   0
    ## NWHI,_FFS                       0.32                0.26                   0
    ## NWHI,_Kure                      0.21                0.50                   0
    ## NWHI,_Laysan                    0.22                0.33                   0
    ## NWHI,_Midway                    0.15                0.38                   0
    ##                  scaffold_246_1244055 scaffold_246_31712 scaffold_247_950885
    ## Japan,_Torishima                 0.40               0.21                0.58
    ## NWHI,_FFS                        0.32               0.71                0.53
    ## NWHI,_Kure                       0.50               0.38                0.54
    ## NWHI,_Laysan                     0.33               0.56                0.44
    ## NWHI,_Midway                     0.69               0.54                0.46
    ##                  scaffold_249_340732 scaffold_24_1490422 scaffold_253_684619
    ## Japan,_Torishima                0.09                0.30                0.23
    ## NWHI,_FFS                       0.21                0.12                0.44
    ## NWHI,_Kure                      0.12                0.29                0.33
    ## NWHI,_Laysan                    0.11                0.22                0.22
    ## NWHI,_Midway                    0.15                0.31                0.62
    ##                  scaffold_256_454799 scaffold_275_1006390 scaffold_27_1646617
    ## Japan,_Torishima                0.40                 0.51                0.65
    ## NWHI,_FFS                       0.47                 0.44                0.50
    ## NWHI,_Kure                      0.25                 0.42                0.46
    ## NWHI,_Laysan                    0.78                 0.56                0.56
    ## NWHI,_Midway                    0.62                 0.54                0.69
    ##                  scaffold_284_808709 scaffold_286_314167 scaffold_287_423519
    ## Japan,_Torishima                0.21                0.26                0.28
    ## NWHI,_FFS                       0.15                0.53                0.50
    ## NWHI,_Kure                      0.12                0.46                0.25
    ## NWHI,_Laysan                    0.11                0.44                0.44
    ## NWHI,_Midway                    0.31                0.38                0.46
    ##                  scaffold_28_4260845 scaffold_28_4262705 scaffold_298_359540
    ## Japan,_Torishima                0.44                0.51                0.26
    ## NWHI,_FFS                       0.32                0.29                0.26
    ## NWHI,_Kure                      0.50                0.54                0.12
    ## NWHI,_Laysan                    0.44                0.44                0.22
    ## NWHI,_Midway                    0.54                0.46                0.08
    ##                  scaffold_298_460712 scaffold_29_3347800 scaffold_306_328368
    ## Japan,_Torishima                0.05                0.44                0.63
    ## NWHI,_FFS                       0.24                0.53                0.50
    ## NWHI,_Kure                      0.17                0.50                0.67
    ## NWHI,_Laysan                    0.00                0.22                0.22
    ## NWHI,_Midway                    0.15                0.23                0.38
    ##                  scaffold_306_944429 scaffold_310_880276 scaffold_316_437138
    ## Japan,_Torishima                0.30                0.19                0.67
    ## NWHI,_FFS                       0.38                0.24                0.62
    ## NWHI,_Kure                      0.29                0.00                0.33
    ## NWHI,_Laysan                    0.22                0.11                0.56
    ## NWHI,_Midway                    0.54                0.08                0.69
    ##                  scaffold_31_750839 scaffold_324_135439 scaffold_32_3673235
    ## Japan,_Torishima               0.47                0.63                0.35
    ## NWHI,_FFS                      0.44                0.44                0.38
    ## NWHI,_Kure                     0.33                0.42                0.38
    ## NWHI,_Laysan                   0.33                0.00                0.33
    ## NWHI,_Midway                   0.38                0.46                0.31
    ##                  scaffold_32_803437 scaffold_32_811415 scaffold_335_662765
    ## Japan,_Torishima               0.05               0.05                0.23
    ## NWHI,_FFS                      0.47               0.41                0.26
    ## NWHI,_Kure                     0.42               0.42                0.17
    ## NWHI,_Laysan                   0.22               0.33                0.00
    ## NWHI,_Midway                   0.46               0.23                0.08
    ##                  scaffold_336_143440 scaffold_33_2797456 scaffold_33_2811656
    ## Japan,_Torishima                0.42                0.40                0.28
    ## NWHI,_FFS                       0.29                0.29                0.41
    ## NWHI,_Kure                      0.29                0.67                0.21
    ## NWHI,_Laysan                    0.56                0.33                0.44
    ## NWHI,_Midway                    0.08                0.46                0.62
    ##                  scaffold_342_417561 scaffold_343_778064 scaffold_343_863766
    ## Japan,_Torishima                0.09                0.30                0.00
    ## NWHI,_FFS                       0.35                0.41                0.00
    ## NWHI,_Kure                      0.12                0.33                0.00
    ## NWHI,_Laysan                    0.11                0.56                0.11
    ## NWHI,_Midway                    0.23                0.15                0.00
    ##                  scaffold_345_694649 scaffold_34_2385714 scaffold_351_703875
    ## Japan,_Torishima                0.16                0.42                0.40
    ## NWHI,_FFS                       0.38                0.35                0.53
    ## NWHI,_Kure                      0.42                0.33                0.54
    ## NWHI,_Laysan                    0.22                0.44                0.44
    ## NWHI,_Midway                    0.54                0.23                0.38
    ##                  scaffold_356_86112 scaffold_360_269450 scaffold_369_134745
    ## Japan,_Torishima               0.33                   0                0.23
    ## NWHI,_FFS                      0.18                   0                0.21
    ## NWHI,_Kure                     0.29                   0                0.12
    ## NWHI,_Laysan                   0.22                   0                0.44
    ## NWHI,_Midway                   0.00                   0                0.15
    ##                  scaffold_381_698673 scaffold_397_13764 scaffold_40_612943
    ## Japan,_Torishima                0.35               0.28               0.63
    ## NWHI,_FFS                       0.18               0.38               0.50
    ## NWHI,_Kure                      0.25               0.46               0.67
    ## NWHI,_Laysan                    0.00               0.33               0.56
    ## NWHI,_Midway                    0.15               0.62               0.69
    ##                  scaffold_411_530837 scaffold_417_428513 scaffold_41_907611
    ## Japan,_Torishima                0.51                0.00               0.26
    ## NWHI,_FFS                       0.44                0.32               0.47
    ## NWHI,_Kure                      0.29                0.17               0.62
    ## NWHI,_Laysan                    0.67                0.44               0.44
    ## NWHI,_Midway                    0.38                0.31               0.38
    ##                  scaffold_429_778790 scaffold_437_192045 scaffold_45_3079339
    ## Japan,_Torishima                0.16                0.07                0.40
    ## NWHI,_FFS                       0.44                0.35                0.71
    ## NWHI,_Kure                      0.42                0.54                0.62
    ## NWHI,_Laysan                    0.44                0.22                0.44
    ## NWHI,_Midway                    0.23                0.38                0.38
    ##                  scaffold_460_237213 scaffold_461_674999 scaffold_472_10182
    ## Japan,_Torishima                0.37                0.58               0.51
    ## NWHI,_FFS                       0.32                0.50               0.82
    ## NWHI,_Kure                      0.54                0.62               0.62
    ## NWHI,_Laysan                    0.33                0.67               0.56
    ## NWHI,_Midway                    0.23                0.31               0.54
    ##                  scaffold_47_2495309 scaffold_487_133339 scaffold_48_3128896
    ## Japan,_Torishima                0.47                0.70                0.47
    ## NWHI,_FFS                       0.26                0.50                0.62
    ## NWHI,_Kure                      0.17                0.38                0.54
    ## NWHI,_Laysan                    0.11                0.44                0.56
    ## NWHI,_Midway                    0.46                0.38                0.69
    ##                  scaffold_491_362387 scaffold_491_382261 scaffold_4_4005861
    ## Japan,_Torishima                0.40                0.00               0.44
    ## NWHI,_FFS                       0.59                0.18               0.24
    ## NWHI,_Kure                      0.46                0.08               0.17
    ## NWHI,_Laysan                    0.44                0.00               0.22
    ## NWHI,_Midway                    0.46                0.31               0.23
    ##                  scaffold_500_214884 scaffold_519_190199 scaffold_51_1220763
    ## Japan,_Torishima                0.16                0.44                0.28
    ## NWHI,_FFS                       0.38                0.29                0.12
    ## NWHI,_Kure                      0.29                0.54                0.17
    ## NWHI,_Laysan                    0.78                0.56                0.78
    ## NWHI,_Midway                    0.62                0.38                0.31
    ##                  scaffold_526_334093 scaffold_52_723164 scaffold_532_590089
    ## Japan,_Torishima                0.21               0.37                0.44
    ## NWHI,_FFS                       0.53               0.44                0.44
    ## NWHI,_Kure                      0.38               0.25                0.50
    ## NWHI,_Laysan                    0.56               0.44                0.33
    ## NWHI,_Midway                    0.31               0.31                0.38
    ##                  scaffold_547_47246 scaffold_54_1630779 scaffold_552_154281
    ## Japan,_Torishima               0.47                0.16                0.42
    ## NWHI,_FFS                      0.47                0.24                0.53
    ## NWHI,_Kure                     0.38                0.42                0.33
    ## NWHI,_Laysan                   0.33                0.33                0.11
    ## NWHI,_Midway                   0.46                0.31                0.38
    ##                  scaffold_555_302525 scaffold_557_489027 scaffold_55_2977720
    ## Japan,_Torishima                0.40                0.33                0.42
    ## NWHI,_FFS                       0.41                0.29                0.29
    ## NWHI,_Kure                      0.42                0.38                0.58
    ## NWHI,_Laysan                    0.22                0.00                0.44
    ## NWHI,_Midway                    0.08                0.15                0.23
    ##                  scaffold_565_253439 scaffold_56_1290372 scaffold_572_19499
    ## Japan,_Torishima                0.47                0.23                  0
    ## NWHI,_FFS                       0.41                0.68                  0
    ## NWHI,_Kure                      0.58                0.33                  0
    ## NWHI,_Laysan                    0.44                0.78                  0
    ## NWHI,_Midway                    0.31                0.31                  0
    ##                  scaffold_57_1788671 scaffold_582_107987 scaffold_582_222480
    ## Japan,_Torishima                0.26                0.09                0.07
    ## NWHI,_FFS                       0.26                0.35                0.32
    ## NWHI,_Kure                      0.38                0.46                0.29
    ## NWHI,_Laysan                    0.44                0.11                0.78
    ## NWHI,_Midway                    0.31                0.38                0.77
    ##                  scaffold_5_697809 scaffold_605_114686 scaffold_608_211809
    ## Japan,_Torishima              0.07                   0                   1
    ## NWHI,_FFS                     0.35                   0                   1
    ## NWHI,_Kure                    0.29                   0                   1
    ## NWHI,_Laysan                  0.67                   0                   1
    ## NWHI,_Midway                  0.15                   0                   1
    ##                  scaffold_60_341016 scaffold_612_363793 scaffold_621_290581
    ## Japan,_Torishima               0.42                0.30                0.33
    ## NWHI,_FFS                      0.53                0.44                0.44
    ## NWHI,_Kure                     0.42                0.21                0.33
    ## NWHI,_Laysan                   0.33                0.11                0.44
    ## NWHI,_Midway                   0.62                0.23                0.38
    ##                  scaffold_62_2806526 scaffold_633_500454 scaffold_64_2598599
    ## Japan,_Torishima                0.40                0.72                0.30
    ## NWHI,_FFS                       0.62                0.65                0.38
    ## NWHI,_Kure                      0.25                0.58                0.50
    ## NWHI,_Laysan                    0.33                0.33                0.33
    ## NWHI,_Midway                    0.46                0.69                0.62
    ##                  scaffold_65_986791 scaffold_670_51777 scaffold_67_2416699
    ## Japan,_Torishima               0.16               0.07                0.49
    ## NWHI,_FFS                      0.56               0.47                0.50
    ## NWHI,_Kure                     0.29               0.29                0.54
    ## NWHI,_Laysan                   0.56               0.67                0.33
    ## NWHI,_Midway                   0.38               0.46                0.69
    ##                  scaffold_684_229342 scaffold_691_412074 scaffold_694_285663
    ## Japan,_Torishima                0.12                0.56                0.33
    ## NWHI,_FFS                       0.26                0.53                0.29
    ## NWHI,_Kure                      0.42                0.46                0.33
    ## NWHI,_Laysan                    0.11                0.44                0.78
    ## NWHI,_Midway                    0.46                0.15                0.46
    ##                  scaffold_698_186739 scaffold_700_166185 scaffold_71_2209803
    ## Japan,_Torishima                0.21                   1                0.14
    ## NWHI,_FFS                       0.24                   1                0.18
    ## NWHI,_Kure                      0.25                   1                0.25
    ## NWHI,_Laysan                    0.22                   1                0.00
    ## NWHI,_Midway                    0.54                   1                0.31
    ##                  scaffold_72_1463271 scaffold_738_79656 scaffold_756_67966
    ## Japan,_Torishima                0.23                  0               0.37
    ## NWHI,_FFS                       0.12                  0               0.62
    ## NWHI,_Kure                      0.25                  0               0.50
    ## NWHI,_Laysan                    0.11                  0               0.11
    ## NWHI,_Midway                    0.15                  0               0.38
    ##                  scaffold_75_1184494 scaffold_760_180618 scaffold_768_353388
    ## Japan,_Torishima                   1                0.12                0.30
    ## NWHI,_FFS                          1                0.09                0.44
    ## NWHI,_Kure                         1                0.46                0.33
    ## NWHI,_Laysan                       1                0.11                0.11
    ## NWHI,_Midway                       1                0.15                0.38
    ##                  scaffold_76_2683423 scaffold_774_111773 scaffold_786_264628
    ## Japan,_Torishima                0.33                0.42                0.12
    ## NWHI,_FFS                       0.41                0.44                0.59
    ## NWHI,_Kure                      0.21                0.33                0.42
    ## NWHI,_Laysan                    0.67                0.44                0.33
    ## NWHI,_Midway                    0.31                0.31                0.62
    ##                  scaffold_78_1735985 scaffold_793_76520 scaffold_7_15896
    ## Japan,_Torishima                0.19               0.56             0.77
    ## NWHI,_FFS                       0.18               0.56             0.97
    ## NWHI,_Kure                      0.21               0.29             0.92
    ## NWHI,_Laysan                    0.33               0.67             1.00
    ## NWHI,_Midway                    0.15               0.46             0.62
    ##                  scaffold_7_4109482 scaffold_820_286874 scaffold_822_90287
    ## Japan,_Torishima               0.63                0.47               0.21
    ## NWHI,_FFS                      0.38                0.65               0.26
    ## NWHI,_Kure                     0.38                0.58               0.29
    ## NWHI,_Laysan                   0.11                0.44               0.67
    ## NWHI,_Midway                   0.15                0.46               0.23
    ##                  scaffold_82_879210 scaffold_834_252344 scaffold_84_2655661
    ## Japan,_Torishima               0.63                0.19                0.44
    ## NWHI,_FFS                      0.38                0.18                0.44
    ## NWHI,_Kure                     0.25                0.21                0.42
    ## NWHI,_Laysan                   0.44                0.11                0.67
    ## NWHI,_Midway                   0.46                0.15                0.23
    ##                  scaffold_854_86476 scaffold_857_326525 scaffold_869_275845
    ## Japan,_Torishima               0.00                0.16                0.60
    ## NWHI,_FFS                      0.24                0.53                0.56
    ## NWHI,_Kure                     0.25                0.58                0.46
    ## NWHI,_Laysan                   0.22                0.56                0.22
    ## NWHI,_Midway                   0.23                0.46                0.46
    ##                  scaffold_88_684287 scaffold_8_4227715 scaffold_91_2132606
    ## Japan,_Torishima               0.44               0.37                0.40
    ## NWHI,_FFS                      0.44               0.15                0.21
    ## NWHI,_Kure                     0.50               0.17                0.12
    ## NWHI,_Laysan                   0.33               0.11                0.22
    ## NWHI,_Midway                   0.31               0.15                0.23
    ##                  scaffold_91_555426 scaffold_925_195188 scaffold_94_1302246
    ## Japan,_Torishima                  0                0.53                0.30
    ## NWHI,_FFS                         0                0.76                0.21
    ## NWHI,_Kure                        0                0.92                0.42
    ## NWHI,_Laysan                      0                0.67                0.00
    ## NWHI,_Midway                      0                0.85                0.46
    ##                  scaffold_958_223818 scaffold_95_2059721 scaffold_97_1719991
    ## Japan,_Torishima                0.44                0.14                0.49
    ## NWHI,_FFS                       0.41                0.38                0.50
    ## NWHI,_Kure                      0.46                0.42                0.46
    ## NWHI,_Laysan                    0.56                0.33                0.67
    ## NWHI,_Midway                    0.31                0.31                0.38
    ##                  scaffold_984_180851 scaffold_990_59711 scaffold_9_4398937
    ## Japan,_Torishima                0.26                  0               0.28
    ## NWHI,_FFS                       0.15                  0               0.56
    ## NWHI,_Kure                      0.21                  0               0.42
    ## NWHI,_Laysan                    0.56                  0               0.33
    ## NWHI,_Midway                    0.31                  0               0.31

``` r
  # make dataframes to compare  
Hob_df <- rownames_to_column(as.data.frame(Hobs), var = "population") %>%
  pivot_longer(cols = 2:187, names_to = "locus", values_to = "H_obs")

Hexp_df <- rownames_to_column(as.data.frame(Hexp), var = "population") %>%
  pivot_longer(cols = 2:187, names_to = "locus", values_to = "H_exp")
  
Hob_df %>%
  left_join(., Hexp_df) %>%
  ggplot(aes(x = H_exp, y = H_obs, color = population)) +
  geom_point() +
  geom_abline(slope = 1)
```

    ## Joining with `by = join_by(population, locus)`

![](02-locus-HWE-evaluation_files/figure-gfm/unnamed-chunk-16-1.png)<!-- -->
Per population is pretty messy… probably because of small sample sizes.

``` r
Hob_df %>%
  left_join(., Hexp_df) %>%
  filter(population == "NWHI,_FFS") %>%
  ggplot(aes(x = H_exp, y = H_obs, color = population)) +
  geom_point() +
  geom_abline(slope = 1)
```

    ## Joining with `by = join_by(population, locus)`

![](02-locus-HWE-evaluation_files/figure-gfm/unnamed-chunk-17-1.png)<!-- -->
Again, kind of a mess, but basically only over 0.4, which is
interesting.

``` r
par(mar=c(5.5, 4.5, 1, 1))
  Hobs.pop <- apply(Hobs, MARGIN = 1, FUN = mean)
  Hexp.pop <- apply(Hexp, 1, mean) 
  barplot(Hexp.pop, ylim=c(0,1), las=3, ylab="Expected heterozygosity")
```

![](02-locus-HWE-evaluation_files/figure-gfm/unnamed-chunk-18-1.png)<!-- -->

``` r
  barplot(Hobs.pop, ylim=c(0,1), las=3, ylab="Observed heterozygosity")
```

![](02-locus-HWE-evaluation_files/figure-gfm/unnamed-chunk-18-2.png)<!-- -->

``` r
bfal_diversity <- data.frame(Pop = names(Hobs.pop),
                              n_samples = sum_output$n.by.pop,
                              Hobs = Hobs.pop,
                              Hexp = Hexp.pop)
                              #Ar = Richness$mean.richness)
as.tibble(bfal_diversity) %>%
  rename(Population = Pop, H_obs = Hobs, H_exp = Hexp) #%>%
```

    ## Warning: `as.tibble()` was deprecated in tibble 2.0.0.
    ## ℹ Please use `as_tibble()` instead.
    ## ℹ The signature and semantics have changed, see `?as_tibble`.
    ## This warning is displayed once every 8 hours.
    ## Call `lifecycle::last_lifecycle_warnings()` to see where this warning was
    ## generated.

    ## # A tibble: 5 × 4
    ##   Population       n_samples H_obs H_exp
    ##   <chr>                <dbl> <dbl> <dbl>
    ## 1 Japan,_Torishima        43 0.315 0.316
    ## 2 NWHI,_FFS               34 0.364 0.362
    ## 3 NWHI,_Kure              24 0.358 0.353
    ## 4 NWHI,_Laysan             9 0.349 0.373
    ## 5 NWHI,_Midway            13 0.356 0.344

``` r
  #write_csv("csv_outputs/reference_pop_hz_summary.csv")
```

chi^2: value of the classical chi-squared test statistic df: degrees of
freedom of the chi-squared test Pr(chi^2 \>): p-value of the chi-squared
test (‘\>’ indicates that the alternative is ‘greater,’ which is always
the case for a chi-squared test) Pr.exact: p-value from an exact test
based on Monte Carlo permutation of alleles (for diploids only). The
default is B = 1000 permutations (set B = 0 to skip this test). Here we
use the function ‘round’ with argument ‘digits = 3’ to round all values
to 3 decimals.

<https://rdrr.io/cran/pegas/man/hw.test.html>

In this case, the degrees of freedom depend on the number alleles for a
given locus.

``` r
# HWE test for all loci - but really what we want is all loci in each population (farther below)
hwe_output <- round(pegas::hw.test(genind_df, B = 1000), digits = 3)
```

    ## Registered S3 method overwritten by 'pegas':
    ##   method      from
    ##   print.amova ade4

    ## Warning in hw.test.loci(x = x, B = B, ...): The following loci were ignored: scaffold_1164_119437
    ## (not the same ploidy for all individuals, or too many missing data)The following loci were ignored: scaffold_1226_115783
    ## (not the same ploidy for all individuals, or too many missing data)The following loci were ignored: scaffold_143_1756997
    ## (not the same ploidy for all individuals, or too many missing data)The following loci were ignored: scaffold_148_727111
    ## (not the same ploidy for all individuals, or too many missing data)The following loci were ignored: scaffold_14_2621198
    ## (not the same ploidy for all individuals, or too many missing data)The following loci were ignored: scaffold_166_410622
    ## (not the same ploidy for all individuals, or too many missing data)The following loci were ignored: scaffold_204_1432955
    ## (not the same ploidy for all individuals, or too many missing data)The following loci were ignored: scaffold_20_1133858
    ## (not the same ploidy for all individuals, or too many missing data)The following loci were ignored: scaffold_227_1219626
    ## (not the same ploidy for all individuals, or too many missing data)The following loci were ignored: scaffold_229_770334
    ## (not the same ploidy for all individuals, or too many missing data)The following loci were ignored: scaffold_238_668548
    ## (not the same ploidy for all individuals, or too many missing data)The following loci were ignored: scaffold_336_143440
    ## (not the same ploidy for all individuals, or too many missing data)The following loci were ignored: scaffold_33_2811656
    ## (not the same ploidy for all individuals, or too many missing data)The following loci were ignored: scaffold_360_269450
    ## (not the same ploidy for all individuals, or too many missing data)The following loci were ignored: scaffold_369_134745
    ## (not the same ploidy for all individuals, or too many missing data)The following loci were ignored: scaffold_411_530837
    ## (not the same ploidy for all individuals, or too many missing data)The following loci were ignored: scaffold_472_10182
    ## (not the same ploidy for all individuals, or too many missing data)The following loci were ignored: scaffold_4_4005861
    ## (not the same ploidy for all individuals, or too many missing data)The following loci were ignored: scaffold_555_302525
    ## (not the same ploidy for all individuals, or too many missing data)The following loci were ignored: scaffold_64_2598599
    ## (not the same ploidy for all individuals, or too many missing data)The following loci were ignored: scaffold_72_1463271
    ## (not the same ploidy for all individuals, or too many missing data)The following loci were ignored: scaffold_738_79656
    ## (not the same ploidy for all individuals, or too many missing data)The following loci were ignored: scaffold_7_15896
    ## (not the same ploidy for all individuals, or too many missing data)The following loci were ignored: scaffold_84_2655661
    ## (not the same ploidy for all individuals, or too many missing data)The following loci were ignored: scaffold_91_2132606
    ## (not the same ploidy for all individuals, or too many missing data)The following loci were ignored: scaffold_94_1302246
    ## (not the same ploidy for all individuals, or too many missing data)The following loci were ignored: scaffold_95_2059721
    ## (not the same ploidy for all individuals, or too many missing data)

    ## Warning in (function (x, ploidy) : Monte Carlo test available only if all
    ## individuals are diploid

    ## Warning in (function (x, ploidy) : Monte Carlo test available only if all
    ## individuals are diploid

    ## Warning in (function (x, ploidy) : Monte Carlo test available only if all
    ## individuals are diploid

    ## Warning in (function (x, ploidy) : Monte Carlo test available only if all
    ## individuals are diploid

    ## Warning in (function (x, ploidy) : Monte Carlo test available only if all
    ## individuals are diploid

    ## Warning in (function (x, ploidy) : Monte Carlo test available only if all
    ## individuals are diploid

    ## Warning in (function (x, ploidy) : Monte Carlo test available only if all
    ## individuals are diploid

    ## Warning in (function (x, ploidy) : Monte Carlo test available only if all
    ## individuals are diploid

    ## Warning in (function (x, ploidy) : Monte Carlo test available only if all
    ## individuals are diploid

    ## Warning in (function (x, ploidy) : Monte Carlo test available only if all
    ## individuals are diploid

    ## Warning in (function (x, ploidy) : Monte Carlo test available only if all
    ## individuals are diploid

    ## Warning in (function (x, ploidy) : Monte Carlo test available only if all
    ## individuals are diploid

    ## Warning in (function (x, ploidy) : Monte Carlo test available only if all
    ## individuals are diploid

    ## Warning in (function (x, ploidy) : Monte Carlo test available only if all
    ## individuals are diploid

    ## Warning in (function (x, ploidy) : Monte Carlo test available only if all
    ## individuals are diploid

    ## Warning in (function (x, ploidy) : Monte Carlo test available only if all
    ## individuals are diploid

    ## Warning in (function (x, ploidy) : Monte Carlo test available only if all
    ## individuals are diploid

    ## Warning in (function (x, ploidy) : Monte Carlo test available only if all
    ## individuals are diploid

    ## Warning in (function (x, ploidy) : Monte Carlo test available only if all
    ## individuals are diploid

    ## Warning in (function (x, ploidy) : Monte Carlo test available only if all
    ## individuals are diploid

    ## Warning in (function (x, ploidy) : Monte Carlo test available only if all
    ## individuals are diploid

    ## Warning in (function (x, ploidy) : Monte Carlo test available only if all
    ## individuals are diploid

    ## Warning in (function (x, ploidy) : Monte Carlo test available only if all
    ## individuals are diploid

    ## Warning in (function (x, ploidy) : Monte Carlo test available only if all
    ## individuals are diploid

    ## Warning in (function (x, ploidy) : Monte Carlo test available only if all
    ## individuals are diploid

    ## Warning in (function (x, ploidy) : Monte Carlo test available only if all
    ## individuals are diploid

    ## Warning in (function (x, ploidy) : Monte Carlo test available only if all
    ## individuals are diploid

``` r
# which ones are statistically sign.
rownames_to_column(data.frame(hwe_output)) %>%
  filter(Pr.exact < 0.01)
```

    ##                 rowname   chi.2 df Pr.chi.2... Pr.exact
    ## 1    scaffold_0_4543438  26.026  6       0.000    0.000
    ## 2  scaffold_116_1932131  10.461  1       0.001    0.001
    ## 3  scaffold_123_1115977  13.937  3       0.003    0.001
    ## 4  scaffold_127_1105814  19.592  6       0.003    0.001
    ## 5  scaffold_1327_103421 138.216  3       0.000    0.000
    ## 6  scaffold_164_1134575  39.857  1       0.000    0.000
    ## 7  scaffold_185_1133045  11.771  1       0.001    0.001
    ## 8  scaffold_185_1154507   7.736  1       0.005    0.006
    ## 9   scaffold_306_944429   9.816  1       0.002    0.002
    ## 10  scaffold_345_694649  16.403  1       0.000    0.000
    ## 11  scaffold_429_778790   7.807  1       0.005    0.006
    ## 12  scaffold_500_214884   7.129  1       0.008    0.002
    ## 13  scaffold_51_1220763 247.303  6       0.000    0.004
    ## 14   scaffold_572_19499 123.000  1       0.000    0.003
    ## 15  scaffold_582_107987  26.155  1       0.000    0.000
    ## 16  scaffold_582_222480  17.309  1       0.000    0.000
    ## 17  scaffold_608_211809 123.000  1       0.000    0.000
    ## 18  scaffold_612_363793 123.008  3       0.000    0.003
    ## 19  scaffold_621_290581 123.034  3       0.000    0.003
    ## 20   scaffold_65_986791  10.461  1       0.001    0.002
    ## 21  scaffold_684_229342 123.019  3       0.000    0.008
    ## 22  scaffold_698_186739  24.802  1       0.000    0.000
    ## 23  scaffold_700_166185 123.000  3       0.000    0.000
    ## 24  scaffold_75_1184494 123.000 10       0.000    0.000
    ## 25  scaffold_760_180618  29.696  1       0.000    0.000
    ## 26  scaffold_78_1735985  45.417  1       0.000    0.000
    ## 27   scaffold_7_4109482 124.577  3       0.000    0.006
    ## 28  scaffold_925_195188  18.875  6       0.004    0.000

28 loci that are out of HWE globally.

After looking at the loci x population HWE info, I can remove any that
are out of HWE in the majority of populations.

``` r
# Chi-squared test: p-value
HWE.test <- data.frame(sapply(seppop(genind_df), 
                              function(ls) pegas::hw.test(ls, B=0)[,3]))
```

    ## Warning in hw.test.loci(x = x, B = B, ...): The following loci were ignored: scaffold_1164_119437
    ## (not the same ploidy for all individuals, or too many missing data)The following loci were ignored: scaffold_1226_115783
    ## (not the same ploidy for all individuals, or too many missing data)The following loci were ignored: scaffold_143_1756997
    ## (not the same ploidy for all individuals, or too many missing data)The following loci were ignored: scaffold_166_410622
    ## (not the same ploidy for all individuals, or too many missing data)The following loci were ignored: scaffold_20_1133858
    ## (not the same ploidy for all individuals, or too many missing data)The following loci were ignored: scaffold_336_143440
    ## (not the same ploidy for all individuals, or too many missing data)The following loci were ignored: scaffold_33_2811656
    ## (not the same ploidy for all individuals, or too many missing data)The following loci were ignored: scaffold_360_269450
    ## (not the same ploidy for all individuals, or too many missing data)The following loci were ignored: scaffold_369_134745
    ## (not the same ploidy for all individuals, or too many missing data)The following loci were ignored: scaffold_411_530837
    ## (not the same ploidy for all individuals, or too many missing data)The following loci were ignored: scaffold_555_302525
    ## (not the same ploidy for all individuals, or too many missing data)The following loci were ignored: scaffold_738_79656
    ## (not the same ploidy for all individuals, or too many missing data)The following loci were ignored: scaffold_84_2655661
    ## (not the same ploidy for all individuals, or too many missing data)The following loci were ignored: scaffold_94_1302246
    ## (not the same ploidy for all individuals, or too many missing data)The following loci were ignored: scaffold_95_2059721
    ## (not the same ploidy for all individuals, or too many missing data)

    ## Warning in hw.test.loci(x = x, B = B, ...): The following loci were ignored: scaffold_1164_119437
    ## (not the same ploidy for all individuals, or too many missing data)The following loci were ignored: scaffold_143_1756997
    ## (not the same ploidy for all individuals, or too many missing data)The following loci were ignored: scaffold_166_410622
    ## (not the same ploidy for all individuals, or too many missing data)The following loci were ignored: scaffold_20_1133858
    ## (not the same ploidy for all individuals, or too many missing data)The following loci were ignored: scaffold_360_269450
    ## (not the same ploidy for all individuals, or too many missing data)The following loci were ignored: scaffold_738_79656
    ## (not the same ploidy for all individuals, or too many missing data)The following loci were ignored: scaffold_84_2655661
    ## (not the same ploidy for all individuals, or too many missing data)The following loci were ignored: scaffold_95_2059721
    ## (not the same ploidy for all individuals, or too many missing data)

    ## Warning in hw.test.loci(x = x, B = B, ...): The following loci were ignored: scaffold_1164_119437
    ## (not the same ploidy for all individuals, or too many missing data)The following loci were ignored: scaffold_143_1756997
    ## (not the same ploidy for all individuals, or too many missing data)The following loci were ignored: scaffold_14_2621198
    ## (not the same ploidy for all individuals, or too many missing data)The following loci were ignored: scaffold_166_410622
    ## (not the same ploidy for all individuals, or too many missing data)The following loci were ignored: scaffold_20_1133858
    ## (not the same ploidy for all individuals, or too many missing data)The following loci were ignored: scaffold_238_668548
    ## (not the same ploidy for all individuals, or too many missing data)The following loci were ignored: scaffold_336_143440
    ## (not the same ploidy for all individuals, or too many missing data)The following loci were ignored: scaffold_360_269450
    ## (not the same ploidy for all individuals, or too many missing data)The following loci were ignored: scaffold_555_302525
    ## (not the same ploidy for all individuals, or too many missing data)The following loci were ignored: scaffold_738_79656
    ## (not the same ploidy for all individuals, or too many missing data)The following loci were ignored: scaffold_84_2655661
    ## (not the same ploidy for all individuals, or too many missing data)The following loci were ignored: scaffold_95_2059721
    ## (not the same ploidy for all individuals, or too many missing data)

    ## Warning in hw.test.loci(x = x, B = B, ...): The following loci were ignored: scaffold_143_1756997
    ## (not the same ploidy for all individuals, or too many missing data)The following loci were ignored: scaffold_148_727111
    ## (not the same ploidy for all individuals, or too many missing data)The following loci were ignored: scaffold_20_1133858
    ## (not the same ploidy for all individuals, or too many missing data)The following loci were ignored: scaffold_227_1219626
    ## (not the same ploidy for all individuals, or too many missing data)The following loci were ignored: scaffold_360_269450
    ## (not the same ploidy for all individuals, or too many missing data)The following loci were ignored: scaffold_472_10182
    ## (not the same ploidy for all individuals, or too many missing data)The following loci were ignored: scaffold_738_79656
    ## (not the same ploidy for all individuals, or too many missing data)The following loci were ignored: scaffold_84_2655661
    ## (not the same ploidy for all individuals, or too many missing data)The following loci were ignored: scaffold_95_2059721
    ## (not the same ploidy for all individuals, or too many missing data)

    ## Warning in hw.test.loci(x = x, B = B, ...): The following loci were ignored: scaffold_1164_119437
    ## (not the same ploidy for all individuals, or too many missing data)The following loci were ignored: scaffold_1226_115783
    ## (not the same ploidy for all individuals, or too many missing data)The following loci were ignored: scaffold_143_1756997
    ## (not the same ploidy for all individuals, or too many missing data)The following loci were ignored: scaffold_166_410622
    ## (not the same ploidy for all individuals, or too many missing data)The following loci were ignored: scaffold_204_1432955
    ## (not the same ploidy for all individuals, or too many missing data)The following loci were ignored: scaffold_20_1133858
    ## (not the same ploidy for all individuals, or too many missing data)The following loci were ignored: scaffold_229_770334
    ## (not the same ploidy for all individuals, or too many missing data)The following loci were ignored: scaffold_238_668548
    ## (not the same ploidy for all individuals, or too many missing data)The following loci were ignored: scaffold_336_143440
    ## (not the same ploidy for all individuals, or too many missing data)The following loci were ignored: scaffold_360_269450
    ## (not the same ploidy for all individuals, or too many missing data)The following loci were ignored: scaffold_369_134745
    ## (not the same ploidy for all individuals, or too many missing data)The following loci were ignored: scaffold_4_4005861
    ## (not the same ploidy for all individuals, or too many missing data)The following loci were ignored: scaffold_555_302525
    ## (not the same ploidy for all individuals, or too many missing data)The following loci were ignored: scaffold_64_2598599
    ## (not the same ploidy for all individuals, or too many missing data)The following loci were ignored: scaffold_72_1463271
    ## (not the same ploidy for all individuals, or too many missing data)The following loci were ignored: scaffold_738_79656
    ## (not the same ploidy for all individuals, or too many missing data)The following loci were ignored: scaffold_7_15896
    ## (not the same ploidy for all individuals, or too many missing data)The following loci were ignored: scaffold_84_2655661
    ## (not the same ploidy for all individuals, or too many missing data)The following loci were ignored: scaffold_91_2132606
    ## (not the same ploidy for all individuals, or too many missing data)The following loci were ignored: scaffold_95_2059721
    ## (not the same ploidy for all individuals, or too many missing data)

``` r
HWE.test.chisq <- t(data.matrix(HWE.test))
{cat("Chi-squared test (p-values):", "\n")
round(HWE.test.chisq,3)}
```

    ## Chi-squared test (p-values):

    ##                  scaffold_0_4543438 scaffold_102_694013 scaffold_1050_36334
    ## Japan._Torishima              0.026               0.501               0.035
    ## NWHI._FFS                     0.432               0.109               0.717
    ## NWHI._Kure                    0.650               0.145               0.432
    ## NWHI._Laysan                  0.853               0.764               0.764
    ## NWHI._Midway                  0.357               0.512               0.512
    ##                  scaffold_1053_177405 scaffold_1078_116737 scaffold_10_1333563
    ## Japan._Torishima                0.781                0.939               0.879
    ## NWHI._FFS                       0.260                0.489               0.040
    ## NWHI._Kure                      0.540                0.689               0.973
    ## NWHI._Laysan                    0.317                0.249               0.764
    ## NWHI._Midway                    0.494                0.764               0.764
    ##                  scaffold_115_2200290 scaffold_1164_119437 scaffold_116_1932131
    ## Japan._Torishima                    1                   NA                0.022
    ## NWHI._FFS                           1                   NA                0.102
    ## NWHI._Kure                          1                   NA                0.744
    ## NWHI._Laysan                        1                0.038                0.099
    ## NWHI._Midway                        1                   NA                0.710
    ##                  scaffold_116_386790 scaffold_11_3057945 scaffold_121_1829963
    ## Japan._Torishima               0.876                   1                    1
    ## NWHI._FFS                      0.437                   1                    1
    ## NWHI._Kure                     0.744                   1                    1
    ## NWHI._Laysan                   0.860                   1                    1
    ## NWHI._Midway                   0.764                   1                    1
    ##                  scaffold_121_565678 scaffold_1226_115783 scaffold_123_1115977
    ## Japan._Torishima               0.771                   NA                0.089
    ## NWHI._FFS                      0.732                1.000                0.104
    ## NWHI._Kure                     0.202                0.000                0.333
    ## NWHI._Laysan                   0.612                0.029                0.172
    ## NWHI._Midway                   0.170                   NA                0.617
    ##                  scaffold_1273_1937 scaffold_127_1105814 scaffold_127_901882
    ## Japan._Torishima              0.992                0.584                   1
    ## NWHI._FFS                     0.860                0.034                   1
    ## NWHI._Kure                    1.000                0.000                   1
    ## NWHI._Laysan                  1.000                0.317                   1
    ## NWHI._Midway                  1.000                0.782                   1
    ##                  scaffold_12_5002902 scaffold_1327_103421 scaffold_140_259881
    ## Japan._Torishima               0.708                0.063               0.243
    ## NWHI._FFS                      0.546                0.106               0.768
    ## NWHI._Kure                     0.432                0.200               0.959
    ## NWHI._Laysan                   0.317                0.019               0.099
    ## NWHI._Midway                   0.631                0.052               0.391
    ##                  scaffold_143_1756997 scaffold_146_16700 scaffold_148_727111
    ## Japan._Torishima                   NA              0.102               0.408
    ## NWHI._FFS                          NA              0.645               0.868
    ## NWHI._Kure                         NA              0.532               0.973
    ## NWHI._Laysan                       NA              0.391                  NA
    ## NWHI._Midway                       NA              0.609               0.279
    ##                  scaffold_14_2621198 scaffold_14_5075430 scaffold_155_1707144
    ## Japan._Torishima               0.443               0.953                0.953
    ## NWHI._FFS                      0.012               0.661                0.180
    ## NWHI._Kure                        NA               0.367                0.569
    ## NWHI._Laysan                   0.764               0.764                0.708
    ## NWHI._Midway                   0.279               0.391                0.317
    ##                  scaffold_155_468225 scaffold_157_815403 scaffold_1589_18141
    ## Japan._Torishima               0.153               0.606               0.705
    ## NWHI._FFS                      0.287               0.437               0.519
    ## NWHI._Kure                     0.236               0.763               0.501
    ## NWHI._Laysan                   0.284               0.370               0.878
    ## NWHI._Midway                   0.317               0.136               0.928
    ##                  scaffold_15_2923659 scaffold_15_5063205 scaffold_164_1134575
    ## Japan._Torishima               0.574               0.051                0.000
    ## NWHI._FFS                      0.503               0.017                0.006
    ## NWHI._Kure                     0.026               0.759                0.121
    ## NWHI._Laysan                   0.022               0.249                0.134
    ## NWHI._Midway                   0.233               0.494                0.797
    ##                  scaffold_166_410622 scaffold_16_1050955 scaffold_16_1254862
    ## Japan._Torishima                  NA               0.011               0.134
    ## NWHI._FFS                         NA               0.374               0.252
    ## NWHI._Kure                        NA               0.403               0.258
    ## NWHI._Laysan                   0.549               1.000               1.000
    ## NWHI._Midway                      NA               0.512               0.136
    ##                  scaffold_16_31673 scaffold_16_32172 scaffold_177_673090
    ## Japan._Torishima             0.453             0.453               0.686
    ## NWHI._FFS                    0.437             0.382               0.933
    ## NWHI._Kure                   0.586             0.586               0.959
    ## NWHI._Laysan                 0.764             0.613               0.284
    ## NWHI._Midway                 0.391             0.782               0.279
    ##                  scaffold_184_724429 scaffold_184_734991 scaffold_185_1133045
    ## Japan._Torishima               0.478               0.426                0.166
    ## NWHI._FFS                      0.407               0.833                0.546
    ## NWHI._Kure                     0.014               0.921                0.838
    ## NWHI._Laysan                   0.764               1.000                0.612
    ## NWHI._Midway                   0.279               0.975                0.279
    ##                  scaffold_185_1154507 scaffold_190_1605668 scaffold_199_875998
    ## Japan._Torishima                0.166                0.623               0.257
    ## NWHI._FFS                       0.889                0.503               0.436
    ## NWHI._Kure                      0.736                0.327               0.290
    ## NWHI._Laysan                    0.612                0.317               0.134
    ## NWHI._Midway                    0.279                0.764               0.710
    ##                  scaffold_1_5556176 scaffold_204_1432955 scaffold_204_239685
    ## Japan._Torishima              0.051                0.501               0.388
    ## NWHI._FFS                     0.588                0.086               0.577
    ## NWHI._Kure                    0.586                0.744               0.327
    ## NWHI._Laysan                  0.739                0.099               0.708
    ## NWHI._Midway                  0.631                   NA               0.494
    ##                  scaffold_209_721065 scaffold_20_1133858 scaffold_210_1478805
    ## Japan._Torishima               0.288                  NA                0.749
    ## NWHI._FFS                      0.223                  NA                0.400
    ## NWHI._Kure                     0.959                  NA                0.624
    ## NWHI._Laysan                   0.708                  NA                0.072
    ## NWHI._Midway                   0.354                  NA                0.638
    ##                  scaffold_214_606303 scaffold_223_94277 scaffold_224_319624
    ## Japan._Torishima               0.813                  1               0.534
    ## NWHI._FFS                      0.350                  1               0.909
    ## NWHI._Kure                     0.955                  1               0.569
    ## NWHI._Laysan                   0.317                  1               0.764
    ## NWHI._Midway                   0.760                  1               0.764
    ##                  scaffold_227_1219626 scaffold_229_770334 scaffold_234_1146621
    ## Japan._Torishima                0.443               0.744                1.000
    ## NWHI._FFS                       0.117               0.898                0.573
    ## NWHI._Kure                      0.831               0.103                0.744
    ## NWHI._Laysan                       NA               0.878                0.860
    ## NWHI._Midway                    1.000                  NA                0.638
    ##                  scaffold_237_341989 scaffold_238_668548 scaffold_245_674441
    ## Japan._Torishima               0.534               0.278                   1
    ## NWHI._FFS                      0.730               0.006                   1
    ## NWHI._Kure                     0.569                  NA                   1
    ## NWHI._Laysan                   0.708               0.549                   1
    ## NWHI._Midway                   0.764                  NA                   1
    ##                  scaffold_246_1244055 scaffold_246_31712 scaffold_247_950885
    ## Japan._Torishima                0.456              0.899               0.278
    ## NWHI._FFS                       0.788              0.616               0.354
    ## NWHI._Kure                      0.102              0.734               0.202
    ## NWHI._Laysan                    0.549              0.318               0.391
    ## NWHI._Midway                    0.056              0.825               0.797
    ##                  scaffold_249_340732 scaffold_24_1490422 scaffold_253_684619
    ## Japan._Torishima               0.749               0.758               0.836
    ## NWHI._FFS                      0.028               0.012               0.647
    ## NWHI._Kure                     0.744               0.403               0.327
    ## NWHI._Laysan                   0.860               0.708               0.134
    ## NWHI._Midway                   0.993               0.512               0.279
    ##                  scaffold_256_454799 scaffold_275_1006390 scaffold_27_1646617
    ## Japan._Torishima               0.282                0.577               0.028
    ## NWHI._FFS                      0.983                0.297               0.007
    ## NWHI._Kure                     0.624                0.532               0.242
    ## NWHI._Laysan                   0.302                0.970               0.515
    ## NWHI._Midway                   0.391                0.940               0.592
    ##                  scaffold_284_808709 scaffold_286_314167 scaffold_287_423519
    ## Japan._Torishima               0.443               0.336               0.153
    ## NWHI._FFS                      0.036               0.732               0.796
    ## NWHI._Kure                     0.744               0.838               0.624
    ## NWHI._Laysan                   0.072               0.391               0.764
    ## NWHI._Midway                   0.512               0.588               0.764
    ##                  scaffold_28_4260845 scaffold_28_4262705 scaffold_298_359540
    ## Japan._Torishima               0.453               0.738               0.336
    ## NWHI._FFS                      0.252               0.089               0.374
    ## NWHI._Kure                     0.973               0.622               0.744
    ## NWHI._Laysan                   0.391               0.391               0.708
    ## NWHI._Midway                   0.184               0.764               0.025
    ##                  scaffold_298_460712 scaffold_29_3347800 scaffold_306_328368
    ## Japan._Torishima               0.876               0.517               0.511
    ## NWHI._FFS                      0.437               0.222               0.718
    ## NWHI._Kure                     0.656               0.102               0.578
    ## NWHI._Laysan                   1.000               0.134               0.284
    ## NWHI._Midway                   0.764               0.077               0.426
    ##                  scaffold_306_944429 scaffold_310_880276 scaffold_316_437138
    ## Japan._Torishima               0.014               0.534               0.560
    ## NWHI._FFS                      0.240               0.717               0.980
    ## NWHI._Kure                     0.042               1.000               0.987
    ## NWHI._Laysan                   0.284               0.860               0.795
    ## NWHI._Midway                   0.184               0.885               0.303
    ##                  scaffold_31_750839 scaffold_324_135439 scaffold_32_3673235
    ## Japan._Torishima              0.976               0.093               0.930
    ## NWHI._FFS                     0.577               0.765               0.168
    ## NWHI._Kure                    0.102               0.431               0.763
    ## NWHI._Laysan                  0.370               0.003               0.612
    ## NWHI._Midway                  0.935               0.279               0.512
    ##                  scaffold_32_803437 scaffold_32_811415 scaffold_335_662765
    ## Japan._Torishima              0.876              0.876               0.388
    ## NWHI._FFS                     0.223              0.400               0.374
    ## NWHI._Kure                    0.197              0.197               0.243
    ## NWHI._Laysan                  0.708              0.948               1.000
    ## NWHI._Midway                  0.279              0.638               0.885
    ##                  scaffold_336_143440 scaffold_33_2797456 scaffold_33_2811656
    ## Japan._Torishima                  NA               0.282                  NA
    ## NWHI._FFS                      0.558               0.558               0.155
    ## NWHI._Kure                        NA               0.102               0.129
    ## NWHI._Laysan                   0.613               0.549               0.261
    ## NWHI._Midway                      NA               0.764               0.391
    ##                  scaffold_342_417561 scaffold_343_778064 scaffold_343_863766
    ## Japan._Torishima               0.063               0.758                1.00
    ## NWHI._FFS                      0.211               0.594                1.00
    ## NWHI._Kure                     0.744               0.811                1.00
    ## NWHI._Laysan                   0.072               0.739                0.86
    ## NWHI._Midway                   0.354               0.041                1.00
    ##                  scaffold_345_694649 scaffold_34_2385714 scaffold_351_703875
    ## Japan._Torishima               0.389               0.606               0.184
    ## NWHI._FFS                      0.768               0.382               0.354
    ## NWHI._Kure                     0.967               0.959               0.676
    ## NWHI._Laysan                   0.099               0.764               0.764
    ## NWHI._Midway                   0.494               0.638               0.935
    ##                  scaffold_356_86112 scaffold_360_269450 scaffold_369_134745
    ## Japan._Torishima              0.623                  NA                  NA
    ## NWHI._FFS                     0.573                  NA               0.546
    ## NWHI._Kure                    0.403                  NA               0.106
    ## NWHI._Laysan                  0.099                  NA               0.764
    ## NWHI._Midway                  1.000                  NA                  NA
    ##                  scaffold_381_698673 scaffold_397_13764 scaffold_40_612943
    ## Japan._Torishima               0.166              0.811              0.050
    ## NWHI._FFS                      0.382              0.594              0.096
    ## NWHI._Kure                     0.053              0.548              0.394
    ## NWHI._Laysan                   1.000              0.478              0.878
    ## NWHI._Midway                   0.764              0.391              0.802
    ##                  scaffold_411_530837 scaffold_417_428513 scaffold_41_907611
    ## Japan._Torishima                  NA               1.000              0.464
    ## NWHI._FFS                      0.567               0.788              0.907
    ## NWHI._Kure                     0.056               0.656              0.073
    ## NWHI._Laysan                   0.134               0.391              1.000
    ## NWHI._Midway                   0.588               0.512              0.588
    ##                  scaffold_429_778790 scaffold_437_192045 scaffold_45_3079339
    ## Japan._Torishima               0.389               0.813               0.590
    ## NWHI._FFS                      0.765               0.911               0.012
    ## NWHI._Kure                     0.431               0.676               0.186
    ## NWHI._Laysan                   0.391               0.708               0.764
    ## NWHI._Midway                   0.354               0.935               0.391
    ##                  scaffold_460_237213 scaffold_461_674999 scaffold_472_10182
    ## Japan._Torishima               0.182               0.286              0.632
    ## NWHI._FFS                      0.040               0.495              0.110
    ## NWHI._Kure                     0.367               0.186              0.065
    ## NWHI._Laysan                   0.370               0.294                 NA
    ## NWHI._Midway                   0.638               0.317              0.197
    ##                  scaffold_47_2495309 scaffold_487_133339 scaffold_48_3128896
    ## Japan._Torishima               0.501               0.007               0.175
    ## NWHI._FFS                      0.374               0.828               0.187
    ## NWHI._Kure                     0.656               0.400               0.194
    ## NWHI._Laysan                   0.860               0.731               0.562
    ## NWHI._Midway                   0.279               0.864               0.559
    ##                  scaffold_491_362387 scaffold_491_382261 scaffold_4_4005861
    ## Japan._Torishima               0.590               1.000              0.438
    ## NWHI._FFS                      0.985               0.573              0.413
    ## NWHI._Kure                     0.685               0.831              0.978
    ## NWHI._Laysan                   0.764               1.000              0.708
    ## NWHI._Midway                   0.928               0.631                 NA
    ##                  scaffold_500_214884 scaffold_519_190199 scaffold_51_1220763
    ## Japan._Torishima               0.561               0.868               0.606
    ## NWHI._FFS                      0.171               0.164               0.716
    ## NWHI._Kure                     0.056               0.625               0.050
    ## NWHI._Laysan                   0.056               0.613               0.002
    ## NWHI._Midway                   0.391               0.444               0.512
    ##                  scaffold_526_334093 scaffold_52_723164 scaffold_532_590089
    ## Japan._Torishima               0.443              0.882               0.785
    ## NWHI._FFS                      0.036              0.765               0.647
    ## NWHI._Kure                     0.377              0.484               0.889
    ## NWHI._Laysan                   0.249              0.391               0.612
    ## NWHI._Midway                   0.631              0.631               0.588
    ##                  scaffold_547_47246 scaffold_54_1630779 scaffold_552_154281
    ## Japan._Torishima              0.698               0.953               0.960
    ## NWHI._FFS                     0.860               0.717               0.588
    ## NWHI._Kure                    0.804               0.586               0.327
    ## NWHI._Laysan                  0.370               0.023               0.860
    ## NWHI._Midway                  0.764               0.512               0.391
    ##                  scaffold_555_302525 scaffold_557_489027 scaffold_55_2977720
    ## Japan._Torishima                  NA               0.623               0.606
    ## NWHI._FFS                      0.634               0.945               0.038
    ## NWHI._Kure                        NA               0.804               0.044
    ## NWHI._Laysan                   0.284               0.003               0.764
    ## NWHI._Midway                      NA               0.764               0.077
    ##                  scaffold_565_253439 scaffold_56_1290372 scaffold_572_19499
    ## Japan._Torishima               0.324               0.388                  1
    ## NWHI._FFS                      0.670               0.039                  0
    ## NWHI._Kure                     0.536               0.157                  1
    ## NWHI._Laysan                   1.000               0.056                  1
    ## NWHI._Midway                   0.481               0.170                  1
    ##                  scaffold_57_1788671 scaffold_582_107987 scaffold_582_222480
    ## Japan._Torishima               0.983               0.749               0.017
    ## NWHI._FFS                      0.401               0.211               0.158
    ## NWHI._Kure                     0.533               0.432               0.116
    ## NWHI._Laysan                   1.000               0.072               0.056
    ## NWHI._Midway                   0.631               0.391               0.048
    ##                  scaffold_5_697809 scaffold_605_114686 scaffold_608_211809
    ## Japan._Torishima             0.813                   1               0.000
    ## NWHI._FFS                    0.586                   1               0.000
    ## NWHI._Kure                   0.393                   1               0.000
    ## NWHI._Laysan                 0.134                   1               0.003
    ## NWHI._Midway                 0.764                   1               0.000
    ##                  scaffold_60_341016 scaffold_612_363793 scaffold_621_290581
    ## Japan._Torishima              0.400               0.744               0.564
    ## NWHI._FFS                     0.688               0.577               0.303
    ## NWHI._Kure                    0.033               0.422               0.327
    ## NWHI._Laysan                  0.612               0.029               0.029
    ## NWHI._Midway                  0.740               0.638               0.935
    ##                  scaffold_62_2806526 scaffold_633_500454 scaffold_64_2598599
    ## Japan._Torishima               0.048               0.348               0.080
    ## NWHI._FFS                      0.467               0.915               0.514
    ## NWHI._Kure                     0.229               0.745               0.303
    ## NWHI._Laysan                   0.382               0.317               0.317
    ## NWHI._Midway                   0.681               0.624                  NA
    ##                  scaffold_65_986791 scaffold_670_51777 scaffold_67_2416699
    ## Japan._Torishima              0.561              0.813               0.106
    ## NWHI._FFS                     0.148              0.661               0.897
    ## NWHI._Kure                    0.116              0.393               0.069
    ## NWHI._Laysan                  0.739              0.294               0.612
    ## NWHI._Midway                  0.444              0.279               0.132
    ##                  scaffold_684_229342 scaffold_691_412074 scaffold_694_285663
    ## Japan._Torishima               0.145               0.272               0.915
    ## NWHI._FFS                      0.401               0.036               0.315
    ## NWHI._Kure                     0.197               0.736               0.327
    ## NWHI._Laysan                   0.029               1.000               0.096
    ## NWHI._Midway                   0.279               0.764               0.279
    ##                  scaffold_698_186739 scaffold_700_166185 scaffold_71_2209803
    ## Japan._Torishima               0.000               0.000               0.623
    ## NWHI._FFS                      0.007               0.000               0.382
    ## NWHI._Kure                     0.102               0.000               0.484
    ## NWHI._Laysan                   0.134               0.029               0.003
    ## NWHI._Midway                   0.782               0.000               0.512
    ##                  scaffold_72_1463271 scaffold_738_79656 scaffold_756_67966
    ## Japan._Torishima               0.836                 NA              0.317
    ## NWHI._FFS                      0.117                 NA              0.091
    ## NWHI._Kure                     0.236                 NA              0.540
    ## NWHI._Laysan                   0.072                 NA              0.022
    ## NWHI._Midway                      NA                 NA              0.391
    ##                  scaffold_75_1184494 scaffold_760_180618 scaffold_768_353388
    ## Japan._Torishima               0.000               0.000               0.243
    ## NWHI._FFS                      0.000               0.002               0.845
    ## NWHI._Kure                     0.000               0.689               0.221
    ## NWHI._Laysan                   0.174               0.860               0.860
    ## NWHI._Midway                   0.005               0.015               0.935
    ##                  scaffold_76_2683423 scaffold_774_111773 scaffold_786_264628
    ## Japan._Torishima               0.564               0.321               0.686
    ## NWHI._FFS                      0.730               0.567               0.152
    ## NWHI._Kure                     0.422               0.102               0.586
    ## NWHI._Laysan                   0.321               0.764               0.370
    ## NWHI._Midway                   0.631               0.170               0.109
    ##                  scaffold_78_1735985 scaffold_793_76520 scaffold_7_15896
    ## Japan._Torishima               0.000              0.245            0.011
    ## NWHI._FFS                      0.000              0.976            0.000
    ## NWHI._Kure                     0.004              0.615            0.007
    ## NWHI._Laysan                   0.370              0.558            0.212
    ## NWHI._Midway                   0.041              0.978               NA
    ##                  scaffold_7_4109482 scaffold_820_286874 scaffold_822_90287
    ## Japan._Torishima              0.080               0.698              0.443
    ## NWHI._FFS                     0.543               0.086              0.889
    ## NWHI._Kure                    0.533               0.414              0.403
    ## NWHI._Laysan                  0.004               0.391              0.294
    ## NWHI._Midway                  0.015               0.928              0.638
    ##                  scaffold_82_879210 scaffold_834_252344 scaffold_84_2655661
    ## Japan._Torishima              0.093               0.501                  NA
    ## NWHI._FFS                     0.514               0.573                  NA
    ## NWHI._Kure                    0.484               0.422                  NA
    ## NWHI._Laysan                  0.391               0.860                  NA
    ## NWHI._Midway                  0.928               0.764                  NA
    ##                  scaffold_854_86476 scaffold_857_326525 scaffold_869_275845
    ## Japan._Torishima              1.000               0.561               0.112
    ## NWHI._FFS                     0.717               0.588               0.489
    ## NWHI._Kure                    0.484               0.414               0.145
    ## NWHI._Laysan                  0.708               0.613               0.284
    ## NWHI._Midway                  0.638               0.928               0.279
    ##                  scaffold_88_684287 scaffold_8_4227715 scaffold_91_2132606
    ## Japan._Torishima              0.639              0.134               0.330
    ## NWHI._FFS                     0.303              0.644               0.503
    ## NWHI._Kure                    0.744              0.243               0.744
    ## NWHI._Laysan                  0.612              0.860               0.099
    ## NWHI._Midway                  0.631              0.764                  NA
    ##                  scaffold_91_555426 scaffold_925_195188 scaffold_94_1302246
    ## Japan._Torishima                  1               0.268                  NA
    ## NWHI._FFS                         1               0.187               0.546
    ## NWHI._Kure                        1               0.025               0.197
    ## NWHI._Laysan                      1               0.558               0.003
    ## NWHI._Midway                      1               0.149               0.797
    ##                  scaffold_958_223818 scaffold_95_2059721 scaffold_97_1719991
    ## Japan._Torishima               0.639                  NA               0.605
    ## NWHI._FFS                      0.303                  NA               0.707
    ## NWHI._Kure                     0.993                  NA               0.672
    ## NWHI._Laysan                   0.613                  NA               0.946
    ## NWHI._Midway                   0.631                  NA               0.766
    ##                  scaffold_984_180851 scaffold_990_59711 scaffold_9_4398937
    ## Japan._Torishima               0.464                  1              0.267
    ## NWHI._FFS                      0.644                  1              0.795
    ## NWHI._Kure                     0.422                  1              0.149
    ## NWHI._Laysan                   0.739                  1              0.261
    ## NWHI._Midway                   0.512                  1              0.143

``` r
# these are the p-values for the per-pop calcs
loci_to_toss_hwe <- rownames_to_column(as.data.frame(HWE.test.chisq), var = "population") %>%
  pivot_longer(cols= 2:length(rownames_to_column(as.data.frame(HWE.test.chisq))), names_to = "locus", values_to = "p_val") %>%
  filter(p_val < 0.05) %>% # which p-values are significant?
  group_by( locus) %>%
  tally() %>% # how many loci are out of HWE in how many populations?
  arrange(desc(n)) %>%
  filter(n>3)
```

Remove loci that are statistically out of HWE in \>3 populations (the
majority).

That’s only 4 loci.

``` r
loci_to_toss_hwe %>%
  write_csv("csv_outputs/loci_to_toss_hwe.csv")

loci_to_keep <- rownames_to_column(as.data.frame(HWE.test.chisq), var = "population") %>%
  pivot_longer(cols= 2:length(rownames_to_column(as.data.frame(HWE.test.chisq))), names_to = "locus", values_to = "p_val") %>%
  select(locus) %>%
  unique() %>%
  anti_join(., loci_to_toss_hwe)
```

    ## Joining with `by = join_by(locus)`

``` r
loci_to_keep %>%
  write_csv("csv_outputs/loci_to_keep_hwe.csv")
# from here, I can generate a "final" dataset.
```

## Additional analyses

25 March 2024

In prep for the manuscript… Looking at LD, primarily Consider
implementation of a Bonferroni correction…

``` r
library(poppr)
```

    ## This is poppr version 2.9.4. To get started, type package?poppr
    ## OMP parallel support: unavailable

    ## 
    ## Attaching package: 'poppr'

    ## The following object is masked from 'package:radiator':
    ## 
    ##     private_alleles

``` r
poppr(genind_df)
```

    ##                Pop   N MLG  eMLG       SE    H     G lambda   E.5  Hexp     Ia
    ## 1 Japan,_Torishima  43  43 10.00 9.46e-07 3.76  43.0  0.977 1.000 0.326 0.0265
    ## 2        NWHI,_FFS  34  34 10.00 1.00e-06 3.53  34.0  0.971 1.000 0.371 1.7239
    ## 3       NWHI,_Kure  24  23  9.84 3.69e-01 3.12  22.2  0.955 0.977 0.365 1.4276
    ## 4     NWHI,_Laysan   9   9  9.00 0.00e+00 2.20   9.0  0.889 1.000 0.400 5.7740
    ## 5     NWHI,_Midway  13  13 10.00 7.30e-08 2.56  13.0  0.923 1.000 0.374 0.1004
    ## 6            Total 123 122  9.99 7.72e-02 4.80 121.0  0.992 0.995 0.380 2.4540
    ##      rbarD      File
    ## 1 0.000164 genind_df
    ## 2 0.010239 genind_df
    ## 3 0.008633 genind_df
    ## 4 0.035777 genind_df
    ## 5 0.000612 genind_df
    ## 6 0.014397 genind_df

``` r
locus_table(genind_df)
```

    ## 
    ## allele = Number of observed alleles

    ## 
    ## 1-D = Simpson index
    ## Hexp = Nei's 1978 gene diversity
    ## ------------------------------------------

    ##                       summary
    ## locus                  allele    1-D   Hexp Evenness
    ##   scaffold_0_4543438   4.0000 0.6151 0.6177   0.8217
    ##   scaffold_102_694013  2.0000 0.3381 0.3394   0.7468
    ##   scaffold_1050_36334  2.0000 0.3091 0.3104   0.7117
    ##   scaffold_1053_177405 2.0000 0.3888 0.3904   0.8141
    ##   scaffold_1078_116737 2.0000 0.4973 0.4994   0.9947
    ##   scaffold_10_1333563  2.0000 0.4997 0.5017   0.9994
    ##   scaffold_115_2200290 1.0000      .      .         
    ##   scaffold_1164_119437 5.0000 0.4750 0.4771   0.6774
    ##   scaffold_116_1932131 2.0000 0.4935 0.4955   0.9872
    ##   scaffold_116_386790  2.0000 0.1216 0.1221   0.5090
    ##   scaffold_11_3057945  1.0000      .      .         
    ##   scaffold_121_1829963 1.0000      .      .         
    ##   scaffold_121_565678  2.0000 0.4960 0.4980   0.9921
    ##   scaffold_1226_115783 7.0000 0.5679 0.5704   0.7697
    ##   scaffold_123_1115977 3.0000 0.6410 0.6436   0.9448
    ##   scaffold_1273_1937   3.0000 0.0479 0.0481   0.3577
    ##   scaffold_127_1105814 4.0000 0.6093 0.6117   0.8147
    ##   scaffold_127_901882  1.0000      .      .         
    ##   scaffold_12_5002902  2.0000 0.4239 0.4256   0.8660
    ##   scaffold_1327_103421 3.0000 0.3769 0.3784   0.7424
    ##   scaffold_140_259881  2.0000 0.3560 0.3575   0.7697
    ##   scaffold_143_1756997 4.0000 0.4277 0.4298   0.6493
    ##   scaffold_146_16700   3.0000 0.4733 0.4752   0.7730
    ##   scaffold_148_727111  3.0000 0.5039 0.5060   0.9696
    ##   scaffold_14_2621198  3.0000 0.4273 0.4290   0.8364
    ##   scaffold_14_5075430  2.0000 0.4301 0.4318   0.8757
    ##   scaffold_155_1707144 2.0000 0.3381 0.3394   0.7468
    ##   scaffold_155_468225  2.0000 0.3604 0.3618   0.7754
    ##   scaffold_157_815403  2.0000 0.3190 0.3203   0.7235
    ##   scaffold_1589_18141  5.0000 0.5091 0.5112   0.9095
    ##   scaffold_15_2923659  4.0000 0.3648 0.3663   0.6631
    ##   scaffold_15_5063205  2.0000 0.4106 0.4123   0.8458
    ##   scaffold_164_1134575 2.0000 0.3964 0.3980   0.8248
    ##   scaffold_166_410622  3.0000 0.2290 0.2300   0.5005
    ##   scaffold_16_1050955  2.0000 0.2499 0.2509   0.6450
    ##   scaffold_16_1254862  2.0000 0.3810 0.3826   0.8032
    ##   scaffold_16_31673    2.0000 0.4810 0.4829   0.9632
    ##   scaffold_16_32172    2.0000 0.4840 0.4860   0.9689
    ##   scaffold_177_673090  2.0000 0.3141 0.3154   0.7176
    ##   scaffold_184_724429  2.0000 0.4973 0.4994   0.9947
    ##   scaffold_184_734991  4.0000 0.6330 0.6356   0.9048
    ##   scaffold_185_1133045 2.0000 0.4944 0.4964   0.9889
    ##   scaffold_185_1154507 2.0000 0.4992 0.5012   0.9984
    ##   scaffold_190_1605668 2.0000 0.2263 0.2272   0.6199
    ##   scaffold_199_875998  3.0000 0.4790 0.4810   0.9227
    ##   scaffold_1_5556176   2.0000 0.4444 0.4463   0.8990
    ##   scaffold_204_1432955 3.0000 0.2854 0.2866   0.6528
    ##   scaffold_204_239685  2.0000 0.3141 0.3154   0.7176
    ##   scaffold_209_721065  2.0000 0.3091 0.3104   0.7117
    ##   scaffold_20_1133858  3.0000 0.5622 0.5646   0.8837
    ##   scaffold_210_1478805 2.0000 0.2382 0.2392   0.6325
    ##   scaffold_214_606303  3.0000 0.2369 0.2379   0.5445
    ##   scaffold_223_94277   1.0000      .      .         
    ##   scaffold_224_319624  2.0000 0.2937 0.2949   0.6938
    ##   scaffold_227_1219626 3.0000 0.1367 0.1373   0.4935
    ##   scaffold_229_770334  4.0000 0.6375 0.6401   0.9197
    ##   scaffold_234_1146621 2.0000 0.1001 0.1005   0.4842
    ##   scaffold_237_341989  2.0000 0.2382 0.2392   0.6325
    ##   scaffold_238_668548  3.0000 0.5005 0.5026   0.9369
    ##   scaffold_245_674441  1.0000      .      .         
    ##   scaffold_246_1244055 3.0000 0.3444 0.3458   0.7224
    ##   scaffold_246_31712   4.0000 0.4227 0.4244   0.6306
    ##   scaffold_247_950885  2.0000 0.4997 0.5017   0.9994
    ##   scaffold_249_340732  3.0000 0.1371 0.1377   0.4636
    ##   scaffold_24_1490422  2.0000 0.2612 0.2623   0.6574
    ##   scaffold_253_684619  2.0000 0.3810 0.3826   0.8032
    ##   scaffold_256_454799  3.0000 0.4691 0.4710   0.9050
    ##   scaffold_275_1006390 4.0000 0.5003 0.5023   0.6214
    ##   scaffold_27_1646617  3.0000 0.6117 0.6142   0.9120
    ##   scaffold_284_808709  2.0000 0.2017 0.2026   0.5941
    ##   scaffold_286_314167  2.0000 0.4741 0.4760   0.9505
    ##   scaffold_287_423519  2.0000 0.4106 0.4123   0.8458
    ##   scaffold_28_4260845  2.0000 0.4777 0.4796   0.9571
    ##   scaffold_28_4262705  2.0000 0.4935 0.4955   0.9872
    ##   scaffold_298_359540  2.0000 0.2017 0.2026   0.5941
    ##   scaffold_298_460712  2.0000 0.1216 0.1221   0.5090
    ##   scaffold_29_3347800  2.0000 0.4682 0.4702   0.9399
    ##   scaffold_306_328368  3.0000 0.5729 0.5753   0.8426
    ##   scaffold_306_944429  2.0000 0.4759 0.4778   0.9538
    ##   scaffold_310_880276  2.0000 0.1629 0.1635   0.5535
    ##   scaffold_316_437138  5.0000 0.6001 0.6025   0.7785
    ##   scaffold_31_750839   2.0000 0.4952 0.4973   0.9906
    ##   scaffold_324_135439  2.0000 0.4793 0.4813   0.9602
    ##   scaffold_32_3673235  2.0000 0.3334 0.3348   0.7410
    ##   scaffold_32_803437   2.0000 0.2612 0.2623   0.6574
    ##   scaffold_32_811415   3.0000 0.2393 0.2403   0.6026
    ##   scaffold_335_662765  2.0000 0.1890 0.1898   0.5809
    ##   scaffold_336_143440  3.0000 0.4432 0.4451   0.7425
    ##   scaffold_33_2797456  2.0000 0.4854 0.4874   0.9716
    ##   scaffold_33_2811656  6.0000 0.4771 0.4791   0.5877
    ##   scaffold_342_417561  2.0000 0.2080 0.2088   0.6006
    ##   scaffold_343_778064  4.0000 0.3648 0.3663   0.6680
    ##   scaffold_343_863766  2.0000 0.0081 0.0081   0.3047
    ##   scaffold_345_694649  2.0000 0.4995 0.5015   0.9989
    ##   scaffold_34_2385714  2.0000 0.4106 0.4123   0.8458
    ##   scaffold_351_703875  2.0000 0.4967 0.4987   0.9934
    ##   scaffold_356_86112   2.0000 0.2556 0.2566   0.6512
    ##   scaffold_360_269450  2.0000 0.1374 0.1381   0.5264
    ##   scaffold_369_134745  3.0000 0.2953 0.2966   0.6276
    ##   scaffold_381_698673  2.0000 0.2668 0.2679   0.6635
    ##   scaffold_397_13764   3.0000 0.4039 0.4055   0.7014
    ##   scaffold_40_612943   3.0000 0.6657 0.6685   0.9979
    ##   scaffold_411_530837  3.0000 0.4966 0.4986   0.9552
    ##   scaffold_417_428513  2.0000 0.1826 0.1833   0.5741
    ##   scaffold_41_907611   3.0000 0.4263 0.4280   0.8350
    ##   scaffold_429_778790  2.0000 0.4239 0.4256   0.8660
    ##   scaffold_437_192045  2.0000 0.3287 0.3300   0.7352
    ##   scaffold_45_3079339  2.0000 0.4662 0.4681   0.9362
    ##   scaffold_460_237213  2.0000 0.4759 0.4778   0.9538
    ##   scaffold_461_674999  2.0000 0.4935 0.4955   0.9872
    ##   scaffold_472_10182   4.0000 0.6216 0.6242   0.8803
    ##   scaffold_47_2495309  2.0000 0.3041 0.3053   0.7058
    ##   scaffold_487_133339  3.0000 0.4781 0.4800   0.8434
    ##   scaffold_48_3128896  5.0000 0.6731 0.6759   0.8634
    ##   scaffold_491_362387  4.0000 0.5250 0.5271   0.8562
    ##   scaffold_491_382261  2.0000 0.1073 0.1078   0.4927
    ##   scaffold_4_4005861   4.0000 0.3174 0.3187   0.5405
    ##   scaffold_500_214884  2.0000 0.4497 0.4516   0.9078
    ##   scaffold_519_190199  4.0000 0.4833 0.4852   0.8150
    ##   scaffold_51_1220763  4.0000 0.2858 0.2870   0.6243
    ##   scaffold_526_334093  2.0000 0.3810 0.3826   0.8032
    ##   scaffold_52_723164   2.0000 0.3810 0.3826   0.8032
    ##   scaffold_532_590089  2.0000 0.4926 0.4946   0.9853
    ##   scaffold_547_47246   2.0000 0.4973 0.4994   0.9947
    ##   scaffold_54_1630779  3.0000 0.2870 0.2882   0.6190
    ##   scaffold_552_154281  2.0000 0.4000 0.4017   0.8301
    ##   scaffold_555_302525  4.0000 0.5260 0.5283   0.6818
    ##   scaffold_557_489027  2.0000 0.2989 0.3001   0.6998
    ##   scaffold_55_2977720  2.0000 0.4572 0.4590   0.9205
    ##   scaffold_565_253439  3.0000 0.4960 0.4980   0.7277
    ##   scaffold_56_1290372  2.0000 0.4548 0.4566   0.9164
    ##   scaffold_572_19499   2.0000 0.0161 0.0162   0.3390
    ##   scaffold_57_1788671  2.0000 0.3646 0.3661   0.7810
    ##   scaffold_582_107987  2.0000 0.4979 0.4999   0.9958
    ##   scaffold_582_222480  2.0000 0.4944 0.4964   0.9889
    ##   scaffold_5_697809    2.0000 0.2723 0.2734   0.6696
    ##   scaffold_605_114686  1.0000      .      .         
    ##   scaffold_608_211809  2.0000 0.5000 0.5020   1.0000
    ##   scaffold_60_341016   4.0000 0.5507 0.5530   0.8472
    ##   scaffold_612_363793  3.0000 0.3121 0.3134   0.6629
    ##   scaffold_621_290581  3.0000 0.3810 0.3825   0.7479
    ##   scaffold_62_2806526  3.0000 0.5435 0.5457   0.8766
    ##   scaffold_633_500454  4.0000 0.6467 0.6493   0.8714
    ##   scaffold_64_2598599  3.0000 0.4174 0.4191   0.8214
    ##   scaffold_65_986791   2.0000 0.4935 0.4955   0.9872
    ##   scaffold_670_51777   2.0000 0.3334 0.3348   0.7410
    ##   scaffold_67_2416699  2.0000 0.4618 0.4637   0.9285
    ##   scaffold_684_229342  3.0000 0.2693 0.2704   0.6152
    ##   scaffold_691_412074  2.0000 0.4417 0.4435   0.8945
    ##   scaffold_694_285663  2.0000 0.3287 0.3300   0.7352
    ##   scaffold_698_186739  2.0000 0.4722 0.4741   0.9471
    ##   scaffold_700_166185  3.0000 0.5040 0.5061   0.9699
    ##   scaffold_71_2209803  2.0000 0.1890 0.1898   0.5809
    ##   scaffold_72_1463271  3.0000 0.2401 0.2411   0.6034
    ##   scaffold_738_79656   2.0000 0.1295 0.1300   0.5177
    ##   scaffold_756_67966   2.0000 0.4572 0.4590   0.9205
    ##   scaffold_75_1184494  5.0000 0.6189 0.6214   0.8603
    ##   scaffold_760_180618  2.0000 0.3516 0.3531   0.7640
    ##   scaffold_768_353388  2.0000 0.3604 0.3618   0.7754
    ##   scaffold_76_2683423  3.0000 0.3790 0.3806   0.7672
    ##   scaffold_774_111773  2.0000 0.4967 0.4987   0.9934
    ##   scaffold_786_264628  2.0000 0.4904 0.4925   0.9812
    ##   scaffold_78_1735985  2.0000 0.4973 0.4994   0.9947
    ##   scaffold_793_76520   4.0000 0.4969 0.4990   0.7000
    ##   scaffold_7_15896     5.0000 0.5892 0.5916   0.8158
    ##   scaffold_7_4109482   3.0000 0.4890 0.4910   0.9157
    ##   scaffold_820_286874  2.0000 0.4825 0.4845   0.9661
    ##   scaffold_822_90287   2.0000 0.2612 0.2623   0.6574
    ##   scaffold_82_879210   2.0000 0.4331 0.4348   0.8805
    ##   scaffold_834_252344  2.0000 0.1761 0.1768   0.5673
    ##   scaffold_84_2655661  4.0000 0.5245 0.5268   0.7313
    ##   scaffold_854_86476   2.0000 0.1562 0.1568   0.5464
    ##   scaffold_857_326525  2.0000 0.4793 0.4813   0.9602
    ##   scaffold_869_275845  2.0000 0.4682 0.4702   0.9399
    ##   scaffold_88_684287   2.0000 0.4106 0.4123   0.8458
    ##   scaffold_8_4227715   2.0000 0.2142 0.2150   0.6071
    ##   scaffold_91_2132606  3.0000 0.2908 0.2920   0.6588
    ##   scaffold_91_555426   1.0000      .      .         
    ##   scaffold_925_195188  4.0000 0.5992 0.6017   0.8291
    ##   scaffold_94_1302246  3.0000 0.3920 0.3936   0.7848
    ##   scaffold_958_223818  2.0000 0.4868 0.4888   0.9742
    ##   scaffold_95_2059721  4.0000 0.3909 0.3926   0.5601
    ##   scaffold_97_1719991  4.0000 0.4913 0.4933   0.6156
    ##   scaffold_984_180851  2.0000 0.2723 0.2734   0.6696
    ##   scaffold_990_59711   1.0000      .      .         
    ##   scaffold_9_4398937   3.0000 0.4486 0.4505   0.6785
    ##   mean                 2.5484 0.3788 0.3804   0.7820

### LD calcs

``` r
JP <- popsub(genind_df, "Japan,_Torishima")

set.seed(335)
ia(JP, sample = 999)
```

![](02-locus-HWE-evaluation_files/figure-gfm/subset-for-just-Torishima-LD-analysis-1.png)<!-- -->

    ##           Ia         p.Ia        rbarD         p.rD 
    ## 0.0265119470 0.4690000000 0.0001637092 0.4670000000

Overall?

``` r
set.seed(335)
ia(genind_df, sample = 999)
```

![](02-locus-HWE-evaluation_files/figure-gfm/check-LD-over-all-pops-1.png)<!-- -->

    ##         Ia       p.Ia      rbarD       p.rD 
    ## 2.45398138 0.00100000 0.01439691 0.00100000

``` r
pairwise_ld_by_locus <- pair.ia(JP, sample = 10, plot = F)

jp_pvals <- rownames_to_column(as.data.frame(pairwise_ld_by_locus)) #%>%
  #filter(rbarD > 0.8)

# test using a
# bonferroni-adjustment
# for multiple tests
as.data.frame(round(p.adjust(jp_pvals$p.rD, "bonferroni"), 3)) %>%
  rename(corrected_pval = `round(p.adjust(jp_pvals$p.rD, "bonferroni"), 3)`) %>%
  filter(corrected_pval < 0.05)
```

    ## [1] corrected_pval
    ## <0 rows> (or 0-length row.names)

LD present in 3 pairs (6 loci) in the Torishima subset. Let’s check the
rest of the dataset.

``` r
#LD_by_locus_all <- pair.ia(genind_df, plot = F, sample = 50)

#ld_for_bonfer <- rownames_to_column(as.data.frame(LD_by_locus_all))

## save that R data frame because it takes a long time to run the pair.ia with the sample call
# ld_for_bonfer %>%
#   write_rds("csv_outputs/LD_for_bonferroni.rds")

ld_for_bonfer <- read_rds("csv_outputs/LD_for_bonferroni.rds")

# bonferroni-adjustment
# for multiple tests
as.data.frame(round(p.adjust(ld_for_bonfer$p.rD, "bonferroni"), 3)) %>%
  rename(corrected_pval = `round(p.adjust(ld_for_bonfer$p.rD, "bonferroni"), 3)`) %>%
  filter(corrected_pval < 0.05)
```

    ## [1] corrected_pval
    ## <0 rows> (or 0-length row.names)

The below analysis of LD in each of the pops is kind of interesting, but
it’s not necessary to remove any of the loci because of the Bonferroni
correction.

``` r
# linked_loci <- rownames_to_column(as.data.frame(LD_by_locus_all)) %>%
#   filter(rbarD > 0.8)
# 
# linked_loci
```

If I need to remove 4 loci because of linkage, can I be strategic about
minimizing missing data? Also, very interesting that the locus on
scaffold 343 is linked with a locus on scaffold 700. \*\*Most likely
unnecessary because of multiple comparisons (Bonferroni correction)

``` r
# locs_to_remove_for_LD <- linked_loci %>%
#   separate(rowname, into = c("loc1", "loc2"), sep = ":") %>%
#   pivot_longer(1:2, names_to = "loc", values_to = "locus") %>%
#   select(locus) %>%
#   left_join(., genos_locs_ind_filtered) %>%
#   group_by(locus) %>%
#   summarise(median(depth)) %>%
#   filter(locus %in% c("scaffold_16_31673", "scaffold_185_1133045", "scaffold_28_4260845", "scaffold_343_863766")) %>%
#   select(locus)
# 
# locs_to_remove_for_LD %>%
#   write_csv("csv_outputs/locs_to_remove_for_LD.csv")
```

Ok, pretty clear which loci generally have higher read depths for each
of those pairs.

Out of curiosity, what about any of the other subpops?

``` r
# What is going on with Laysan??
pop1 <- popsub(genind_df, "NWHI,_Laysan")
pairwise_ld_by_locus <- pair.ia(pop1)
```

![](02-locus-HWE-evaluation_files/figure-gfm/other-pops-LD-1.png)<!-- -->

``` r
rownames_to_column(as.data.frame(pairwise_ld_by_locus)) %>%
  filter(rbarD > 0.8) %>%
  separate(rowname, into = c("loc1", "loc2"), sep = ":")
```

    ##                    loc1                 loc2        Ia     rbarD
    ## 1  scaffold_1053_177405  scaffold_140_259881 0.8022284 0.8060684
    ## 2   scaffold_10_1333563  scaffold_491_362387 1.0000000 1.0000000
    ## 3  scaffold_1164_119437  scaffold_116_386790 0.7027027 0.8534918
    ## 4  scaffold_1226_115783 scaffold_227_1219626 0.7133758 0.8854377
    ## 5  scaffold_1226_115783  scaffold_324_135439 0.8853755 0.8854377
    ## 6  scaffold_1226_115783  scaffold_343_863766 0.7133758 0.8854377
    ## 7  scaffold_1226_115783  scaffold_557_489027 0.8853755 0.8854377
    ## 8  scaffold_1226_115783  scaffold_582_107987 0.9792531 0.9799367
    ## 9  scaffold_1226_115783  scaffold_700_166185 0.7133758 0.8854377
    ## 10 scaffold_1226_115783  scaffold_71_2209803 0.8853755 0.8854377
    ## 11 scaffold_1327_103421 scaffold_227_1219626 0.6896552 0.8058230
    ## 12 scaffold_1327_103421  scaffold_324_135439 0.8000000 0.8058230
    ## 13 scaffold_1327_103421  scaffold_343_863766 0.6896552 0.8058230
    ## 14 scaffold_1327_103421  scaffold_557_489027 0.8000000 0.8058230
    ## 15 scaffold_1327_103421  scaffold_700_166185 0.6896552 0.8058230
    ## 16 scaffold_1327_103421  scaffold_71_2209803 0.8000000 0.8058230
    ## 17 scaffold_155_1707144   scaffold_854_86476 1.0000000 1.0000000
    ## 18  scaffold_15_2923659 scaffold_210_1478805 0.8091895 0.8097254
    ## 19  scaffold_15_2923659  scaffold_342_417561 0.8091895 0.8097254
    ## 20 scaffold_185_1133045 scaffold_185_1154507 1.0000000 1.0000000
    ## 21 scaffold_210_1478805 scaffold_227_1219626 0.7027027 0.8534918
    ## 22 scaffold_210_1478805  scaffold_324_135439 0.8524590 0.8534918
    ## 23 scaffold_210_1478805  scaffold_342_417561 1.0000000 1.0000000
    ## 24 scaffold_210_1478805  scaffold_343_863766 0.7027027 0.8534918
    ## 25 scaffold_210_1478805  scaffold_557_489027 0.8524590 0.8534918
    ## 26 scaffold_210_1478805  scaffold_700_166185 0.7027027 0.8534918
    ## 27 scaffold_210_1478805  scaffold_71_2209803 0.8524590 0.8534918
    ## 28 scaffold_210_1478805  scaffold_75_1184494 0.8198758 0.9135003
    ## 29 scaffold_227_1219626  scaffold_324_135439 0.8000000 1.0000000
    ## 30 scaffold_227_1219626  scaffold_342_417561 0.7027027 0.8534918
    ## 31 scaffold_227_1219626  scaffold_343_863766 1.0000000 1.0000000
    ## 32 scaffold_227_1219626  scaffold_557_489027 0.8000000 1.0000000
    ## 33 scaffold_227_1219626  scaffold_582_107987 0.7027027 0.8534918
    ## 34 scaffold_227_1219626  scaffold_612_363793 0.7133758 0.8854377
    ## 35 scaffold_227_1219626  scaffold_684_229342 0.7133758 0.8854377
    ## 36 scaffold_227_1219626  scaffold_700_166185 1.0000000 1.0000000
    ## 37 scaffold_227_1219626  scaffold_71_2209803 0.8000000 1.0000000
    ## 38  scaffold_237_341989   scaffold_4_4005861 1.0000000 1.0000000
    ## 39  scaffold_249_340732  scaffold_47_2495309 1.0000000 1.0000000
    ## 40  scaffold_27_1646617  scaffold_316_437138 0.8000000 0.8000000
    ## 41  scaffold_284_808709  scaffold_437_192045 0.8000000 0.8958280
    ## 42  scaffold_284_808709  scaffold_72_1463271 1.0000000 1.0000000
    ## 43  scaffold_284_808709   scaffold_8_4227715 0.7027027 0.8534918
    ## 44  scaffold_28_4260845  scaffold_28_4262705 1.0000000 1.0000000
    ## 45  scaffold_306_328368  scaffold_612_363793 0.7989789 0.8044837
    ## 46  scaffold_306_328368  scaffold_684_229342 0.7989789 0.8044837
    ## 47  scaffold_306_944429  scaffold_94_1302246 0.7967807 0.8363145
    ## 48  scaffold_310_880276  scaffold_834_252344 1.0000000 1.0000000
    ## 49  scaffold_324_135439  scaffold_342_417561 0.8524590 0.8534918
    ## 50  scaffold_324_135439  scaffold_343_863766 0.8000000 1.0000000
    ## 51  scaffold_324_135439  scaffold_557_489027 1.0000000 1.0000000
    ## 52  scaffold_324_135439  scaffold_582_107987 0.8524590 0.8534918
    ## 53  scaffold_324_135439  scaffold_612_363793 0.8853755 0.8854377
    ## 54  scaffold_324_135439  scaffold_684_229342 0.8853755 0.8854377
    ## 55  scaffold_324_135439  scaffold_700_166185 0.8000000 1.0000000
    ## 56  scaffold_324_135439  scaffold_71_2209803 1.0000000 1.0000000
    ## 57  scaffold_342_417561  scaffold_343_863766 0.7027027 0.8534918
    ## 58  scaffold_342_417561  scaffold_557_489027 0.8524590 0.8534918
    ## 59  scaffold_342_417561  scaffold_700_166185 0.7027027 0.8534918
    ## 60  scaffold_342_417561  scaffold_71_2209803 0.8524590 0.8534918
    ## 61  scaffold_342_417561  scaffold_75_1184494 0.8198758 0.9135003
    ## 62  scaffold_343_863766  scaffold_557_489027 0.8000000 1.0000000
    ## 63  scaffold_343_863766  scaffold_582_107987 0.7027027 0.8534918
    ## 64  scaffold_343_863766  scaffold_612_363793 0.7133758 0.8854377
    ## 65  scaffold_343_863766  scaffold_684_229342 0.7133758 0.8854377
    ## 66  scaffold_343_863766  scaffold_700_166185 1.0000000 1.0000000
    ## 67  scaffold_343_863766  scaffold_71_2209803 0.8000000 1.0000000
    ## 68  scaffold_360_269450  scaffold_552_154281 1.0000000 1.0000000
    ## 69  scaffold_437_192045  scaffold_72_1463271 0.8000000 0.8958280
    ## 70  scaffold_51_1220763  scaffold_684_229342 0.8807986 0.8811342
    ## 71  scaffold_54_1630779  scaffold_621_290581 0.8277817 0.8277836
    ## 72  scaffold_557_489027  scaffold_582_107987 0.8524590 0.8534918
    ## 73  scaffold_557_489027  scaffold_612_363793 0.8853755 0.8854377
    ## 74  scaffold_557_489027  scaffold_684_229342 0.8853755 0.8854377
    ## 75  scaffold_557_489027  scaffold_700_166185 0.8000000 1.0000000
    ## 76  scaffold_557_489027  scaffold_71_2209803 1.0000000 1.0000000
    ## 77  scaffold_582_107987  scaffold_700_166185 0.7027027 0.8534918
    ## 78  scaffold_582_107987  scaffold_71_2209803 0.8524590 0.8534918
    ## 79  scaffold_582_222480  scaffold_612_363793 0.7337278 0.8360078
    ## 80  scaffold_612_363793  scaffold_700_166185 0.7133758 0.8854377
    ## 81  scaffold_612_363793  scaffold_71_2209803 0.8853755 0.8854377
    ## 82  scaffold_684_229342  scaffold_700_166185 0.7133758 0.8854377
    ## 83  scaffold_684_229342  scaffold_71_2209803 0.8853755 0.8854377
    ## 84  scaffold_700_166185  scaffold_71_2209803 0.8000000 1.0000000
    ## 85  scaffold_72_1463271   scaffold_8_4227715 0.7027027 0.8534918
    ## 86   scaffold_82_879210  scaffold_95_2059721 0.8400646 0.8405485

50 unique loci involved in LD in the Laysan population.

A couple of options: there’s something funky with the evolutionary
history of hybridization that’s present in the Laysan birds, or there
are functionally duplicate samples… somehow?

``` r
# are there cloned multilocus genotypes, somehow?
Laysan_pair <- clonecorrect(pop1) %>% pair.ia
```

![](02-locus-HWE-evaluation_files/figure-gfm/unnamed-chunk-23-1.png)<!-- -->

``` r
rownames_to_column(as.data.frame(Laysan_pair)) %>%
  filter(rbarD > 0.8)
```

    ##                                      rowname        Ia     rbarD
    ## 1   scaffold_1053_177405:scaffold_140_259881 0.8022284 0.8060684
    ## 2    scaffold_10_1333563:scaffold_491_362387 1.0000000 1.0000000
    ## 3   scaffold_1164_119437:scaffold_116_386790 0.7027027 0.8534918
    ## 4  scaffold_1226_115783:scaffold_227_1219626 0.7133758 0.8854377
    ## 5   scaffold_1226_115783:scaffold_324_135439 0.8853755 0.8854377
    ## 6   scaffold_1226_115783:scaffold_343_863766 0.7133758 0.8854377
    ## 7   scaffold_1226_115783:scaffold_557_489027 0.8853755 0.8854377
    ## 8   scaffold_1226_115783:scaffold_582_107987 0.9792531 0.9799367
    ## 9   scaffold_1226_115783:scaffold_700_166185 0.7133758 0.8854377
    ## 10  scaffold_1226_115783:scaffold_71_2209803 0.8853755 0.8854377
    ## 11 scaffold_1327_103421:scaffold_227_1219626 0.6896552 0.8058230
    ## 12  scaffold_1327_103421:scaffold_324_135439 0.8000000 0.8058230
    ## 13  scaffold_1327_103421:scaffold_343_863766 0.6896552 0.8058230
    ## 14  scaffold_1327_103421:scaffold_557_489027 0.8000000 0.8058230
    ## 15  scaffold_1327_103421:scaffold_700_166185 0.6896552 0.8058230
    ## 16  scaffold_1327_103421:scaffold_71_2209803 0.8000000 0.8058230
    ## 17   scaffold_155_1707144:scaffold_854_86476 1.0000000 1.0000000
    ## 18  scaffold_15_2923659:scaffold_210_1478805 0.8091895 0.8097254
    ## 19   scaffold_15_2923659:scaffold_342_417561 0.8091895 0.8097254
    ## 20 scaffold_185_1133045:scaffold_185_1154507 1.0000000 1.0000000
    ## 21 scaffold_210_1478805:scaffold_227_1219626 0.7027027 0.8534918
    ## 22  scaffold_210_1478805:scaffold_324_135439 0.8524590 0.8534918
    ## 23  scaffold_210_1478805:scaffold_342_417561 1.0000000 1.0000000
    ## 24  scaffold_210_1478805:scaffold_343_863766 0.7027027 0.8534918
    ## 25  scaffold_210_1478805:scaffold_557_489027 0.8524590 0.8534918
    ## 26  scaffold_210_1478805:scaffold_700_166185 0.7027027 0.8534918
    ## 27  scaffold_210_1478805:scaffold_71_2209803 0.8524590 0.8534918
    ## 28  scaffold_210_1478805:scaffold_75_1184494 0.8198758 0.9135003
    ## 29  scaffold_227_1219626:scaffold_324_135439 0.8000000 1.0000000
    ## 30  scaffold_227_1219626:scaffold_342_417561 0.7027027 0.8534918
    ## 31  scaffold_227_1219626:scaffold_343_863766 1.0000000 1.0000000
    ## 32  scaffold_227_1219626:scaffold_557_489027 0.8000000 1.0000000
    ## 33  scaffold_227_1219626:scaffold_582_107987 0.7027027 0.8534918
    ## 34  scaffold_227_1219626:scaffold_612_363793 0.7133758 0.8854377
    ## 35  scaffold_227_1219626:scaffold_684_229342 0.7133758 0.8854377
    ## 36  scaffold_227_1219626:scaffold_700_166185 1.0000000 1.0000000
    ## 37  scaffold_227_1219626:scaffold_71_2209803 0.8000000 1.0000000
    ## 38    scaffold_237_341989:scaffold_4_4005861 1.0000000 1.0000000
    ## 39   scaffold_249_340732:scaffold_47_2495309 1.0000000 1.0000000
    ## 40   scaffold_27_1646617:scaffold_316_437138 0.8000000 0.8000000
    ## 41   scaffold_284_808709:scaffold_437_192045 0.8000000 0.8958280
    ## 42   scaffold_284_808709:scaffold_72_1463271 1.0000000 1.0000000
    ## 43    scaffold_284_808709:scaffold_8_4227715 0.7027027 0.8534918
    ## 44   scaffold_28_4260845:scaffold_28_4262705 1.0000000 1.0000000
    ## 45   scaffold_306_328368:scaffold_612_363793 0.7989789 0.8044837
    ## 46   scaffold_306_328368:scaffold_684_229342 0.7989789 0.8044837
    ## 47   scaffold_306_944429:scaffold_94_1302246 0.7967807 0.8363145
    ## 48   scaffold_310_880276:scaffold_834_252344 1.0000000 1.0000000
    ## 49   scaffold_324_135439:scaffold_342_417561 0.8524590 0.8534918
    ## 50   scaffold_324_135439:scaffold_343_863766 0.8000000 1.0000000
    ## 51   scaffold_324_135439:scaffold_557_489027 1.0000000 1.0000000
    ## 52   scaffold_324_135439:scaffold_582_107987 0.8524590 0.8534918
    ## 53   scaffold_324_135439:scaffold_612_363793 0.8853755 0.8854377
    ## 54   scaffold_324_135439:scaffold_684_229342 0.8853755 0.8854377
    ## 55   scaffold_324_135439:scaffold_700_166185 0.8000000 1.0000000
    ## 56   scaffold_324_135439:scaffold_71_2209803 1.0000000 1.0000000
    ## 57   scaffold_342_417561:scaffold_343_863766 0.7027027 0.8534918
    ## 58   scaffold_342_417561:scaffold_557_489027 0.8524590 0.8534918
    ## 59   scaffold_342_417561:scaffold_700_166185 0.7027027 0.8534918
    ## 60   scaffold_342_417561:scaffold_71_2209803 0.8524590 0.8534918
    ## 61   scaffold_342_417561:scaffold_75_1184494 0.8198758 0.9135003
    ## 62   scaffold_343_863766:scaffold_557_489027 0.8000000 1.0000000
    ## 63   scaffold_343_863766:scaffold_582_107987 0.7027027 0.8534918
    ## 64   scaffold_343_863766:scaffold_612_363793 0.7133758 0.8854377
    ## 65   scaffold_343_863766:scaffold_684_229342 0.7133758 0.8854377
    ## 66   scaffold_343_863766:scaffold_700_166185 1.0000000 1.0000000
    ## 67   scaffold_343_863766:scaffold_71_2209803 0.8000000 1.0000000
    ## 68   scaffold_360_269450:scaffold_552_154281 1.0000000 1.0000000
    ## 69   scaffold_437_192045:scaffold_72_1463271 0.8000000 0.8958280
    ## 70   scaffold_51_1220763:scaffold_684_229342 0.8807986 0.8811342
    ## 71   scaffold_54_1630779:scaffold_621_290581 0.8277817 0.8277836
    ## 72   scaffold_557_489027:scaffold_582_107987 0.8524590 0.8534918
    ## 73   scaffold_557_489027:scaffold_612_363793 0.8853755 0.8854377
    ## 74   scaffold_557_489027:scaffold_684_229342 0.8853755 0.8854377
    ## 75   scaffold_557_489027:scaffold_700_166185 0.8000000 1.0000000
    ## 76   scaffold_557_489027:scaffold_71_2209803 1.0000000 1.0000000
    ## 77   scaffold_582_107987:scaffold_700_166185 0.7027027 0.8534918
    ## 78   scaffold_582_107987:scaffold_71_2209803 0.8524590 0.8534918
    ## 79   scaffold_582_222480:scaffold_612_363793 0.7337278 0.8360078
    ## 80   scaffold_612_363793:scaffold_700_166185 0.7133758 0.8854377
    ## 81   scaffold_612_363793:scaffold_71_2209803 0.8853755 0.8854377
    ## 82   scaffold_684_229342:scaffold_700_166185 0.7133758 0.8854377
    ## 83   scaffold_684_229342:scaffold_71_2209803 0.8853755 0.8854377
    ## 84   scaffold_700_166185:scaffold_71_2209803 0.8000000 1.0000000
    ## 85    scaffold_72_1463271:scaffold_8_4227715 0.7027027 0.8534918
    ## 86    scaffold_82_879210:scaffold_95_2059721 0.8400646 0.8405485

``` r
pop2 <- popsub(genind_df, "NWHI,_Midway")
pairwise_ld_by_locus <- pair.ia(pop2)
```

![](02-locus-HWE-evaluation_files/figure-gfm/other-pops-LD2-1.png)<!-- -->

``` r
rownames_to_column(as.data.frame(pairwise_ld_by_locus)) %>%
  filter(rbarD > 0.8)
```

    ##                                   rowname        Ia     rbarD
    ## 1 scaffold_116_386790:scaffold_834_252344 1.0000000 1.0000000
    ## 2     scaffold_16_31673:scaffold_16_32172 0.8172674 0.8179576
    ## 3 scaffold_28_4260845:scaffold_28_4262705 0.7884615 0.8034583
    ## 4 scaffold_298_359540:scaffold_298_460712 0.8217822 0.9179260
    ## 5 scaffold_298_359540:scaffold_335_662765 0.6923077 0.8689862

Okay, well, Laysan Island aside, we have the four loci to remove based
on LD from the full dataset.
