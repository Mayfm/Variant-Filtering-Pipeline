# Family-Based Rare Variant Filtering Pipeline


# Load libraries

library(tidyverse)
library(readxl)
library(tibble)


# Import annotated variant table

variants <- read.csv("FAM7.csv")


#Initial data exploration####

dim(variants)
colnames(variants)

# Check genotype encoding
unique(variants$Proband.ZYG)
unique(variants$Father.ZYG)
unique(variants$Mother.ZYG)
unique(variants$AfSib.ZYG)
unique(variants$UASib.ZYG)

# Check available annotations and filters
unique(variants$ANNOTATION)
unique(variants$FILTER)
unique(variants$X.CHR)

# Quality filtering####
rf1 <- variants %>%
  filter(
    QUAL > 30,
    FILTER == "PASS",
    GQ > 20
  )

cat("After quality filter:", nrow(rf1), "variants\n")

# Population frequency filtering####

rf2 <- rf1 %>%
  filter(
    is.na(TGP_FREQ)     | TGP_FREQ < 0.01,
    is.na(ESP_FREQ)     | ESP_FREQ < 0.01,
    is.na(EVE_ALT_FREQ) | EVE_ALT_FREQ < 0.01
  )

cat("After frequency filter:", nrow(rf2), "variants\n")


# Functional consequence filtering####

rf3 <- rf2 %>%
  filter(
    ANNOTATION %in% c(
      "splicing",
      "nonsynonymous SNV",
      "stopgain SNV",
      "frameshift insertion",
      "frameshift deletion",
      "exonic;splicing",
      "stoploss SNV"
    )
  )

cat("After functional filter:", nrow(rf3), "variants\n")


# Autosomal recessive segregation filter####

rf4 <- rf3 %>%
  filter(
    Proband.ZYG == "Proband:hom",
    Father.ZYG  == "Father:het",
    Mother.ZYG  == "Mother:het",
    AfSib.ZYG   == "AfSib:hom",
    UASib.ZYG %in% c(
      "UASib:het",
      "UASib:na"
    )
  )

cat("After segregation filter:", nrow(rf4), "variants\n")


# Candidate variant summary#####

candidate_variants <- rf4 %>%
  select(
    X.CHR,
    START,
    END,
    REF,
    ALT,
    GENE,
    GENE_NAME,
    ANNOTATION,
    NT_CHANGE,
    AA_CHANGE,
    
    Proband.ZYG,
    Father.ZYG,
    Mother.ZYG,
    AfSib.ZYG,
    UASib.ZYG,
    
    TGP_FREQ,
    ESP_FREQ,
    EVE_ALT_FREQ,
    
    SIFT_PRED,
    PPH2_PRED,
    MTT_PRED,
    
    GERP,
    PHASTCONS,
    GRANTHAM_SC
  ) %>%
  arrange(GENE)

candidate_variants


# Export final candidates#####

write_csv(
  candidate_variants,
  "candidate_variants.csv"
)


# Filtering summary####

cat("Initial variants:      ", nrow(variants), "\n")
cat("After quality filter:  ", nrow(rf1), "\n")
cat("After frequency filter:", nrow(rf2), "\n")
cat("After functional filter:", nrow(rf3), "\n")
cat("After segregation filter:", nrow(rf4), "\n")


# Alternative inheritance models####

# X-linked model

x_linked <- rf2 %>%
  filter(
    X.CHR == "X",
    ANNOTATION %in% c(
      "splicing",
      "nonsynonymous SNV",
      "stopgain SNV",
      "frameshift insertion",
      "frameshift deletion",
      "exonic;splicing",
      "stoploss SNV"
    ),
    Proband.ZYG == "Proband:hom",
    Father.ZYG %in% c(
      "Father:het",
      "Father:hom"
    ),
    Mother.ZYG == "Mother:het"
  )

nrow(x_linked)
unique(x_linked$GENE)


# Relaxed autosomal recessive model

ar_relaxed <- rf3 %>%
  filter(
    Proband.ZYG == "Proband:hom",
    Father.ZYG  == "Father:het",
    Mother.ZYG  == "Mother:het"
  ) %>%
  select(
    GENE,
    ANNOTATION,
    NT_CHANGE,
    AA_CHANGE,
    AfSib.ZYG,
    UASib.ZYG,
    TGP_FREQ,
    ESP_FREQ
  )

ar_relaxed


# Exploratory de novo model

de_novo <- rf3 %>%
  filter(
    Proband.ZYG == "Proband:het",
    Father.ZYG  %in% c("Father:na"),
    Mother.ZYG  %in% c("Mother:na")
  )

de_novo %>%
  count(GENE, sort = TRUE)