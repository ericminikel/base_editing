
# Import reference
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")

#BiocManager::install(c("rtracklayer", "GenomicRanges"),force =TRUE)


#BiocManager::install("BSgenome.Hsapiens.UCSC.hg19", force = TRUE)
#BiocManager::install("BSgenome.Mmusculus.UCSC.mm10", force = TRUE)

if (interactive()) setwd('~/src/be')
library(BSgenome.Mmusculus.UCSC.mm10)
library(BSgenome.Hsapiens.UCSC.hg19)
library(rtracklayer)
library(GenomicRanges)
library(Biostrings)
library(tidyverse)
library(openxlsx)
library(dplyr)

# mouse ####
# Amplicon coordinates####
amplicon_coords <- read.table("data/amplicon/amplicon_coordinates_mouse.tsv", header = FALSE, col.names = c("chr", "start", "end"))
#amplicon_coordinates_mouse.tsv

granges_coords <- GRanges(
  seqnames = amplicon_coords$chr,
  ranges = IRanges(start = amplicon_coords$start, end = amplicon_coords$end)
)

# Matching sequence####
sequences <- getSeq(BSgenome.Mmusculus.UCSC.mm10, granges_coords)
#BSgenome.Mmusculus.UCSC.mm10


# Save fasta file
#writeXStringSet(sequences, filepath = "extracted_sequences.fasta")

# Output table to Crispresso analysis

# name
amplicon_coords$AMPLICON_NAME <- paste(amplicon_coords$chr, amplicon_coords$start, amplicon_coords$end, sep = ":")
amplicon_coords$AMPLICON_NAME <- gsub(":", "-", amplicon_coords$AMPLICON_NAME, fixed = TRUE)

# combine

amplicon_sequence <- data.frame(
  AMPLICON_NAME = amplicon_coords$AMPLICON_NAME,
  AMPLICON_SEQUENCE = as.character(sequences),
  stringsAsFactors = FALSE
)




# dedup
#amplicons %>% filter(!duplicated(AMPLICON_SEQUENCE)) %>% mutate(gRNA="GGCAGCCGATACCCGGGGCAGGG") -> dedup
#write_tsv(dedup, "dedup_amplicon_sequence.tsv")



# 37Xmouse_top100_BED.xlsx ####
amplicons=read.xlsx("data/rhampseq/37Xmouse_top100_BED.xlsx", sheet="amplicon_pool", colNames =FALSE )

#37Xmouse_top100_BED.xlsx"


colnames(amplicons) = c("chrom","start","end","id","zero","strand")

guides=read.xlsx("data/rhampseq/37Xmouse_top100_BED.xlsx", sheet="gRNA_pool",colNames = FALSE)
#37Xmouse_top100_BED.xlsx

colnames(guides) = c("chrom","start","end","id","zero","strand")



# matching_gRNA ####
matching_gRNA <- amplicons %>% 
  inner_join(guides, by= c("chrom","id"), suffix = c("_a","_g")) %>% 
  select(-zero_a,-zero_g) %>%
  mutate(amplicon_sequence= as.character(getSeq(BSgenome.Mmusculus.UCSC.mm10, chrom, start_a, end_a))) %>%
  mutate(input_gRNA=as.character(getSeq(BSgenome.Mmusculus.UCSC.mm10, chrom, start_g+1, end_g,strand= strand_g)))%>%
   mutate(input_gRNA=substr(input_gRNA,1,20))%>%
   mutate(base1pos=case_when(strand_g =="+" ~ start_g, 
                             strand_g =="-" ~ end_g))
   


#output strand

matching_gRNA$AMPLICON_NAME <- paste(matching_gRNA$start_a, matching_gRNA$end_a, sep = "-")
matching_gRNA$AMPLICON_NAME <- paste(matching_gRNA$chr, matching_gRNA$AMPLICON_NAME, sep = "_")




matching_gRNA %>% 
  select(AMPLICON_NAME,amplicon_sequence,input_gRNA,strand_g)%>%
  rename(strand=strand_g) -> ampli_gRNA

write_tsv(ampli_gRNA, "output/ampli_gRNA_mouse_strand.tsv")

ampli_gRNA %>% select(-strand) ->ampli_gRNA

write_tsv(ampli_gRNA, "output/ampli_gRNA_mouse.tsv")








# human ####

## Amplicon coordinates####
amplicon_coords <- read.table("data/amplicon/amplicon_coordinates_human.tsv", header = FALSE, col.names = c("chr", "start", "end"))

granges_coords <- GRanges(
  seqnames = amplicon_coords$chr,
  ranges = IRanges(start = amplicon_coords$start, end = amplicon_coords$end)
)


## Matching sequence ####
sequences <- getSeq(BSgenome.Hsapiens.UCSC.hg19, granges_coords)


## Save fasta file
writeXStringSet(sequences, filepath = "extracted_sequences.fasta")

## Output table to Crispresso analysis

amplicon_coords$AMPLICON_NAME <- paste(amplicon_coords$chr, amplicon_coords$start, amplicon_coords$end, sep = ":")
amplicon_coords$AMPLICON_NAME <- gsub(":", "-", amplicon_coords$AMPLICON_NAME, fixed = TRUE)

# combine
amplicon_sequence <- data.frame(
  AMPLICON_NAME = amplicon_coords$AMPLICON_NAME,
  AMPLICON_SEQUENCE = as.character(sequences),
  stringsAsFactors = FALSE
)




## 37Xhouse_top100_BED.xlsx ####
amplicons=read.xlsx("data/rhampseq/37Xhuman_top100_BED.xlsx", sheet="amplicon_pool", colNames =FALSE )
colnames(amplicons) = c("chrom","start","end","id","zero","strand")

guides=read.xlsx("data/rhampseq/37Xhuman_top100_BED.xlsx", sheet="gRNA_pool",colNames = FALSE)
colnames(guides) = c("chrom","start","end","id","zero","strand")



## matching_gRNA ####
matching_gRNA <- amplicons %>% 
  inner_join(guides, by= c("chrom","id"), suffix = c("_a","_g")) %>% 
  select(-zero_a,-zero_g) %>%
  mutate(amplicon_sequence= as.character(getSeq(BSgenome.Hsapiens.UCSC.hg19, chrom, start_a, end_a))) %>%
  mutate(input_gRNA=as.character(getSeq(BSgenome.Hsapiens.UCSC.hg19, chrom, start_g+1, end_g,strand= strand_g)))%>%
  mutate(input_gRNA=substr(input_gRNA,1,20))%>%
  mutate(base1pos=case_when(strand_g =="+" ~ start_g, 
                            strand_g =="-" ~ end_g))



#output strand

matching_gRNA$AMPLICON_NAME <- paste(matching_gRNA$start_a, matching_gRNA$end_a, sep = "-")
matching_gRNA$AMPLICON_NAME <- paste(matching_gRNA$chr, matching_gRNA$AMPLICON_NAME, sep = "_")




matching_gRNA %>% 
  select(AMPLICON_NAME,amplicon_sequence,input_gRNA,strand_g)%>%
  rename(strand=strand_g)-> ampli_gRNA

write_tsv(ampli_gRNA, "output/ampli_gRNA_human_strand.tsv")

ampli_gRNA %>% select(-strand) ->ampli_gRNA

write_tsv(ampli_gRNA, "output/ampli_gRNA_human.tsv")







