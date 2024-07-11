  if (interactive()) {
    setwd('~/src/be/')
  }
  suppressPackageStartupMessages({
    library(tidyverse)
    library(janitor)
    library(dplyr)
    library(readxl)
    rename =dplyr::rename
    
  })
  
# MOUSE VEP ####
  data_m <- suppressMessages(read_tsv( "data/rhampseq/37XMouse_identified_matched.tsv", col_types=cols(),show_col_types = FALSE) %>% clean_names())
  data_m$genomic_coordinate <- gsub(":", "_", data_m$genomic_coordinate)
  data_m <- data_m %>%
    mutate(grna = genomic_coordinate) %>%
    select(grna,chromosome,site_sequence,site_sequence_gaps_allowed,strand)%>%
    separate(grna, into = c("chromosome", "positions"), sep = "_") %>%
    separate(positions, into = c("start", "end"), sep = "-") %>%
    mutate(chromosome = str_replace(chromosome, "chr", ""),
           start = as.integer(start),
           end = as.integer(end))
  
  
  data_m_na<- data_m %>%
    filter(is.na(site_sequence))%>%
    mutate(site_sequence=site_sequence_gaps_allowed)
  
  
  data_m_not_na<- data_m %>%
    filter(!is.na(site_sequence))
  
  data_m <- data_m_na %>%
    bind_rows(data_m_not_na) %>%
    select(-site_sequence_gaps_allowed)
  
  data_m$site_sequence=gsub("-","",data_m$site_sequence)
  
  
  
  
  data_m_check <- data_m %>%
    mutate(
      relative_position = pmap(list(site_sequence, strand, start, end), function(seq, strand, start, end) {
        gregexpr("C", seq)[[1]] - 1
      })
    ) %>%
    unnest(relative_position) %>%
    filter(relative_position != -1) %>% # filter out non-matches
    mutate(pos = case_when(strand=='+' ~ start + relative_position+1,
                           strand=='-' ~ end - relative_position)) %>%
    mutate(REF = case_when(strand=='+' ~ 'C',
                           strand=='-' ~ 'G')) %>%
    mutate(ALT = case_when(strand=='+' ~ 'T',
                           strand=='-' ~ 'A')) %>%
    mutate(locus = paste0(chromosome, '_', start, '_', end)) %>%
    mutate(ALT = ifelse(strand == "+", "T", "A"),
           ID = paste(chromosome, pos, REF, ALT, sep = "_"))
  
  
  data_m_out <- data_m_check %>%
    rename(CHROM=chromosome, POS=pos) %>%
    select(CHROM, POS, ID, REF, ALT) %>%
    mutate(
      QUAL = ".",
      FILTER = ".",
      INFO = "."
    ) 
  
  
  
  
  ## export vcf mouse ####
  #write_tsv(data_m_out,"output/vcf_m.tsv")
  
  
  
  
  # mm10 to mm39 ####
  
  # get the data from converting tools
  new_vep_m = read_tsv('data/vep/vcf_mm39_coordinate.tsv',show_col_types = FALSE) %>%
    select(-x6,-x7,-x8) %>%
    mutate(ID = paste(chromosome, pos_grcm39, REF, ALT, sep = "_"))%>%
    rename(CHROM=chromosome, POS=pos_grcm39) %>%
    select(CHROM, POS, ID, REF, ALT,pos_id_mm10) %>%
    mutate(
      QUAL = ".",
      FILTER = ".",
      INFO = "."
    ) 
  
  
  write_tsv(new_vep_m,'output/vcf_m_GRCm39.tsv') # export to vep tools
  
  ## match mouse pos ####
  match_gRNA = suppressMessages(read_excel('data/rhampseq/37Xmouse_top100_BED.xlsx', sheet=1,col_names = FALSE))
  match_amplicon = suppressMessages(read_excel('data/rhampseq/37Xmouse_top100_BED.xlsx', sheet=2, col_names = FALSE))
  
  match_gRNA$Genomic_Coordinate <- paste(paste(match_gRNA$...1, match_gRNA$...2, sep = "_"),match_gRNA$...3, sep = "-")  
  match_gRNA <- match_gRNA %>% 
    mutate(id =...4) %>%
    select(Genomic_Coordinate, id)
  
  
  match_amplicon$Genomic_Coordinate <- paste(paste(match_amplicon$...1, match_amplicon$...2, sep = "_"),match_amplicon$...3, sep = "-")  
  match_amplicon_m <- match_amplicon %>% 
    mutate(id =...4) %>%
    select(Genomic_Coordinate, id) %>%
    inner_join(match_gRNA, by= 'id')  %>%
    mutate(amplicon =Genomic_Coordinate.x, gRNA =Genomic_Coordinate.y) %>%
    select(amplicon,gRNA) 
  
  
  match_amplicon_m <- match_amplicon_m %>% 
    mutate(gRNA_pos=gRNA)%>%
    separate(gRNA, into = c("match", "match2"), sep = "-") %>%
    select(-match2)
  
  

  
  
  
  
  
  ### vep result ####
  vep_m_GRCm39 = read_tsv('data/vep/vep_result_m.txt',show_col_types = FALSE) %>% clean_names() %>%
    mutate(
      impactvalue = case_when(
        impact == "HIGH" ~ 4,
        impact == "MODERATE" ~ 3,
        impact == "LOW" ~ 2,
        impact == "MODIFIER" ~ 1,
        TRUE ~ 0
      )
    ) %>%
    mutate(gRNA_pos=location)%>%
    separate(gRNA_pos, into = c("match", "match2"), sep = "-") %>%
    select(-match2) %>%
    mutate(match = paste0("chr", match))
  
  
  vep_m_GRCm39$match <- gsub(":", "_", vep_m_GRCm39$match)
  vep_m_GRCm39 <- vep_m_GRCm39 %>%
    separate(match, into = c("chr", "pos"), sep = "_") %>%
    mutate(pos=as.numeric(pos)-7) %>%
    mutate(match = paste(chr, pos,sep = "_"))
  
  data_m_check %>%
    inner_join(new_vep_m, by=c('ID'='pos_id_mm10')) %>% # match to original locus position
    inner_join(vep_m_GRCm39, by=c('ID.y'='number_uploaded_variation')) %>%
    arrange(desc(impactvalue)) %>%
    group_by(locus) %>%
    slice(1) %>%
    select(orig_locus=locus, variant_GRCm39=ID.y, impact, consequence, symbol) -> worst_vep_per_locus_m_GRCm39 
  
  
  
  write_tsv(worst_vep_per_locus_m_GRCm39, 'output/worst_vep_per_locus_m_GRCm39.tsv', na='')
  
  
  
  worst_impact_m_GRCm39 <- vep_m_GRCm39 %>%
    group_by(number_uploaded_variation) %>%
    filter(impactvalue == max(impactvalue)) %>%
    ungroup() 
  
  vep_mons_smry_mouse_GRCm39 <- worst_impact_m_GRCm39 %>% 
    distinct(number_uploaded_variation, .keep_all = TRUE) %>%
    mutate(consequence= sub(",.*", "", consequence)) %>%
    rename(pos_id =number_uploaded_variation) %>%
    group_by(pos_id)  %>%
    group_by(consequence,impact) %>%
    summarise(count=n(), 
              .groups ="keep" ) %>%
    mutate(sample ='mouse')
  
  
  worst_vep_per_variant_m_GRCm39 <- worst_impact_m_GRCm39 %>%
    distinct(number_uploaded_variation, .keep_all = TRUE) %>%
    rename(pos_id =number_uploaded_variation) %>%
    group_by(pos_id) %>%
    select(pos_id_GRCm39=pos_id, consequence, impact, symbol, biotype)%>%
    mutate()
  
  
  
  
  write_tsv(vep_mons_smry_mouse_GRCm39, "output/impact_cons_smry_m_GRCm39.tsv")
  write_tsv(worst_vep_per_variant_m_GRCm39, "output/worst_vep_per_variant_m_GRCm39.tsv")
  
  
  
  
  
  
  
# HUMAN VEP ####
  data_h <- suppressMessages(read_tsv( "data/rhampseq/37XHuman_identified_matched.tsv",show_col_types = FALSE) %>% clean_names())
  data_h$genomic_coordinate <- gsub(":", "_", data_h$genomic_coordinate)
  data_h <- data_h %>%
    mutate(grna = genomic_coordinate) %>%
    select(grna,chromosome,site_sequence,site_sequence_gaps_allowed,strand)%>%
    separate(grna, into = c("chromosome", "positions"), sep = "_") %>%
    separate(positions, into = c("start", "end"), sep = "-") %>%
    mutate(chromosome = str_replace(chromosome, "chr", ""),
           start = as.integer(start),
           end = as.integer(end))
  
  
  
  
  data_h_na<- data_h %>%
    filter(is.na(site_sequence))%>%
    mutate(site_sequence=site_sequence_gaps_allowed)
  
  
  data_h_not_na<- data_h %>%
    filter(!is.na(site_sequence))
  
  data_h <- data_h_na %>%
    bind_rows(data_h_not_na) %>%
    select(-site_sequence_gaps_allowed)
  
  data_h$site_sequence=gsub("-","",data_h$site_sequence)
  
  
  
  
  
  data_h_check <- data_h %>%
    mutate(
      relative_position = pmap(list(site_sequence, strand, start, end), function(seq, strand, start, end) {
        gregexpr("C", seq)[[1]] - 1
      })
    ) %>%
    unnest(relative_position) %>%
    filter(relative_position != -1) %>% # filter out non-matches
    mutate(pos = case_when(strand=='+' ~ start + relative_position+1,
                           strand=='-' ~ end - relative_position)) %>%
    mutate(REF = case_when(strand=='+' ~ 'C',
                           strand=='-' ~ 'G')) %>%
    mutate(ALT = case_when(strand=='+' ~ 'T',
                           strand=='-' ~ 'A')) %>%
    mutate(locus = paste0(chromosome, '_', start, '_', end)) %>%
    mutate(ALT = ifelse(strand == "+", "T", "A"),
           ID = paste(chromosome, pos, REF, ALT, sep = "_"))
  
  
  
  
  data_h_out <- data_h_check %>%
    rename(CHROM=chromosome, POS=pos) %>%
    select(CHROM, POS, ID, REF, ALT) %>%
    mutate(
      QUAL = ".",
      FILTER = ".",
      INFO = "."
    ) 
  
  
  
  
  ## export vcf human ####
  write_tsv(data_h_out,"output/vcf_h_GRCh37.tsv")
  
  
  
  
  
  ## match human pos ####
  match_gRNA = suppressMessages(read_excel('data/rhampseq/37Xhuman_top100_BED.xlsx', sheet=1,col_names = FALSE))
  match_amplicon = suppressMessages(read_excel('data/rhampseq/37Xhuman_top100_BED.xlsx', sheet=2, col_names = FALSE))
  
  match_gRNA$Genomic_Coordinate <- paste(paste(match_gRNA$...1, match_gRNA$...2, sep = "_"),match_gRNA$...3, sep = "-")  
  match_gRNA <- match_gRNA %>% 
    mutate(id =...4) %>%
    select(Genomic_Coordinate, id)
  
  
  match_amplicon$Genomic_Coordinate <- paste(paste(match_amplicon$...1, match_amplicon$...2, sep = "_"),match_amplicon$...3, sep = "-")  
  match_amplicon_h <- match_amplicon %>% 
    mutate(id =...4) %>%
    select(Genomic_Coordinate, id) %>%
    inner_join(match_gRNA, by= 'id')  %>%
    mutate(amplicon =Genomic_Coordinate.x, gRNA =Genomic_Coordinate.y) %>%
    select(amplicon,gRNA)
  
  
  
  match_amplicon_h <- match_amplicon_h %>% 
    mutate(gRNA_pos=gRNA)%>%
    separate(gRNA, into = c("match", "match2"), sep = "-") %>%
    select(-match2)
  
  
  
  
  
  all_amplicon_h = read_tsv('data/rhampseq/HEK/HEK.all.allele.edit.tsv',show_col_types = FALSE) %>% clean_names()%>%
    distinct(sample)%>%
    rename(amplicon=sample)
  
  ### vep result ####
  vep_h= read_tsv('data/vep/vep_result_h.txt',show_col_types = FALSE) %>% clean_names() %>%
    mutate(
      impactvalue = case_when(
        impact == "HIGH" ~ 4,
        impact == "MODERATE" ~ 3,
        impact == "LOW" ~ 2,
        impact == "MODIFIER" ~ 1,
        TRUE ~ 0
      )
    ) %>%
    mutate(gRNA_pos=location)%>%
    separate(gRNA_pos, into = c("match", "match2"), sep = "-") %>%
    select(-match2) %>%
    mutate(match = paste0("chr", match))
  

  
  vep_h$match <- gsub(":", "_", vep_h$match)
  vep_h <- vep_h %>%
    separate(match, into = c("chr", "pos"), sep = "_") %>%
    mutate(pos=as.numeric(pos)-7) %>%
    mutate(match = paste(chr, pos,sep = "_"))
  
  
  data_h_check %>%
    inner_join(vep_h, by=c('ID'='number_uploaded_variation')) %>%
    arrange(desc(impactvalue)) %>%
    group_by(locus) %>%
    slice(1) %>%
    select(locus, variant=ID, impact, consequence, symbol) -> worst_vep_per_locus_h
  
  write_tsv(worst_vep_per_locus_h, 'output/worst_vep_per_locus_h_GRCh37.tsv', na='') 
  
  worst_impact_h <- vep_h %>%
    group_by(number_uploaded_variation) %>%
    filter(impactvalue == max(impactvalue)) %>%
    ungroup() 
  
  
  vep_mons_smry_human <- worst_impact_h %>%
    distinct(number_uploaded_variation, .keep_all = TRUE) %>%
    mutate(consequence= sub(",.*", "", consequence)) %>%
    rename(pos_id =number_uploaded_variation) %>%
    group_by(pos_id)  %>%
    group_by(consequence,impact) %>%
    summarise(count=n(),.groups = "keep") %>%
    mutate(sample ='human')
  
  worst_vep_per_variant_h <- worst_impact_h %>%
    distinct(number_uploaded_variation, .keep_all = TRUE) %>%
    rename(pos_id =number_uploaded_variation) %>%
    group_by(pos_id) %>%
    select(pos_id, consequence, impact, symbol, biotype)
  
  
  write_tsv(vep_mons_smry_human, "output/impact_cons_smry_h_GRCh37.tsv")
  write_tsv(worst_vep_per_variant_h, "output/worst_vep_per_variant_h_GRCh37.tsv")
  

  
  
  
  
  
