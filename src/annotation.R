setwd('~/src/be')
library(tidyverse)
library(janitor)
library(dplyr)
rename =dplyr::rename


data <- read_tsv("data/result-human-gRNA-mouse-genome.txt", col_types=cols())
#data <- read_tsv("data/results_vep.tsv", col_types=cols())


# Transform the data
vcf_data <- data %>%
  mutate(
    POS = Position + 6,
    REF = substr(DNA, 7, 7),
    ALT = case_when(
      REF == "T" ~ "C",
      REF == "C" ~ "T",
      REF == "A" ~ "G",
      REF == "G" ~ "A",
      TRUE ~ REF
    ),
    ID = paste(Chromosome, POS, REF, ALT, sep = "_")
  ) %>%
  select(Chromosome, POS, ID, REF, ALT) %>%
  mutate(
    QUAL = ".",
    FILTER = ".",
    INFO = "."
  )


write.table(vcf_data,file = "human-gRNA-mouse-genome_variants.vcf",quote = FALSE,sep = "\t",col.names = FALSE,row.names = FALSE)



#  Impact annotations data - VEP ####
## human ####
annotations <- read_tsv("data/circleseq/vep/URaejT6W3owC6nfq.txt",col_types=cols()) %>%clean_names() #human
#x726ggX3l1eY5Aps.txt human-mouse


annotations <- annotations %>%
  mutate(
    impactvalue = case_when(
      impact == "HIGH" ~ 4,
      impact == "MODERATE" ~ 3,
      impact == "LOW" ~ 2,
      impact == "MODIFIER" ~ 1,
      TRUE ~ 0
    )
  )

worst_impact_annotations <- annotations %>%
  group_by(number_uploaded_variation) %>%
  filter(impactvalue == max(impactvalue)) %>%
  ungroup() 

vep_cons_smry_human <- worst_impact_annotations %>%
  mutate(consequence= sub(",.*", "", consequence)) %>%
  group_by(consequence) %>%
  summarise(count=n()) %>%
  mutate(sample ='human')

unique_worst_impact_annotations_human <- worst_impact_annotations %>%
  distinct(number_uploaded_variation, .keep_all = TRUE) %>%
  rename(pos_id =number_uploaded_variation) %>%
  group_by(pos_id) %>%
  select(pos_id, consequence, impact, symbol, biotype)%>%
  mutate(sample ='human')
  

write_csv(unique_worst_impact_annotations_human, "output/impact_annot_human.csv")
write_csv(vep_cons_smry_human, "output/impact_cons_smry_human.csv")




## mouse ####

annotations <- read_tsv("data/circleseq/vep/x726ggX3l1eY5Aps.txt",col_types=cols()) %>%clean_names() #human-mouse

annotations <- annotations %>%
  mutate(
    impactvalue = case_when(
      impact == "HIGH" ~ 4,
      impact == "MODERATE" ~ 3,
      impact == "LOW" ~ 2,
      impact == "MODIFIER" ~ 1,
      TRUE ~ 0
    )
  )

worst_impact_annotations <- annotations %>%
  group_by(number_uploaded_variation) %>%
  filter(impactvalue == max(impactvalue)) %>%
  ungroup() 

vep_cons_smry_mouse <- worst_impact_annotations %>%
  mutate(consequence= sub(",.*", "", consequence)) %>%
  group_by(consequence) %>%
  summarise(count=n()) %>%
  mutate(sample ='mouse')

unique_worst_impact_annotations_mouse <- worst_impact_annotations %>%
  distinct(number_uploaded_variation, .keep_all = TRUE) %>%
  rename(pos_id =number_uploaded_variation) %>%
  group_by(pos_id) %>%
  select(pos_id, consequence, impact, symbol, biotype) %>%
  mutate(sample ='mouse')


write_csv(unique_worst_impact_annotations_mouse, "output/impact_annot_mouse.csv")
write_csv(vep_cons_smry_mouse, "output/impact_cons_smry_mouse.csv")







mouse_human_match <- read_tsv("data/jax_HOM_MouseHumanSequence.tsv",col_types=cols()) %>% clean_names()

mouse_human_match <- mouse_human_match %>%
  select(db_class_key, common_organism_name,symbol) %>%
  pivot_wider(names_from = common_organism_name, values_from = symbol)%>%
  clean_names()%>%
  unnest(c(mouse_laboratory),keep_empty = TRUE) %>%
  unnest(c(human),keep_empty = TRUE)%>%
  rename(mouse= mouse_laboratory)


human_tumor_suppressors <- read_tsv("data/circleseq/intogen_tumor_suppressors.tsv",col_types=cols())

tumor_suppressors <- human_tumor_suppressors %>%
  rename(human= gene) %>%
  inner_join(mouse_human_match, by= 'human') %>%
  select(-db_class_key)

human_tumor_suppressor  <- unique_worst_impact_annotations_human %>%
  clean_names() %>%
  mutate(match = symbol %in% tumor_suppressors$human)


write_csv(human_tumor_suppressor, "match_suppressor_human.csv")










