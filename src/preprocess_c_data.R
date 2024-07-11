options(stringsAsFactors=F)
setwd('~/d/sci/src/be')
library(reshape2)
library(survival)
library(tidyverse)

expts = c('C')

for (expt in expts) {
  
  #master = read.table(paste('../data/animals/raw/',expt,'_master.tsv',sep=''),sep='\t',header=T,quote='',comment.char='')
  #blind = read.table(paste('../data/animals/raw/blind.tsv',sep=''),sep='\t',header=T,quote='',comment.char='')

  raw_weights = read.table(paste('../data/animals/raw/',expt,'_weights.tsv',sep=''),sep='\t',header=T,quote='',comment.char='')
  raw_weights = raw_weights[raw_weights$animal != '',]
  surgery_notes = read.table(paste('../data/animals/raw/',expt,'_surgery_notes.tsv',sep=''),sep='\t',header=T,quote='',comment.char='')
  master = read.table(paste('../data/animals/raw/',expt,'_master.tsv',sep=''),sep='\t',header=T,quote='',comment.char='')

  melted_weights = melt(raw_weights,id.vars=c('animal'))
  colnames(melted_weights) = c('animal','datetext','weight')
  melted_weights$weight = suppressWarnings(as.numeric(melted_weights$weight)) # this will force text, empty, and NA fields to NA
  melted_weights = melted_weights[!is.na(melted_weights$weight),]
  melted_weights$date = as.Date(gsub('X','',melted_weights$datetext),format='%Y.%m.%d')
  
  if ('baseline_date' %in% colnames(raw_weights)) {
    baselines = data.frame(animal=raw_weights$animal, date=as.Date(raw_weights$baseline_date,format='%Y-%m-%d'), weight=raw_weights$baseline_weight)
    melted_weights = rbind(baselines[,c('animal','date','weight')], melted_weights[,c('animal','date','weight')])
  } 
  
  melted_weights$inoculation_date = as.Date(master$inoculation_date[match(melted_weights$animal, master$animal)],format='%Y-%m-%d')
  melted_weights$dpi = suppressMessages(as.integer(melted_weights$date - melted_weights$inoculation_date))

  melted_weights = melted_weights[!is.na(melted_weights$dpi),]
  melted_weights$weight = suppressMessages(as.numeric(melted_weights$weight))
  melted_weights = melted_weights[!is.na(melted_weights$weight),]
  
  write.table(melted_weights[,c('animal','dpi','weight')], paste('../data/animals/processed/',expt,'_weights.tsv',sep=''),sep='\t',row.names=F,col.names=T,quote=F)
  
  nests = read_tsv(paste0('../data/animals/raw/',expt,'_nests.tsv'))
  #surgery_notes$cage = master$cage[match(surgery_notes$animal, master$animal)]
  nests$inoculation_date = as.Date(master$inoculation_date[match(nests$cage, master$cage)],format='%Y-%m-%d')
  nests %>%
    pivot_longer(cols=-c(cage,inoculation_date)) %>%
    select(cage, inoculation_date, date=name, score=value) %>%
    mutate(edry = suppressMessages(as.numeric(gsub('\\/.*','',score))),
           nslt = suppressMessages(as.numeric(gsub('.*\\/','',score)))) %>%
    mutate(comb = (edry + nslt)/2,
           dpi = suppressMessages(as.integer(as.Date(date) - inoculation_date))) %>%
    filter(!is.na(comb)) %>%
    select(cage, dpi, comb, edry, nslt) -> nests_processed
  
  write.table(nests_processed,paste('../data/animals/processed/',expt,'_nests.tsv',sep=''),sep='\t',row.names=F,col.names=T,quote=F)
  
  

}





