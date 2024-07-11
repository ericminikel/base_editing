# cd /broad/prions/U19_proj1-2/output/2024-03-21-1000/allele_edit_summary/
# awk '{print FILENAME "\t" $0}' *.tsv > Next_seq.ALL.allele.edit.tsv

# merge all allele_edit_summary to one file


suppressPackageStartupMessages({
  
if (interactive()) setwd('~/src/be') # setwd('~/d/sci/src/be')
library(readxl)
library(tidyverse)
library(janitor)
library(dplyr)
})





# USEFUL FUNCTIONS ####

percent = function(x, digits=0, signed=F) gsub(' ','',paste0(ifelse(x > 0 & signed, '+', ''),formatC(100*x,format='f',digits=digits),'%'))

upper = function(x, ci=0.95) { 
  alpha = 1 - ci
  sds = qnorm(1-alpha/2)
  mean(x, na.rm=T) + sds*sd(x, na.rm=T)/sqrt(sum(!is.na(x)))
}
lower = function(x, ci=0.95) { 
  alpha = 1 - ci
  sds = qnorm(1-alpha/2)
  mean(x, na.rm=T) - sds*sd(x, na.rm=T)/sqrt(sum(!is.na(x)))
}

alpha = function(rgb_hexcolor, proportion) {
  hex_proportion = sprintf("%02x",round(proportion*255))
  rgba = paste(rgb_hexcolor,hex_proportion,sep='')
  return (rgba)
}

ci_alpha = 0.35 # degree of transparency for shading confidence intervals in plot

clipdist = function(x, minx, maxx) {
  return (pmin(maxx,pmax(minx,x)))
}

cast_quietly = function(x) suppressWarnings(as.numeric(x))

clipcopy = function(tbl) {
  clip = pipe("pbcopy", "w")  
  write.table(tbl, file=clip, sep = '\t', quote=F, row.names = F, na='')
  close(clip)
}

# READ IN DATA #####

datasets = tribble(
  ~edits_file, ~samples_file, ~meta_file, ~output_suffix,
  'Nextseq_allele_35dpi.tsv', 'Nextseq_sample_35dpi.tsv', 'Nextseq_meta.tsv', '35dpi',
  'Nextseq_allele_100dpi.tsv', 'Nextseq_sample_100dpi.tsv', 'Nextseq_meta.tsv', '100dpi',
  
  )


for (i in 1:nrow(datasets)) {



edits_raw = read_tsv(paste0('data/rhampseq/Next_seq/',datasets$edits_file[i]), col_types=cols(),show_col_types = FALSE) %>%
  clean_names() %>%
  filter(!(sample %in% c('AMPLICON_NAME','Sample'))) %>%
  dplyr::rename(amplicon=sample,
         grna=g_rna,
         numerator=number_reads) %>%
  mutate(numerator = as.numeric(numerator),
         denominator = as.numeric(denominator)) %>%
  mutate(sample = gsub('\\.allele\\.edit\\.tsv','',filename)) %>%
  select(-filename)


edits_raw %>% group_by(amplicon) %>% summarize(.groups='keep', mean_denom = mean(denominator)) -> mean_denoms

edits_raw %>%
inner_join(filter(mean_denoms, mean_denom >= 1000), by='amplicon') -> edits



# JOIN & PROCESS DATA #####

samples = read_tsv(paste0('data/rhampseq/Next_seq/',datasets$samples_file[i]), col_types=cols(),show_col_types = FALSE)
meta = read_tsv(paste0('data/rhampseq/Next_seq/',datasets$meta_file[i]), col_types=cols(),show_col_types = FALSE)



#sort
edits$amplicon_sort <- edits$amplicon
suppressWarnings({
  
edits <- edits %>%
  separate(amplicon_sort, into = c("chromosome", "positions"), sep = "_") %>%
  separate(positions, into = c("start", "end"), sep = "-") %>%
  mutate(chromosome = str_replace(chromosome, "chr", ""),
         chromosome = as.integer(chromosome),
         start = as.integer(start),
         end = as.integer(end))
})
edits <- edits %>% arrange(chromosome, start) %>%
  select(-chromosome,-start,-end)


edits %>%
  filter(strand != 'File not exist') %>%
  distinct(amplicon) %>%
  mutate(x = row_number()) %>%
  mutate(y = max(x) - x + 1) -> amplicons



edits %>%
  inner_join(amplicons, by='amplicon') %>%
  inner_join(samples, by='sample') %>%
  inner_join(meta, by='treatment') %>%
  mutate(p_edited = numerator/denominator) -> plottable




# T test 
library(dplyr)

# GFP-sorted TadCBEd-treated
# TadCBEd TadCBEd-treated

BE39max <- plottable %>%
  filter(treatment == 'BE39max' | treatment == 'untreated')
TadCBEd <- plottable %>%
  filter(treatment == 'TadCBEd' | treatment == 'untreated') 



## T TEST ####
BE39max_p <- BE39max %>%
  group_by(amplicon) %>%
  mutate(p_edited = replace_na(p_edited,0)) %>%
  summarize(.groups = 'keep',
            ttest_p = t.test(p_edited[treatment=='untreated'],p_edited[treatment=='BE39max'],alternative='less')$p.value) %>%
  mutate(treatment="BE39max")

BE39max_untreated_p <- BE39max_p %>%
  mutate(treatment="untreated")

TadCBEd_p <- TadCBEd %>%
  group_by(amplicon) %>%
  mutate(p_edited = replace_na(p_edited,0)) %>%
  summarize(.groups = 'keep',
            ttest_p = t.test(p_edited[treatment=='untreated'],p_edited[treatment=='TadCBEd'],alternative='less')$p.value) %>%
  mutate(treatment="TadCBEd")

TadCBEd_untreated_p <- TadCBEd_p %>%
  mutate(treatment="untreated", ttest_p = 1.0)



merged_pvalue <- bind_rows(BE39max_p,BE39max_untreated_p,TadCBEd_p) 

BE39max_p <- bind_rows(BE39max_p, BE39max_untreated_p)
TadCBEd_p <- bind_rows(TadCBEd_p, TadCBEd_untreated_p)




leg = tibble(
  treatment = c('BE39max','BE39max','TadCBEd','TadCBEd'),
  significant = c(TRUE, FALSE, TRUE, FALSE),
  tx_disp = c('BE39max', 'BE39max', 'TadCBEd', 'TadCBEd'),
  p_disp = c('P < 0.01', 'P ≥ 0.01','P < 0.01', 'P ≥ 0.01'),
  pch = c(19, 1, 19, 1),
  color = c('#fc6767', '#fc6767','#1080fe', '#1080fe')
)



BE39max_p<- plottable %>%
  merge(BE39max_p, by = c("amplicon", "treatment")) %>%
  mutate(significant = ttest_p <=0.01) %>%
  select(-color,-x,-y ,-offset) 

TadCBEd_p<- plottable %>%
  merge(TadCBEd_p, by = c("amplicon", "treatment")) %>%
  mutate(significant = ttest_p <=0.01) %>%
  select(-color,-x,-y ,-offset) 


plottable <- merged_pvalue %>%
  inner_join(plottable, by = c("amplicon", "treatment"))





## plots for read counts
match_gRNA = suppressMessages(read_excel('data/rhampseq/37Xmouse_top100_BED.xlsx', sheet=1,col_names = FALSE))
match_amplicon = suppressMessages(read_excel('data/rhampseq/37Xmouse_top100_BED.xlsx', sheet=2, col_names = FALSE))

match_gRNA$Genomic_Coordinate <- paste(paste(match_gRNA$...1, match_gRNA$...2, sep = "_"),match_gRNA$...3, sep = "-")  
match_gRNA <- match_gRNA %>% 
  mutate(id =...4) %>%
  select(Genomic_Coordinate, id)


match_amplicon$Genomic_Coordinate <- paste(paste(match_amplicon$...1, match_amplicon$...2, sep = "_"),match_amplicon$...3, sep = "-")  
match_amplicon <- match_amplicon %>% 
  mutate(id =...4) %>%
  select(Genomic_Coordinate, id) %>%
  inner_join(match_gRNA, by= 'id')  %>%
  mutate(amplicon =Genomic_Coordinate.x, grna =Genomic_Coordinate.y) %>%
  select(amplicon,grna)

read_count = read_tsv('data/rhampseq/37XMouse_identified_matched.tsv', col_types=cols(),show_col_types = FALSE)%>%clean_names()
read_count$genomic_coordinate <- paste(gsub(":", "_", read_count$genomic_coordinate), sep='')
read_count <- read_count %>%
  mutate(grna = genomic_coordinate) %>%
  select(grna,nuclease_read_count)




BE39max_plot_reads <- BE39max_p %>% 
  group_by(amplicon,strand,significant) %>%
  summarize(.groups = "keep", 
            p_BE39max = mean(ttest_p[treatment == "BE39max"]),
            mean_editing_BE39max =mean(p_edited[treatment == 'BE39max'], na.rm = TRUE) - mean(p_edited[treatment == "untreated"], na.rm = TRUE),
  ) %>%
  ungroup() %>%
  mutate(mean_editing_BE39max=pmax(mean_editing_BE39max,0)) %>%
  inner_join(match_amplicon, by= 'amplicon') %>%
  inner_join(read_count, by='grna') %>%
  mutate(read_count =nuclease_read_count, treatment="BE39max") %>%
  select(-grna,-nuclease_read_count) %>%
  merge(leg, by=c("significant", "treatment"))
  

TadCBEd_plot_reads <- TadCBEd_p %>% 
  group_by(amplicon,strand,significant) %>%
  summarize(.groups = "keep", 
            p_TadCBEd = mean(ttest_p[treatment == "TadCBEd"]),
            mean_editing_TadCBEd = mean(p_edited[treatment == 'TadCBEd'], na.rm = TRUE) - mean(p_edited[treatment == "untreated"],na.rm = TRUE),
             ) %>%
  ungroup() %>%
  mutate(mean_editing_sorted=pmax(mean_editing_TadCBEd,0)) %>%
  inner_join(match_amplicon, by= 'amplicon') %>%
  inner_join(read_count, by='grna') %>%
  mutate(read_count =nuclease_read_count, treatment="TadCBEd") %>%
  select(-grna,-nuclease_read_count) %>%
  merge(leg, by=c("significant", "treatment"))




# PLOT #######


resx=300
png(paste0('display_items/Nextseq_rhampseq_summary_',datasets$output_suffix[i],'.png'),width=3.25*resx,height=8.5*resx, res=resx)
par(mar=c(4,6,1,1))
ylims = range(plottable$y) + c(-0.5, 0.5)
#xlims = c(0, max(plottable$p_edited, na.rm=T)*1.1)
xlims = c(0, 0.3)

xats = 0:100/20
xbigs = 0:10/10
plot(NA, NA, xlim=xlims, ylim=ylims, axes=F, ann=F, xaxs='i', yaxs='i')
axis(side=1, at=xats, tck=-0.025, labels=NA)
axis(side=1, at=xbigs, tck=-0.05, labels=percent(xbigs))
mtext(side=1, line=2.5, text='% edited')
axis(side=2, at=ylims, lwd.ticks=0, labels=NA)
mtext(side=2, at=amplicons$y, text=amplicons$amplicon, las=2, line=0.25, cex=.4)
points(x=plottable$p_edited, y=plottable$y + plottable$offset, col=plottable$color, pch=1, cex=0.4)



# LABEL P-VALUE 0.01
significant_points <- plottable %>% 
  filter(ttest_p < 0.01) %>%
  distinct(amplicon, .keep_all = TRUE) 
 

#points(x = (significant_points$ttest_p)*50, y = significant_points$y, col = "black", pch = 4, cex = 0.5)

normal_points <- amplicons %>%
  filter(!(amplicon %in% significant_points$amplicon))


#points(x = (significant_points$ttest_p)*10, y = significant_points$y, col = "black", pch = 4, cex = 0.3)
#mtext(side=4, line=0.25, at=significant_points$y, text=rep("*",nrow(significant_points)), las=2)
#mtext(side=2, at=normal_points$y, text=normal_points$amplicon, las=2, line=0.25, cex=.4)
#mtext(side=2, at=significant_points$y, text=significant_points$amplicon, las=2, line=0.25, cex=.4, font=2,col = 'red')

plottable %>%
  group_by(amplicon, y, offset, color, treatment) %>%
  summarize(.groups='keep',
            mean= mean(p_edited),
            l95 = lower(p_edited),
            u95 = upper(p_edited)) %>%
  ungroup() -> smry






barwidth= 0.15
arrowwidth = 0.015
rect(xleft=rep(0, nrow(smry)), xright=smry$mean, ybottom=smry$y + smry$offset - barwidth, ytop=smry$y + smry$offset+ barwidth, col=alpha(smry$color,ci_alpha), lwd=1.5, border=NA)
suppressWarnings({
  
arrows(x0=smry$l95, x1=smry$u95, y0=smry$y + smry$offset, angle=90, code=3, length=arrowwidth, col=smry$color, lwd=1.5)
})
abline(v=0.01, lty=1, lwd=0.1)
#legend('topright', legend=meta$treatment, col=meta$color, text.col=meta$color, pch=1, bty='n', cex=0.5)
dev.off()




## PLOT READ COUNTS

resx=300
png(paste0('display_items/Nextseq_read_count_summary_',datasets$output_suffix[i] ,'.png'),width=6*resx,height=5*resx, res=resx)


ylims = range(0,0.25)
xlims = range(min(BE39max_plot_reads$read_count, TadCBEd_plot_reads$read_count, na.rm = TRUE),max(BE39max_plot_reads$read_count,TadCBEd_plot_reads$read_count), na.rm = T)*c(0.9,1.2)
plot(NA, NA, xlim=xlims, ylim=ylims, axes=F, ann=F, xaxs='i', yaxs='i', log ="x")

xats = rep(1:9,4) * 10^rep(0:3,each=9)
xbigs = 10^(0:3)
axis(side=1, at=xats, tck=-0.025, lwd.ticks=1, labels=NA)
axis(side=1, at=xbigs, tck=-0.05, lwd.ticks=1, labels=NA)
axis(side=1, at=xbigs, tck=-0.05, lwd=0, line=-0.5, labels=xbigs, cex.axis=0.8)
mtext(side=1, line=2.5, text='read count')

yats = 0:100/100
ybigs = 0:10/10
axis(side=2, at=yats, tck=-0.025, lwd.ticks=1, labels=NA)
axis(side=2, at=ybigs, tck=-0.05, lwd.ticks=1, labels=NA)
axis(side=2, at=ybigs, tck=-0.05, lwd=0, line=-0.25, labels=percent(ybigs), cex.axis=0.8, las=2)
axis(side=2, at=0.01, tck=-0.05, lwd=0, line=-0.25, labels='1%', cex.axis=0.8, las=2)
mtext(side=2, line=2.5, text='edited (%)')

points(x=BE39max_plot_reads$read_count, y=BE39max_plot_reads$mean_editing_BE39max, pch=BE39max_plot_reads$pch, cex=1,col = BE39max_plot_reads$color)
points(x=TadCBEd_plot_reads$read_count, y=TadCBEd_plot_reads$mean_editing_TadCBEd, pch=TadCBEd_plot_reads$pch, cex=1,col = TadCBEd_plot_reads$color)
abline(h=0.01, lty=3)

legend("topright", paste0(leg$tx_disp, ' ', leg$p_disp), col=leg$color, pch=leg$pch,bty = "n", cex=0.8)

#legend(x=1900, y=0.7, pch=19, legend = c("p<=0.01", "p>0.01"), col = c("black","grey"),title ="TadCBEd_treated_GFP-sorted")

dev.off()



# SUMMARY TABLE

amplicon_smmary <- plottable[!is.na(plottable$p_edited), ] %>%
  group_by(amplicon, strand, treatment,sample,numerator,denominator) %>%
  mutate(ttest_p = ifelse(treatment == "untreated", 1, ttest_p)) %>%
  summarize(.groups='keep',
            edited= mean(p_edited)) %>%
  ungroup()

write_tsv(amplicon_smmary, paste0('output/NextSeq_',datasets$output_suffix[i],'_indiv.tsv' ) )

BE39max_reads <- BE39max_p%>% 
  group_by(amplicon,strand) %>%
  summarize(.groups = "keep", 
            p_value = mean(ttest_p[treatment == "BE39max"]),
            editing_BE39max =mean(p_edited[treatment == 'BE39max'], na.rm = TRUE),
            editing_untreated =mean(p_edited[treatment == "untreated"], na.rm = TRUE),
            difference= editing_BE39max-editing_untreated,
  ) %>%
  ungroup() %>%
  inner_join(match_amplicon, by= 'amplicon') %>%
  select(-grna) %>%
  filter(!is.na(editing_untreated))

TadCBEd_reads <- TadCBEd_p%>% 
  group_by(amplicon,strand) %>%
  summarize(.groups = "keep", 
            p_value = mean(ttest_p[treatment == "TadCBEd"]),
            editing_TadCBEd =mean(p_edited[treatment == 'TadCBEd'], na.rm = TRUE),
            editing_untreated =mean(p_edited[treatment == "untreated"], na.rm = TRUE),
            difference= editing_TadCBEd-editing_untreated
  ) %>%
  ungroup() %>%
  inner_join(match_amplicon, by= 'amplicon') %>%
  select(-grna) %>%
  filter(!is.na(editing_untreated))


write_tsv(BE39max_reads, paste0('output/NextSeq_',datasets$output_suffix[i],'_smry_BE39max.tsv' ) )
write_tsv(TadCBEd_reads, paste0('output/NextSeq_',datasets$output_suffix[i],'_smry_TadCBEd.tsv' ) )


}




