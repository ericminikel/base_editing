# cd /broad/prions/U19_proj1-2/output/2024-02-05-1509/allele_edit_summary/
# cd /broad/prions/U19_proj1-2/output/2024-02-15-102809/allele_edit_summary/

# awk '{print FILENAME "\t" $0}' *.tsv > HEK_TadCBEd.ALL.allele.edit.tsv
# awk '{print FILENAME "\t" $0}' *.tsv > HEK_TadCBEd.ALL.allele.edit.tsv

# merge all allele_edit_summary to one file

suppressPackageStartupMessages({
if (interactive()) setwd('~/src/be') 
library(tidyverse)
library(janitor)
library(readxl)
  library(openxlsx)
  
rename =dplyr::rename
})

# untreated (NT)
# unsorted TadCBEd-treated (US)
# GFP-sorted TadCBEd-treated (GFP)





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

edits_raw = read_tsv('data/rhampseq/HEK/HEK_TadCBEd.ALL.allele.edit.tsv', col_types=cols(),show_col_types = FALSE)  %>%
  clean_names() %>%
  filter(!(sample %in% c('AMPLICON_NAME','Sample'))) %>%
  rename(gRNA =g_rna,
         amplicon=sample,
         numerator=number_reads)  %>%
  mutate(numerator = as.numeric(numerator),
         denominator = as.numeric(denominator)) %>%
  mutate(sample = gsub('\\.allele\\.edit\\.tsv','',filename)) %>%
  select(-filename)


edits_raw %>% 
  group_by(amplicon) %>% 
  summarize(.groups='keep', mean_denom = mean(denominator)) -> mean_denoms

edits_raw %>%
  inner_join(filter(mean_denoms, mean_denom >= 1000), by='amplicon') -> edits_HEK



# JOIN & PROCESS DATA #####

HEK_samples = read_tsv('data/rhampseq/HEK/HEK_TadCBEd_samples.tsv', col_types=cols(),show_col_types = FALSE)
HEK_meta = read_tsv('data/rhampseq/HEK/HEK_TadCBEd_meta.tsv', col_types=cols(),show_col_types = FALSE)

edits_HEK$amplicon_sort <- edits_HEK$amplicon

suppressWarnings({
  
edits_HEK <- edits_HEK %>%
  separate(amplicon_sort, into = c("chromosome", "positions"), sep = "_") %>%
  separate(positions, into = c("start", "end"), sep = "-") %>%
  mutate(chromosome = str_replace(chromosome, "chr", ""),
         chromosome = as.integer(chromosome),
         start = as.integer(start),
         end = as.integer(end))

})
edits_HEK <- edits_HEK %>% arrange(chromosome, start) %>%
  select(-chromosome,-start,-end)


edits_HEK %>%
  filter(strand != 'File not exist') %>%
  distinct(amplicon) %>%
  mutate(x = row_number()) %>%
  mutate(y = max(x) - x + 1) -> amplicons



edits_HEK %>%
  inner_join(amplicons, by='amplicon') %>%
  inner_join(HEK_samples, by='sample') %>%
  inner_join(HEK_meta, by='treatment') %>%
  mutate(p_edited = numerator/denominator) -> plottable



## T TEST #### 



GFP_sorted <- plottable %>%
  filter(treatment == 'GFP-sorted TadCBEd-treated' | treatment == 'untreated')
unsorted <- plottable %>%
  filter(treatment == 'unsorted TadCBEd-treated' | treatment == 'untreated') 


GFP_sorted_p <- GFP_sorted %>%
  group_by(amplicon) %>%
  mutate(p_edited = replace_na(p_edited,0)) %>%
  summarize(.groups = 'keep',
            ttest_p = t.test(p_edited[treatment=='untreated'],p_edited[treatment=='GFP-sorted TadCBEd-treated'],alternative='less')$p.value) %>%
    mutate(treatment="GFP-sorted TadCBEd-treated")

sorted_untreated_p <- GFP_sorted_p %>%
  mutate(treatment="untreated")

unsorted_p <- unsorted %>%
  group_by(amplicon) %>%
  mutate(p_edited = replace_na(p_edited,0)) %>%
  summarize(.groups = 'keep',
            ttest_p = t.test(p_edited[treatment=='untreated'],p_edited[treatment=='unsorted TadCBEd-treated'],alternative='less')$p.value) %>%
    mutate(treatment="unsorted TadCBEd-treated")

unsorted_untreated_p <- unsorted_p %>%
  mutate(treatment="untreated")



merged_pvalue <- bind_rows(GFP_sorted_p,sorted_untreated_p,unsorted_p) 

GFP_sorted_p <- bind_rows(GFP_sorted_p, sorted_untreated_p)
unsorted_p <- bind_rows(unsorted_p, unsorted_untreated_p)




leg = tibble(
  treatment = c('GFP-sorted TadCBEd-treated','GFP-sorted TadCBEd-treated','unsorted TadCBEd-treated','unsorted TadCBEd-treated'),
  significant = c(TRUE, FALSE, TRUE, FALSE),
  tx_disp = c('sorted', 'sorted', 'unsorted', 'unsorted'),
  p_disp = c('P < 0.01', 'P ≥ 0.01','P < 0.01', 'P ≥ 0.01'),
  pch = c(19, 19, 1, 1),
  color = c('#000000', '#CCCCCC','#000000', '#CCCCCC')
)



GFP_sorted_p<- plottable %>%
  merge(GFP_sorted_p, by = c("amplicon", "treatment")) %>%
  mutate(significant = ttest_p <=0.01) %>%
  select(-color,-x,-y ,-offset) 

unsorted_p<- plottable %>%
  merge(unsorted_p, by = c("amplicon", "treatment")) %>%
  mutate(significant = ttest_p <=0.01) %>%
  select(-color,-x,-y ,-offset) 


plottable <- merged_pvalue %>%
  inner_join(plottable, by = c("amplicon", "treatment"))



#### plots for read counts
match_gRNA = suppressMessages(read_excel('data/rhampseq/37Xhuman_top100_BED.xlsx', sheet=1,col_names = FALSE))
match_amplicon = suppressMessages(read_excel('data/rhampseq/37Xhuman_top100_BED.xlsx', sheet=2, col_names = FALSE))

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

read_count = read_tsv('data/rhampseq/37XHuman_identified_matched.tsv', col_types=cols(),show_col_types = FALSE)%>% clean_names()
read_count$genomic_coordinate <- paste("chr",gsub(":", "_", read_count$genomic_coordinate), sep='')
read_count <- read_count %>%
  mutate(grna = genomic_coordinate) %>%
  select(grna,nuclease_read_count)





####
sorted_plot_reads <- GFP_sorted_p %>% 
  group_by(amplicon,gRNA,strand,significant) %>%
  summarize(.groups = "keep", 
            #p_unsorted = mean(ttest_p[treatment == "unsorted TadCBEd-treated"]),
            p_sorted = mean(ttest_p[treatment == "GFP-sorted TadCBEd-treated"]),
            #mean_editing_unsorted = mean(p_edited[treatment == 'unsorted TadCBEd-treated'], na.rm = TRUE) - mean(p_edited[treatment == "untreated"],na.rm = TRUE),
            mean_editing_sorted =mean(p_edited[treatment == 'GFP-sorted TadCBEd-treated'], na.rm = TRUE) - mean(p_edited[treatment == "untreated"], na.rm = TRUE),
  ) %>%
  ungroup() %>%
  mutate(mean_editing_sorted=pmax(mean_editing_sorted,0)) %>%
  inner_join(match_amplicon, by= 'amplicon') %>%
  inner_join(read_count, by='grna') %>%
  mutate(read_count =nuclease_read_count, treatment="GFP-sorted TadCBEd-treated") %>%
  select(-grna,-nuclease_read_count) %>%
  merge(leg, by=c("significant", "treatment"))
  

unsorted_plot_reads <- unsorted_p %>% 
  group_by(amplicon,gRNA,strand,significant) %>%
  summarize(.groups = "keep", 
            p_unsorted = mean(ttest_p[treatment == "unsorted TadCBEd-treated"]),
            mean_editing_unsorted = mean(p_edited[treatment == 'unsorted TadCBEd-treated'], na.rm = TRUE) - mean(p_edited[treatment == "untreated"],na.rm = TRUE),
             ) %>%
  ungroup() %>%
  mutate(mean_editing_sorted=pmax(mean_editing_unsorted,0)) %>%
  inner_join(match_amplicon, by= 'amplicon') %>%
  inner_join(read_count, by='grna') %>%
  mutate(read_count =nuclease_read_count, treatment="unsorted TadCBEd-treated") %>%
  select(-grna,-nuclease_read_count) %>%
  merge(leg, by=c("significant", "treatment"))





# PLOT #######

## HEK_TadCBEd_sorted

resx=300
png('display_items/HEK_TadCBEd_rhampseq_summary_sorted.png',width=4*resx,height=8.5*resx, res=resx)

par(mar=c(4,6,1,2))
ylims = range(plottable$y) + c(-0.5, 0.5)
xlims = c(0, 1)
#xats = 0:100/100
xbigs = 0:10/10
plot(NA, NA, xlim=xlims, ylim=ylims, axes=F, ann=F, xaxs='i', yaxs='i')
axis(side=1, at=seq(0,1,0.2), tck=-0.025, lwd.ticks=1, labels=percent(seq(0,1,0.2)),cex.axis=0.7)
axis(side=1, at=seq(0,1,0.1), tck=-0.015, lwd.ticks=1, labels=NA)

#axis(side=1, at=xats, tck=-0.025, labels=NA)
#axis(side=1, at=xbigs, tck=-0.05, labels=percent(xbigs))
mtext(side=1, line=2.5, text='% edited')
axis(side=2, at=ylims, lwd.ticks=0, labels=NA)

mtext(side=2, at=amplicons$y, text=amplicons$amplicon, las=2, line=0.25, cex=.4)




points(x=plottable$p_edited, y=plottable$y + plottable$offset, col=plottable$color, pch=1, cex=0.4)



# plot p-value 0.01
#significant_points <- plottable %>% 
#  distinct(amplicon, .keep_all = TRUE) %>%
#  filter(ttest_p < 0.01)

significant_points <- plottable %>% 
  filter( amplicon=="chr20_4679914-4680144")

#points(x = (significant_points$ttest_p)*50, y = significant_points$y, col = "black", pch = 4, cex = 0.5)

normal_points <- amplicons %>%
  filter(!(amplicon %in% significant_points$amplicon))


#points(x = (significant_points$ttest_p)*10, y = significant_points$y, col = "black", pch = 4, cex = 0.3)
#mtext(side=4, line=0.25, at=significant_points$y, text=rep("*",nrow(significant_points)), las=2)
#mtext(side=2, at=normal_points$y, text=normal_points$amplicon, las=2, line=0.25, cex=.4)
#mtext(side=2, at=significant_points$y, text=significant_points$amplicon, las=2, line=0.25, cex=.4, font=2,col = 'red')

smry<- plottable %>%
  group_by(amplicon, y, offset, color, treatment) %>%
  summarize(.groups='keep',
            mean= mean(p_edited),
            l95 = lower(p_edited),
            u95 = upper(p_edited)) %>%
  ungroup() 





barwidth= 0.15
arrowwidth = 0.015
segments(x0=smry$mean, y0=smry$y + smry$offset - barwidth, y1=smry$y + smry$offset+ barwidth, col=smry$color, lwd=1.5)
rect(xleft=rep(0, nrow(smry)), xright=smry$mean, ybottom=smry$y + smry$offset - barwidth, ytop=smry$y + smry$offset+ barwidth, col=alpha(smry$color,ci_alpha), lwd=1.5, border=NA)
suppressWarnings({
arrows(x0=smry$l95, x1=smry$u95, y0=smry$y + smry$offset, angle=90, code=3, length=arrowwidth, col=smry$color, lwd=1.5)
})

abline(v=0.01, lty=1, lwd=0.1)
par(xpd=T)
#legend(0.7, 83,legend=HEK_meta$treatment, col=HEK_meta$color, text.col=HEK_meta$color, pch=1, bty='n', cex=0.5)
par(xpd=F)

unnecssesarymessage = dev.off()



## HEK_TadCBEd_unsorted

resx=300
png('display_items/HEK_TadCBEd_rhampseq_summary_unsorted.png',width=4*resx,height=8.5*resx, res=resx)


plottable <- plottable %>%
  filter(treatment != "GFP-sorted TadCBEd-treated")

par(mar=c(4,6,1,2))
ylims = range(plottable$y) + c(-0.5, 0.5)
xlims = c(0, 1)
#xats = 0:100/100
xbigs = 0:10/10
plot(NA, NA, xlim=xlims, ylim=ylims, axes=F, ann=F, xaxs='i', yaxs='i')
axis(side=1, at=seq(0,1,0.2), tck=-0.025, lwd.ticks=1, labels=percent(seq(0,1,0.2)),cex.axis=0.7)
axis(side=1, at=seq(0,1,0.1), tck=-0.015, lwd.ticks=1, labels=NA)

mtext(side=1, line=2.5, text='% edited')
axis(side=2, at=ylims, lwd.ticks=0, labels=NA)
mtext(side=2, at=amplicons$y, text=amplicons$amplicon, las=2, line=0.25, cex=.4)



points(x=plottable$p_edited, y=plottable$y + plottable$offset, col=plottable$color, pch=1, cex=0.4)



## plot p-value 0.01
# significant_points <- plottable %>% 
#  distinct(amplicon, .keep_all = TRUE) %>%
#  filter(ttest_p < 0.01)

significant_points <- plottable %>% 
  filter( amplicon=="chr20_4679914-4680144")

# points(x = (significant_points$ttest_p)*50, y = significant_points$y, col = "black", pch = 4, cex = 0.5)

normal_points <- amplicons %>%
  filter(!(amplicon %in% significant_points$amplicon))


#points(x = (significant_points$ttest_p)*10, y = significant_points$y, col = "black", pch = 4, cex = 0.3)
#mtext(side=4, line=0.25, at=significant_points$y, text=rep("*",nrow(significant_points)), las=2)
#mtext(side=2, at=normal_points$y, text=normal_points$amplicon, las=2, line=0.25, cex=.4)
#mtext(side=2, at=significant_points$y, text=significant_points$amplicon, las=2, line=0.25, cex=.4, font=2,col = 'red')

smry<- plottable %>%
  group_by(amplicon, y, offset, color, treatment) %>%
  summarize(.groups='keep',
            mean= mean(p_edited),
            l95 = lower(p_edited),
            u95 = upper(p_edited)) %>%
  ungroup()  %>% filter(treatment != "GFP-sorted TadCBEd-treated")  






barwidth= 0.15
arrowwidth = 0.015
segments(x0=smry$mean, y0=smry$y + smry$offset - barwidth, y1=smry$y + smry$offset+ barwidth, col=smry$color, lwd=1.5)
rect(xleft=rep(0, nrow(smry)), xright=smry$mean, ybottom=smry$y + smry$offset - barwidth, ytop=smry$y + smry$offset+ barwidth, col=alpha(smry$color,ci_alpha), lwd=1.5, border=NA)
suppressWarnings({
arrows(x0=smry$l95, x1=smry$u95, y0=smry$y + smry$offset, angle=90, code=3, length=arrowwidth, col=smry$color, lwd=1.5)
})
HEK_meta<- HEK_meta %>%
  filter(treatment != "GFP-sorted TadCBEd-treated") 


abline(v=0.01, lty=1, lwd=0.1)
par(xpd=T)
#legend(0.7, 83,legend=HEK_meta$treatment, col=HEK_meta$color, text.col=HEK_meta$color, pch=1, bty='n', cex=0.5)
par(xpd=F)

unnecssesarymessage = dev.off()




## plot read counts for amplicon

resx=300
png('display_items/HEK_TadCBEd_read_count_summary.png',width=6*resx,height=5*resx, res=resx)


ylims = range(0,0.25)
xlims = range(min(sorted_plot_reads$read_count, unsorted_plot_reads$read_count, na.rm = TRUE),max(sorted_plot_reads$read_count,unsorted_plot_reads$read_count), na.rm = T)*c(0.9,1.2)
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

points(x=sorted_plot_reads$read_count, y=sorted_plot_reads$mean_editing_sorted, pch=sorted_plot_reads$pch, cex=1,col = sorted_plot_reads$color)
points(x=unsorted_plot_reads$read_count, y=unsorted_plot_reads$mean_editing_unsorted, pch=unsorted_plot_reads$pch, cex=1,col = unsorted_plot_reads$color)
abline(h=0.01, lty=3)

#legend("topright", paste0(leg$tx_disp, ' ', leg$p_disp), col=leg$color, pch=leg$pch,bty = "n", cex=0.8)
#legend(x=1900, y=0.7, pch=19, legend = c("p<=0.01", "p>0.01"), col = c("black","grey"),title ="TadCBEd_treated_GFP-sorted")

unnecssesarymessage = dev.off()



## SUMMARY

amplicon_smmary <- plottable[!is.na(plottable$p_edited), ] %>%
  group_by(amplicon, strand, treatment,sample,numerator,denominator) %>%
  mutate(ttest_p = ifelse(treatment == "untreated", 1, ttest_p)) %>%
  summarize(.groups='keep',
            edited= mean(p_edited)) %>%
  ungroup()

write_tsv(amplicon_smmary, "output/20231014_HEK_TadCBEd_indiv.tsv") 


GFP_sorted_reads <- GFP_sorted_p %>% 
  group_by(amplicon,strand) %>%
  summarize(.groups = "keep", 
            p_value = mean(ttest_p[treatment == "GFP-sorted TadCBEd-treated"]),
            editing_GFP_sorted =mean(p_edited[treatment == 'GFP-sorted TadCBEd-treated'], na.rm = TRUE),
            editing_untreated =mean(p_edited[treatment == "untreated"], na.rm = TRUE),
            difference= editing_GFP_sorted-editing_untreated,
  ) %>%
  ungroup() %>%
  inner_join(match_amplicon, by= 'amplicon') %>%
  select(-grna) %>%
  filter(!is.na(difference),!is.na(p_value))
  

unsorted_reads <- unsorted_p %>% 
  group_by(amplicon,strand) %>%
  summarize(.groups = "keep", 
            p_value = mean(ttest_p[treatment == "unsorted TadCBEd-treated"]),
            editing_unsorted =mean(p_edited[treatment == 'unsorted TadCBEd-treated'], na.rm = TRUE),
            editing_untreated =mean(p_edited[treatment == "untreated"], na.rm = TRUE),
            difference= editing_unsorted-editing_untreated,
  ) %>%
  ungroup() %>%
  inner_join(match_amplicon, by= 'amplicon') %>%
  select(-grna) %>%
  filter(!is.na(difference),!is.na(p_value))



write_tsv(GFP_sorted_reads, "output/20231014_HEK_sorted_TadCBEd_smry.tsv") 
write_tsv(unsorted_reads, "output/20231014_HEK_unsorted_TadCBEd_smry.tsv") 




# Create a new workbook
wb <- createWorkbook()

# Add sheets and write data to each sheet
addWorksheet(wb, "GFP_sorted_reads")
writeData(wb, "GFP_sorted_reads", GFP_sorted_reads)

addWorksheet(wb, "unsorted_reads")
writeData(wb, "unsorted_reads", unsorted_reads)

# Save the workbook to an Excel file
#saveWorkbook(wb, "output/20231014_HEK_TadCBEd_smry.xlsx", overwrite = TRUE)







