options(stringsAsFactors=FALSE)
  
  
  
  
suppressPackageStartupMessages({
if (interactive()) {
  setwd('~/src/be/')
}
library(tidyverse)
library(survival)
library(dplyr)
library(janitor)

})



percent = function(x, format='f', digits=0, signed=F) { paste0(ifelse(signed & x > 0, '+', ''), formatC(x*100,format=format,digits=digits),'%') }

alpha = function(rgb_hexcolor, proportion) {
    hex_proportion = sprintf("%02x",round(proportion*255))
    rgba = paste(rgb_hexcolor,hex_proportion,sep='')
    return (ifelse(is.na(rgb_hexcolor),NA,rgba))
}
ci_alpha = 0.35 # degree of transparency for shading confidence intervals in plot

meta = read_tsv('data/animals/raw/C_meta.tsv',show_col_types = FALSE)
blind = read_tsv('data/animals/raw/C_blind.tsv',show_col_types = FALSE)
master = read_tsv('data/animals/raw/C_master.tsv',show_col_types = FALSE)
surgery_notes = read_tsv('data/animals/raw/C_surgery_notes.tsv',show_col_types = FALSE)
weights = read_tsv('data/animals/processed/C_weights.tsv',show_col_types = FALSE)
nests = read_tsv('data/animals/processed/C_nests.tsv',show_col_types = FALSE)

meta = meta[with(meta, order(display)),]

master %>% 
  inner_join(blind, by='animal') %>%
  inner_join(meta, by='cohort') %>%
  select(animal, cohort, inoculum, treatment, dpi, acm, prion_endpoint, color, lty, inoculation_date) -> survdata

survdata$censored = is.na(survdata$dpi)
survdata$acm[survdata$censored] = FALSE
survdata$dpi[survdata$censored] = 600



# survdata %>% group_by(cohort) %>% summarize(.groups='keep', n=n(), mean=mean(dpi), sd=sd(dpi))
# survdata %>% filter(inoculum!='none') %>% group_by(treatment) %>% summarize(.groups='keep', n=n(), mean=mean(dpi), sd=sd(dpi))
# survdata %>% filter(acm) %>% group_by(cohort) %>% summarize(.groups='keep', mean=mean(dpi), sd=sd(dpi))
# survdata %>% group_by(cohort) %>% summarize(.groups='keep', mean=mean(dpi), sd=sd(dpi))

sf = survfit(Surv(dpi, acm) ~ cohort, data=survdata)
sf$color = meta$color[match(gsub('cohort=','',names(sf$strata)), meta$cohort)]
sf$lty = meta$lty[match(gsub('cohort=','',names(sf$strata)), meta$cohort)]

#survdiff(Surv(dpi, acm) ~ treatment, data=subset(survdata, inoculum != 'none'))
#survdiff(Surv(dpi, acm) ~ treatment, data=subset(survdata, inoculum %in% 'RES-03'))
#survdiff(Surv(dpi, acm) ~ treatment, data=subset(survdata, inoculum %in% 'RES-07'))

xlims = c(0,600)
ylims = c(0,1.05)


  
# c_surv ####

resx=300 
png('display_items/c_surv.png',width=7.2*resx,height=2.4*resx,res=resx)
par(mar=c(3,4,1,2))
plot(sf, col=sf$color, lwd=2, lty=sf$lty, xlim=xlims, ylim=ylims, axes=F, ann=F, xaxs='i', yaxs='i')

medians <-data.frame(summary(sf)$table[, "median"])%>% clean_names() %>%
  dplyr::rename(median=summary_sf_table_median)

medians$cohort=gsub("cohort=","",row.names(medians))
rownames(medians) <- 1:nrow(medians)

medians = medians %>%
  inner_join(meta, by='cohort')%>%
  select(cohort,median,display,color,lty) %>%
  mutate(plot_text=median,
         pos=median,
         signif_lev=NA)


medians[4, "pos"] <- 280
medians[3, "pos"] <- 495
medians[5, "plot_text"] <- '441'
medians[3, "plot_text"] <- '491'
medians[5, "signif_lev"] <- '**'
medians[3, "signif_lev"] <- '***'


axis(side=1, at=seq(0,600,100), labels=NA, tck=-0.05)
axis(side=1, at=seq(0,600,100), lwd=0, labels=seq(0,600,100), line=-0.5)
axis(side=1, at=seq(0,600,50), labels=NA, tck=-0.025)
mtext(side=1, line=1.5, text='days post-inoculation')
axis(side=2, at=0:4/4, labels=percent(0:4/4), las=2)
mtext(side=2, line=2.75, text='survival')
abline(h=0.5, lty=3)
par(xpd=T)
#text(medians$pos, 1.05, labels=medians$plot_text, cex=0.7, col=medians$color)
#text(medians$pos+17, 1.1, labels=medians$signif_lev, cex=0.7, col=medians$color) ## add sig level
#legend(x='bottomleft', bty='n', cex=.7, lwd=2, gsub('SpCBE3.9','CBE 1.1',meta$display), 
 #      col=meta$color, 
  #     text.col=meta$color, 
   #    lty=meta$lty, title.col='#000000', title='prion strain + treatment')
par(xpd=F)
unnecessary_message = dev.off()       
       



weights %>%
  arrange(dpi) %>%
  group_by(animal) %>%
  dplyr::slice(1) %>%
  select(animal, baseline_weight=weight) -> baselines

weights %>%
  inner_join(baselines, by='animal') %>%
  mutate(rel=weight/baseline_weight-1) %>%
  inner_join(blind, by='animal') %>%
  inner_join(meta, by='cohort') %>%
  group_by(cohort, dpi, color, lty) %>%
  summarize(.groups='keep',
            mean = mean(rel),
            sd = sd(rel),
            n = n(),
            l95 = mean(rel) - 1.96*sd(rel)/sqrt(n()),
            u95 = mean(rel) + 1.96*sd(rel)/sqrt(n())) %>%
  filter(n > 1) -> weights_rel


# c_wts ####
resx=300
png('display_items/c_wts.png',width=6.5*resx,height=3*resx,res=resx)
weights_rel -> weights_curr
par(mar=c(4,4,1,4))
xlims = c(0, 600)
ylims = c(-0.4, 0.8)
plot(NA, NA, xlim=xlims, ylim=ylims, axes=F, ann=F, xaxs='i', yaxs='i')
axis(side=1, at=seq(0,600,100), labels=NA, tck=-0.05)
axis(side=1, at=seq(0,600,100), lwd=0, labels=seq(0,600,100), line=-0.5)
axis(side=1, at=seq(0,600,50), labels=NA, tck=-0.025)
mtext(side=1, line=1.5, text='days post-inoculation')
axis(side=2, at=seq(-0.4,0.8,.2), labels=NA, las=2, tck=-0.05)
axis(side=2, at=seq(-0.4,0.8,.2), labels=percent(seq(-0.4,0.8,.2), signed=T), lwd=0, las=2, line=-0.5)
#axis(side=2, at=seq(-.3,.4,.05), labels=NA, las=2, tck=-0.025)
abline(h=0)
abline(h=-.2, col='red', lty=3)
mtext(side=4, at=-.2, las=2, col='red', text='endpoint', cex=.8)
mtext(side=2, line=2.75, text='weight (% change)')
for (coh in unique(weights_curr$cohort)) {
  subs = weights_curr[weights_curr$cohort==coh,]
  points(subs$dpi, subs$mean, type='l', lwd=2, col=subs$color, lty=subs$lty)
  polygon(x=c(subs$dpi,rev(subs$dpi)),y=c(subs$l95,rev(subs$u95)),col=alpha(subs$color,ci_alpha),border=NA)
}
par(xpd=T)
#legend(x=max(xlims),y=max(ylims),lwd=3,bty='n',cex=.8,
#       legend=meta$display[meta$cohort %in% weights_curr$cohort],
#       col=meta$color[meta$cohort %in% weights_curr$cohort],
 #      text.col=meta$color[meta$cohort %in% weights_curr$cohort],
#       lty=meta$lty[meta$cohort %in% weights_curr$cohort])
par(xpd=F)
unnecessary_message = dev.off()




master %>%
  inner_join(blind, by=c('animal', 'cage')) %>%
  inner_join(meta, by='cohort') %>%
  group_by(cage, cohort, color, lty) %>%
  summarize(.groups='keep', n=n()) -> cages



nests %>%
  inner_join(cages, by='cage') %>%
  group_by(cohort, dpi, color, lty) %>%
  summarize(.groups='keep',
            mean_comb = mean(comb)) -> nest_smry


  
# c_nests ####
resx=300

png('display_items/c_nests.png',width=6.5*resx,height=3*resx,res=resx)
nest_smry -> nest_curr
xlims = c(0, 600)
ylims = c(0, 2.05)
par(mar=c(3,4,1,2))
plot(NA, NA, axes=F, ann=F, xaxs='i', yaxs='i', xlim=xlims, ylim=ylims)
axis(side=1, at=seq(0,600,100), labels=NA, tck=-0.025)
axis(side=1, at=seq(0,600,100), lwd=0, labels=seq(0,600,100), line=-0.5)
axis(side=1, at=seq(0,600,50), labels=NA, tck=-0.015)
mtext(side=1, line=1.5, text='days post-inoculation')
axis(side=2, at=0:2, las=2)
axis(side=2, at=0:4/2, las=2, labels=NA, tck=-0.02)
mtext(side=2, line=2, text='nest score')
for (coh in unique(nest_curr$cohort)) {
  nest_curr %>%
    filter(cohort==coh) %>%
    arrange(dpi) -> subs
  points(x=subs$dpi, y=subs$mean_comb, type='l', lwd=0.25, col=subs$color, lty=subs$lty)
  #points(x=subs$dpi, y=subs$mean_comb, pch=20, cex=.25, col=subs$color)
}

# trend line for each group
for (coh in unique(nest_curr$cohort)) {
  subs = nest_curr %>% filter(cohort==coh) 
  loess_model = loess(mean_comb ~ dpi, data =subs, span = 1.5)
  dpi = seq(50,600,0.1)
  score = predict(loess_model, dpi)
  points(dpi, score, type='l', col = subs$color[1], lwd = 2, lty=subs$lty)
}


unnecssesarymessage = dev.off()

#par(xpd=T)
#legend(x=max(xlims),y=max(ylims),lwd=3,bty='n',cex=.8,
#       legend=meta$display[meta$cohort %in% weights_curr$cohort],
 #      col=meta$color[meta$cohort %in% weights_curr$cohort],
  #     text.col=meta$color[meta$cohort %in% weights_curr$cohort],
   #    lty=meta$lty[meta$cohort %in% weights_curr$cohort])
#par(xpd=F)







