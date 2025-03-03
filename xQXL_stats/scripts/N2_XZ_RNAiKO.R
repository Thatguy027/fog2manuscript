library(xQTLStats)
library(qs)
library(tidyverse)
library(vcfR)
library(AlphaSimR)
library(Rfast)
library(ggpubr)

setwd(glue::glue("{dirname(rstudioapi::getActiveDocumentContext()$path)}/.."))

#pull from github
xQTLSims.dir = '/Users/Stefan/github_repos/xQTLSims/'

source.dir=paste0(xQTLSims.dir, 'R/')
data.dir=paste0(xQTLSims.dir, 'data/')
project.dir=paste0(xQTLSims.dir, 'projects/032924/')
# project.dir=paste0(xQTLSims.dir, 'projects/old_ECA191/')

#function to simulate crosses, treatement of X is incomplete/broken
#and additional helper functions
source(paste0(source.dir, 'simWormCrosses.R'))
source(paste0(source.dir, 'helperFxs.R'))
source(paste0(source.dir, 'makeCountTables.R'))

# modified to deal with old data

makeCountTables_old <- function(sample.key, vcf, gt) {
  countdfs=list() 
  for(i in 1:nrow(sample.key) ) { 
    #nrow(sample.key)) {
    p1=sample.key$'parent 1'[i]
    p2=sample.key$'parent 2'[i]
    
    #in this case, mate herm N2 with male CB
    p.names=c(p1, p2) #B4856') #XZ1516')
    mating.matrix=matrix(c(1,2), nrow=1)
    # in this case, given one-way cross, potentially allow N2 herm to self 
    founderPop = createFounderPop(vcf,gt, p.names, gmap, X.drop=F) #c('N2', 'XZ1516'))
    
    genMap=getGenMap(founderPop)
    
    
    sn=sample.key$'sample name'[i]
    scounts=read_tsv(paste0(sample.dir, sn, sample.suffix)) #'.table'))
    scounts.sub=scounts[paste0(scounts$contig, '_', scounts$position) %in% genMap$id,]
    scounts=data.frame(id=paste0(scounts.sub$contig, '_', scounts.sub$position),ref=scounts.sub$refCount, alt=scounts.sub$altCount)
    scounts=left_join(genMap, scounts, by='id') %>%
      dplyr::distinct(id, .keep_all=T)
    #come on GATK aseReadCounter, wtf is up with different length output given different BAM input
    #for know fill out with 0s 
    scounts$ref[is.na(scounts$ref)]=0
    scounts$alt[is.na(scounts$alt)]=0
    names(scounts)[1]='ID'
    
    #phase it 
    countdf=phaseBiparental(scounts, p.names[1], founderPop, genMap)
    
    #note, we need a better structure for keeping track of which parent is which 
    attr(countdf, 'p1')=p.names[1]
    attr(countdf, 'p2')=p.names[2]
    
    snn=paste(sample.key[i,], collapse='_')
    countdfs[[snn]]=countdf
  }
  
  return(countdfs)
}

#unique chromosomes 
uchr=c(as.character(as.roman(1:5)), 'X') #paste0('chr', as.roman(1:16))

gmap.file=paste0(data.dir, 'geneticMapXQTLsnplist.rds')
gmap=restructureGeneticMap(gmap.file)

#pretty intensive memory usage here
#include web link to vcf file 
#elegans.isotypes.vcf=paste0(data.dir,'WI.20220216.impute.isotype.vcf.gz')

#filtered vcf as qsave objects
# !!!! find the premade objects here folks !!!! :
# /u/project/kruglyak/jsbloom/xQTL/elegans/ref/
#and place in your data.dir
elegans.isotypes.vcf.qs=paste0(data.dir,'WI.20220216.vcf.qs')

#filtered numeric gt calls from vcf  as qsave object
elegans.isotypes.vcf.gt.qs=paste0(data.dir,'WI.20220216.vcf.GT.qs')

#run once, Laura skip this ============================================================
#use premade objects to save yourself the memory related headache, but this is how the objects are made
#preprocessVCF(elegans.isotypes.vcf,elegans.isotypes.vcf.qs,elegans.isotypes.vcf.gt.qs)
#======================================================================================

vcf=qread(elegans.isotypes.vcf.qs)
gt=qread(elegans.isotypes.vcf.gt.qs)

#some samples had low depth from "old" stuff below
#this script combined allele counts from multiple sequencing runs
# ~/UCLA/Projects/bulkGWAS/2024_nic_crosses2/scripts/20240430_process_run1_run2.R

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # rnai - ko 
sample.key.file=paste0('../xQTLSims/projects/032924/samplekey_032924.txt')

sample.key=read_tsv(sample.key.file) %>%
  dplyr::select(`sample name`:condition)
sample.dir='../xQXL_stats/data/2024_combined_xQTL/'
sample.suffix=".table"
# countdfs=makeCountTables(sample.key, vcf,gt)

# sample.key.small <- sample.key %>%
#   dplyr::filter(`sample name` %in% c("S1", "S2", "S3", "S12", "S25", "S26"))
# countdfs=makeCountTables(sample.key.small, vcf, gt)
sample.key.mico <- sample.key %>%
  dplyr::filter(`sample name` %in% c("S2", "S3"))

countdfs=makeCountTables(sample.key.mico, vcf, gt)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # ko - ht115
sample.key.file=paste0('../xQTLSims/projects/070824/2040708_Sample_key.tsv')

sample.key=read_tsv(sample.key.file) %>%
  dplyr::select(`sample`:condition)
sample.dir='../xQXL_stats/data/20240708_HT115_v_KO/'
sample.suffix=".table"
# countdfs=makeCountTables(sample.key, vcf,gt)

sample.key.small <- sample.key %>%
  dplyr::filter(`sample` %in% c("S2", "S6", "S8", "S12")) %>%
  dplyr::rename(`sample name` = sample)
countdfs=makeCountTables(sample.key.small, vcf, gt)


for(cts in 1:length(countdfs)){
  countdfs[[cts]] <- countdfs[[cts]] %>%
    dplyr::filter(!(ref==0 & alt == 0) & ref+alt > 10 & !(ref==0 | alt ==0))
}

#calculate the allele frequencies and SEs, params here need some more dialing in 
#should parallelize this .... yawn 
afds=lapply(names(countdfs), function(snn) {
  calcAFD(countdfs[[snn]], experiment.name=snn,sample.size=1e4, sel.strength=.95, bin.width=9000, eff.length=2000, uchr=uchr)
})
names(afds)=names(countdfs)


plots=lapply(names(afds), function(snn) {
  plotIndividualExperiment(afds[[snn]], snn) 
})
names(plots)=names(countdfs)
#e.g.  to visualize just call
#plots[[1]]

#dump individual plots somewhere ---------------
plot.dir= '../plots/'
for(snn in names(plots)) {
  ggsave(paste0( plot.dir, snn, '_9000_rmv0_dp10.png'), plots[[snn]], width=16)
}


# comp_df <- data.frame(contrast1=c(21,1,2, 17, 18,13,15,23,25,21,24,26,28),
#                       contrast2=c(22,12,3, 19,20, 19,20,12,12,12,1,1,1))
# 
# comp_df <- data.frame(contrast1=c(1,3),
#                       contrast2=c(2,4))
# ko - rnai - need to fix 
comp_df <- data.frame(contrast1=c(1),
                      contrast2=c(2))

# ko - ht115
comp_df <- data.frame(contrast1=c(1,1,2,3),
                      contrast2=c(2,3,4,4))

#I have to restructure this stuff  ---------------------------------------------

results<-list()
for(comp_plots in 1:nrow(comp_df)){
  ind1=comp_df$contrast1[comp_plots]
  ind2=comp_df$contrast2[comp_plots]
  cont1=afds[[ind1]]
  cont2=afds[[ind2]]
  results[[comp_plots]]=calcContrastStats(results=list(cont1, cont2),
                            L=paste0('_', names(afds)[ind1]), R=paste0('_', names(afds)[ind2]) ) #'_high1', R='_unsel1')
  
  results[[comp_plots]] <- na.omit(results[[comp_plots]])
  
  interval_df <- results[[comp_plots]] %>%
    dplyr::group_by(chrom)%>%
    # dplyr::filter(chrom=="II", physical.position < 5e6) %>%
    dplyr::filter(LOD > max(LOD)-2) %>%
    dplyr::mutate(lcon = min(physical.position),
                  rcon = max(physical.position)) %>%
    dplyr::arrange(desc(LOD))%>%
    dplyr::distinct(chrom, lcon, rcon, .keep_all = T) %>%
    dplyr::mutate(marker = paste0(chrom,":",lcon,"-",rcon)) %>%
    dplyr::select(marker, physical.position, LOD)
  
  suffix.1=gsub("ref_","",colnames(results[[comp_plots]])[6])
  suffix.2=gsub("ref_","",colnames(results[[comp_plots]])[17])
  
  sc=plotContrast(results[[comp_plots]], suffix1=suffix.1, suffix2=suffix.2) #S3_N2_XZ1516_7_KO')
  
  s=plotSummary(results[[comp_plots]], effective.n.tests=2000)
  
  pare1=str_split(suffix.1, pattern = "_")[[1]][2]
  pare2=str_split(suffix.2, pattern = "_")[[1]][3]
  cond1=str_split(suffix.1, pattern = "_")[[1]][5]
  cond2=str_split(suffix.2, pattern = "_")[[1]][5]
  gen1=str_split(suffix.1, pattern = "_")[[1]][4]
  gen2=str_split(suffix.2, pattern = "_")[[1]][4]
  
  write.table(interval_df, file = glue::glue('{plot.dir}{pare1}_{pare2}_F{gen1}-{gen2}_contrast_{cond1}-{cond2}_9000_rmv0_dp10.tsv'), col.names = T, row.names = F, quote = F, sep = "\t")
  
  a=ggarrange(plots[[ind1]],plots[[ind2]], s,  nrow=3)
  ggsave(filename = glue::glue('{plot.dir}{pare1}_{pare2}_F{gen1}-{gen2}_contrast_{cond1}-{cond2}_9000_rmv0_dp10.png'),
         plot = a, height = 16, width = 16)
}

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # ht115 - ko 
my_y_title <- expression(paste("-log10(", italic("p"), ")"))

colnames(results[[4]])

ses <- results[[4]] %>%
  dplyr::select(chrom:ID, rnai=afd.se_S8_N2_XZ1516_4_HT115,
                ko=afd.se_S12_N2_XZ1516_8_HT115) %>%
  tidyr::gather(expt,se,-(chrom:ID))

sup35 <- data.frame(chrom = "III",
                    physical.position = 11119226)

mll1 <- data.frame(chrom = "V",
                   physical.position = 20470485)

fqplot <- results[[4]] %>%
  dplyr::select(chrom:ID, ko=afd_S12_N2_XZ1516_8_HT115, rnai=afd_S8_N2_XZ1516_4_HT115) %>%
  tidyr::gather(expt,fq,-(chrom:ID))%>%
  dplyr::left_join(.,ses, by = c("chrom", "physical.position", "ID","expt")) %>%
  ggplot()+
  aes(x = physical.position/1e6, y = fq)+
  geom_vline(aes(xintercept = physical.position/1e6),data = sup35, color = "#D41159")+
  geom_vline(aes(xintercept = physical.position/1e6),data = mll1,color =  "#D41159")+
  geom_line(aes(color = expt))+
  scale_fill_manual(values = c("#FFC20A", "#0C7BDC"))+
  scale_color_manual(values = c("#FFC20A", "#0C7BDC"))+
  geom_ribbon(aes(ymin = fq - 1.96*se, ymax = fq + 1.96*se, fill = expt), alpha = 0.5)+
  facet_grid(.~chrom, space = "free", scales = "free")+
  ylim(0,1)+
  theme_bw(18)+
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        axis.text.x = element_blank())+
  labs(x = "Genomic Position (Mb)", y = "N2 allele frequency")

ggsave(filename = "/Volumes/HDbulkGWA/fog2_microPUB/fog2_v_HT115_f4-8.pdf", height = 4, width = 12)
ggsave(filename = "/Volumes/HDbulkGWA/fog2_microPUB/fog2_v_HT115_f4-8.png", height = 4, width = 12)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # rnai - ko 
my_y_title <- expression(paste("-log10(", italic("p"), ")"))

colnames(results[[1]])

ses <- results[[1]] %>%
  dplyr::select(chrom:ID, rnai=afd.se_S12_N2_XZ1516_6_RNAi,
                ko=afd.se_S1_N2_XZ1516_6_KO) %>%
  tidyr::gather(expt,se,-(chrom:ID))

sup35 <- data.frame(chrom = "III",
                    physical.position = 11119226)

mll1 <- data.frame(chrom = "V",
                   physical.position = 20470485)

fqplot <- results[[1]] %>%
  dplyr::select(chrom:ID, ko=afd_S1_N2_XZ1516_6_KO, rnai=afd_S12_N2_XZ1516_6_RNAi) %>%
  tidyr::gather(expt,fq,-(chrom:ID))%>%
  dplyr::left_join(.,ses, by = c("chrom", "physical.position", "ID","expt")) %>%
  ggplot()+
  aes(x = physical.position/1e6, y = fq)+
  geom_vline(aes(xintercept = physical.position/1e6),data = sup35, color = "#D41159")+
  geom_vline(aes(xintercept = physical.position/1e6),data = mll1,color =  "#D41159")+
  geom_line(aes(color = expt))+
  scale_fill_manual(values = c("#FFC20A", "#0C7BDC"))+
  scale_color_manual(values = c("#FFC20A", "#0C7BDC"))+
  geom_ribbon(aes(ymin = fq - 1.96*se, ymax = fq + 1.96*se, fill = expt), alpha = 0.5)+
  facet_grid(.~chrom, space = "free", scales = "free")+
  ylim(0,1)+
  theme_bw(18)+
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        axis.text.x = element_blank())+
  labs(x = "Genomic Position (Mb)", y = "N2 allele frequency")

ggsave(filename = "/Volumes/HDbulkGWA/fog2_microPUB/fog2_expt_f6.pdf", height = 4, width = 12)

lodplot <- results[[1]] %>%
  dplyr::select(chrom:ID, p=p) %>%
  ggplot()+
  aes(x = physical.position/1e6, y = -log10(p))+
  geom_vline(aes(xintercept = physical.position/1e6),data = sup35, color = "#D41159")+
  geom_vline(aes(xintercept = physical.position/1e6),data = mll1,color =  "#D41159")+
  geom_line()+
  scale_fill_manual(values = c("#FFC20A", "#0C7BDC"))+
  scale_color_manual(values = c("#FFC20A", "#0C7BDC"))+
  # geom_ribbon(aes(ymin = fq - 1.96*se, ymax = fq + 1.96*se, fill = expt), alpha = 0.5)+
  facet_grid(.~chrom, space = "free", scales = "free")+
  # ylim(0,1)+
  theme_bw(18)+
  scale_y_continuous(breaks = c(0,10,20))+
  theme(legend.position = "none",
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        strip.background = element_blank(),
        strip.text.x = element_blank())+
  labs(x = "Genomic Position (Mb)", y = "LOD")

ggsave(filename = "/Volumes/HDbulkGWA/fog2_microPUB/fog2_expt_f6_lod.pdf", height = 2, width = 12)

ggpubr::ggarrange(fqplot, lodplot, nrow = 2, heights = c(2,1), align = "v")

ggsave(filename = "/Volumes/HDbulkGWA/fog2_microPUB/fog2_expt_f6_lod_fq.pdf", height = 6, width = 12)

View(results[[2]])
colnames(results[[2]])

ses <- results[[2]] %>%
  dplyr::select(chrom:ID, rnai=afd.se_S2_N2_XZ1516_7_RNAi,
                ko=afd.se_S3_N2_XZ1516_7_KO) %>%
  tidyr::gather(expt,se,-(chrom:ID))

sup35 <- data.frame(chrom = "III",
                    physical.position = 11119226)

mll1 <- data.frame(chrom = "V",
                   physical.position = 20470485)

fqplot <- results[[2]] %>%
  dplyr::select(chrom:ID, ko=afd_S3_N2_XZ1516_7_KO, rnai=afd_S2_N2_XZ1516_7_RNAi) %>%
  tidyr::gather(expt,fq,-(chrom:ID))%>%
  dplyr::left_join(.,ses, by = c("chrom", "physical.position", "ID","expt")) %>%
  ggplot()+
  aes(x = physical.position/1e6, y = fq)+
  geom_vline(aes(xintercept = physical.position/1e6),data = sup35, color = "#D41159")+
  geom_vline(aes(xintercept = physical.position/1e6),data = mll1,color =  "#D41159")+
  geom_line(aes(color = expt))+
  scale_fill_manual(values = c("#FFC20A", "#0C7BDC"))+
  scale_color_manual(values = c("#FFC20A", "#0C7BDC"))+
  geom_ribbon(aes(ymin = fq - 1.96*se, ymax = fq + 1.96*se, fill = expt), alpha = 0.5)+
  facet_grid(.~chrom, space = "free", scales = "free")+
  ylim(0,1)+
  theme_bw(18)+
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        axis.text.x = element_blank())+
  labs(x = "Genomic Position (Mb)", y = "N2 allele frequency")

ggsave(fqplot, filename = "/Volumes/HDbulkGWA/fog2_microPUB/fog2_expt_f7_dp10.pdf", height = 4, width = 12)

lodplot <- results[[2]] %>%
  dplyr::select(chrom:ID, p=p) %>%
  ggplot()+
  aes(x = physical.position/1e6, y = -log10(p))+
  geom_vline(aes(xintercept = physical.position/1e6),data = sup35, color = "#D41159")+
  geom_vline(aes(xintercept = physical.position/1e6),data = mll1,color =  "#D41159")+
  geom_line()+
  scale_fill_manual(values = c("#FFC20A", "#0C7BDC"))+
  scale_color_manual(values = c("#FFC20A", "#0C7BDC"))+
  geom_hline(aes(yintercept = -log10(.05/2000)), color = "gray50", linetype = 2)+
  facet_grid(.~chrom, space = "free", scales = "free")+
  # ylim(0,1)+
  theme_bw(18)+
  scale_y_continuous(breaks = c(0,10,20))+
  theme(legend.position = "none",
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        strip.background = element_blank(),
        strip.text.x = element_blank())+
  labs(x = "Genomic Position (Mb)", y = my_y_title)

ggsave(lodplot, filename = "/Volumes/HDbulkGWA/fog2_microPUB/fog2_expt_f7_dp10_lod.pdf", height = 2, width = 12)

ggpubr::ggarrange(fqplot, lodplot, nrow = 2, heights = c(2,1), align = "v")

ggsave(filename = "/Volumes/HDbulkGWA/fog2_microPUB/fog2_expt_f7_dp10_lod_fq.pdf", height = 6, width = 12)

# compare earlier timepoints to address leonid feedback 
indx=10

colnames(results[[indx]])[grep("afd.se",colnames(results[[indx]]))[1]]
colnames(results[[indx]])[grep("afd.se",colnames(results[[indx]]))[2]]

ses <- results[[indx]] %>%
  dplyr::select(chrom:ID, early=colnames(results[[indx]])[grep("afd.se",colnames(results[[indx]]))[1]],
                late=colnames(results[[indx]])[grep("afd.se",colnames(results[[indx]]))[2]]) %>%
  tidyr::gather(expt,se,-(chrom:ID))

sup35 <- data.frame(chrom = "III",
                    physical.position = 11119226)

mll1 <- data.frame(chrom = "V",
                   physical.position = 20470485)

fqplot <- results[[indx]] %>%
  dplyr::select(chrom:ID, chrom:ID, early=colnames(results[[indx]])[grep("afd_",colnames(results[[indx]]))[1]],
                late=colnames(results[[indx]])[grep("afd_",colnames(results[[indx]]))[2]]) %>%
  tidyr::gather(expt,fq,-(chrom:ID))%>%
  dplyr::left_join(.,ses, by = c("chrom", "physical.position", "ID","expt")) %>%
  ggplot()+
  aes(x = physical.position/1e6, y = fq)+
  geom_vline(aes(xintercept = physical.position/1e6),data = sup35, color = "#D41159")+
  geom_vline(aes(xintercept = physical.position/1e6),data = mll1,color =  "#D41159")+
  geom_line(aes(color = expt))+
  scale_fill_manual(values = c("#FFC20A", "#0C7BDC"))+
  scale_color_manual(values = c("#FFC20A", "#0C7BDC"))+
  geom_ribbon(aes(ymin = fq - 1.96*se, ymax = fq + 1.96*se, fill = expt), alpha = 0.5)+
  facet_grid(.~chrom, space = "free", scales = "free")+
  ylim(0,1)+
  theme_bw(18)+
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        axis.text.x = element_blank())+
  labs(x = "Genomic Position (Mb)", y = "N2 allele frequency")

ggsave(fqplot, filename = "/Volumes/HDbulkGWA/fog2_microPUB/fog2_expt_f4-f7_dp10_KO.pdf", height = 4, width = 12)

lodplot <- results[[indx]] %>%
  dplyr::select(chrom:ID, p=p) %>%
  ggplot()+
  aes(x = physical.position/1e6, y = -log10(p))+
  geom_vline(aes(xintercept = physical.position/1e6),data = sup35, color = "#D41159")+
  geom_vline(aes(xintercept = physical.position/1e6),data = mll1,color =  "#D41159")+
  geom_line()+
  scale_fill_manual(values = c("#FFC20A", "#0C7BDC"))+
  scale_color_manual(values = c("#FFC20A", "#0C7BDC"))+
  geom_hline(aes(yintercept = -log10(.05/2000)), color = "gray50", linetype = 2)+
  facet_grid(.~chrom, space = "free", scales = "free")+
  # ylim(0,1)+
  theme_bw(18)+
  scale_y_continuous(breaks = c(0,10,20))+
  theme(legend.position = "none",
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        strip.background = element_blank(),
        strip.text.x = element_blank())+
  labs(x = "Genomic Position (Mb)", y = my_y_title)

ggsave(lodplot, filename = "/Volumes/HDbulkGWA/fog2_microPUB/fog2_expt_f4-f7_dp10_lod_ko.pdf", height = 2, width = 12)

combplot <- ggpubr::ggarrange(fqplot, lodplot, nrow = 2, heights = c(2,1), align = "v")

ggsave(combplot, filename = "/Volumes/HDbulkGWA/fog2_microPUB/fog2_expt_f4-f7_dp10_lod_fq_KO.pdf", height = 6, width = 12)

lods_ko <- results[[indx]] %>%
  dplyr::select(chrom:ID, p=p)

indx=9

colnames(results[[indx]])[grep("afd.se",colnames(results[[indx]]))[1]]
colnames(results[[indx]])[grep("afd.se",colnames(results[[indx]]))[2]]

ses <- results[[indx]] %>%
  dplyr::select(chrom:ID, early=colnames(results[[indx]])[grep("afd.se",colnames(results[[indx]]))[1]],
                late=colnames(results[[indx]])[grep("afd.se",colnames(results[[indx]]))[2]]) %>%
  tidyr::gather(expt,se,-(chrom:ID))

sup35 <- data.frame(chrom = "III",
                    physical.position = 11119226)

mll1 <- data.frame(chrom = "V",
                   physical.position = 20470485)

fqplot <- results[[indx]] %>%
  dplyr::select(chrom:ID, chrom:ID, early=colnames(results[[indx]])[grep("afd_",colnames(results[[indx]]))[1]],
                late=colnames(results[[indx]])[grep("afd_",colnames(results[[indx]]))[2]]) %>%
  tidyr::gather(expt,fq,-(chrom:ID))%>%
  dplyr::left_join(.,ses, by = c("chrom", "physical.position", "ID","expt")) %>%
  ggplot()+
  aes(x = physical.position/1e6, y = fq)+
  geom_vline(aes(xintercept = physical.position/1e6),data = sup35, color = "#D41159")+
  geom_vline(aes(xintercept = physical.position/1e6),data = mll1,color =  "#D41159")+
  geom_line(aes(color = expt))+
  scale_fill_manual(values = c("#FFC20A", "#0C7BDC"))+
  scale_color_manual(values = c("#FFC20A", "#0C7BDC"))+
  geom_ribbon(aes(ymin = fq - 1.96*se, ymax = fq + 1.96*se, fill = expt), alpha = 0.5)+
  facet_grid(.~chrom, space = "free", scales = "free")+
  ylim(0,1)+
  theme_bw(18)+
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        axis.text.x = element_blank())+
  labs(x = "Genomic Position (Mb)", y = "N2 allele frequency")

ggsave(fqplot, filename = "/Volumes/HDbulkGWA/fog2_microPUB/fog2_expt_f4-f7_dp10_RNAi.pdf", height = 4, width = 12)

lodplot <- results[[indx]] %>%
  dplyr::select(chrom:ID, p=p) %>%
  ggplot()+
  aes(x = physical.position/1e6, y = -log10(p))+
  geom_vline(aes(xintercept = physical.position/1e6),data = sup35, color = "#D41159")+
  geom_vline(aes(xintercept = physical.position/1e6),data = mll1,color =  "#D41159")+
  geom_line()+
  scale_fill_manual(values = c("#FFC20A", "#0C7BDC"))+
  scale_color_manual(values = c("#FFC20A", "#0C7BDC"))+
  geom_hline(aes(yintercept = -log10(.05/2000)), color = "gray50", linetype = 2)+
  facet_grid(.~chrom, space = "free", scales = "free")+
  # ylim(0,1)+
  theme_bw(18)+
  scale_y_continuous(breaks = c(0,10,20))+
  theme(legend.position = "none",
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        strip.background = element_blank(),
        strip.text.x = element_blank())+
  labs(x = "Genomic Position (Mb)", y = my_y_title)

ggsave(lodplot, filename = "/Volumes/HDbulkGWA/fog2_microPUB/fog2_expt_f4-f7_dp10_lod_RNAi.pdf", height = 2, width = 12)

combplot <- ggpubr::ggarrange(fqplot, lodplot, nrow = 2, heights = c(2,1), align = "v")

ggsave(combplot, filename = "/Volumes/HDbulkGWA/fog2_microPUB/fog2_expt_f4-f7_dp10_lod_fq_RNAi.pdf", height = 6, width = 12)

lods_rnai <- results[[indx]] %>%
  dplyr::select(chrom:ID, p=p)

lod_comp <- ggplot(lods_ko %>% dplyr::filter(ID%in%lods_rnai$ID))+
  aes(x = physical.position/1e6, y = -log10(p))+
  geom_vline(aes(xintercept = physical.position/1e6),data = sup35, color = "#D41159")+
  geom_vline(aes(xintercept = physical.position/1e6),data = mll1,color =  "#D41159")+
  geom_line(color = "#FFC20A")+
  geom_line(data = lods_rnai, color = "#0C7BDC")+
  # scale_fill_manual(values = c("#FFC20A", "#0C7BDC"))+
  # scale_color_manual(values = c("#FFC20A", "#0C7BDC"))+
  geom_hline(aes(yintercept = -log10(.05/2000)), color = "gray50", linetype = 2)+
  facet_grid(.~chrom, space = "free", scales = "free")+
  # ylim(0,1)+
  theme_bw(18)+
  scale_y_continuous(breaks = c(0,10,20))+
  theme(legend.position = "none",
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        strip.background = element_blank(),
        strip.text.x = element_blank())+
  labs(x = "Genomic Position (Mb)", y = my_y_title)

ggsave(lod_comp, filename = "/Volumes/HDbulkGWA/fog2_microPUB/fog2_expt_f4-f7_dp10_lod_KO-RNAi.pdf", height = 2, width = 12)

calcContrastStats_null=function(results, L='meta_high', R='meta_low'){
  
  combined.results = results %>% purrr::reduce(dplyr::left_join) %>% dplyr::ungroup() #, by=c('ID', 'chrom', 'physical.position') #%>% dplyr::ungroup()
  
  combined.results=combined.results %>% 
    dplyr::mutate(p= 2*pnorm(abs(z), lower.tail=F) ) %>%
    dplyr::mutate(LOD=PvalToLOD(p)) %>%suppressWarnings()
  
  return(combined.results)
}

# small
comp_df <- data.frame(contrast1=c(1,2,5,5,5,6,6,7,7,8),
                      contrast2=c(4,3,6,4,2,1,3,8,2,3))

#I have to restructure this stuff  ---------------------------------------------

results<-list()
for(comp_plots in 1:nrow(comp_df)){
ind1=comp_df$contrast1[comp_plots]
ind2=comp_df$contrast2[comp_plots]
cont1=afds[[ind1]]
cont2=afds[[ind2]]
results[[11]]=calcContrastStats(results=list(cont2),
                                        L=paste0('_', names(afds)[ind2]), R=NULL ) #'_high1', R='_unsel1')

combined.results = results %>% purrr::reduce(dplyr::left_join) %>% dplyr::ungroup() #, by=c('ID', 'chrom', 'physical.position') #%>% dplyr::ungroup()

combined.results=combined.results %>% 
  dplyr::mutate(p= 2*pnorm(abs(z), lower.tail=F) ) %>%
  dplyr::mutate(LOD=PvalToLOD(p)) %>%suppressWarnings()



