

#title: temporal dynamics of fecal microbiome and AMR in female pigs (breeding sows)

#load libaries
#library("devtools")
library("devtools")
library("dada2")
library("ShortRead")
library(phyloseq)
library(ggplot2)
library(dplyr)
library(plyr)
library(reshape)
library(reshape2)
library(scales)
library(dunn.test)
library(readr)
library(stringr)
library(phyloseq)
library(vegan)
library(tidyr)
library(data.table)
library(Maaslin2)
library(plotly)
library(synchrony)
library(rlist)
library(DirichletMultinomial)
library(ggpubr)
library(RColorBrewer)
library(jcolors)
library(microbiome)

# read microbiome -asv table, taxa table and metadata file stored in rds file
#phyloseq_age_sow_bacteriome=readRDS("./data/microbiome_16SrRNA/phyloseq_age_sow_bacteriome.rds")

## normalized microbial counts using CSS method
library(metagenomeSeq)
phyloseq_age_sow.metaseq <- phyloseq_to_metagenomeSeq(phyloseq_age_sow_bacteriome)
phyloseq_age_sow.metaseq.norm<- cumNorm(phyloseq_age_sow.metaseq, p=cumNormStat(phyloseq_age_sow.metaseq))

CSS_phyloseq_age_sow <- MRcounts(phyloseq_age_sow.metaseq.norm, norm = TRUE)
CSS_phyloseq_age_sow.ps <- merge_phyloseq(otu_table(CSS_phyloseq_age_sow, taxa_are_rows = T),
                                          sample_data(phyloseq_age_sow_bacteriome),
                                          tax_table(phyloseq_age_sow_bacteriome))

#sample_data_md=read.csv("./data/microbiome_16s/age_cohort_16s_metadata.csv", row.names = 1) #updated sample data
CSS_phyloseq_age_sow.ps <- merge_phyloseq(phyloseq::otu_table(CSS_phyloseq_age_sow.ps, taxa_are_rows = T),
                                          phyloseq::sample_data(sample_data_md),
                                          phyloseq::tax_table(CSS_phyloseq_age_sow.ps))

CSS_phyloseq_age_sow.ps <- readRDS("./data/microbiome_16SrRNA/CSS_phyloseq_age_sow.ps")

#-------------------------------------------------------------------------------------------------
# summarized number of raw reads distribution
#-------------------------------------------------------------------------------------------------
library(lme4)
library(nlme)
##Sample metadata
age_cohort_sow_md= read.csv("./data/microbiome_16SrRNA/age_cohort_16s_metadata.csv")
#track_dada2=read.csv("./data/track_dada2_sow.csv")
#age_cohort_sow_md= left_join(age_cohort_sow_md, track_dada2, by="S.Name")
boxplot(RawInput/1000 ~ factor(Age_days), data=age_cohort_sow_md)
lme_raw_reads16s <- lme(RawInput/1000 ~ factor(Age_week)+factor(Cohort), random = ~ 1|PigID, data=age_cohort_sow_md, method= "REML")
summary(lme_raw_reads16s)
lme_raw_reads16s_anova = anova(lme_raw_reads16s,type="marginal")
car::Anova(lme_raw_reads16s,type="III",test="F")
lme_raw_reads16s_anova
intervals(lme_raw_reads16s)
emmeans::emmeans(lme_raw_reads16s, pairwise ~ Age_week)
## facility/housing
lme_raw_reads16s <- lme(RawInput ~ Housing_facility, random = ~ 1|PigID, data=age_cohort_sow_md, method= "REML")
summary(lme_raw_reads16s)
lme_raw_reads16s_anova = anova(lme_raw_reads16s,type="marginal")
lme_raw_reads16s_anova
intervals(lme_raw_reads16s)

#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------
# DMM models
#-----------------------------------------------------------------------
library(phylosmith)
# filtered low abundant taxa
length(phyloseq::taxa_sums(CSS_phyloseq_age_sow.ps)[phyloseq::taxa_sums(CSS_phyloseq_age_sow.ps) == 0])
pseq_taxa_fil=taxa_filter(CSS_phyloseq_age_sow.ps, frequency = 0.05)
pseq.genus = tax_glom(pseq_taxa_fil, "Genus")

glomTax=tax_table(pseq.genus)[,"Genus"]
glomOTU=otu_table(pseq.genus) 
glomTable= merge(glomOTU, glomTax, by=0,all=TRUE)
rownames(glomTable)=glomTable[,"Genus"]
glomTable$Row.names <- NULL
glomTable$Genus <- NULL
dmm_input_genus_sow= as.matrix(t(glomTable))    

fit <- lapply(1:10, dmn, count = dmm_input_genus_sow , verbose=TRUE) #dmm_input_genus
best <- fit[[which.min(unlist(lplc))]]
mixture=mixturewt(best)
head(mixture)

assignment_all <- apply(mixture(best), 1, which.max)
assignment_all_dmm_df= as.data.frame (assignment_all)
xtabs(~assignment_all, assignment_all_dmm_df)

#write.csv(assignment_all_dmm_df, "assignment_all_dmm_df_sow.csv")

for (k in seq(ncol(fitted(best)))) {
  d <- melt(fitted(best))
  colnames(d) <- c("OTU", "cluster", "value")
  d <- subset(d, cluster == k) %>%
    # Arrange OTUs by assignment strength
    arrange(value) %>%
    mutate(OTU = factor(OTU, levels = unique(OTU)))
}  

for (k in seq(ncol(fitted(best)))) {
  d <- melt(fitted(best))
  colnames(d) <- c("OTU", "cluster", "value")
  d <- subset(d, cluster == k) %>%
    # Arrange OTUs by assignment strength
    arrange(value) %>%
    mutate(OTU = factor(OTU, levels = unique(OTU))) %>%
    # Only show the most important drivers
    filter(abs(value) > quantile(abs(value), 0.93))     
  
  p <- ggplot(d, aes(x = OTU, y = value)) +
    geom_bar(stat = "identity") +
    coord_flip() +
    labs(title = paste("Top drivers: community type", k))
  print(p)
}

#merged into metadata::age_cohort_16s_metadata.csv
#write.csv(assignment_all_dmm_df, "assignment_all_dmm_df_sow.csv")
age_cohort_sow_md= read.csv("./data/microbiome_16s/age_cohort_16s_metadata.csv")
age_cohort_sow_md_dmm=age_cohort_sow_md

# figure 1
#top 10 bacterial genera
ps1.div.fecal.fil=phylosmith::taxa_filter(CSS_phyloseq_age_sow.ps, frequency = 0.10) #taxa that are seen in 25% of the samples.
ps1.div.genus=tax_glom(CSS_phyloseq_age_sow.ps, "Genus")
sample_data(ps1.div.genus)$Age_week= as.factor(sample_data(ps1.div.genus1)$assignment_all_dmm_sow_rearrange)

OTUnames10 = names(sort(taxa_sums(ps1.div.genus), TRUE)[1:10])
GP10  = prune_taxa(OTUnames10,  ps1.div.genus)
mGP10 = prune_taxa(OTUnames10, ps1.div.genus)
sample_data(mGP10)$assignment_all_dmm_sow_rearrange= as.character(sample_data(mGP10)$assignment_all_dmm_sow_rearrange)
str(sample_data(mGP10)$assignment_all_dmm_sow_rearrange)
ps1.div.genus1.dmm<- merge_samples(mGP10, "assignment_all_dmm_sow_rearrange")
ps_rel_abund = phyloseq::transform_sample_counts(mGP10, function(x){x / sum(x)})
ps_rel_abund=psmelt(ps_rel_abund)

phy_relative_long_bi <- ps_rel_abund %>%
  group_by(Sample,Genus) %>%
  mutate(mean_relative_abund = Abundance/sum(Abundance))
grouped_df(phy_relative_long_bi, "assignment_all_dmm_sow_rearrange")

phy_relative_long_bi$Genus <- as.character(phy_relative_long_bi$Genus)

top_genera_dmm.plt <- 
  phy_relative_long_bi %>%
  ggplot(aes(x = assignment_all_dmm_sow_rearrange, y = Abundance,fill = reorder(Genus, -Abundance))) +
  geom_bar(stat = "identity",position = "fill" ) +
  theme(axis.ticks.x = element_blank())+
  theme(legend.key.size = unit(0.5, "cm"))+
  ggthemes::scale_fill_tableau("Classic 20",direction = -1) +  #Tableau 20
  scale_y_continuous(labels = scales::percent_format())+
  theme(axis.title.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        #axis.text.x = element_text(size = 8),
        #axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y = element_text(size = 14),
        legend.title = element_text(size = 14),
        legend.text = element_text( size = 14),
        plot.title = element_text(size=12, face = "bold", hjust = 0.5),
        panel.background = element_rect(fill = "white"),
        strip.background = element_rect(fill="white"),
        strip.text = element_text(face = "bold"),
        text=element_text(family="Times New Roman",  size=12))+
  guides(fill=guide_legend(title="Genus"))+xlab("Enterotype")
top_genera_dmm.plt

#----------------------------------------------------------------------
## dmm over time code is adapted Stewart, C.J., Ajami, N.J., O’Brien, J.L. et al. 
#Temporal development of the gut microbiome in early childhood from the TEDDY study. Nature 562, 583–588 (2018).
library(reshape2)
library(dplyr)
library(scales)
library(tidyverse)
library(ggplot2)
library(ggthemes)
library(Hmisc)
library(microbiome)
library(DirichletMultinomial)
library(magrittr)
library(qiime2R)
library(igraph)
library(slam)
library(lattice)
library(xtable)
#----------------------------------------------------------------------
sow_pig_dmm_md     <- age_cohort_sow_md_dmm               ## metadata after merging dmm cluster [assignment_dmm]
sow_pig_dmm_md$Age_week
SAMPLING_COLUMN     <- "Age_week"
TIME_POINTS     <- c(3,6,10,12,22,32,49,50,53)  #c(21,42,70,84,154,224,343,350,371) #  converted into week in final manuscript c(3,6,10,12,22,32,49,50,53) 
TIME_FUDGE      <- 0               # exact time points
REQUIRE_ALL_TP  <- FALSE           # required all time points if not set false
CLUSTER_COLUMN  <- "assignment_all_dmm_sow_rearrange"
SUBJECT_COLUMN  <- "PigID"
MINIMUM_PERCENT <- 0.001
NODE_COLOR      <- "#072b5a"
EDGE_COLOR      <- c("white", "orange", "darkorange", "orangered", "darkred")
EDGE_PCT_RANGE  <- c(0, 1)   


#----------------------------------------------------------------------
# check column names
#----------------------------------------------------------------------

vars <- c(SAMPLING_COLUMN, CLUSTER_COLUMN, SUBJECT_COLUMN)
if (!all(vars %in% colnames(sow_pig_dmm_md)))
  stop(sprintf("Column name not found: %s", paste(collapse=", ", setdiff(vars, colnames(sow_pig_dmm_md)))))
if (nrow(sow_pig_dmm_md) == 0)
  stop("Metadata object 'sow_pig_dmm_md' doesn't have any rows.")

clusterNames <- unique(as.character(sort(sow_pig_dmm_md[[CLUSTER_COLUMN]])))
timePtNames  <- unique(as.character(sort(TIME_POINTS)))
stateNames   <- apply(expand.grid(clusterNames, timePtNames), 1L, paste, collapse="@")
nStates      <- length(stateNames)

#----------------------------------------------------------------------
# Filter out rows missing data for time, animal or cluster
#----------------------------------------------------------------------

sow_pig_dmm_md <- sow_pig_dmm_md[!is.na(sow_pig_dmm_md[[SAMPLING_COLUMN]]) & !is.na(sow_pig_dmm_md[[CLUSTER_COLUMN]]) & !is.na(sow_pig_dmm_md[[SUBJECT_COLUMN]]), vars, drop=FALSE]
str(sow_pig_dmm_md)


#----------------------------------------------------------------------
# Map samples to time points
#----------------------------------------------------------------------

closestTP <- sapply(sow_pig_dmm_md[[SAMPLING_COLUMN]], function (x) {
  x <- TIME_POINTS[which.min(abs(x - TIME_POINTS))]
})
residuals <- abs(closestTP - sow_pig_dmm_md[[SAMPLING_COLUMN]])

sow_pig_dmm_md[[SAMPLING_COLUMN]] <- closestTP

if (TIME_FUDGE <  1) inRange <- residuals <= TIME_FUDGE * closestTP
if (TIME_FUDGE >= 1) inRange <- residuals <= TIME_FUDGE

sow_pig_dmm_md        <- sow_pig_dmm_md[inRange,,drop=FALSE]
residuals <- residuals[inRange]

sow_pig_dmm_md <- sow_pig_dmm_md[order(residuals),,drop=FALSE]
sow_pig_dmm_md <- sow_pig_dmm_md[!duplicated(paste(sow_pig_dmm_md[[SAMPLING_COLUMN]], 
                                                   sow_pig_dmm_md[[SUBJECT_COLUMN]])),,drop=FALSE]

#----------------------------------------------------------------------
# require a animal/ssample to have samples from all the time points
# In this case, make sure REQUIRE_ALL_TP is set to FALSE
#----------------------------------------------------------------------
if (REQUIRE_ALL_TP) {
  sow_pig_dmm_md <- plyr::ddply(sow_pig_dmm_md, SUBJECT_COLUMN, function (x) {
    if (nrow(x) == length(TIME_POINTS)) return (x)
    return (NULL)
  })
}
timePtCounts <- table(sow_pig_dmm_md[[SAMPLING_COLUMN]])
nodeCounts   <- unlist(as.list(table(paste(sep="@", sow_pig_dmm_md[[CLUSTER_COLUMN]], sow_pig_dmm_md[[SAMPLING_COLUMN]]))))
nodeCounts   <- setNames(nodeCounts / timePtCounts[sub("^.*@", "", names(nodeCounts))], names(nodeCounts))

#----------------------------------------------------------------------
# number of animal at each time point and/or cluster
#----------------------------------------------------------------------


table(sow_pig_dmm_md$Age_week, dnn = SAMPLING_COLUMN)
table(sow_pig_dmm_md$assignment_all_dmm_sow  ,   dnn = CLUSTER_COLUMN)
table(sow_pig_dmm_md$assignment_all_dmm_sow  , sow_pig_dmm_md$Age_week, dnn = c(CLUSTER_COLUMN, SAMPLING_COLUMN))

#----------------------------------------------------------------------
# number of each transition
#----------------------------------------------------------------------

transitionMatrix_sow <- matrix(0, nrow=nStates, ncol=nStates, dimnames=list(stateNames, stateNames))

sow_pig_dmm_md[[CLUSTER_COLUMN]] <- paste(sep="@", sow_pig_dmm_md[[CLUSTER_COLUMN]], sow_pig_dmm_md[[SAMPLING_COLUMN]])
sow_pig_dmm_md[[SAMPLING_COLUMN]]    <- as.numeric(factor(sow_pig_dmm_md[[SAMPLING_COLUMN]]))

plyr::d_ply(sow_pig_dmm_md[,vars], SUBJECT_COLUMN, function (x) {
  
  if (nrow(x) < 2) return (NULL)
  
  x <- x[order(x[[SAMPLING_COLUMN]]),,drop=FALSE]
  
  for (i in seq_len(nrow(x) - 1)) {
    
    t1 <- x[i,   SAMPLING_COLUMN]
    t2 <- x[i+1, SAMPLING_COLUMN]
    if (t2 - t1 > 1) next
    
    c1 <- as.character(x[i,   CLUSTER_COLUMN])
    c2 <- as.character(x[i+1, CLUSTER_COLUMN])
    transitionMatrix_sow[c1, c2] <<- transitionMatrix_sow[c1, c2] + 1
  }
})

#----------------------------------------------------------------------
# Rescale to percentages
#----------------------------------------------------------------------
indices <- sort(rep(1:length(timePtNames), length(clusterNames)))
for (timePt in 1:(length(timePtNames) - 1)) {
  i <- which(indices == timePt)
  transitionMatrix_sow[i,] <- transitionMatrix_sow[i,] / sum(transitionMatrix_sow[i,])
}

#----------------------------------------------------------------------
# Dropwith zero observations 
#----------------------------------------------------------------------

transitionMatrix_sow[which(transitionMatrix_sow < MINIMUM_PERCENT)] <- 0

stateNames   <- stateNames[rowSums(transitionMatrix_sow) | colSums(transitionMatrix_sow)]
clusterNames <- clusterNames[sapply(clusterNames, function (x) any(grep(sprintf("^%s@", x), stateNames)))]
timePtNames  <- timePtNames[ sapply(timePtNames,  function (x) any(grep(sprintf("@%s$", x), stateNames)))]

transitionMatrix_sow <- transitionMatrix_sow[stateNames, stateNames, drop=FALSE]

#----------------------------------------------------------------------
# fig
#----------------------------------------------------------------------
library(scales)
g <- graph.adjacency(transitionMatrix_sow, mode="upper", weighted=T)
#rescale(E(g)$weight)
layout <- matrix(unlist(strsplit(stateNames, "@", fixed=TRUE)), ncol=2, byrow=TRUE)
layout <- matrix(ncol=2, c(
  as.numeric(factor(layout[,2], levels=timePtNames)), 
  as.numeric(factor(layout[,1], levels=rev(clusterNames))) ))
vertex_weights <- rescale(as.vector(nodeCounts[names(V(g))]))
edge_weights   <- rescale((E(g)$weight), from=EDGE_PCT_RANGE)
edge_weights <- ifelse(edge_weights < 0.0, 0, edge_weights) #  ignore less than 0.1
edgeColorRamp <- colorRamp(EDGE_COLOR)

plot.igraph(
  x                = g,
  xlim             = c(0,1),
  ylim             = c(-1, 1),
  layout           = layout,
  vertex.label     = NA,
  vertex.size      = rescale(sqrt(vertex_weights / 3.14)) * 25,
  vertex.color     = c("#5ab4ac"),
  edge.width       = edge_weights*9 , 
  #vertex.shape="pie", 
  edge.color       = rgb(edgeColorRamp(edge_weights) / 255) )
xLabelPos = (seq_along(timePtNames)  - 1) *  2 / (length(timePtNames)  - 1) - 1
yLabelPos = (seq_along(clusterNames) - 1) * -2 / (length(clusterNames) - 1) + 1

text(-1.3, yLabelPos, clusterNames)
text(xLabelPos, 1.2, timePtNames)
#text(xLabelPos, 1.15, sprintf("n = %i", timePtCounts), cex=.5)

with(
  cbreaks(EDGE_PCT_RANGE, pretty_breaks(10), percent),
  legend(
    x      = 1.5, 
    y      = 1, 
    fill   = rgb(edgeColorRamp(rescale(breaks)) /255), 
    legend = labels, 
    title  = "Transition\nFrequency", 
    bty    = "n" ))
#eoffice::topptx(filename = "./dmm.pptx", width=8, height=6.5) 
#figure exported to powerpoint--for manual annotation time-point/age-week
#---------------------------------------------------------------------------------------------------
# microbiome maturation estimated as described by  https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8686833/
microbiome_maturation_results <- readRDS("./data/microbiome_16SrRNA/microbiome_maturation_results.rds")
microbiome_maturation_dmm.plt <- 
  ggplot(microbiome_maturation_results, aes(x=Age_week, y=weightedaverage, color=assignment_all_dmm_sow_rearrange))+ 
  geom_point(size = 2.5)+
  xlab("Age, week")+
  scale_color_distiller(type = "div")+ylab("Maturation Score")+
  stat_smooth(data=wiehtedaverage_results, 
              aes(x=Age_week, y = weightedaverage),
              method = "lm", 
              formula = y ~ x + I(x^2), 
              size = 1, 
              color="black", 
              level=0.95)+ #95% Confidence interval
  theme_classic()+ylim(1, NA)+guides(color=guide_legend(title="Enterotype"))

microbiome_maturation_dmm.plt

# upset plot common ASVs
library(UpSetR)
library("MicrobiotaProcess")
library(phyloseq)
phyloseq_upset_CSS_phyloseq_age_sow.ps=CSS_phyloseq_age_sow.ps
sample_data(phyloseq_upset_CSS_phyloseq_age_sow.ps)$sampling1 =as.factor(sample_data(phyloseq_upset_CSS_phyloseq_age_sow.ps)$sampling1)

sample_data(phyloseq_upset_CSS_phyloseq_age_sow.ps)$sampling1 <- revalue(sample_data(phyloseq_upset_CSS_phyloseq_age_sow.ps)$sampling1,
                                                                         c("Growing" = "3-22 week of age",
                                                                           "Breeding_Gestation"= "32 and 49 week of age",
                                                                           "Parturition" ="50 week of age",
                                                                           "Postpartum"  ="53 week of age"))
sampling1=unique(data.frame (sample_data(phyloseq_upset_CSS_phyloseq_age_sow.ps)$sampling1))
#library(ComplexUpset)
vennlist <- get_vennlist(obj=phyloseq_upset_CSS_phyloseq_age_sow.ps, factorNames="sampling1")
upsetda <- get_upset(obj=phyloseq_upset_CSS_phyloseq_age_sow.ps, factorNames="sampling1")
upset(upsetda, sets=c("53 week of age","50 week of age", "32 and 49 week of age","3-22 week of age"), keep.order = T, 
      set_size.show =FALSE ) # manually annotate y label 
#upset(upsetda,t(sampling1),width_ratio=0.1,set_sizes=FALSE)

#-----------------------------------------------------------
# figure 2- alpha and beta diversity
#-----------------------------------------------------------
CSS_phyloseq_age_sow.ps <- readRDS("./data/microbiome_16SrRNA/CSS_sow_merged.ps.rds")
head(sample_data(CSS_phyloseq_age_sow.ps))
phylum_sow.ps = tax_glom(CSS_phyloseq_age_sow.ps, "Phylum")
class_sow.ps = tax_glom(CSS_phyloseq_age_sow.ps, "Class")
genus_sow.ps = tax_glom(CSS_phyloseq_age_sow.ps, "Genus")
ASV_sow.ps1 = CSS_phyloseq_age_sow.ps

div_phy <- microbiome::alpha(phylum_sow.ps, index = c("observed", "diversity_shannon", "evenness_pielou"))
div_class <- microbiome::alpha(class_sow.ps, index = c("observed", "diversity_shannon", "evenness_pielou"))
div_genus <- microbiome::alpha(genus_sow.ps, index = c("observed", "diversity_shannon", "evenness_pielou"))
div_asv1 <- microbiome::alpha(ASV_sow.ps1, index = c("observed", "diversity_shannon", "evenness_pielou"))

ps1.meta.sow <- meta(CSS_phyloseq_age_sow.ps)
ps1.meta.sow$observed_phy <- div_phy$observed 
ps1.meta.sow$shannon_phy <- div_phy$diversity_shannon 
ps1.meta.sow$evenness_phy <- div_phy$evenness_pielou 

ps1.meta.sow$observed_class <- div_class$observed 
ps1.meta.sow$shannon_class <- div_class$diversity_shannon 
ps1.meta.sow$evenness_class <- div_class$evenness_pielou 

ps1.meta.sow$observed_genus <- div_genus$observed 
ps1.meta.sow$shannon_genus <- div_genus$diversity_shannon 
ps1.meta.sow$evenness_genus <- div_genus$evenness_pielou 

ps1.meta.sow$observed_asv1 <- div_asv1$observed 
ps1.meta.sow$shannon_asv1 <- div_asv1$diversity_shannon 
ps1.meta.sow$evenness_asv1 <- div_asv1$evenness_pielou 
#saveRDS(ps1.meta.sow, "ps1.meta.sow.rds")
div_md=readRDS("./data/microbiome_16SrRNA/ps1.meta.sow.rds")


richenss_div.plt <- 
  ggplot(div_md,aes(as.factor(Age_week),observed_asv1, group=PigID))+
  geom_line(color="grey") +
  geom_point(size=1.3)+
  ylim(0, NA)+
  #annotate("rect", xmin = c(0.9, 5.2, 7.8, 8.5), xmax = c(5, 7.7, 8.4, 9.2), 
  #         ymin = -Inf, ymax = Inf, alpha = .18, fill = c("#4E79A7", "#00A087B2", "red", "#B07AA1"))+
  theme_classic()+
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        text=element_text(family="Times New Roman",  size=14))+
  theme(legend.position = 'none')+ xlab("Age, week")+ylab("Richness")
richenss_div.plt

shannon_div.plt <- 
  ggplot(div_md,aes(as.factor(Age_week),shannon_asv1, group=PigID))+
  geom_line(color="grey") +
  geom_point(size=1.3)+
  ylim(2, NA)+
  #annotate("rect", xmin = c(0.9, 5.2, 7.8, 8.5), xmax = c(5, 7.7, 8.4, 9.2), 
  #         ymin = -Inf, ymax = Inf, alpha = .18, fill = c("#4E79A7", "#00A087B2", "red", "#B07AA1"))+
  theme_classic()+
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        text=element_text(family="Times New Roman",  size=14))+
  theme(legend.position = 'none')+ xlab("Age, week")+ylab("Shannon diversity")

#---------------------------------------------------------------------------------------------------
#beta-diversity plots-- figures 2A/B
CSS_phyloseq_age_sow.ps
metadata_asv_sow <- data.frame(sample_data(CSS_phyloseq_age_sow.ps))
#---------------------------------------------------------------------------------------------------
#NMDS
set.seed(1)
otu_2veg_asv_sow <- t(otu_table(CSS_phyloseq_age_sow.ps))
veg.dist.asv.sow <- vegdist(otu_2veg_asv_sow, method="bray")
ord_asv_sow <- metaMDS(veg.dist.asv.sow, k=2, trymax = 10000)
ord_asv_sow$stress  ## 0.10

data.scores = as.data.frame(scores(ord_asv_sow))
data.scores$Age_week = metadata_asv_sow$Age_week
data.scores$Age_week=as.factor(data.scores$Age_week)
cent_asv_sow <- aggregate(cbind(NMDS1, NMDS2) ~ Age_week,data = data.scores, FUN = mean)


segs_asv_sow <- merge(data.scores, setNames(cent_asv_sow, c('Age_week', 'oNMDS1','oNMDS2')),by = c('Age_week'),sort = FALSE)

nmds_all_samples.plt
ggplot(data.scores, aes(x = NMDS1, y = NMDS2)) +
  geom_point(size=2.5, aes(color=factor(Age_week))) + #pch = 21, alpha = 0.9
  theme_bw()+
  theme_minimal(base_size = 12)+
  geom_segment(data = segs_asv_sow,mapping = aes(xend = oNMDS1, yend = oNMDS2, color=factor(Age_week)), size = 0.08)+
  #stat_ellipse(aes(linetype = as.factor(Housing_facility.1), fill=as.factor(Housing_facility.1)),size = 0.5)+
  #geom_polygon(data = hulls, aes(color = sampling1), size=0.8, alpha=0.1) +
  scale_fill_manual(values = c("red", "blue", "black", "seagreen", "#e377c2","mediumpurple", "tan1","turquoise1", "brown")) +
  scale_color_manual(values = c("red", "blue", "black", "seagreen", "#e377c2","mediumpurple", "tan1","turquoise1", "brown")) +
  scale_colour_brewer(palette = "Set2")+
  scale_color_manual(values=c("blue","red","black", "#fc7d0b", "cyan4", "brown","mediumpurple", "darkgrey","#e377c2"))+ #"#00a2b3" #7b848f
  #ggsci::scale_color_rickandmorty()+
  #scale_color_manual(values = c("#4E79A7", "#00A087B2", "red", "#B07AA1"))+
  #ggthemes::scale_color_tableau("Tableau 20",direction = 1)+
  ggtitle("ASVs, stress=0.07")+
  # scale_fill_brewer(type="qual", palette="Set2")+
  theme(axis.text.x = element_text(size = 9),
        axis.text.y = element_text(size = 9),
        axis.text = element_text(size = 9),
        #text=element_text(family="Times New Roman",  size=11),
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 12))+
  theme(legend.key.height = unit(0.5, "cm"))+
  guides(colour = guide_legend(override.aes = list(size=4)))
nmds_all_samples.plt

#---------------------------------------------------------------------------------------------------
# only estrus to weaning/lacation 
#---------------------------------------------------------------------------------------------------
CSS_phyloseq_age_sow_estrus.ps=subset_samples(CSS_phyloseq_age_sow.ps, Age_week%in%c("32","49", "50", "53"))
metadata_asv_sow <- data.frame(sample_data(CSS_phyloseq_age_sow_estrus.ps))

metadata_asv_sow <- data.frame(sample_data(CSS_phyloseq_age_sow_estrus.ps))
otu_2veg_asv <- t(otu_table(CSS_phyloseq_age_sow_estrus.ps))
veg.dist.asv <- vegdist(otu_2veg_asv, method="bray")
metadata_asv_sow$assignment_all_dmm_sow=as.factor(metadata_asv_sow$assignment_all_dmm_sow)

# with adonis2
perm <- how(nperm = 999)
setBlocks(perm) <- with(metadata_asv_sow, PigID)
pair_wise.asv.age = adonis2(veg.dist.asv ~ as.factor(Age_week)+as.factor(Cohort),  data = metadata_asv_sow, permutations = perm)
pair_wise.asv.age

#NMDS
set.seed(1)
otu_2veg_asv_sow <- t(otu_table(CSS_phyloseq_age_sow_estrus.ps))
veg.dist.asv.sow <- vegdist(otu_2veg_asv_sow, method="bray")
ord_asv_sow <- metaMDS(veg.dist.asv.sow, k=2, trymax = 10000)
ord_asv_sow$stress  ## 0.10


data.scores = as.data.frame(scores(ord_asv_sow))
data.scores$Age_week = metadata_asv_sow$Age_week
data.scores$Age_week=as.factor(data.scores$Age_week)

cent_asv_sow <- aggregate(cbind(NMDS1, NMDS2) ~ Age_week,  data = data.scores, FUN = mean)
segs_asv_sow <- merge(data.scores, setNames(cent_asv_sow, c('Age_week','oNMDS1','oNMDS2')),by = c('Age_week' ),sort = FALSE)

nmds_estrus_wean<-                      
  ggplot(data.scores, aes(x = NMDS1, y = NMDS2)) +
  geom_point(size=3, aes(color=factor(Age_week))) + #pch = 21, alpha = 0.9
  theme_bw()+
  theme_minimal(base_size = 12)+
  stat_ellipse(geom = "polygon",aes(color = factor(Age_week)), alpha=0.0,level = 0.90, type = 't',size = 1.1)+
  scale_fill_manual(values = c("red", "blue", "black", "seagreen", "#e377c2","mediumpurple", "tan1","turquoise1", "brown")) +
  scale_color_manual(values = c("red", "blue", "black", "seagreen", "#e377c2","mediumpurple", "tan1","turquoise1", "brown")) +
  scale_colour_brewer(palette = "Set2")+
  scale_color_manual(values=c("blue","red","black", "#fc7d0b", "cyan4", "brown","mediumpurple", "grey60","#e377c2"))+ #"#00a2b3" #7b848f
  scale_color_manual(values=c("brown","mediumpurple", "darkgrey","#e377c2"))+ #"#00a2b3" #7b848f
  ggtitle("ASVs, stress=0.07")+
  # scale_fill_brewer(type="qual", palette="Set2")+
  theme(axis.text.x = element_text(size = 9),
        axis.text.y = element_text(size = 9),
        axis.text = element_text(size = 9),
        #text=element_text(family="Times New Roman",  size=11),
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 12))+
  theme(legend.key.height = unit(0.5, "cm"))+
  guides(colour = guide_legend(override.aes = list(size=4)))+
  guides(fill=guide_legend(title="Age, week"))
nmds_estrus_wean
#eoffice::topptx(filename = "./fig/beta_age_sow_asv6_edit.pptx", width=5.2, height=4.0) 


#---------------------------------------------------------------------------------------------------
#BC-distance boxplots
library("metagMisc")
#---------------------------------------------------------------------------------------------------
sample_data(CSS_phyloseq_age_sow.ps)$Age_week <- as.factor(sample_data(CSS_phyloseq_age_sow.ps)$Age_week)
sampling_bray_age = phyloseq_group_dissimilarity(CSS_phyloseq_age_sow.ps, group =c("Age_week"), between_groups =F, notch = T)
sampling_bray_age$data
sampling_bray_age_df= as.data.frame(sampling_bray_age$data)

# difference BC distances
summ_bc=glm(Dist~Group, data=sampling_bray_arg_df)
summary(summ_bc)
emmeans::emmeans(summ_bc, pairwise~Group, adjust="tukey")

#library(PupillometryR)
#devtools::install_github('erocoar/gghalves')
#library(gghalves)

bray_age_dist.plt <- ggplot(sampling_bray_age_df,aes(Group, y=Dist, fill=Group))+
  geom_boxplot(position= position_dodge(width = 0.5),  width=0.4, fatten = 2, alpha=0.5, notch =F)+
  ylim(0,NA)+
  geom_flat_violin(position = position_nudge(x = .25), alpha = .5)+
  geom_point(aes(color = Group), 
             position = position_jitter(width = .15),
             size = .3, 
             alpha = .5,
             show.legend = F)+
  scale_fill_manual(values=c("blue","red","black", "#fc7d0b", "cyan4", "brown","mediumpurple", "darkgrey","#e377c2"))+ #"#00a2b3" #7b848f
  scale_color_manual(values=c("blue","red","black", "#fc7d0b", "cyan4", "brown","mediumpurple", "darkgrey","#e377c2"))+ #"#00a2b3" #7b848f
  
  theme_minimal(base_size = 12)+ylim(0,1)+
  #theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  guides(fill=guide_legend(title="Sampling Point"))+
  theme(legend.position = 'none')+ ylab("Pairwise dissimiarity, Bray-Curtis")+xlab("Age, week")

bray_age_dist.plt  



#### fig3
#---------------------------------------------------------------------------------------------------
#RA plot --phylum level
#---------------------------------------------------------------------------------------------------
#ps1.div.phylum= metagMisc::phyloseq_filter_prevalence(CSS_phyloseq_age_sow.ps, 
#                                                      prev.trh = 0.2, abund.trh = 5, 
#                                                      abund.type = "mean", threshold_condition = "AND")  
ps1.div.phylum=tax_glom(CSS_phyloseq_age_sow.ps, "Phylum")
phylum_relative_sow <- transform_sample_counts(ps1.div.phylum, function(x) x / sum(x) )
phylum_relative_long_sow <- psmelt(phylum_relative_sow)
phylum_relative_long_sow <- phylum_relative_long_sow %>%
  group_by(Phylum) %>%
  mutate(mean_relative_abund = mean(Abundance))

unique(phylum_relative_long_sow$Sample)
sample=unique(sort(phylum_relative_long_sow$Sample))
sample=rev(sample)
phylum_relative_long_sow$Sample <- factor(phylum_relative_long_sow$Sample, levels=c(sample))

phylum_relative_long_sow$Phylum[phylum_relative_long_sow$Abundance < 0.02] <- "Other(< 2%)"  ## mean_relative_abund < 0.005
#fam_relative_long_sow$Sample=factor(fam_relative_long_sow$Sample, levels=c("3", "6", "10", "12", "22", "32", "49", "50", "53"))
phylum_relative_long_sow.plt<-phylum_relative_long_sow %>%
  ggplot(aes(x = Sample, y = Abundance, fill = reorder(Phylum, -Abundance))) +
  geom_bar(stat = "identity",width = 1) +
  #facet_wrap(~Age_week, scales = "free_x", nrow=1) +
  geom_bar(stat = "identity", alpha=0.9)+
  scale_y_continuous(expand = c(0, 0))+
  geom_vline(data = lines, aes(xintercept = vlines), color="white")+
  #geom_point(aes(x=Age_week, y=0))+
  theme(legend.key.size = unit(0.45, 'cm'))+
  ggthemes::scale_fill_tableau("Classic 20",direction = -1)+  #Tableau 20
  # scale_fill_manual(values = friendly_pal("bright_seven"))+
  theme(axis.title.x = element_text(size = 14),
        #axis.text.y = element_text(size = 14),
        #axis.text.x = element_text(size = 8),
        axis.text.y=element_blank(),
        #axis.ticks.y=element_blank(),
        axis.title.y = element_text(size = 14),
        legend.title = element_text(size = 14),
        legend.text = element_text( size = 15),
        plot.title = element_text(size=12, face = "bold", hjust = 0.5),
        panel.background = element_rect(fill = "white"),
        strip.background = element_rect(fill="white"),
        strip.text = element_text(face = "bold"),
        text=element_text(family="Times New Roman",  size=15))+
  guides(fill=guide_legend(title="Phylum"))+xlab("")+ coord_flip() 
phylum_relative_long_sow.plt # figure further manually annoted
#---------------------------------------------------------------------------------------------------
#microbial_genera_traject.plt
#---------------------------------------------------------------------------------------------------
maaslin_result_nonsign=read.csv("./data/microbiome_16SrRNA/maaslin_test_results_microbial_genera.csv")
maaslin_result_nonsign=subset(maaslin_result_nonsign,!trajec_sig%in%c("nonsign_desc", "nonsign_inc")) #"nonsign_desc", "nonsign_inc"
genera_maasling_nonsign= maaslin_result_nonsign$Genus
genus_timeseries_subset_nonsign=phyloseq::subset_taxa(CSS_phyloseq_age_sow_Genus.ps,Genus%in%genera_maasling_nonsign)
unique(tax_table(genus_timeseries_subset_nonsign)[,"Genus"])
OTUnames10_nonsign = names(sort(taxa_sums(genus_timeseries_subset_nonsign), TRUE)[1:6])

#GP10  = prune_taxa(OTUnames10_nonsign,  genus_timeseries_subset_nonsign)
#mGP10 = subset_taxa(genus_timeseries_subset_nonsign, OTUnames10_nonsign)
melt_genus_nonsign=psmelt(genus_timeseries_subset_nonsign)%>%
  group_by(Genus)

#top 5 phylum
# "Firmicutes"       "Bacteroidota"     "Spirochaetota"    "Proteobacteria"   "Actinobacteriota"
melt_genus_merge_nonsign=left_join(melt_genus_nonsign, maaslin_result_nonsign, by=c("Genus") )
melt_genus_merge_nonsign_subset=subset(melt_genus_merge_nonsign, Phylum%in%c("Firmicutes", "Bacteroidota", "Spirochaetota", "Proteobacteria","Actinobacteriota"))
melt_genus_merge_nonsign_subset$Abundance1=log10(melt_genus_merge_nonsign_subset$Abundance+1)
tgc_nonsign <- Rmisc::summarySE(melt_genus_merge_nonsign_subset, measurevar="Abundance1", groupvars=c("Genus","Age_week", "trajectories", "Phylum"))
tgc_nonsign
tgc_nonsign$Phylum=factor(tgc_nonsign$Phylum, levels=c("Firmicutes", "Bacteroidota",  "Actinobacteriota","Proteobacteria","Spirochaetota"))
tgc_nonsign$trajectories=factor(tgc_nonsign$trajectories, levels=c("increase", "descrease"))

microbial_genera_traject.plt <- 
  ggplot(tgc_nonsign, aes(x=as.factor(Age_week), y=Abundance1, color=trajectories, group=Genus)) + 
  geom_smooth(position = pd, se=F, size=0.5)+
  scale_color_manual(values=c("steelblue2", "darkgrey"))+
  facet_wrap(~Phylum, nrow=1)+ylab("Mean log10 abudnance")+
  theme_classic()+
  theme(axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        text=element_text(family="Times New Roman",  size=10),
        panel.background = element_blank(), panel.border = element_rect(fill = NA, colour = "grey30"))+ xlab("Age, week")
microbial_genera_traject.plt  
# manually annoted the taxa names 

#------------------------------------------  
# microbiome supplementary figures
#--------------------------------------------
#--------------------------------------------
#supp fig- model fit DMM
#--------------------------------------------
lplc <- base::sapply(fit, DirichletMultinomial::laplace) # AIC / BIC / Laplace
plot(lplc, type="l", xlab="Number of Dirichlet Components(DMM)", ylab="Model Fit",bty="n")
points(lplc, pch=21, cex=1.4, bg='grey80', col="black")
#--------------------------------------------
## DMM with envit
#--------------------------------------------
metadata_genus_sows <- data.frame(sample_data(genus_phyloseq_age_sow.ps))
otu_2veg_genus <- t(otu_table(genus_phyloseq_age_sow.ps))
veg.dist.genus <- vegdist(otu_2veg_genus, method="bray")
metadata_genus_sows$assignment_all_dmm_sow_rearrange=as.factor(metadata_genus_sows$assignment_all_dmm_sow_rearrange)

# with adonis2
perm <- how(nperm = 999)
setBlocks(perm) <- with(metadata_genus_sows, PigID)
pair_wise.genus.age = adonis2(veg.dist.genus ~ Age_week+Gender+Housing_facility,  data = metadata_genus_sows, permutations = perm)
pair_wise.genus.age
pairwise.adonis2(veg.dist.genus ~ as.factor(Age_week), data = metadata_genus_sows,p.adjust.m ='BH')

pair_wise.genus.housing = adonis2(veg.dist.genus ~  Housing_facility,  data = metadata_genus_sows, permutations = perm, type='margin')
pair_wise.genus.housing

set.seed(1)
hellin_genus.ps <- microbiome::transform(genus_phyloseq_age_sow.ps, "hellinger") ## hellinger
otu_2veg_genus_sows <- t(otu_table(genus_phyloseq_age_sow.ps))
veg.dist.genus.sows <- vegdist(otu_2veg_genus_sows, method="bray")
ord_genus_sows <- metaMDS(veg.dist.genus, k=2, trymax = 1000)
ord_genus_sows$stress  ## 0.09


data.scores = as.data.frame(scores(ord_genus_sows))
data.scores$Age_week = metadata_genus_sows$Age_week
data.scores$assignment_all_dmm_sow_rearrange = metadata_genus_sows$assignment_all_dmm_sow_rearrange


cent_genus_sows <- aggregate(cbind(NMDS1, NMDS2) ~ Age_week+assignment_all_dmm_sow_rearrange, data = data.scores, FUN = mean)
segs_genus_sows <- merge(data.scores, setNames(cent_genus_sows, c('Age_days',"assignment_all_dmm_sow_rearrange",
                                                                  'oNMDS1','oNMDS2')),
                         by = c('Age_days', "assignment_all_dmm_sow_rearrange"), sort = FALSE)


# ENV WITH dmm
sd2veg_dmm_class1= select(metadata_genus_sows, c("assignment_all_dmm_sow_rearrange"))
dune.spp.fit <- envfit(ord_genus_sows, otu_2veg_genus, permutations = 999)
head(dune.spp.fit)
dune.envfit <- envfit(ord_genus_sows, sd2veg_dmm_class1, permutations = 999) # this fits environmental vectors

site.scrs <- as.data.frame(scores(ord_genus_sows, display = "sites")) #save NMDS results into dataframe
site.scrs <- cbind(site.scrs, Age_days = metadata_genus_sows$Age_week) #add grouping variable "Management" to dataframe
head(site.scrs)
spp.scrs <- as.data.frame(scores(dune.spp.fit, display = "vectors")) #save species intrinsic values into dataframe
spp.scrs <- cbind(spp.scrs, class = rownames(spp.scrs)) #add species names to dataframe
spp.scrs <- cbind(spp.scrs, pval = dune.spp.fit$vectors$pvals) #add pvalues to dataframe so you can select species which are significant
spp.scrs<- cbind(spp.scrs, abrev = abbreviate(spp.scrs$class, minlength = 15)) #abbreviate species names
sig.spp.scrs <- subset(spp.scrs, pval<=0.05) #subset data to show species significant at 0.05
head(spp.scrs)


env.scores.dune <- as.data.frame(scores(dune.envfit, display = "factors")) #extracts relevant scores from envifit
env.scores.dune <- cbind(env.scores.dune, env.variables = rownames(env.scores.dune)) #and then gives them their names
env.scores.dune <- cbind(env.scores.dune, abrev = abbreviate(env.scores.dune$env.variables, minlength = 2)) #and then gives them their names

env.scores.dune <- cbind(env.scores.dune, pval = dune.envfit$factors$pvals) # add pvalues to dataframe
sig.env.scrs <- subset(env.scores.dune, pval<=0.05) #subset data to show variables significant at 0.05
head(env.scores.dune)

sig.env.scrs$env.variables <- revalue(sig.env.scrs$env.variables,
                                      c("assignment_all_dmm_sow_rearrange1"="DMM1",
                                        "assignment_all_dmm_sow_rearrange2"="DMM2",
                                        "assignment_all_dmm_sow_rearrange3"="DMM3",
                                        "assignment_all_dmm_sow_rearrange4"="DMM4"))
dmm_envit_age.plt<-
  nmds.plot.dune <- ggplot(data.scores, aes(x = NMDS1, y = NMDS2))+
  geom_point(aes(color=factor(Age_week)),size = 3)+ #adds site points 
  #ggtitle("Microbiome genus,stress=0.09")+
  scale_color_manual(values = c("red", "blue", "black", "seagreen", "#e377c2","mediumpurple", "tan1","turquoise1", "brown")) +
  scale_fill_manual(values = c("red", "blue", "black", "seagreen", "#e377c2","mediumpurple", "tan1","turquoise1", "brown")) +
  guides(fill = guide_legend(override.aes = list(size=3)))+
  theme_classic()+
  geom_segment(aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2), 
               data = env.scores.dune, size =0.5, alpha = 0.5, colour = "grey30")+
  guides(color=guide_legend(title="Age, weeks"))

nmds.plot.dune+geom_segment(data = sig.env.scrs, aes(x = 0, xend=NMDS1, y=0, yend=NMDS2), 
                            arrow = arrow(length = unit(0.1, "cm")), colour = "grey10", lwd=0.1) +
  ggrepel::geom_text_repel(data = sig.env.scrs, aes(x=NMDS1, y=NMDS2, label = env.variables), 
                           cex = 3, direction = "both", segment.size = 0.1)

dmm_envit_age.plt
#------------------------------------------------------------
#alpha diversity DMM
#---------------------------------------------------------
richness_asv_bydmm.plt= ggplot(div_md,aes(as.factor(assignment_all_dmm_sow_rearrange),observed_asv1))+
  geom_boxplot(fill='#3C5488B2', fatten = 1.4, outlier.color = NA,notch = F)+
  stat_boxplot(geom="errorbar", width=0.2) +
  geom_jitter(shape=21, fill='black', width=0.1, size=1.3, alpha=0.70)+
  theme_classic()+
  ylim(0,NA)+
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        text=element_text(family="Times New Roman",  size=14))+
  theme(legend.position = 'none')+ggtitle("Richness")+ylab("")+xlab("Enterotype")

shannon_asv_bydmm.plt= ggplot(div_md,aes(as.factor(assignment_all_dmm_sow_rearrange),shannon_asv1, fill=as.factor(assignment_all_dmm_sow_rearrange)))+
  stat_boxplot(geom="errorbar", width=0.2) +
  geom_boxplot(fatten = 1.4, outlier.color = NA,notch = F)+
  geom_jitter(shape=21, fill='black', width=0.1, size=1.3, alpha=0.70)+
  scale_fill_manual(values=c('#3C5488B2','#3C5488B2','#3C5488B2','#3C5488B2'))+
  theme_classic()+
  ylim(1,5)+
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        text=element_text(family="Times New Roman",  size=14))+
  theme(legend.position = 'none')+ggtitle("Shannon diversity")+ylab("")+xlab("Enterotype")




#------------------------------------------ 
# NMDS phylum level
#----------------------------------------
CSS_phyloseq_age_sow.ps <- readRDS("./data/microbiome_16SrRNA/CSS_phyloseq_age_sow.ps")
CSS_phyloseq_age_sow.fil= subset_samples(CSS_phyloseq_age_sow.ps, !SampleID%in%c("VV078.A6"))
phylum_phyloseq_age_sow.ps=tax_glom(CSS_phyloseq_age_sow.fil, "Phylum")
metadata_asv_sow_phylum <- data.frame(sample_data(phylum_phyloseq_age_sow.ps))
otu_2veg_asv_sow_phylum <- t(otu_table(phylum_phyloseq_age_sow.ps))
veg.dist.asv.sow.phylum <- vegdist(otu_2veg_asv_sow_phylum, method="bray")

perm <- how(nperm = 999)
setBlocks(perm) <- with(metadata_asv_sow_phylum, PigID)

pair_wise.asv.sampling1 = adonis2(veg.dist.asv.sow.phylum ~  sampling1+as.factor(Cohort),  data = metadata_asv_sow_phylum, permutations = perm, type='margin')
pair_wise.asv.sampling1

set.seed(1)
#hellin_asv.ps <- microbiome::transform(CSS_phyloseq_age_sow.ps, "hellinger") ## hellinger
otu_2veg_asv_sow_phylum <- t(otu_table(phylum_phyloseq_age_sow.ps))
veg.dist.asv.sow.phylum <- vegdist(otu_2veg_asv_sow_phylum, method="bray")
ord_asv_sow_phylum <- metaMDS(veg.dist.asv.sow.phylum, k=2, trymax = 10000)
ord_asv_sow_phylum$stress  ## 0.09

data.scores.phylum = as.data.frame(scores(ord_asv_sow_phylum))
data.scores.phylum$Age_week = metadata_asv_sow_phylum$Age_week
data.scores.phylum$sampling1 = metadata_asv_sow_phylum$sampling1


cent_asv_sow_phylum <- aggregate(cbind(NMDS1, NMDS2) ~ Age_week+sampling1, data = data.scores.phylum, FUN = mean)
segs_asv_sow_phylum <- merge(data.scores.phylum, setNames(cent_asv_sow_phylum, c('Age_week',"sampling1", 
                                                                                 'oNMDS1','oNMDS2')),by = c('Age_week',"sampling1"), sort = FALSE) 


phylum_nmds.plt <- 
  ggplot(segs_asv_sow_phylum, aes(x = NMDS1, y = NMDS2)) +
  geom_point(size=2.5, aes(color=factor(Age_week))) +
  theme_bw()+
  scale_fill_manual(values = c("red", "blue", "black", "seagreen", "#e377c2","mediumpurple", "tan1","turquoise1", "brown")) +
  scale_color_manual(values = c("#4E79A7", "#00A087B2", "red", "#B07AA1"))+
  #ggtitle("Phylum, stress=0.09")+
  #scale_fill_brewer(type="qual", palette="Set2")+
  theme(legend.key.height = unit(0.4, "cm"))+
  scale_color_manual(values=c("blue","red","black", "#fc7d0b", "cyan4", "brown","mediumpurple", "darkgrey","#e377c2"))+ #"#00a2b3" #7b848f
  theme(legend.key.height = unit(0.4, "cm"))+
  theme(axis.text.x = element_text(size = 9),
        axis.text.y = element_text(size = 9),
        axis.text = element_text(size = 9),
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 12))+
  theme(legend.key.height = unit(0.5, "cm"))+
  guides(colour = guide_legend(override.aes = list(size=4)))+
  guides(color=guide_legend(title="Age, week"))
phylum_nmds.plt  


# beta-dispersion
phylum_sampling_betadisper <- betadisper(veg.dist.asv.sow.phylum, metadata_asv_sow_phylum$Age_week,
                                         type = c("centroid"),bias.adjust=TRUE)
phylum_sampling_betadisper

phylum_sampling_betadisper_anova <- anova(betadisper(veg.dist.asv.sow.phylum, 
                                                     metadata_asv_sow_phylum$Age_week,
                                                     type = c("centroid"),bias.adjust=TRUE))
phylum_sampling_betadisper_anova 

## Permutation test for F
permutest(phylum_sampling_betadisper, permutations = 999)
TukeyHSD(phylum_sampling_betadisper)
beta_disper_phylum_df <- data.frame(phylum_sampling_betadisper$distances,phylum_sampling_betadisper$group)

distance_tocentroid_phylum.plt <- 
  ggboxplot(beta_disper_phylum_df, x = "phylum_sampling_betadisper.group", y = "phylum_sampling_betadisper.distances", 
            bxp.errorbar=T, width = 0.75, 
            color = "phylum_sampling_betadisper.group", notch = F)+
  scale_color_manual(values = c("#4E79A7", "#00A087B2", "red", "#B07AA1"))+
  scale_color_manual(values = c("red", "blue", "black", "seagreen", "#e377c2","mediumpurple", "tan1","turquoise1", "brown")) +
  scale_color_manual(values=c("blue","red","black", "#fc7d0b", "cyan4", "brown","mediumpurple", "darkgrey","#e377c2"))+ #"#00a2b3" #7b848f
  
  theme(axis.title.x = element_blank())+ theme_minimal()+
  theme(legend.position = "none")+
  theme(axis.text.x = element_blank())+
  ##theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))+
  ylab("Distance to centrioid")+xlab("")
distance_tocentroid_phylum.plt

#-------------------------------
# NMDS-class level
#-------------------------------
class_phyloseq_age_sow.ps=tax_glom(CSS_phyloseq_age_sow.fil, "Class")
metadata_asv_sow_class <- data.frame(sample_data(class_phyloseq_age_sow.ps))
otu_2veg_asv_sow_class <- t(otu_table(class_phyloseq_age_sow.ps))
veg.dist.asv.sow.class <- vegdist(otu_2veg_asv_sow_class, method="bray")

perm <- how(nperm = 999)
setBlocks(perm) <- with(metadata_asv_sow_class, PigID)

pair_wise.asv.sampling1 = adonis2(veg.dist.asv.sow.class ~  sampling1+as.factor(Cohort),  data = metadata_asv_sow_class, permutations = perm, type='margin')
pair_wise.asv.sampling1

set.seed(1)
#hellin_asv.ps <- microbiome::transform(CSS_phyloseq_age_sow.ps, "hellinger") ## hellinger
otu_2veg_asv_sow_class <- t(otu_table(class_phyloseq_age_sow.ps))
veg.dist.asv.sow.class <- vegdist(otu_2veg_asv_sow_class, method="bray")
ord_asv_sow_class <- metaMDS(veg.dist.asv.sow.class, k=2, trymax = 10000)
ord_asv_sow_class$stress  ## 0.09

data.scores.class = as.data.frame(scores(ord_asv_sow_class))
data.scores.class$Age_week = metadata_asv_sow_class$Age_week
data.scores.class$sampling1 = metadata_asv_sow_class$sampling1


cent_asv_sow_class <- aggregate(cbind(NMDS1, NMDS2) ~ Age_week+sampling1, data = data.scores.class, FUN = mean)
segs_asv_sow_class <- merge(data.scores.class, setNames(cent_asv_sow_class, c('Age_week',"sampling1", 
                                                                              'oNMDS1','oNMDS2')),by = c('Age_week',"sampling1"), sort = FALSE) 


class_nmds.plt <- 
  ggplot(segs_asv_sow_class, aes(x = NMDS1, y = NMDS2)) +
  geom_point(size=2.5, aes(color=factor(Age_week))) +
  theme_bw()+
  scale_fill_manual(values = c("red", "blue", "black", "seagreen", "#e377c2","mediumpurple", "tan1","turquoise1", "brown")) +
  scale_color_manual(values = c("#4E79A7", "#00A087B2", "red", "#B07AA1"))+
  #ggtitle("class, stress=0.09")+
  #scale_fill_brewer(type="qual", palette="Set2")+
  theme(legend.key.height = unit(0.4, "cm"))+
  scale_color_manual(values=c("blue","red","black", "#fc7d0b", "cyan4", "brown","mediumpurple", "darkgrey","#e377c2"))+ #"#00a2b3" #7b848f
  theme(legend.key.height = unit(0.4, "cm"))+
  theme(axis.text.x = element_text(size = 9),
        axis.text.y = element_text(size = 9),
        axis.text = element_text(size = 9),
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 12))+
  theme(legend.key.height = unit(0.5, "cm"))+
  guides(colour = guide_legend(override.aes = list(size=4)))+
  guides(color=guide_legend(title="Age, week"))
class_nmds.plt  


# beta-dispersion
class_sampling_betadisper <- betadisper(veg.dist.asv.sow.class, metadata_asv_sow_class$Age_week,
                                        type = c("centroid"),bias.adjust=TRUE)
class_sampling_betadisper

class_sampling_betadisper_anova <- anova(betadisper(veg.dist.asv.sow.class, 
                                                    metadata_asv_sow_class$Age_week,
                                                    type = c("centroid"),bias.adjust=TRUE))
class_sampling_betadisper_anova 

## Permutation test for F
permutest(class_sampling_betadisper, permutations = 999)
TukeyHSD(class_sampling_betadisper)
beta_disper_class_df <- data.frame(class_sampling_betadisper$distances,class_sampling_betadisper$group)

distance_tocentroid_class.plt <- 
  ggboxplot(beta_disper_class_df, x = "class_sampling_betadisper.group", y = "class_sampling_betadisper.distances", 
            bxp.errorbar=T, width = 0.75, 
            color = "class_sampling_betadisper.group", notch = F)+
  scale_color_manual(values = c("#4E79A7", "#00A087B2", "red", "#B07AA1"))+
  scale_color_manual(values = c("red", "blue", "black", "seagreen", "#e377c2","mediumpurple", "tan1","turquoise1", "brown")) +
  scale_color_manual(values=c("blue","red","black", "#fc7d0b", "cyan4", "brown","mediumpurple", "darkgrey","#e377c2"))+ #"#00a2b3" #7b848f
  
  theme(axis.title.x = element_blank())+ theme_minimal()+
  theme(legend.position = "none")+
  theme(axis.text.x = element_blank())+
  ##theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))+
  ylab("Distance to centrioid")+xlab("")
distance_tocentroid_class.plt

#-------------------------------
# NMDS-- Family level
#-------------------------------
family_phyloseq_age_sow.ps=tax_glom(CSS_phyloseq_age_sow.fil, "Family")
metadata_asv_sow_family <- data.frame(sample_data(family_phyloseq_age_sow.ps))
otu_2veg_asv_sow_family <- t(otu_table(family_phyloseq_age_sow.ps))
veg.dist.asv.sow.family <- vegdist(otu_2veg_asv_sow_family, method="bray")

perm <- how(nperm = 999)
setBlocks(perm) <- with(metadata_asv_sow_family, PigID)

pair_wise.asv.sampling1 = adonis2(veg.dist.asv.sow.family ~  sampling1+as.factor(Cohort),  data = metadata_asv_sow_family, permutations = perm, type='margin')
pair_wise.asv.sampling1

set.seed(1)
#hellin_asv.ps <- microbiome::transform(CSS_phyloseq_age_sow.ps, "hellinger") ## hellinger
otu_2veg_asv_sow_family <- t(otu_table(family_phyloseq_age_sow.ps))
veg.dist.asv.sow.family <- vegdist(otu_2veg_asv_sow_family, method="bray")
ord_asv_sow_family <- metaMDS(veg.dist.asv.sow.family, k=2, trymax = 10000)
ord_asv_sow_family$stress  ## 0.09

data.scores.family = as.data.frame(scores(ord_asv_sow_family))
data.scores.family$Age_week = metadata_asv_sow_family$Age_week
data.scores.family$sampling1 = metadata_asv_sow_family$sampling1


cent_asv_sow_family <- aggregate(cbind(NMDS1, NMDS2) ~ Age_week+sampling1, data = data.scores.family, FUN = mean)
segs_asv_sow_family <- merge(data.scores.family, setNames(cent_asv_sow_family, c('Age_week',"sampling1", 
                                                                                 'oNMDS1','oNMDS2')),by = c('Age_week',"sampling1"), sort = FALSE) 


family_nmds.plt <- 
  ggplot(segs_asv_sow_family, aes(x = NMDS1, y = NMDS2)) +
  geom_point(size=2.5, aes(color=factor(Age_week))) +
  theme_bw()+
  scale_fill_manual(values = c("red", "blue", "black", "seagreen", "#e377c2","mediumpurple", "tan1","turquoise1", "brown")) +
  scale_color_manual(values = c("#4E79A7", "#00A087B2", "red", "#B07AA1"))+
  #ggtitle("family, stress=0.09")+
  #scale_fill_brewer(type="qual", palette="Set2")+
  theme(legend.key.height = unit(0.4, "cm"))+
  scale_color_manual(values=c("blue","red","black", "#fc7d0b", "cyan4", "brown","mediumpurple", "darkgrey","#e377c2"))+ #"#00a2b3" #7b848f
  theme(legend.key.height = unit(0.4, "cm"))+
  theme(axis.text.x = element_text(size = 9),
        axis.text.y = element_text(size = 9),
        axis.text = element_text(size = 9),
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 12))+
  theme(legend.key.height = unit(0.5, "cm"))+
  guides(colour = guide_legend(override.aes = list(size=4)))+
  guides(color=guide_legend(title="Age, week"))
family_nmds.plt  


# beta-dispersion
family_sampling_betadisper <- betadisper(veg.dist.asv.sow.family, metadata_asv_sow_family$Age_week,
                                         type = c("centroid"),bias.adjust=TRUE)
family_sampling_betadisper

family_sampling_betadisper_anova <- anova(betadisper(veg.dist.asv.sow.family, 
                                                     metadata_asv_sow_family$Age_week,
                                                     type = c("centroid"),bias.adjust=TRUE))
family_sampling_betadisper_anova 

## Permutation test for F
permutest(family_sampling_betadisper, permutations = 999)
TukeyHSD(family_sampling_betadisper)
beta_disper_family_df <- data.frame(family_sampling_betadisper$distances,family_sampling_betadisper$group)

distance_tocentroid_family.plt <- 
  ggboxplot(beta_disper_family_df, x = "family_sampling_betadisper.group", y = "family_sampling_betadisper.distances", 
            bxp.errorbar=T, width = 0.75, 
            color = "family_sampling_betadisper.group", notch = F)+
  scale_color_manual(values = c("#4E79A7", "#00A087B2", "red", "#B07AA1"))+
  scale_color_manual(values = c("red", "blue", "black", "seagreen", "#e377c2","mediumpurple", "tan1","turquoise1", "brown")) +
  scale_color_manual(values=c("blue","red","black", "#fc7d0b", "cyan4", "brown","mediumpurple", "darkgrey","#e377c2"))+ #"#00a2b3" #7b848f
  
  theme(axis.title.x = element_blank())+ theme_minimal()+
  theme(legend.position = "none")+
  theme(axis.text.x = element_blank())+
  ##theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))+
  ylab("Distance to centrioid")+xlab("")
distance_tocentroid_family.plt
#-------------------------------  
#NMDS--Genus level
#-------------------------------
genus_phyloseq_age_sow.ps=tax_glom(CSS_phyloseq_age_sow.fil, "Genus")
metadata_asv_sow_genus <- data.frame(sample_data(genus_phyloseq_age_sow.ps))
otu_2veg_asv_sow_genus <- t(otu_table(genus_phyloseq_age_sow.ps))
veg.dist.asv.sow.genus <- vegdist(otu_2veg_asv_sow_genus, method="bray")

perm <- how(nperm = 999)
setBlocks(perm) <- with(metadata_asv_sow_genus, PigID)

pair_wise.asv.sampling1 = adonis2(veg.dist.asv.sow.genus ~  sampling1+as.factor(Cohort),  data = metadata_asv_sow_genus, permutations = perm, type='margin')
pair_wise.asv.sampling1

set.seed(1)
#hellin_asv.ps <- microbiome::transform(CSS_phyloseq_age_sow.ps, "hellinger") ## hellinger
otu_2veg_asv_sow_genus <- t(otu_table(genus_phyloseq_age_sow.ps))
veg.dist.asv.sow.genus <- vegdist(otu_2veg_asv_sow_genus, method="bray")
ord_asv_sow_genus <- metaMDS(veg.dist.asv.sow.genus, k=2, trymax = 10000)
ord_asv_sow_genus$stress  ## 0.09

data.scores.genus = as.data.frame(scores(ord_asv_sow_genus))
data.scores.genus$Age_week = metadata_asv_sow_genus$Age_week
data.scores.genus$sampling1 = metadata_asv_sow_genus$sampling1


cent_asv_sow_genus <- aggregate(cbind(NMDS1, NMDS2) ~ Age_week+sampling1, data = data.scores.genus, FUN = mean)
segs_asv_sow_genus <- merge(data.scores.genus, setNames(cent_asv_sow_genus, c('Age_week',"sampling1", 
                                                                              'oNMDS1','oNMDS2')),by = c('Age_week',"sampling1"), sort = FALSE) 


genus_nmds.plt <- 
  ggplot(segs_asv_sow_genus, aes(x = NMDS1, y = NMDS2)) +
  geom_point(size=2.5, aes(color=factor(Age_week))) +
  theme_bw()+
  scale_fill_manual(values = c("red", "blue", "black", "seagreen", "#e377c2","mediumpurple", "tan1","turquoise1", "brown")) +
  scale_color_manual(values = c("#4E79A7", "#00A087B2", "red", "#B07AA1"))+
  #ggtitle("genus, stress=0.09")+
  #scale_fill_brewer(type="qual", palette="Set2")+
  theme(legend.key.height = unit(0.4, "cm"))+
  scale_color_manual(values=c("blue","red","black", "#fc7d0b", "cyan4", "brown","mediumpurple", "darkgrey","#e377c2"))+ #"#00a2b3" #7b848f
  theme(legend.key.height = unit(0.4, "cm"))+
  theme(axis.text.x = element_text(size = 9),
        axis.text.y = element_text(size = 9),
        axis.text = element_text(size = 9),
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 12))+
  theme(legend.key.height = unit(0.5, "cm"))+
  guides(colour = guide_legend(override.aes = list(size=4)))+
  guides(color=guide_legend(title="Age, week"))
genus_nmds.plt  


# beta-dispersion
genus_sampling_betadisper <- betadisper(veg.dist.asv.sow.genus, metadata_asv_sow_genus$Age_week,
                                        type = c("centroid"),bias.adjust=TRUE)
genus_sampling_betadisper

genus_sampling_betadisper_anova <- anova(betadisper(veg.dist.asv.sow.genus, 
                                                    metadata_asv_sow_genus$Age_week,
                                                    type = c("centroid"),bias.adjust=TRUE))
genus_sampling_betadisper_anova 

## Permutation test for F
permutest(genus_sampling_betadisper, permutations = 999)
TukeyHSD(genus_sampling_betadisper)
beta_disper_genus_df <- data.frame(genus_sampling_betadisper$distances,genus_sampling_betadisper$group)

distance_tocentroid_genus.plt <- 
  ggboxplot(beta_disper_genus_df, x = "genus_sampling_betadisper.group", y = "genus_sampling_betadisper.distances", 
            bxp.errorbar=T, width = 0.75, 
            color = "genus_sampling_betadisper.group", notch = F)+
  scale_color_manual(values = c("#4E79A7", "#00A087B2", "red", "#B07AA1"))+
  scale_color_manual(values = c("red", "blue", "black", "seagreen", "#e377c2","mediumpurple", "tan1","turquoise1", "brown")) +
  scale_color_manual(values=c("blue","red","black", "#fc7d0b", "cyan4", "brown","mediumpurple", "darkgrey","#e377c2"))+ #"#00a2b3" #7b848f
  
  theme(axis.title.x = element_blank())+ theme_minimal()+
  theme(legend.position = "none")+
  theme(axis.text.x = element_blank())+
  ##theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))+
  ylab("Distance to centrioid")+xlab("")
distance_tocentroid_genus.plt

#--------------------------------------------
#genus level relative abundance plot 
#--------------------------------------------
CSS_phyloseq_age_sow_Genus.ps= tax_glom(CSS_phyloseq_age_sow.ps, "Genus")
table  <-CSS_phyloseq_age_sow_Genus.ps 
TRESHOLD = 0.005 # abundance treshold to say when a Genus is present
sow_ps_abundance_melt <- psmelt(table)
test <-  aggregate(Abundance~Genus, sow_ps_abundance_melt, sum)
test <- test[rev(order(test$Abundance)),]
Genera_tokeep <- test$Genus[test$Genus %in% head(test$Genus,15)]

CSS_phyloseq_age_sow_genus.fil=merge_samples(CSS_phyloseq_age_sow_Genus.ps, "Age_week")
OTUnames10 = names(sort(taxa_sums(CSS_phyloseq_age_sow_Genus.ps), TRUE)[1:15])
GP10  = prune_taxa(OTUnames10,  CSS_phyloseq_age_sow_Genus.ps)
mGP10 = subset_taxa(CSS_phyloseq_age_sow_Genus.ps, Genus%in%Genera_tokeep)
ps_rel_abund = phyloseq::transform_sample_counts(mGP10, function(x){x / sum(x)})
ps_rel_abund=psmelt(ps_rel_abund)

ps_rel_abund <- ps_rel_abund %>% group_by(Genus)
ps_rel_abund$Abundance[is.na(ps_rel_abund$Abundance)]<-0
ps_rel_abund <- aggregate(Abundance ~Genus+Age_week,FUN=sum,  data=ps_rel_abund)


lines <- data.frame(vlines = c(6,10,12,22, 32,49,50))
lines <- data.frame(vlines = c(32))
scale_x_continuous( expand = c(0, 0))+
  scale_y_continuous( expand = c(0, 0))

PAL_color1 <- rev(c("#4DBBD5B2", "#00A087B2", "#3C5488B2", "#F39B7FB2", "#8491B4B2", "#91D1C2B2", "#DC0000B2", "#7E6148B2",
                    '#DFE07C','#7F8E39','#3F4521','#42858C','#205062',
                    '#1D3554','#513D30','#70513A','#AF8A5F'))
show_col(PAL_color1)

scale_fill_color1 <- function(number_of_colors){
  structure(list(
    scale_fill_manual(values = get_palette(palette = PAL_color1 ,k=number_of_colors))   
  ))
}
ps_rel_genus_abund.plt <- 
  ggplot(ps_rel_abund, 
         aes(x=Age_week,y=as.numeric(Abundance),fill=reorder(Genus, -Abundance))) + 
  geom_area(position = 'fill')+
  scale_x_continuous( expand = c(0, 0))+
  scale_y_continuous( expand = c(0, 0))+
  ggthemes::scale_fill_tableau("Classic 20",direction = 1)+  #Tableau 20
  scale_fill_color1(length(unique(ps_rel_abund$Genus)))+
  #ggsci::scale_fill_material("deep-orange")+
  scale_x_continuous(breaks = seq(from = 3, to = 58, by = 10))+
  geom_vline(data = lines, aes(xintercept = vlines), color="grey")+
  #geom_point(aes(x=Age_week, y=0))+
  theme_minimal()+
  theme(legend.key.size = unit(0.45, 'cm'))+ ylab("Relative abundance")
theme(text=element_text(family="Times New Roman",  size=16))+xlab("Pig age in week")
ps_rel_genus_abund.plt

#-----------------------------------------------------------------------  
# mycobiome composition
#-----------------------------------------------------------------------    

fungal_CSS_phyloseq_age_sow.ps <- readRDS("./data/mycobiome/fungal_CSS_phyloseq_age_sow.ps")
head(sample_data(fungal_CSS_phyloseq_age_sow.ps)) 


# taxonomy -- RA plot at fungal genera
fungal_CSS_phyloseq_age_sow_genus.ps=tax_glom(fungal_CSS_phyloseq_age_sow.ps, "Genus")
fungal_CSS_phyloseq_age_sow_genus_merge=merge_samples(fungal_CSS_phyloseq_age_sow_genus.ps, "Age_week")
fungal_CSS_phyloseq_age_sow_genus_merge = phyloseq::transform_sample_counts(fungal_CSS_phyloseq_age_sow_genus_merge, function(x){x / sum(x)})
fungal_CSS_phyloseq_age_sow_genus_merge=psmelt(fungal_CSS_phyloseq_age_sow_genus_merge)
fungal_CSS_phyloseq_age_sow_genus_merge=fungal_CSS_phyloseq_age_sow_genus_merge%>%group_by(Genus)

fungal_CSS_phyloseq_age_sow_genus_merge$Sample=factor(fungal_CSS_phyloseq_age_sow_genus_merge$Sample, 
                                                      levels=c("3", "6", "10", "12", "22", "32", "49", "50", "53"))

fungal_taxa <- ggplot(fungal_CSS_phyloseq_age_sow_genus_merge, aes(Sample, Abundance, fill=Genus))+
  geom_bar(stat='identity', color='black', width=0.97)+
  #scale_fill_manual(values = PAL_LB1)+
  ggthemes::scale_fill_tableau("Tableau 20",direction = -2)+
  theme_minimal()+ ylab("Relative abundance")+
  theme(legend.text = element_text( size = 14))+ xlab("Age, week")
fungal_taxa  

#--------------------------------------------------------------------------------------------------- 
# fungal alpha diversity
library(knitr)
#---------------------------------------------------------------------------------------------------
tab_mycobiome <- microbiome::alpha(fungal_CSS_phyloseq_age_sow.ps, 
                                   index = c("observed","diversity_inverse_simpson", "diversity_shannon", "dominance_gini"))
kable(head(tab_mycobiome))
head(tab_mycobiome)

ps1_meta_mycobiome <- meta(fungal_CSS_phyloseq_age_sow.ps)
kable(head(ps1_meta_mycobiome))

ps1_meta_mycobiome$observed <- tab_mycobiome$observed 
ps1_meta_mycobiome$diversity_inverse_simpson <- tab_mycobiome$diversity_inverse_simpson
ps1_meta_mycobiome$diversity_shannon <- tab_mycobiome$diversity_shannon 
ps1_meta_mycobiome$dominance_gini <- tab_mycobiome$dominance_gini 

ps1_meta_mycobiome$diversity_inverse_simpson <- tab_mycobiome$diversity_inverse_simpson

boxplot(diversity_shannon ~ Age_week, data = ps1_meta_mycobiome, frame = FALSE, ylim=c(0,3), xlab = "Age, week", ylab = "")
title(main = "Shannon diversity")

#---------------------------------------------------------------------------------------------------
# mycobiome--supplementary figures
#---------------------------------------------------------------------------------------------------

# family level
fungal_CSS_phyloseq_age_sow_family.ps=tax_glom(fungal_CSS_phyloseq_age_sow.ps, "Family")
fungal_CSS_phyloseq_age_sow_family_merge=merge_samples(fungal_CSS_phyloseq_age_sow_family.ps, "Age_week")
fungal_CSS_phyloseq_age_sow_family_merge = phyloseq::transform_sample_counts(fungal_CSS_phyloseq_age_sow_family_merge, function(x){x / sum(x)})
fungal_CSS_phyloseq_age_sow_family_merge=psmelt(fungal_CSS_phyloseq_age_sow_family_merge)
fungal_CSS_phyloseq_age_sow_family_merge=fungal_CSS_phyloseq_age_sow_family_merge%>% group_by(Family)


fungal_CSS_phyloseq_age_sow_family_merge$Sample=factor(fungal_CSS_phyloseq_age_sow_family_merge$Sample, levels=
                                                         rev(c("3", "6", "10", "12", "22", "32", "49", "50", "53")))

pal_manual <- c("#4DBBD5B2", "#00A087B2", "#3C5488B2", "#F39B7FB2", "#8491B4B2", "#91D1C2B2", "#DC0000B2", "#7E6148B2",
                '#DFE07C','#7F8E39','#3F4521','#42858C','#205062',
                '#1D3554','#513D30','#70513A','#AF8A5F')
show_col(pal_manual)

fungal_family_taxa <- ggplot(fungal_CSS_phyloseq_age_sow_family_merge, aes(Sample, Abundance, fill=Family))+
  geom_bar(stat='identity', color='black', width=0.97)+
  scale_fill_manual(values = rev(pal_manual))+
  #ggthemes::scale_fill_tableau("Tableau 20",direction = -2)+
  theme_minimal()+xlab("Age, week")+
  theme(legend.text = element_text( size = 13))+ coord_flip()
fungal_family_taxa

# phylum
pal_manual1 <- rev(c('#DFE07C','#7F8E39','#3F4521','#42858C','#205062',
                     '#1D3554','#513D30','#70513A','#AF8A5F'))
show_col(pal_manual)
fungal_CSS_phyloseq_age_sow_phylum.ps=tax_glom(fungal_CSS_phyloseq_age_sow.ps, "Genus", NArm = F)
fungal_CSS_phyloseq_age_sow_phylum_merge = phyloseq::transform_sample_counts(fungal_CSS_phyloseq_age_sow_phylum.ps, function(x){x / sum(x)})
fungal_CSS_phyloseq_age_sow_phylum_merge=psmelt(fungal_CSS_phyloseq_age_sow_phylum_merge)
fungal_CSS_phyloseq_age_sow_phylum_merge=fungal_CSS_phyloseq_age_sow_phylum_merge%>%group_by(Phylum)

#fungal_CSS_phyloseq_age_sow_phylum_merge$Phylum

fungal_phylum.plt <- ggplot(fungal_CSS_phyloseq_age_sow_phylum_merge, aes(Sample, Abundance, fill=Phylum))+
  geom_bar(stat='identity', width=0.98)+
  #scale_fill_manual(values = pal_manual1)+
  theme_minimal()+
  facet_wrap(~Age_week, scales='free_x', nrow=1)+
  ggthemes::scale_fill_tableau("Color Blind",direction = -2)+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank())
#fungal_phylum.plt<- fungal_phylum.plt+ scale_fill_discrete(labels=c("Unclassified", "Ascomycota", "Basidiomycota"))

# fungal bc distance
fungal_CSS_phyloseq_age_sow_genus.ps=tax_glom(fungal_CSS_phyloseq_age_sow.ps, "Genus")
sample_data(fungal_CSS_phyloseq_age_sow_genus.ps)$Age_week=as.factor(sample_data(fungal_CSS_phyloseq_age_sow_genus.ps)$Age_week)

hellin_fungal_genus.ps <- microbiome::transform(fungal_CSS_phyloseq_age_sow_genus.ps, "hellinger") ## hellinger
hellin_fungal_genus.ps <- prune_samples(sample_sums(hellin_fungal_genus.ps)>0, hellin_fungal_genus.ps)  ## number
#fix(hellin_fungal_genus_melt)
hellin_fungal_genus_melt=psmelt(hellin_fungal_genus.ps)
#write.csv(hellin_fungal_genus_melt, "./data/hellin_fungal_genus_melt.csv")
#VV131
hellin_fungal_genus.ps=subset_samples(hellin_fungal_genus.ps,!SampleID%in%c("VV131.E11"))
metadata_fungal_genus_sows= data.frame(sample_data(hellin_fungal_genus.ps))
otu_2veg_fungal_genus_sows <- t(otu_table(hellin_fungal_genus.ps))
veg.dist.fungal_genus.sows <- vegdist(otu_2veg_fungal_genus_sows, method="bray")


fungal_genus_sampling_betadisper <- betadisper(veg.dist.fungal_genus.sows, 
                                               metadata_fungal_genus_sows$Age_week,type = c("centroid"),bias.adjust=TRUE)
fungal_genus_sampling_betadisper

fungal_genus_sampling_betadisper_anova <- anova(betadisper(veg.dist.fungal_genus.sows, 
                                                           metadata_fungal_genus_sows$Age_week,
                                                           type = c("centroid"),bias.adjust=TRUE))
fungal_genus_sampling_betadisper_anova # 

df_fungal_betadisp <- data.frame(fungal_genus_sampling_betadisper$distances,fungal_genus_sampling_betadisper$group)

df_fungal_betadisp <-
  ggboxplot(df_fungal_betadisp, x = "fungal_genus_sampling_betadisper.group", y = "fungal_genus_sampling_betadisper.distances", bxp.errorbar=T, width = 0.75, 
            color = "fungal_genus_sampling_betadisper.group", notch = F)+
  scale_color_manual(values = c("#4E79A7", "#00A087B2", "red", "#B07AA1"))+
  scale_color_manual(values=c("blue","red","black", "#fc7d0b", "cyan4", "brown","mediumpurple", "darkgrey","#e377c2"))+ #"#00a2b3" #7b848f
  theme(axis.title.x = element_blank())+ theme_minimal()+
  theme(legend.position = "none")+
  ##theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))+
  ylab("Distance to centriod, fungal genera composition")+xlab("Age, week")
df_fungal_betadisp

#-----------------------------------------------------------------------------
# AMR data
#-----------------------------------------------------------------------------
library(lme4)
library(lmerTest)
#-----------------------------------------------------------------------------
target_amr_sow= read.csv("./data/phenotypic_genotypic_amr/targeted_amrgenes_breeder.csv")
target_amr_sow_melt=melt(target_amr_sow, id=c("PigID", "cohort",  "Age_weeks", "sampling1"))
target_amr_sow_melt$sampling1=factor(target_amr_sow_melt$sampling1, levels=c("Growing","Breeding_Gestation","Farrowing", "Weaning"))

#------------------------------
# abundance of tetA gene 
#--------------------------------
target_amr_sow_melt_non_stand_teta=subset(target_amr_sow_melt, variable%in%c("log10tet"))
model_lmer_amr_genes= lmer(as.numeric(value)~as.factor(Age_weeks)+(cohort)+(1| PigID), data=target_amr_sow_melt_non_stand_teta)
summary(model_lmer_amr_genes)
anova(model_lmer_amr_genes)
performance::icc(model_lmer_amr_genes)
emmeans::emmeans(model_lmer_amr_genes,  pairwise ~ Age_weeks)
d.tet1 <- summary(emmeans::emmeans(model_lmer_amr_genes, ~Age_weeks))
ggplot(d.tet1, aes(as.factor(Age_weeks))) + 
  geom_point(aes(y = emmean), size = 1.5, shape = 21, fill='black',position = position_dodge(width = 0.5)) +
  geom_line(aes(y = emmean,group=1)) + #linetype = "twodash"
  geom_jitter(data=target_amr_sow_melt_non_stand_teta, aes(y=as.numeric(value)),size = 1, alpha=0.4)+
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), linetype = "solid", width = 0.19, position = position_dodge(width = 0.5)) +
  scale_color_brewer( palette="Set1")+
  theme_classic()+ylab("Log10 gene copy/gram")+ xlab("Age, week")
#-----------------------------------
# abundance of standardize tetA gene
#-----------------------------------
target_amr_sow_melt_standteta=subset(target_amr_sow_melt, variable%in%c("standardize.tetA..log.10.")) 
model_lmer_amr_genes= lmer(as.numeric(value)~as.factor(Age_weeks)+(cohort)+(1|PigID), data=target_amr_sow_melt_standteta)
summary(model_lmer_amr_genes)
anova(model_lmer_amr_genes)
performance::icc(model_lmer_amr_genes)
emmeans::emmeans(model_lmer_amr_genes,pairwise ~ Age_weeks)
d.tet2 <- summary(emmeans::emmeans(model_lmer_amr_genes, ~Age_weeks))  
ggplot(d.tet2, aes(as.factor(Age_weeks))) + 
  geom_point(aes(y = emmean), size = 1.5, shape = 21, fill='black',position = position_dodge(width = 0.5)) +
  geom_line(aes(y = emmean,group=1),linetype = "dashed" ) + #linetype = "twodash"
  geom_jitter(data=target_amr_sow_melt_standteta, aes(y=as.numeric(value)),size = 1, alpha=0.4)+
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), linetype = "solid", width = 0.19, position = position_dodge(width = 0.5)) +
  scale_color_brewer( palette="Set1")+
  ylim(-6,NA)+
  theme_classic()+ylab("Log10 gene copy/gram (standardized")+ xlab("Age, week")

#------------------------------
# abundance of blactxm gene
#--------------------------------
target_amr_sow_melt_non_stand_ctxm=subset(target_amr_sow_melt, variable%in%c("log10ctxm"))
model_lmer_amr_genes= lmer(as.numeric(value)~as.factor(Age_weeks)+(cohort)+(1| PigID), data=target_amr_sow_melt_non_stand_ctxm)
summary(model_lmer_amr_genes)
anova(model_lmer_amr_genes)
performance::icc(model_lmer_amr_genes)
emmeans::emmeans(model_lmer_amr_genes, pairwise~ Age_weeks)
d.blactxm <- summary(emmeans::emmeans(model_lmer_amr_genes, ~Age_weeks))   

ggplot(d.blactxm, aes(as.factor(Age_weeks))) + 
  geom_point(aes(y = emmean), size = 1.5, shape = 21, fill='black',position = position_dodge(width = 0.5)) +
  geom_line(aes(y = emmean,group=1),linetype = "solid" ) + #linetype = "twodash"
  geom_jitter(data=target_amr_sow_melt_non_stand_ctxm, aes(y=as.numeric(value)),size = 1, alpha=0.4)+
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), linetype = "solid", width = 0.19, position = position_dodge(width = 0.5)) +
  scale_color_brewer( palette="Set1")+
  ylim(4,6)+
  theme_classic()+ylab("Log10 gene copy/gram")+ xlab("Age, week")

#-------------------------------------
# abundance of standardized blactxm gene
#-------------------------------------
target_amr_sow_melt_stand_ctxm=subset(target_amr_sow_melt, variable%in%c("stand.ctxm.log.10."))
model_lmer_amr_genes= lmer(as.numeric(value)~as.factor(Age_weeks)+(cohort)+(1| PigID), data=target_amr_sow_melt_stand_ctxm)
summary(model_lmer_amr_genes)
anova(model_lmer_amr_genes)
emmeans::emmeans(model_lmer_amr_genes, pairwise~ Age_weeks)
d.standblactxm <- summary(emmeans::emmeans(model_lmer_amr_genes, ~Age_weeks))   

ggplot(d.standblactxm, aes(as.factor(Age_weeks))) + 
  geom_point(aes(y = emmean), size = 1.5, shape = 21, fill='black',position = position_dodge(width = 0.5)) +
  geom_line(aes(y = emmean,group=1),linetype = "dashed" ) + #linetype = "twodash"
  geom_jitter(data=target_amr_sow_melt_stand_ctxm, aes(y=as.numeric(value)),size = 1, alpha=0.4)+
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), linetype = "solid", width = 0.19, position = position_dodge(width = 0.5)) +
  scale_color_brewer( palette="Set1")+ylab("blactxm_standarizw")+
  theme_classic()+ylab("Log10 gene copy/gram (standardized)")+ xlab("Age, week")
#-----------------------
#16s rRNA gene copies
#-----------------------
ggplot(target_amr_sow, aes(as.factor(Age_weeks), y=as.numeric(log10.16s)))+
  #geom_point(size=1)+
  geom_line(aes(group=PigID), color='grey')+
  geom_point(size=1, alpha=0.55)+ #color='#3C5488B2',
  scale_color_manual(values=c("black"))+
  geom_vline(xintercept = 31, colour = "grey", size=0.5, linetype='dashed')+
  ylab("Log10 16S rRNA gene copies per gram")+ xlab("Age, week")+
  ylim(1,NA)+
  theme_classic()+theme(legend.position = 'none')
#--------------------------------------------------------------------------------------------------------------------------
#phenotypic AMR data
pheno_all_fil <- read.csv("./data/phenotypic_amr_genes_data/breeder_sow_amr_data.csv")

pheno_allbreed_ent<- subset(pheno_all_fil, group%in%c("ENT"))
pheno_allbreed_ec<- subset(pheno_all_fil, group%in%c("MAC"))
#------------------------------------
#E.coli-Tetracyclines
#------------------------------------
Tetracyclines_EC=subset(pheno_allbreed_ec, variable%in%c("Tetracyclines"))
str(Tetracyclines_EC)
Tetracyclines_EC$pigID=as.factor(Tetracyclines_EC$pigID)
Tetracyclines_EC$gender_cohort=as.factor(Tetracyclines_EC$gender_cohort)

Tetracyclines_EC$group
Tetracyclines_EC$variable
basic_model_tet <- mgcv::gam(as.numeric(value) ~ s(age_week,bs='cr', k=5)+
                               as.factor(gender_cohort),
                             correlation = corARMA(form = ~ 1|age_week, p = 1),
                             data = Tetracyclines_EC, 
                             method = "REML")
basic_summary_tet <- summary(basic_model_tet)
basic_summary_tet
#------------------------------------
# E.coli-AminoPenicillins
#------------------------------------
ec_pen_md=subset(pheno_allbreed_ec, variable%in%c("Penicillins"))
str(ec_pen_md)
ec_pen_md$pigID=as.factor(ec_pen_md$pigID)
ec_pen_md$gender_cohort=as.factor(ec_pen_md$gender_cohort)
ec_pen_model <- mgcv::gam(as.numeric(value) ~ s(age_week,bs='cr', k=5)+
                            as.factor(gender_cohort),
                          correlation = corARMA(form = ~ 1|age_week, p = 1),
                          data = ec_pen_md, 
                          method = "REML")
ec_pen_model_summary <- summary(ec_pen_model)
#------------------------------------
# E.coli- Beta-lactam/3rd gen cephalosporins
#------------------------------------
ec_betalactam_md=subset(pheno_allbreed_ec, variable%in%c("Beta-lactam"))
str(ec_betalactam_md)
ec_betalactam_md$pigID=as.factor(ec_betalactam_md$pigID)
ec_betalactam_md$gender_cohort=as.factor(ec_betalactam_md$gender_cohort)
ec_betalactam_model <- mgcv::gam(as.numeric(value) ~ s(age_week,bs='cr', k=5)+
                                   as.factor(gender_cohort),
                                 correlation = corARMA(form = ~ 1|age_week, p = 1),
                                 data = ec_betalactam_md, 
                                 method = "REML")
ec_betalactam_model_summary <- summary(ec_betalactam_model)
#------------------------------------
#E.coli-Macrolides
#------------------------------------
ec_macrolids_md=subset(pheno_allbreed_ec, variable%in%c("Macrolids"))
str(ec_macrolids_md)
ec_macrolids_md$pigID=as.factor(ec_macrolids_md$pigID)
ec_macrolids_md$gender_cohort=as.factor(ec_macrolids_md$gender_cohort)
ec_macrolids_model <- mgcv::gam(as.numeric(value) ~ s(age_week,bs='cr', k=5)+
                                  as.factor(gender_cohort),
                                correlation = corARMA(form = ~ 1|age_week, p = 1),
                                data = ec_macrolids_md, 
                                method = "REML")
ec_macrolids_model_summary <- summary(ec_macrolids_model)
ec_macrolids_model_summary
#------------------------------------
# E.coli-Aminoglycosides 
#------------------------------------
ec_aminoglycoside_md=subset(pheno_allbreed_ec, variable%in%c("Aminoglycosides"))
str(ec_aminoglycoside_md)
ec_aminoglycoside_md$pigID=as.factor(ec_aminoglycoside_md$pigID)
ec_aminoglycoside_md$gender_cohort=as.factor(ec_aminoglycoside_md$gender_cohort)
ec_aminoglycoside_model <- mgcv::gam(as.numeric(value) ~ s(age_week,bs='cr', k=5)+
                                       as.factor(gender_cohort),
                                     correlation = corARMA(form = ~ 1|age_week, p = 1),
                                     data = ec_aminoglycoside_md, 
                                     method = "REML")
ec_aminoglycoside_model_summary <- summary(ec_aminoglycoside_model)
ec_aminoglycoside_model_summary

#------------------------------------
#E.coli - Sulphonamides
#------------------------------------
ec_sulphonamide_md=subset(pheno_allbreed_ec, variable%in%c("Sulphonamides"))
str(ec_sulphonamide_md)
ec_sulphonamide_md$pigID=as.factor(ec_sulphonamide_md$pigID)
ec_sulphonamide_md$gender_cohort=as.factor(ec_sulphonamide_md$gender_cohort)
ec_sulphonamide_model <- mgcv::gam(as.numeric(value) ~ s(age_week,bs='cr', k=5)+
                                     as.factor(gender_cohort),
                                   correlation = corARMA(form = ~ 1|age_week, p = 1),
                                   data = ec_sulphonamide_md, 
                                   method = "REML")
ec_sulphonamide_model_summary <- summary(ec_sulphonamide_model)
ec_sulphonamide_model_summary
#------------------------------------
#E.coli --Phenicols
#------------------------------------
ec_phenicol_md=subset(pheno_allbreed_ec, variable%in%c("Phenicols"))
ec_phenicol_md$pigID=as.factor(ec_phenicol_md$pigID)
ec_phenicol_md$gender_cohort=as.factor(ec_phenicol_md$gender_cohort)
ec_phenicol_model <- mgcv::gam(as.numeric(value) ~ s(age_week,bs='cr', k=5)+
                                 as.factor(gender_cohort),
                               correlation = corARMA(form = ~ 1|age_week, p = 1),
                               data = ec_phenicol_md, 
                               method = "REML")
ec_phenicol_model_summary <- summary(ec_phenicol_model)
ec_phenicol_model_summary

#--------------------------
#-----------------------------
# enterococcus - Penicillins
#------------------------------------
ent_pen_md=subset(pheno_allbreed_ent, variable%in%c("Penicillins"))
ent_pen_md$pigID=as.factor(ent_pen_md$pigID)
ent_pen_md$gender_cohort=as.factor(ent_pen_md$gender_cohort)
ent_pen_model <- mgcv::gam(as.numeric(value) ~ s(age_week,bs='cr', k=5),
                           correlation = corARMA(form = ~ 1|age_week, p = 1),
                           data = ent_pen_md, 
                           method = "REML")
ent_pen_model_summary <- summary(ent_pen_model)
ent_pen_model_summary

#-------------------------------------
#enterococcus-Quinolones/Fluoroquinolone
#--------------------------------------
ent_quino_md=subset(pheno_allbreed_ent, variable%in%c("Quinolones/Fluoroquinolone"))
str(ent_quino_md)
ent_quino_md$pigID=as.factor(ent_quino_md$pigID)
ent_quino_md$gender_cohort=as.factor(ent_quino_md$gender_cohort)
ent_quino_model <- mgcv::gam(as.numeric(value) ~ s(age_week,bs='cr', k=5)+
                               gender_cohort,
                             correlation = corARMA(form = ~ 1|age_week, p = 1),
                             data = ent_quino_md, 
                             method = "REML")
ent_quino_model_summary <- summary(ent_quino_model)
ent_quino_model_summary

#-------------------------------------
#enterococcus-Tetracycline
#--------------------------------------
ent_tet_md=subset(pheno_allbreed_ent, variable%in%c("Tetracyclines"))
str(ent_tet_md)
ent_tet_md$pigID=as.factor(ent_tet_md$pigID)
ent_tet_md$gender_cohort=as.factor(ent_tet_md$gender_cohort)
ent_tet_model <- mgcv::gam(as.numeric(value) ~ s(age_week,bs='cr', k=5)+
                             gender_cohort,
                           correlation = corARMA(form = ~ 1|age_week, p = 1),
                           data = ent_tet_md, 
                           method = "REML")
ent_tet_model_summary <- summary(ent_tet_model)
#------------------------------------
# enterococcus -- Macrolides
#------------------------------------
ent_macrolide_md=subset(pheno_allbreed_ent, variable%in%c("Macrolids"))
ent_macrolide_md$pigID=as.factor(ent_macrolide_md$pigID)
ent_macrolide_md$gender_cohort=as.factor(ent_macrolide_md$gender_cohort)
ent_macrolide_model <- mgcv::gam(as.numeric(value) ~ s(age_week,bs='cr', k=5)+
                                   gender_cohort,
                                 correlation = corARMA(form = ~ 1|age_week, p = 1),
                                 data = ent_macrolide_md, 
                                 method = "REML")
ent_macrolide_model_summary <- summary(ent_macrolide_model)
ent_macrolide_model_summary

#-------------------------------------
#enterococcus-aminoglycosides
#--------------------------------------
ent_aminoglycoside_md=subset(pheno_allbreed_ent, variable%in%c("Aminoglycosides"))
str(ent_aminoglycoside_md)
ent_aminoglycoside_md$pigID=as.factor(ent_aminoglycoside_md$pigID)
ent_aminoglycoside_md$gender_cohort=as.factor(ent_aminoglycoside_md$gender_cohort)
ent_aminoglycoside_model <- mgcv::gam(as.numeric(value) ~ s(age_week,bs='cr', k=5)+
                                        gender_cohort,
                                      correlation = corARMA(form = ~ 1|age_week, p = 1),
                                      data = ent_aminoglycoside_md, 
                                      method = "REML")
ent_aminoglycoside_model_summary <- summary(ent_aminoglycoside_model)
ent_aminoglycoside_model_summary
#------------------------------------
#enterococcus -Lincosamides
#------------------------------------
ent_lincosamide_md=subset(pheno_allbreed_ent, variable%in%c("Lincosamides"))
ent_lincosamide_md$pigID=as.factor(ent_lincosamide_md$pigID)
ent_lincosamide_md$gender_cohort=as.factor(ent_lincosamide_md$gender_cohort)
ent_lincosamide_model <- mgcv::gam(as.numeric(value) ~ s(age_week,bs='cr', k=5)+
                                     gender_cohort,
                                   correlation = corARMA(form = ~ 1|age_week, p = 1),
                                   data = ent_lincosamide_md, 
                                   method = "REML")
ent_lincosamide_model_summary <- summary(ent_lincosamide_model)
ent_lincosamide_model_summary
#------------------------------------
#enterococcus--Nitrofurans
#------------------------------------
ent_nitrofuran_md=subset(pheno_allbreed_ent, variable%in%c("Nitrofurans"))
ent_nitrofuran_md$pigID=as.factor(ent_nitrofuran_md$pigID)
ent_nitrofuran_md$gender_cohort=as.factor(ent_nitrofuran_md$gender_cohort)
ent_nitrofuran_model <- mgcv::gam(as.numeric(value) ~ s(age_week,bs='cr', k=5)+
                                    gender_cohort,
                                  correlation = corARMA(form = ~ 1|age_week, p = 1),
                                  data = ent_nitrofuran_md, 
                                  method = "REML")
ent_nitrofuran_model_summary <- summary(ent_nitrofuran_model)
ent_nitrofuran_model_summary
#-----------------------------------------------------------------------------------------------------