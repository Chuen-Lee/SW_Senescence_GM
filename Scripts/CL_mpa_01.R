# MPA script ----
setwd("/Users/senmon/Documents/PhD/R_analysis/Metagenomics/")
## libraries ----
library(sjPlot)
library(DHARMa)
library(car)
library(microbiome)
library(tidyverse)
library(lme4)
library(speedyseq)
library(ggthemes)
library(mgcv)
library(ggeffects)
library(tidyverse)
library(hms)
library(gllvm)
library(vegan)
library(ANCOMBC)
library(MuMIn)
library(lmerTest)
library(ggpubr)
## functions ----
vegan_otu <- function(physeq) {
  OTU <- otu_table(physeq)
  if (taxa_are_rows(OTU)) {
    OTU <- t(OTU)
  }
  return(as(OTU, "matrix"))
}

alpha_estimate <- function(physeq) {
  ttt <- transform_sample_counts(physeq, function(x) round(x, 0))
  richnessEstRare<-estimate_richness(ttt, split=TRUE, measures= c("Chao1", "Shannon", "Observed"))
  RareMeta <- as.data.frame(sample_data(ttt))
  RareMeta$Chao1 <- richnessEstRare$Chao1
  RareMeta$Shannon <- richnessEstRare$Shannon
  RareMeta$Observed <- richnessEstRare$Observed
  AlphaMeta <- as_tibble(RareMeta) %>% mutate(survive = as.factor(survive)) %>% mutate(TerminalYear = as.factor(TerminalYear)) %>% mutate(TQscale = scales::rescale(TQ, to =c(0,1))) %>% mutate(Timeinfridge = as.numeric(Timeinfridge)) %>% 
    mutate(CatchTime = tidyr::replace_na(CatchTime, 43292)) %>% mutate(CatchTime = as.numeric(CatchTime))
  return(AlphaMeta)
}

## Loading raw files ----
# load functions and metadata

# metaphlan4
mpamat <- read.csv("InputTables/merged_reads_abundance.txt", sep = "\t") %>% pivot_wider(names_from = TubeNo, values_from = Reads) %>% mutate(OTU = paste0("OTU", 1:nrow(.)))
# taxa table
mpataxa <- mpamat %>% dplyr::select(OTU,OUT) %>% separate_wider_delim(OUT, delim = "|", names = c("Kingdom","Phylum","Order","Class","Family","Genus","Species","SGB"), too_few = "align_start") %>% column_to_rownames(., var="OTU")
mpataxo7 <- mpataxa %>% filter(!is.na(SGB))
# otu table
mpaL7 <- mpamat  %>% dplyr::select(-OUT) %>% column_to_rownames(., var="OTU") %>% replace(is.na(.), 0)
mpaL7 <- mpaL7[rownames(mpaL7) %in% row.names(mpataxo7), ]

# sample data
st <- read.csv("InputTables/st.csv")
st4 <- st %>% mutate(SexEstimate = as.factor(SexEstimate),SampleYear=as.factor(SampleYear),season=as.factor(season), TerminalYear=as.factor(TerminalYear)) %>% column_to_rownames("TubeNumber") 
stdata <- sample_data(st4)

mpaOTUL7 <- otu_table(as.matrix(mpaL7), taxa_are_rows = TRUE)

mpaTAXL7 = tax_table(as.matrix(mpataxo7))

mpapL7 <- phyloseq(mpaOTUL7, mpaTAXL7, stdata) %>% filter_sample_data(SamplingAge > 0.5) %>% filter_sample_data(Type %in% "F") %>% mutate_sample_data(Age_scaled = arm::rescale(SamplingAge,"center"),Timeinfridge_scaled = arm::rescale(Timeinfridge, "center"),CatchTime_scaled = arm::rescale(CatchTime, "center"),TQcorrected_scaled = arm::rescale(TQcorrected, "center")) 

mpapL7.bac <- subset_taxa(mpapL7, Kingdom %in% "k__Bacteria") %>% filter_sample_data(Type%in%"F") %>% filter_sample_data(!(TubeNo %in% c("SW1168","SW1183")))
mpapL7samplesums <- data.frame(sample_sums(mpapL7.bac))


## How many MPA per sample?
mpapL7.bacdf <- psmelt(mpapL7.bac) %>% filter(Abundance > 0) #

mpapersamp <- mpapL7.bacdf %>% group_by(Sample) %>% summarise(count=n())

quantile(mpapersamp$count)
mean(mpapersamp$count)
sd(mpapersamp$count)/sqrt(length(mpapersamp$count))
ggplot(mpapersamp,aes(count))+geom_boxplot()

mpapersamp2 <- sd(mpapersamp$count)/sqrt(length(mpapersamp$count))

mpaunique <- mpapL7.bacdf %>% distinct(OTU,.keep_all = T) # 1025 different ASVs

#### How many terminal year?
mpapL7.bacst <- data.frame(sample_data(mpapL7.bac)) %>% group_by(TerminalYear) %>% distinct(BirdID,.keep_all = T) %>% summarise(count=n())

#### Terminal year mean stats
data.frame(sample_data(mpapL7.bac)) %>% group_by(TerminalYear) %>% summarise(mean=mean(SamplingAge),sd=sd(SamplingAge))

#### How age range and mean and se
mpal7age <- mpapL7.bacdf %>% distinct(Sample,.keep_all = T)

quantile(mpal7age$SamplingAge)
mean(mpal7age$SamplingAge)
sd(mpal7age$SamplingAge)/sqrt(length(mpal7age$SamplingAge))

#### main phylums
mpapL7.bacdf %>% group_by(Phylum) %>% summarise(n = n()) %>% arrange(n)
mpapL7.bacdf %>% summarise(n = n())

# plot of age distribution of samples


#tiff("Output/samplingageplot.tiff", res =200, units = "in", width = 9, height = 12, bg = "white")
data.frame(sample_data(mpapL7.bac)) %>%
  ggplot(., aes(x=SamplingAge,y=as.factor(BirdID), colour=SexEstimate, shape = TerminalYear)) + 
  geom_line(aes(group=BirdID)) + 
  geom_point() +
  ylab("Bird ID") +
  xlab("Sampling Age") + 
  scale_colour_colorblind(name = "Sex") +
  theme_tufte(base_size = 10, base_family = "Arial") + 
  theme(axis.line = element_line(colour = "black", linetype=1))  
#dev.off()

#### storage at 4C
data.frame(sample_data(mpapL7.bac)) %>% summarise(min4 = min(Timeinfridge),max4 = max(Timeinfridge), median(Timeinfridge),mean(Timeinfridge), sd(Timeinfridge),se=sd(Timeinfridge)/sqrt(n()))

## Alpha diversity ----
# ggiNEXT - takes a long time to run
# library(iNEXT)
# library(readr)
### rarefaction curve ###
#### using iNEXT ####

# #make the ASV abundance table into a matrix and check by printing first 2 rows/columns
# abund <- as(otu_table(mpapL7), "matrix")
# abund[1:2,1:2] #check
#
# #convert to a dataframe
# abund2 <- as.data.frame(abund)
# str(abund2)
#
# #iNEXT only takes numeric values, so change all values in the dataframe to numeric values instead of integers.
# df2 <- mutate_all(abund2, function(x) as.numeric(x))  %>% mutate(OTU = rownames(.)) %>% dplyr::select(OTU, everything())
# str(df2)
# write_csv(df2,"Output/MFFiNEXTmpa.csv")
# 
# df2 <- read.csv("Output/MFFiNEXTmpa.csv", row.names = 1, sep = ",")
# df3 <- df2 %>% dplyr::select(-SW1168,-SW1183) %>% .[,1:20]
# inext_test<-iNEXT(df3, q=0, datatype="abundance", endpoint=20000)
# dput(inext, "inextTrial.txt")
# 
# inext7 <- dget("InputTables/inextTrial.txt")
# 
# inext7check <- inext7[["iNextEst"]][["size_based"]]
# 
# #plot rarefaction curve
# rarefaction<- ggiNEXT(inext7, type=1, se=FALSE, grey= TRUE) + theme(legend.position = "none")+ xlab("Sequencing depth") + ylab("Observed Taxa") 
# rf7 <- rarefaction + geom_vline(xintercept=5500, alpha=0.5, linetype=2) +scale_shape_manual(values=rep(20,164))
# 
# #Plot sample completeness curve
# completeness<-ggiNEXT(inext7, type=2, se=TRUE, grey=TRUE)+scale_shape_manual(values=rep(20,164))+ theme(legend.position = "none") +xlab("Read count") +ylab("Sample completeness")
# c7 <-  completeness + geom_vline(xintercept=5500, alpha=0.5, linetype=2) + geom_hline(yintercept = 0.95, linetype=2)
# 
# ggiNEXT(inext7, type=3) + theme(legend.position = "none")
# 
# inextarrange7 <- ggarrange(rf7,c7)

#rarefy
mpapL7.bacrare <- rarefy_even_depth(mpapL7.bac, rngseed=88,sample.size = 5500, replace=F)
mpal7alpha <- alpha_estimate(mpapL7.bacrare)

mpaalphameta4 <- mpal7alpha %>% group_by(BirdID) %>% mutate(capage=pmin(SamplingAge,12),deltaAge = capage - mean(capage), meanAge = mean(capage), SamplingAge2 = round(SamplingAge),sampling_id = 1:n(), termageint = interaction(Age_scaled, TerminalYear, sep = ":"), TerminalYearbird = case_when(any(TerminalYear %in% "1") ~ "1", TRUE ~ "0"))


### capped aged (12 years) ----
#### species richeness and capped aged ----
capageglm2 <- glmer.nb(Observed ~ capage +TerminalYear + season + SexEstimate+ Timeinfridge_scaled + CatchTime_scaled + TQcorrected_scaled + SampleYear + (1|BirdID), data=mpaalphameta4,nAGQ=0 )
summary(capageglm2)
car::Anova(capageglm2,3)
#tab_model(capageglm2, string.est	= "Estimate", show.se=TRUE,show.ci=FALSE, show.fstat = TRUE,show.re.var= TRUE, show.stat = TRUE,show.loglik = TRUE)
vif(capageglm2)
simulateResiduals(capageglm2,plot=T)

capageglmdata1 <- ggpredict(capageglm2, terms = c("capage"))
capageglmdatap1 <- plot(capageglmdata1)  + geom_line(data=mpaalphameta4, aes(x=capage,y=Observed, group=BirdID), linetype = 5, inherit.aes = FALSE, colour="gray") + geom_point(data=mpaalphameta4, aes(x=capage,y=Observed), inherit.aes = FALSE) +
  labs(title="",
       x ="Age (years)", y = "Species richness (observed)") +
  theme_tufte(base_size = 28, base_family = "Arial")+ theme(axis.line = element_line(colour = "black", linetype=1))

capageglmdata2 <- ggpredict(capageglm2, terms = c("TerminalYear"))
capageglmdatap2 <- plot(capageglmdata2, show_data = T, alpha = 0.2, jitter=T) +
  ylim(0,30) +
  labs(title="",
       x ="Terminal Year", y = "Species richness (observed)")  +
  theme_tufte(base_size = 28, base_family = "Arial")+ theme(axis.line = element_line(colour = "black", linetype=1))

capageplots <- ggarrange(capageglmdatap1,capageglmdatap2, labels = c("A","B"))

tiff("Output/capageobsage.tiff", res =400, units = "in", width = 12, height = 9, bg = "transparent")
capageplots
dev.off()

#### species diversity and capped age ----
capshannonglm2 <- lmer(Shannon ~ capage + TerminalYear + season + as.numeric(SampleYear) + SexEstimate+ Timeinfridge_scaled + CatchTime_scaled + TQcorrected_scaled  + (1|BirdID), data=mpaalphameta4 )
summary(capshannonglm2)
#tab_model(capshannonglm2)
vif(capshannonglm2)
r.squaredGLMM(capshannonglm2)

capshannonglmdata <- data.frame(ggpredict(capshannonglm2, terms=c("capage [all]","TerminalYear"))) %>% mutate(capage = x) 

capshanplot <-ggplot(capshannonglmdata, aes(x=capage,y=predicted, colour=group, fill=group)) +geom_line() + geom_ribbon(aes(ymin=conf.low,ymax=conf.high), alpha=0.2, colour=NA) + geom_point(data=mpaalphameta4, aes(x=capage,y=Shannon, colour = TerminalYear), inherit.aes = FALSE)  + geom_line(data=mpaalphameta4, aes(x=capage,y=Shannon, group=BirdID), linetype = 5, inherit.aes = FALSE, colour="gray") +
  ylab("Species diversity (shannon)") +
  xlab("Age") +
  scale_colour_colorblind(name="Terminal Year",labels=c("No","Yes")) +
  scale_fill_colorblind(name="Terminal Year",labels=c("No","Yes")) +
  theme_tufte(base_size = 15, base_family = "Arial")

### delta age - mean age ----
#### species richness delta ----
mpaalphameta4 %>% filter(!(deltaAge %in% "0")) %>% distinct(BirdID,.keep_all = T) %>% ungroup %>% reframe(mean(meanAge), sd(meanAge))

meancenteredglm5 <- glmer.nb(Observed ~ deltaAge*meanAge+ TerminalYearbird + meanAge + season + SexEstimate+ Timeinfridge_scaled + CatchTime_scaled + TQcorrected_scaled+ as.numeric(SampleYear) + (1|BirdID), data = mpaalphameta4 ,nAGQ=0)
summary(meancenteredglm5)
car::Anova(meancenteredglm5,type="III",test.statistic=c("Chisq"))
#tab_model(meancenteredglm5, string.est	= "Estimate", show.se=TRUE,show.ci=FALSE, show.fstat = TRUE,show.re.var= TRUE, show.stat = TRUE)
r.squaredGLMM(meancenteredglm5)
vif(meancenteredglm5)
simulateResiduals(meancenteredglm5,plot=T)

plot(ggpredict(meancenteredglm5, terms=c("deltaAge [all]","meanAge [3,7]"), back_transform=TRUE, type = "fixed"), show_data =T, dot_alpha = 1)

#delta age 
meancenteredglm5data <- data.frame(ggpredict(meancenteredglm5, terms=c("deltaAge [all]","meanAge [3,7]"), back_transform=TRUE, type = "fixed")) #%>% filter((group %in% "3" & x >= -2.25 & x <= 2.25) | (group %in% "5" & x >= -3.63 & x <= 2.75) | (group %in% "7" & x >= -2.67 & x <= 2.07)  )

mpaalphameta4369 <- mpaalphameta4 %>% mutate(age369 = case_when(meanAge <= 6 ~ "<6", meanAge >=6 ~ ">=6"))

mpaalphameta4 %>% ungroup %>% distinct(BirdID,.keep_all = T) %>% summarise(quantile(meanAge))

mpaalphameta4369 %>% group_by(age369) %>% summarise(min=min(deltaAge),max=max(deltaAge), n=n())

deltaobsplot <- ggplot(meancenteredglm5data, aes(x=x,y=predicted,colour=group))   + geom_line(data=mpaalphameta4369, aes(x=deltaAge,y=Observed, group=BirdID), linetype = 5, inherit.aes = FALSE, colour="gray") +  
  geom_point(data=mpaalphameta4369, aes(x=deltaAge,y=Observed,colour=age369), inherit.aes = FALSE) +
  geom_line() + 
  geom_ribbon(aes(ymin=conf.low,ymax=conf.high,fill=group),alpha=0.2,colour=NA) +
  ylab("Species richness (observed)") +
  xlab("Delta Age (years)") +
  scale_colour_manual(name="Mean Age",values = c("black","#E69F00","black","#E69F00")) +
  scale_fill_colorblind(guide = "none") +
  theme_tufte(base_size = 15, base_family = "Arial") + coord_cartesian(ylim = c(0, 120))+ theme(axis.line = element_line(colour = "black", linetype=1))  

tiff("Output/deltaobsplot.tiff", res =400, units = "in", width = 10, height = 8, bg = "white")
deltaobsplot
dev.off()

#### species diversity delta ----
shancenteredglm1 <- lmer(Shannon ~ deltaAge + meanAge + TerminalYear+ SampleYear + season + SexEstimate+ Timeinfridge_scaled + CatchTime_scaled + TQcorrected_scaled + (1|BirdID), data = mpaalphameta4)
summary(shancenteredglm1)
vif(shancenteredglm1)
simulateResiduals(shancenteredglm1,plot=T)

shancenteredglm5 <- lmer(Shannon ~ deltaAge  + TerminalYearbird + meanAge   + SampleYear  + season + SexEstimate+ Timeinfridge_scaled + CatchTime_scaled + TQcorrected_scaled + (1|BirdID), data = mpaalphameta4 )
summary(shancenteredglm5)
car::Anova(shancenteredglm5,type="III")
#tab_model(shancenteredglm5, string.est	= "Estimate", show.se=TRUE,show.ci=FALSE, show.fstat = TRUE,show.re.var= TRUE, show.stat = TRUE)
r.squaredGLMM(shancenteredglm5)



##Beta diversity ----
mpapl7bacmembers <- transform(mpapL7.bac, "compositional") %>% core_members(.,detection = 0.001, prevalence = 0.001)
mpapL7clr <- prune_taxa(mpapl7bacmembers,mpapL7.bac) %>%
  transform(., "clr") %>% 
  mutate_sample_data(AgeClass = case_when(SamplingAge <= 3 ~ "1-3",SamplingAge > 3 & SamplingAge <= 7 ~ "4-6",SamplingAge > 7 & SamplingAge <= 10 ~ "6-9", SamplingAge > 10 ~ "9+"), meancenteredage = SamplingAge - mean(SamplingAge), meanSamplingAge = mean(SamplingAge), SamplingAge2 = round(SamplingAge),sampling_id = 1:n(), cappedAge=pmin(SamplingAge,10), meancenteredagegroup = case_when(meancenteredage < -3 ~ "A", meancenteredage >= -3 & meancenteredage <= 2 ~"B", meancenteredage > 2 ~ "C" ),capage=pmin(SamplingAge,12)) %>% mutate_sample_data(agess = case_when(capage <= 1 ~ "1",capage <= 3 ~ "3",capage <= 5 ~ "5",capage <= 7 ~ "7",capage <= 9 ~ "9",capage <= 11 ~ "11", capage <= 12 ~ "12"), agess = factor(agess, levels = c("1", "3","5","7","9","11","12"))) 
mpapL7clrr <- data.frame(sample_sums(mpapL7clr))
mpa_clrmat<-vegan_otu(mpapL7clr)
mpaclrpca_st <- as(sample_data(mpapL7clr),"data.frame") %>% 
  mutate(survive = as.factor(survive),SampleYear=as.factor(SampleYear)) %>% group_by(BirdID) %>%
  mutate(deltaAge = capage - mean(capage), meanAge = mean(capage),TerminalYearbird = case_when(any(TerminalYear %in% "1") ~ "1", TRUE ~ "0")) %>% dplyr::select(capage,TerminalYear, season, SampleYear, SexEstimate, Timeinfridge_scaled,CatchTime_scaled,TQcorrected_scaled,BirdID,deltaAge,meanAge,TerminalYearbird)

#%>% merge(.,OffspringParent, by= "BirdID",all.x=T) %>% mutate(FatherAge = replace_na(FatherAge,1785), MotherAge = replace_na(MotherAge,1451))
#quantile(mpaclrpca_st$MotherAge,na.rm=T)

perm <- how(nperm = 9999, blocks = mpaclrpca_st$BirdID)
set.seed(888)
mpaclrpermsurv<- adonis2(mpa_clrmat ~ capage +  capage*TerminalYear + TerminalYear + season + SampleYear + SexEstimate+ Timeinfridge_scaled + CatchTime_scaled + TQcorrected_scaled, data=mpaclrpca_st, permutations = perm, method = "euclidean", by= "margin")

mpaclrpermsurv2<- adonis2(mpa_clrmat ~ capage + I(capage^2) +TerminalYear + season + SampleYear + SexEstimate+ Timeinfridge_scaled + CatchTime_scaled + TQcorrected_scaled, data=mpaclrpca_st, permutations = perm, method = "euclidean", by= "margin")

mpaclrpermsurv3<- adonis2(mpa_clrmat ~ capage +TerminalYear + season + SampleYear + SexEstimate+ Timeinfridge_scaled + CatchTime_scaled + TQcorrected_scaled, data=mpaclrpca_st, permutations = perm, method = "euclidean", by= "margin")

# pca plot
mpapL7ord <- phyloseq::ordinate(mpapL7clr, method = "RDA", distance = "euclidean")

orddfage <- plot_ordination(mpapL7clr, ordination= mpapL7ord, color="capage", type = "samples", shape = "TerminalYear",justDF = TRUE) 

orddfage2 <- orddfage %>% mutate(agess = case_when(capage <= 1 ~ "1",capage <= 3 ~ "3",capage <= 5 ~ "5",capage <= 7 ~ "7",capage <= 9 ~ "9",capage <= 11 ~ "11", capage <= 12 ~ "12"), agess = factor(agess, levels = c("1", "3","5","7","9","11","12"))) %>% group_by(agess) %>% mutate(agessmeanAxis1 = mean(PC1),agessmeanAxis2 = mean(PC2),agessseAxis1 = sd(PC1)/sqrt(n()),agessseAxis2 = sd(PC2)/sqrt(n()), yminimum = agessmeanAxis2 - agessseAxis2,ymaximum = agessmeanAxis2 + agessseAxis2,xminimum = agessmeanAxis1 - agessseAxis1,xmaximum = agessmeanAxis1 + agessseAxis1 ) 

orddfage3 <- orddfage2 %>% distinct(agess, .keep_all = T) %>% mutate(agess = as.numeric(levels(agess))[agess])

betagroupplot2 <- plot_ordination(mpapL7clr, ordination= mpapL7ord, color="capage", type = "samples", shape = "TerminalYear") +
  geom_point(size=3,aes(fill=capage))+
  xlab("PC1 (6.4%)") +
  ylab("PC2 (5.2%)") +
  geom_point(data=orddfage3, aes(x=agessmeanAxis1,y=agessmeanAxis2, fill=as.numeric(agess)),colour="black", inherit.aes = F,size=6, pch=23) +
  geom_errorbar(data=orddfage3,aes(x =agessmeanAxis1, ymin = yminimum, ymax= ymaximum, colour=as.numeric(agess))) +
  geom_errorbarh(data=orddfage3,aes(xmin =xminimum - agessseAxis1, xmax= xmaximum, y= agessmeanAxis2,colour=as.numeric(agess))) +
  scale_color_gradientn(colors = c("darkolivegreen3","dodgerblue","gold","red"), name = "Age (years)") +
  scale_fill_gradientn(colors = c("darkolivegreen3","dodgerblue","gold","red"), name = "Age (years)")+
  theme_tufte(base_size = 20, base_family = "Arial") + theme(axis.line = element_line(colour = "black", linetype=1))

tiff("Output/capagebeta.tiff", res =400, units = "in", width = 12, height = 8, bg = "transparent")
betagroupplot2
dev.off()

## Differential abundance analysis ----
mpacoremembers <- microbiome::transform(mpapL7.bac, "compositional") %>%
  core_members(., detection = 0.0001, prevalence = 0.2)

### Ancombc2 ----
mpapL7filteredd <- prune_taxa( mpacoremembers,mpapL7.bac) %>% mutate_sample_data(AgeClass = case_when(SamplingAge <= 5 ~ "1-5",SamplingAge > 5 & SamplingAge <= 10 ~ "6-10", SamplingAge > 10 ~ "10+"), meancenteredage = SamplingAge - mean(SamplingAge), meanSamplingAge = mean(SamplingAge), SamplingAge2 = round(SamplingAge),sampling_id = 1:n(),capage=pmin(SamplingAge,12), capage2 = I(capage^2))

mpacatancomid <- ancombc2(mpapL7filteredd,fix_formula = "capage +TerminalYear + season + SampleYear + SexEstimate+ Timeinfridge_scaled + CatchTime_scaled + TQcorrected_scaled",rand_formula = "(1|BirdID)" , p_adj_method = "BH")

mpacatancomidres <- mpacatancomid$res 

#ancom age
mpacatancomidresage <- mpacatancomidres %>% dplyr::select(taxon, ends_with("capage")) %>% mutate(direct = case_when(lfc_capage < 0 & p_capage < 0.05 ~ "p < 0.05", lfc_capage > 0 & p_capage < 0.05 ~ "p < 0.05", TRUE ~ "p > 0.05")) %>% merge(.,mpataxaspecies,by="taxon",all.x=T)

mpacatancomidresageplot <- ggplot(mpacatancomidresage, aes(x=lfc_capage,y= reorder(Species, +lfc_capage),color=direct)) +
  geom_point() +
  geom_errorbar(aes(xmin = lfc_capage - 1.96*se_capage, xmax = lfc_capage + 1.96*se_capage)) +
  geom_vline(xintercept = 0, linetype = 3) +
  ylab("") +
  xlab("Log fold change with age")+
  ggtitle("A. ANCOMBC2 Age")+
  scale_color_manual(values = c("p < 0.05"="olivedrab4","p < 0.05"="orangered3","p > 0.05"="seashell3"), name="Significance") +
  theme_tufte(base_size = 10, base_family = "Arial")+ theme(axis.line = element_line(colour = "black", linetype=1))

mpacatancomidresterm <- mpacatancomidres %>% dplyr::select(taxon, ends_with("TerminalYear1")) %>% mutate(direct = case_when(lfc_TerminalYear1 < 0 & p_TerminalYear1 < 0.05 ~ "p < 0.05", lfc_TerminalYear1 > 0 & p_TerminalYear1 < 0.05 ~ "p < 0.05", TRUE ~ "p > 0.05")) %>% merge(.,mpataxaspecies,by="taxon",all.x=T)

# ancom terminal year
mpacatancomidrestermplot <- ggplot(mpacatancomidresterm, aes(x=lfc_TerminalYear1,y= reorder(Species, +lfc_TerminalYear1), colour=direct)) +
  geom_point() +
  geom_errorbar(aes(xmin = lfc_TerminalYear1 - 1.96*se_TerminalYear1, xmax = lfc_TerminalYear1 + 1.96*se_TerminalYear1)) +
  geom_vline(xintercept = 0, linetype = 3) +
  ylab("") +
  xlab("Log fold change with Terminal Year") +
  ggtitle("C. ANCOMBC2 Terminal Year")+
  scale_color_manual(values = c("p < 0.05"="olivedrab4","p < 0.05"="orangered3","p > 0.05"="seashell3"),name="Significance") +
  theme_tufte(base_size = 10, base_family = "Arial")+ theme(axis.line = element_line(colour = "black", linetype=1))

### gllvm ----
mpaclr <- mpapL7.bac %>% microbiome::transform("clr")
mpacoreclr <- prune_taxa(mpacoremembers,mpaclr)
mpacoreclrmat <- vegan_otu(mpacoreclr)

mpacorestclr <- as(sample_data(mpacoreclr),"data.frame") %>% 
  mutate(survive = as.factor(survive)) %>% 
  mutate(capage=pmin(SamplingAge,12),SampleYear=as.numeric(SampleYear))%>%     
  dplyr::select(TerminalYear,BirdID,season,SampleYear,capage, SexEstimate,CatchTime_scaled,Timeinfridge_scaled,TQcorrected_scaled) 
mpagllvm <- gllvm(mpacoreclrmat, mpacorestclr, formula = ~ capage+TerminalYear + season + SampleYear + SexEstimate+ Timeinfridge_scaled + CatchTime_scaled + TQcorrected_scaled , family = "gaussian",num.lv = 1,sd.errors = TRUE, row.eff = ~(1|BirdID), seed = 88)

mpavgllvm1clrsig <- data.frame(summary(mpagllvm)[["Coef.tableX"]]) %>% mutate(terms= row.names(.)) %>%
  separate(terms, into = c("terms", "taxon"), sep = "\\:(?!.*:)", remove = FALSE) 

mpataxaspecies <- mpataxa %>% rownames_to_column("taxon") %>% separate_wider_delim(Species, delim = "__", names = c("s","Species")) %>% dplyr::select(taxon,Species)


# gllvm age
mpagllvmcoefage <- mpavgllvm1clrsig %>% filter(terms %in% "capage") %>% rename(pval = Pr...z.., se = Std..Error) %>% mutate(direct = case_when(Estimate < 0 & pval < 0.05 ~ "p < 0.05", Estimate > 0 & pval < 0.05 ~ "p < 0.05", TRUE ~ "p > 0.05")) %>% merge(.,mpataxaspecies,by="taxon",all.x=T)

mpagllvmresageplot <- ggplot(mpagllvmcoefage, aes(x=Estimate,y= reorder(Species, +Estimate),color=direct)) +
  geom_point() +
  geom_errorbar(aes(xmin = Estimate - 1.96*se, xmax = Estimate + 1.96*se)) +
  geom_vline(xintercept = 0, linetype = 3) +
  ggtitle("B. GLLVM Age") + 
  ylab("") +
  xlab("Log fold change with age") +
  scale_color_manual(values = c("p < 0.05"="olivedrab4","p < 0.05"="orangered3","p > 0.05"="seashell3"), name="Significance") +
  theme_tufte(base_size = 10, base_family = "Arial")+ theme(axis.line = element_line(colour = "black", linetype=1))

ggarrange(mpacatancomidresageplot,mpagllvmresageplot)


# gllvm terminal year
mpagllvmcoefterm <- mpavgllvm1clrsig %>% filter(terms %in% "TerminalYear1") %>% rename(pval = Pr...z.., se = Std..Error) %>% mutate(direct = case_when(Estimate < 0 & pval < 0.05 ~ "p < 0.05", Estimate > 0 & pval < 0.05 ~ "p < 0.05", TRUE ~ "p > 0.05")) %>% merge(.,mpataxaspecies,by="taxon",all.x=T)

mpagllvmrestermplot <- ggplot(mpagllvmcoefterm, aes(x=Estimate,y= reorder(Species, +Estimate),color=direct)) +
  geom_point() +
  geom_errorbar(aes(xmin = Estimate - 1.96*se, xmax = Estimate + 1.96*se)) +
  geom_vline(xintercept = 0, linetype = 3) +
  ggtitle("D. GLLVM Terminal Year") + 
  ylab("") +
  xlab("Log fold change with Terminal Year") +
  scale_color_manual(values = c("p < 0.05"="olivedrab4","p < 0.05"="orangered3","p > 0.05"="seashell3"), name="Significance") +
  theme_tufte(base_size = 10, base_family = "Arial")+ theme(axis.line = element_line(colour = "black", linetype=1))

ggarrange(mpacatancomidresageplot,mpagllvmresageplot,mpacatancomidrestermplot,mpagllvmrestermplot)


