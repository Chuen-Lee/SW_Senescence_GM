# eggNOG script ----
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


## phyloseq objects ----
### eggNOG phyloseq ----
#load otu table
otuNOG <- read.delim("InputTables/EM.NOGL0.txt", row.names=1, sep = "\t")
#load tax table
taxNOG <- read.csv("InputTables/NOG.Chuen.annotations.csv",row.names=1)
#convert otu and tax tables to matrix
NOGotumat <- as.matrix(otuNOG)
NOGtaxmat <- as.matrix(taxNOG)

# sample data
st <- read.csv("InputTables/st.csv")
st4 <- st %>% mutate(SexEstimate = as.factor(SexEstimate),SampleYear=as.factor(SampleYear),season=as.factor(season), TerminalYear=as.factor(TerminalYear)) %>% column_to_rownames("TubeNumber") 
stdata <- sample_data(st4)

#nog phyloseq object
NOGOTU <- otu_table(NOGotumat, taxa_are_rows = TRUE)
NOGTAX = tax_table(NOGtaxmat)
stdata <- sample_data(st4)
NOGp = phyloseq(NOGOTU, NOGTAX, stdata)

NOGpo <- NOGp %>% mutate_sample_data(Type = case_when(SampleType %in% c("F","f") ~ "F", TRUE ~ "C"))
NOGpf <- NOGpo %>% filter_sample_data(Type=="F") %>% filter_sample_data(SamplingAge > 0.5) %>% mutate_sample_data(SamplingAge2 = round(SamplingAge)) %>% mutate_sample_data(AgeClass = case_when(SamplingAge2 <= 5 ~ "A", SamplingAge2 <= 10 ~ "B", SamplingAge2 > 10 ~ "C"), capage=pmin(SamplingAge,12)) %>% mutate_sample_data(Age_scaled = arm::rescale(SamplingAge,"center"),Timeinfridge_scaled = arm::rescale(Timeinfridge, "center"),CatchTime_scaled = arm::rescale(CatchTime, "center"),TQcorrected_scaled = arm::rescale(TQcorrected, "center")) 

#number of NOGs in total
psmelt(NOGpf) %>%
  {n_distinct(.$OTU)}


## Alpha diversity ----
# ggiNEXT - takes a long time to run
# library(iNEXT)
# library(readr)
### rarefaction curve ###
#### using iNEXT ####

# #make the ASV abundance table into a matrix and check by printing first 2 rows/columns
# abund <- as(otu_table(NOGpf), "matrix")
# abund[1:2,1:2] #check
#
# #convert to a dataframe
# abund2 <- as.data.frame(abund)
# str(abund2)
#
# #iNEXT only takes numeric values, so change all values in the dataframe to numeric values instead of integers.
# df2 <- mutate_all(abund2, function(x) as.numeric(x))  %>% mutate(OTU = rownames(.)) %>% dplyr::select(OTU, everything())
# str(df2)
# write_csv(df2,"Output/MFFiNEXTNOG.csv")
# 
# df2 <- read.csv("Output/MFFiNEXTNOG.csv", row.names = 1, sep = ",")
# df3 <- df2 %>% dplyr::select(-SW1168,-SW1183) %>% .[,1:20]
# inext_test<-iNEXT(df3, q=0, datatype="abundance", endpoint=500000)
# dput(inext, "inextNOG.txt")
# 
# inext7 <- dget("InputTables/inextNOG.txt")
# 
# inext7check <- inext7[["iNextEst"]][["size_based"]]
# 
# #plot rarefaction curve
# rarefaction<- ggiNEXT(inext7, type=1, se=FALSE, grey= TRUE) + theme(legend.position = "none")+ xlab("Sequencing depth") + ylab("Observed Taxa") 
# rf7 <- rarefaction + geom_vline(xintercept=10000, alpha=0.5, linetype=2) +scale_shape_manual(values=rep(20,164))
# 
# #Plot sample completeness curve
# completeness<-ggiNEXT(inext7, type=2, se=TRUE, grey=TRUE)+scale_shape_manual(values=rep(20,164))+ theme(legend.position = "none") +xlab("Read count") +ylab("Sample completeness")
# c7 <-  completeness + geom_vline(xintercept=10000, alpha=0.5, linetype=2) + geom_hline(yintercept = 0.95, linetype=2)
# 
# ggiNEXT(inext7, type=3) + theme(legend.position = "none")
# 
# inextarrange7 <- ggarrange(rf7,c7)

#rarefy
NOGrare <- rarefy_even_depth(NOGpf, rngseed=88,sample.size = 100000, replace=F)
NOGalpha <- alpha_estimate(NOGrare) %>% mutate(transobs = exp(arm::rescale(Observed)),transshan = exp(arm::rescale(Shannon))) %>% group_by(BirdID) %>% mutate(deltaAge = capage - mean(capage), meanAge = mean(capage),TerminalYearbird = case_when(any(TerminalYear %in% "1") ~ "1", TRUE ~ "0"))


NOGglmobs <- lmer(transobs ~ capage+TerminalYear + season + SexEstimate+ Timeinfridge_scaled + CatchTime_scaled + TQcorrected_scaled + SampleYear+ (1|BirdID), data=NOGalpha)
summary(NOGglmobs)
simulateResiduals(NOGglmobs,plot=T)
NOGglmobsdata <- data.frame(ggpredict(NOGglmobs, terms="capage")) %>% mutate(Observed = (log(predicted)*2*576.4921)+3053.719, logconf.low = (log(conf.low)*2*576.4921)+3053.719,logconf.high = (log(conf.high)*2*576.4921)+3053.719)
nogobsplot <- ggplot(NOGglmobsdata,aes(x=x,y=Observed)) + geom_line() + geom_ribbon(aes(ymin=logconf.low,ymax=logconf.high), alpha = 0.3) +geom_point(data=NOGalpha,aes(x=capage,y=Observed),inherit.aes = F)+ geom_line(data=NOGalpha, aes(x=capage,y=Observed, group=BirdID), linetype = 5, inherit.aes = FALSE, colour="gray")+ xlab("Age (years)") +
  theme_tufte(base_size = 15, base_family = "Arial")+ theme(axis.line = element_line(colour = "black", linetype=1))
tab_model(NOGglmobs)
r.squaredGLMM(NOGglmobs)


NOGlmshan <- lmer(exp(Shannon) ~ capage+TerminalYear + season + SexEstimate+ Timeinfridge_scaled +CatchTime_scaled  + TQcorrected_scaled + (SampleYear) + (1|BirdID), data=NOGalpha)
summary(NOGlmshan)
r.squaredGLMM(NOGlmshan)
vif(NOGlmshan)
simulateResiduals(NOGlmshan,plot=T)
NOGlmshandata <- data.frame(ggpredict(NOGlmshan, terms="capage")) %>% mutate(Shannon = log(predicted), logconf.low = log(conf.low),logconf.high = log(conf.high))
nogshanplot <- ggplot(NOGlmshandata,aes(x=x,y=Shannon)) + geom_line() + geom_ribbon(aes(ymin=logconf.low,ymax=logconf.high), alpha = 0.3) +geom_point(data=NOGalpha,aes(x=capage,y=Shannon),inherit.aes = F)+ geom_line(data=NOGalpha, aes(x=capage,y=Shannon, group=BirdID), linetype = 5, inherit.aes = FALSE, colour="gray")+ xlab("Age (years)") +
  theme_tufte(base_size = 15, base_family = "Arial")+ theme(axis.line = element_line(colour = "black", linetype=1))

ggarrange(nogobsplot,nogshanplot,labels = c("A","B"))

tiff("Output/capagenog.tiff", res =400, units = "in", width = 12, height = 9, bg = "transparent")
ggarrange(nogobsplot,nogshanplot,labels = c("A","B"))
dev.off()

NOGdelobs <- lmer(transobs ~ deltaAge + meanAge +TerminalYearbird + season + SexEstimate+ Timeinfridge_scaled + CatchTime_scaled + TQcorrected_scaled + as.numeric(SampleYear) + (1|BirdID), data=NOGalpha)
summary(NOGdelobs)
car::Anova(NOGdelobs,type="III")
simulateResiduals(NOGdelobs,plot=T)
NOGdelobsdata <- data.frame(ggpredict(NOGdelobs, terms="deltaAge [all]", back_transform = F)) 
nogobsdeltaplot <- ggplot(NOGdelobsdata,aes(x=x,y=predicted)) + geom_line() + geom_ribbon(aes(ymin=conf.low,ymax=conf.high), alpha = 0.3) + geom_line(data=NOGalpha, aes(x=deltaAge,y=transobs, group=BirdID), linetype = 5, inherit.aes = FALSE, colour="gray")+ 
  geom_point(data=NOGalpha,aes(x=deltaAge,y=transobs),inherit.aes = F) + xlab("Delta Age (years)")+
  ylab("Transformed functional richness (observed)")  +
  theme_tufte(base_size = 15, base_family = "Arial")+ theme(axis.line = element_line(colour = "black", linetype=1))
r.squaredGLMM(NOGdelobs)


NOGdelshan <- lmer(exp(Shannon) ~ deltaAge + meanAge +TerminalYearbird + season + SexEstimate+ Timeinfridge_scaled + CatchTime_scaled + TQcorrected_scaled + as.numeric(SampleYear) + (1|BirdID), data=NOGalpha)
summary(NOGdelshan)
simulateResiduals(NOGdelshan,plot=T)
NOGdelshandata <- data.frame(ggpredict(NOGdelshan, terms="deltaAge [all]",back_transform = F)) #%>% mutate(Shannon = log(predicted), logconf.low = log(abs(conf.low)),logconf.high = log(abs(conf.high)))
plot(NOGdelshandata)
nogshandeltaplot <- ggplot(NOGdelshandata,aes(x=x,y=predicted)) + geom_line() + geom_ribbon(aes(ymin=conf.low,ymax=conf.high), alpha = 0.3)  + geom_line(data=NOGalpha, aes(x=deltaAge,y=exp(Shannon), group=BirdID), linetype = 5, inherit.aes = FALSE, colour="gray")+ 
  geom_point(data=NOGalpha,aes(x=deltaAge,y=exp(Shannon)),inherit.aes = F)+ xlab("Delta Age (years)") + ylab("Transformed functional Shannon diversity") +
  theme_tufte(base_size = 15, base_family = "Arial")+ theme(axis.line = element_line(colour = "black", linetype=1))
r.squaredGLMM(NOGdelshan)

tiff("Output/capagenogdelta.tiff", res =400, units = "in", width = 12, height = 9, bg = "transparent")
ggarrange(nogobsdeltaplot,nogshandeltaplot,labels = c("A","B"))
dev.off()

# Beta Diversity Analysis eggNOG ----

NOGpf5p <- transform(NOGpf, "compositional") %>%core_members(., detection = 0.0000001, prevalence = 0.05)

NOGpfclr <- prune_taxa(NOGpf5p,NOGpf) %>% transform(., "clr") 
NOG_clrmat<-vegan_otu(NOGpfclr)
NOGpfclrpca_st <- as(sample_data(NOGpfclr),"data.frame") %>% 
  mutate(survive = as.factor(survive)) %>% group_by(BirdID) %>%
  mutate(deltaAge = capage - mean(capage), meanAge = mean(capage),TerminalYearbird = case_when(any(TerminalYear %in% "1") ~ "1", TRUE ~ "0"))

n_distinct(NOGpfclrpca_st$BirdID)

perms <- how(nperm = 9999, blocks = NOGpfclrpca_st$BirdID)
set.seed(888)
NOGclrperm<- adonis2(NOG_clrmat ~ capage + I(capage^2) +  capage*TerminalYear + I(capage^2)*TerminalYear + season + SampleYear + SexEstimate+ Timeinfridge_scaled + CatchTime_scaled + TQcorrected_scaled, data=NOGpfclrpca_st, permutations = perms, method = "euclidean", by= "margin")
NOGclrperm2<- adonis2(NOG_clrmat ~ capage + I(capage^2) + TerminalYear + season + SampleYear + SexEstimate+ Timeinfridge_scaled + CatchTime_scaled + TQcorrected_scaled, data=NOGpfclrpca_st, permutations = perms, method = "euclidean", by= "margin")

# final adonis2
NOGclrperm3 <- adonis2(NOG_clrmat ~ capage + TerminalYear + season + SampleYear + SexEstimate+ Timeinfridge_scaled + CatchTime_scaled + TQcorrected_scaled, data=NOGpfclrpca_st, permutations = perms, method = "euclidean", by= "margin")


NOGord <- phyloseq::ordinate(NOGpfclr, method = "PCoA", distance = "euclidean")

betaNOGplot <- plot_ordination(NOGpfclr, ordination= NOGord, color="capage", type = "samples") +
  geom_point(size=1,aes(fill=capage))+
  xlab("PC1 (14%)") +
  ylab("PC2 (8.2%)") +
  scale_color_gradientn(colors = c("darkolivegreen3","dodgerblue","gold","red"), name = "Age (years)") +
  scale_fill_gradientn(colors = c("darkolivegreen3","dodgerblue","gold","red"), name = "Age (years)")+
  theme_tufte(base_size = 15, base_family = "Arial") + theme(axis.line = element_line(colour = "black", linetype=1))

NOGorddfage <- plot_ordination(NOGpfclr, ordination= NOGord, color="capage", type = "samples",justDF = TRUE) 
NOGorddfage2 <- NOGorddfage %>% mutate(agess = case_when(capage <= 1 ~ "1",capage <= 3 ~ "3",capage <= 5 ~ "5",capage <= 7 ~ "7",capage <= 9 ~ "9",capage <= 11 ~ "11", capage <= 12 ~ "12"), agess = factor(agess, levels = c("1", "3","5","7","9","11","12"))) %>% group_by(agess) %>% mutate(agessmeanAxis1 = mean(Axis.1),agessmeanAxis2 = mean(Axis.2),agessseAxis1 = sd(Axis.1)/sqrt(n()),agessseAxis2 = sd(Axis.2)/sqrt(n()), yminimum = agessmeanAxis2 - agessseAxis2,ymaximum = agessmeanAxis2 + agessseAxis2,xminimum = agessmeanAxis1 - agessseAxis1,xmaximum = agessmeanAxis1 + agessseAxis1 ) 

NOGorddfage3 <- NOGorddfage2 %>% distinct(agess, .keep_all = T) %>% mutate(agess = as.numeric(levels(agess))[agess])

betaNOGplot2 <- plot_ordination(NOGpfclr, ordination= NOGord, color="capage", type = "samples") +
  geom_point(size=1,aes(fill=capage))+
  xlab("PC1 (14%)") +
  ylab("PC2 (8.2%)") +
  geom_point(data=NOGorddfage3, aes(x=agessmeanAxis1,y=agessmeanAxis2, fill=as.numeric(agess)),colour="black", inherit.aes = F,size=3, pch=23) +
  geom_errorbar(data=NOGorddfage3,aes(x =agessmeanAxis1, ymin = yminimum, ymax= ymaximum, colour=as.numeric(agess))) +
  geom_errorbarh(data=NOGorddfage3,aes(xmin =xminimum - agessseAxis1, xmax= xmaximum, y= agessmeanAxis2,colour=as.numeric(agess))) +
  scale_color_gradientn(colors = c("darkolivegreen3","dodgerblue","gold","red"), name = "Age (years)") +
  scale_fill_gradientn(colors = c("darkolivegreen3","dodgerblue","gold","red"), name = "Age (years)")+
  theme_tufte(base_size = 15, base_family = "Arial") + theme(axis.line = element_line(colour = "black", linetype=1))


# Differential Abundance Analysis eggNOG ----
# NEW CATEGORIES
Chuen_new_COG_tax <- read.csv("~/Downloads/cog_all_results.tsv",sep="\t")

NOGpfmelt <- psmelt(NOGpf) %>% dplyr::select(OTU,Abundance,Sample) %>% merge(.,Chuen_new_COG_tax,by.x="OTU",by.y="COG",all.x=T ) 

NOGpfmelted <- NOGpfmelt %>% mutate(Cat = case_when(is.na(Cat)~"`",Cat == "" ~ "`" ,TRUE ~ Cat)) %>% separate_longer_position(Cat,1) %>% group_by(Sample,OTU) %>% arrange(desc(Cat)) %>% distinct(Sample,.keep_all = T)
NOGpfcogcat <- NOGpfmelted %>% group_by(Sample,Cat) %>% summarise(COGabundance = sum(Abundance)) %>% pivot_wider(names_from = Sample,values_from = COGabundance) %>%  mutate(OTU = paste0("Cat_", Cat)) %>% column_to_rownames(var="OTU") 


#matrix and phyloseq
NOGpfcogcatotu <- NOGpfcogcat %>% dplyr::select(-Cat)
NOGpfcogcatotumat <- otu_table(as.matrix(NOGpfcogcatotu), taxa_are_rows = TRUE)

eggNOG_Cat_description <- read_csv("InputTables/eggNOG_Cat_description.csv")
NOGpfcogcattax <- NOGpfcogcat %>% dplyr::select(Cat) %>% rownames_to_column("rownames") %>% merge(.,eggNOG_Cat_description,by="Cat",all.x=T ) %>% column_to_rownames("rownames")
NOGpfcogcattaxmat = tax_table(as.matrix(NOGpfcogcattax))

NOGcat <- phyloseq(NOGpfcogcatotumat, NOGpfcogcattaxmat, stdata) %>% filter_sample_data(SamplingAge > 0.5) %>% filter_sample_data(Type %in% "F") %>% mutate_sample_data(Age_scaled = arm::rescale(SamplingAge,"center"),Timeinfridge_scaled = arm::rescale(Timeinfridge, "center"),CatchTime_scaled = arm::rescale(CatchTime, "center"),TQcorrected_scaled = arm::rescale(TQcorrected, "center"), capage = pmin(SamplingAge, 12)) 
NOGcat.sam <- NOGcat %>% filter_sample_data(Type%in%"F") 



# new cat ancombc2 ----
NOGcatancomid <- ancombc2(NOGcat.sam,fix_formula = "capage+ TerminalYear + season + SampleYear + SexEstimate+ Timeinfridge_scaled + CatchTime_scaled + TQcorrected_scaled",rand_formula = "(1|BirdID)" , p_adj_method = "BH")
NOGcatancomidres <- NOGcatancomid$res

NOGpfcogcattaxrows <- NOGpfcogcattax %>% rownames_to_column("taxon")

NOGcatancomidresage <- NOGcatancomidres %>% dplyr::select(taxon, ends_with("capage")) %>% mutate(direct = case_when(lfc_capage< 0 & p_capage < 0.05 ~ "Negative_LFC", lfc_capage > 0 & p_capage < 0.05 ~ "p < 0.05", TRUE ~ "p > 0.05")) %>% arrange(desc(lfc_capage)) %>% mutate(order = 1:n()) %>% merge(.,NOGpfcogcattaxrows,by="taxon",all=T)

NOGorder <- NOGcatancomidresage %>% dplyr::select(taxon,order) 

NOGcatancomidresageplot <- ggplot(NOGcatancomidresage, aes(x=lfc_capage,y= reorder(COG_Description, +lfc_capage),color=direct)) +
  geom_point(size=3) +
  geom_errorbar(aes(xmin = lfc_capage - 1.96*se_capage, xmax = lfc_capage + 1.96*se_capage),linewidth =1.25) +
  geom_vline(xintercept = 0, linetype = 3) +
  ggtitle("A. ANCOMBC2")+
  ylab("") +
  xlab("Log fold change with age") +
  scale_color_manual(values = c("p < 0.05"="olivedrab4","p < 0.05"="orangered3","p > 0.05"="seashell3"), name="Significance") +
  theme_tufte(base_size = 15, base_family = "Arial")+ theme(axis.line = element_line(colour = "black", linetype=1))

tiff("Output/NOG_cat_ANCOM.tiff", res =400, units = "in", width = 15, height = 9, bg = "transparent")
NOGcatancomidresageplot
dev.off()

## new cat gllvm ----
NOGcatcaomp <- transform(NOGcat, "compositional") %>%
  core_members(.,detection = 0.001, prevalence = 0.2)

NOGcatclr <- transform(NOGcat,"clr")

NOGcatmat <- vegan_otu(NOGcatclr)
NOGcatstclr <- as(sample_data(NOGcatclr),"data.frame") %>% 
  mutate(TerminalYear = as.factor(TerminalYear)) %>% 
  mutate(capage=pmin(SamplingAge,12))%>%     
  dplyr::select(TerminalYear,BirdID,season,SampleYear,capage, SexEstimate,CatchTime_scaled,Timeinfridge_scaled,TQcorrected_scaled) 
NOGcatgllvm <- gllvm(NOGcatmat, NOGcatstclr, formula = ~ capage+TerminalYear + season + SampleYear + SexEstimate+ Timeinfridge_scaled + CatchTime_scaled + TQcorrected_scaled , family = "gaussian",num.lv = 2,sd.errors = TRUE, row.eff = ~(1|BirdID), seed = 8,starting.val="random")


summary(NOGcatgllvm)
coefplot.gllvm(NOGcatgllvm, which.Xcoef = 1)
coefplot.gllvm(NOGcatgllvm, which.Xcoef = 10)
coefplot.gllvm(NOGcatgllvm, which.Xcoef = 3:6)
coefplot.gllvm(NOGcatgllvm, which.Xcoef = 7:9)

coefplot.gllvm(NOGcatgllvm, which.Xcoef = 3) 
coefplot.gllvm(NOGcatgllvm, which.Xcoef = 4)
coefplot.gllvm(NOGcatgllvm, which.Xcoef = 5)
coefplot.gllvm(NOGcatgllvm, which.Xcoef = 6)
coefplot.gllvm(NOGcatgllvm, which.Xcoef = 7)
coefplot.gllvm(NOGcatgllvm, which.Xcoef = 8)
coefplot.gllvm(NOGcatgllvm, which.Xcoef = 9)

NOGsummarygllvm <- summary(NOGcatgllvm)
NOGcatgllvmcoef <- as.data.frame(NOGsummarygllvm[["Coef.tableX"]]) %>% mutate(terms= row.names(.)) %>% separate_wider_regex(cols = terms, c(terms = ".*", ":", OTU = ".*"))

NOGcatgllvmcoefage <- NOGcatgllvmcoef %>% filter(terms %in% "capage") %>% rename(pval = 'Pr(>|z|)', se = 'Std. Error') %>% mutate(direct = case_when(Estimate < 0 & pval < 0.05 ~ "p < 0.05", Estimate > 0 & pval < 0.05 ~ "p < 0.05", TRUE ~ "p > 0.05")) %>% merge(.,NOGorder,by.x="OTU", by.y="taxon",all=T) %>% merge(.,NOGpfcogcattaxrows,by.x="OTU",by.y="taxon",all=T)

NOGcatgllvmresageplot <- ggplot(NOGcatgllvmcoefage, aes(x=Estimate,y= reorder(COG_Description, -order),color=direct)) +
  geom_point(size = 3) +
  geom_errorbar(aes(xmin = Estimate - 1.96*se, xmax = Estimate + 1.96*se),linewidth=1.25) +
  geom_vline(xintercept = 0, linetype = 3) +
  ggtitle("B. GLLVM") + 
  ylab("") +
  xlab("Log fold change with age") +
  scale_color_manual(values = c("p < 0.05"="olivedrab4","p < 0.05"="orangered3","p > 0.05"="seashell3"), name="Significance") +
  theme_tufte(base_size = 15, base_family = "Arial")+ theme(axis.line = element_line(colour = "black", linetype=1)) + 
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank())

ggarrange(NOGcatancomidresageplot,NOGcatgllvmresageplot, labels = c("A. ANCOMBC2","B. GLLVM"),common.legend = T, vjust = 0.5, legend = "right")

tiff("Output/NOG_cat_DAA.tiff", res =400, units = "in", width = 15, height = 9, bg = "transparent")
ggarrange(NOGcatancomidresageplot,NOGcatgllvmresageplot,common.legend = T, vjust = 0.5, legend = "right", widths = c(3,1))
dev.off()

## cat X only ----
Chuen_new_COG_taxX <- Chuen_new_COG_tax %>% filter(Cat %in% "X")

NOGpfX <- prune_taxa(Chuen_new_COG_taxX$COG, NOGpf) %>% mutate_sample_data(Age_scaled = arm::rescale(SamplingAge,"center"),Timeinfridge_scaled = arm::rescale(Timeinfridge, "center"),CatchTime_scaled = arm::rescale(CatchTime, "center"),TQcorrected_scaled = arm::rescale(TQcorrected, "center"), capage = pmin(SamplingAge, 12)) 
NOGpfXdf <- psmelt(NOGpfX) 
n_distinct(NOGpfXdf$OTU)

NOGpfXcoremembers <- transform(NOGpfX,"compositional") %>%
  core_members(., detection = 0.001, prevalence =0.2)
NOGpfXcore <- NOGpfX %>% mutate_sample_data(capage = SamplingAge)

NOGpfXclr <- transform(NOGpfXcore,"clr")
NOGXtmat <- vegan_otu(NOGpfXclr)
NOGXstclr <- as(sample_data(NOGpfXclr),"data.frame") %>% 
  mutate(TerminalYear = as.factor(TerminalYear)) %>% 
  mutate(capage=pmin(SamplingAge,12))%>%     
  dplyr::select(TerminalYear,BirdID,season,SampleYear,capage, SexEstimate,CatchTime_scaled,Timeinfridge_scaled,TQcorrected_scaled) 
NOGXgllvm <- gllvm(NOGXtmat, NOGXstclr, formula = ~ capage+TerminalYear + season + SampleYear + SexEstimate+ Timeinfridge_scaled + CatchTime_scaled + TQcorrected_scaled , family = "gaussian",num.lv = 1,sd.errors = TRUE, row.eff = ~(1|BirdID), seed = 8,starting.val="random")
coefplot.gllvm(NOGXgllvm,which.Xcoef = 1:2)
summary(NOGXgllvm)

NOGXsummarygllvm <- summary(NOGXgllvm)
NOGXgllvmcoef <- as.data.frame(NOGXsummarygllvm[["Coef.tableX"]]) %>% mutate(terms= row.names(.)) %>% separate_wider_regex(cols = terms, c(terms = ".*", ":", OTU = ".*"))

NOGXgllvmcoefage <- NOGXgllvmcoef %>% filter(terms %in% "capage") %>% rename(pval = 'Pr(>|z|)', se = 'Std. Error') %>% mutate(direct = case_when(Estimate < 0 & pval < 0.05 ~ "p < 0.05", Estimate > 0 & pval < 0.05 ~ "p < 0.05", TRUE ~ "p > 0.05"))

NOGXgllvmresageplot <- ggplot(NOGXgllvmcoefage, aes(x=Estimate,y= reorder(OTU, +Estimate),color=direct)) +
  geom_point() +
  geom_errorbar(aes(xmin = Estimate - 1.96*se, xmax = Estimate + 1.96*se)) +
  geom_vline(xintercept = 0, linetype = 3) +
  ggtitle("B. GLLVM") + 
  ylab("") +
  xlab("Log fold change with age") +
  scale_color_manual(values = c("p < 0.05"="olivedrab4","p < 0.05"="orangered3","p > 0.05"="seashell3"), name="Significance") +
  theme_tufte(base_size = 10, base_family = "Arial")+ theme(axis.line = element_line(colour = "black", linetype=1))


NOGpfXcoremembers <- transform(NOGpfX,"compositional") %>%
  core_members(., detection = 0.001, prevalence =0.2)
NOGpfXcore <- NOGpfX %>% mutate_sample_data(capage = SamplingAge)

NOGXancomid <- ancombc2(NOGpfXcore,fix_formula = "capage+ TerminalYear + season + SampleYear + SexEstimate+ Timeinfridge_scaled + CatchTime_scaled + TQcorrected_scaled",rand_formula = "(1|BirdID)" , p_adj_method = "BH")
NOGXancomidres <- NOGXancomid$res

NOGXancomidresage <- NOGXancomidres %>% dplyr::select(taxon, ends_with("capage")) %>% mutate(direct = case_when(lfc_capage< 0 & p_capage < 0.05 ~ "Negative_LFC", lfc_capage > 0 & p_capage < 0.05 ~ "Positive_LFC", TRUE ~ "Neutral"))

NOGXancomidresageplot <- ggplot(NOGXancomidresage, aes(x=lfc_capage,y= reorder(taxon, +lfc_capage),color=direct)) +
  geom_point() +
  geom_errorbar(aes(xmin = lfc_capage - 1.96*se_capage, xmax = lfc_capage + 1.96*se_capage)) +
  geom_vline(xintercept = 0, linetype = 3) +
  ylab("") +
  xlab("Log fold change with age") +
  ggtitle("A. ANCOMBC2")+
  scale_color_manual("Significance",values = c("Positive_LFC"="olivedrab4","Negative_LFC"="orangered3","Neutral"="seashell3")) +
  theme_tufte(base_size = 10, base_family = "Arial")

tiff("Output/COG_Ancomc2GLLVM.tiff", res =400, units = "in", width = 10, height = 15, bg = "transparent")
ggarrange(NOGXancomidresageplot,NOGXgllvmresageplot, common.legend = T,legend.grob = get_legend(NOGXgllvmresageplot),legend="right" )
dev.off()

## COG2801 DAA ----
NOGpf2801 <- prune_taxa("COG2801",NOGpfXclr)
NOGpf2801df <- psmelt(NOGpf2801) %>% group_by(BirdID) %>% mutate(deltaAge = capage - mean(capage), meanAge = mean(capage),TerminalYearbird = case_when(any(TerminalYear %in% "1") ~ "1", TRUE ~ "0"))

NOG2801gam <- gam(Abundance ~ s(capage,bs="cs")+ TerminalYear + season + SampleYear + SexEstimate+ Timeinfridge_scaled + CatchTime_scaled + TQcorrected_scaled + s(BirdID, bs = "re"), data = NOGpf2801df, family = "gaussian")
summary(NOG2801gam)
NOG2801gamdata <- ggpredict(NOG2801gam, terms="capage", back_transform = TRUE)
plot(NOG2801gamdata) + geom_point(data=NOGpf2801df, aes(x=capage,y=Abundance), inherit.aes = F) 

simulateResiduals(NOG2801gam,plot=T)

ggplot(NOGpf2801df, aes(x=Abundance)) + geom_histogram()

NOG2801glm <- lmer( Abundance ~ capage+ TerminalYear + season  + SexEstimate+ Timeinfridge_scaled + CatchTime_scaled + TQcorrected_scaled+ SampleYear + (1|BirdID), data = NOGpf2801df)
summary(NOG2801glm)
NOG2801glmdata <- ggpredict(NOG2801glm, terms="capage")
COG2801lmplot <- plot(NOG2801glmdata) + geom_point(data=NOGpf2801df, aes(x=capage,y=Abundance), inherit.aes = F) + xlab("Age")
vif(NOG2801glm)

NOG2801gamdm <- gam(Abundance ~ deltaAge+ meanAge +TerminalYearbird + season + SampleYear + SexEstimate+ Timeinfridge_scaled + CatchTime_scaled + TQcorrected_scaled + s(BirdID, bs = "re"), data = NOGpf2801df, family = "gaussian")
summary(NOG2801gamdm)

NOG2801glmdm <- lmer(Abundance ~ deltaAge*meanAge +TerminalYearbird + season  + SexEstimate+ Timeinfridge_scaled + CatchTime_scaled + TQcorrected_scaled + as.numeric(SampleYear) + (1|BirdID), data = NOGpf2801df)
summary(NOG2801glmdm)
car::Anova(NOG2801glmdm,type=3)
r.squaredGLMM(NOG2801glmdm)
NOG2801glmdmdata <- ggpredict(NOG2801glmdm,terms = c("deltaAge [all]"))
COG2801lmwithinplot <- plot(NOG2801glmdmdata) + geom_line(data=NOGpf2801df, aes(x=deltaAge,y=Abundance, group=BirdID), linetype = 5, inherit.aes = FALSE, colour="gray") + geom_point(data=NOGpf2801df, aes(x=deltaAge,y=Abundance), inherit.aes = F)+ ylab("CLR-transformed COG2801 Abundance")+ xlab("Delta Age (years)")  +
  ggtitle("")+
  theme_tufte(base_size = 15, base_family = "Arial")+ theme(axis.line = element_line(colour = "black", linetype=1))
r.squaredGLMM(NOG2801glmdm)

simulateResiduals(NOG2801glmdm, plot=T)

NOG2801glmmmdata <- ggpredict(NOG2801glmdm,terms = c("meanAge [all]"))
COG2801lmwithinmeanplot <- plot(NOG2801glmmmdata) + geom_line(data=NOGpf2801df, aes(x=meanAge,y=Abundance, group=BirdID), linetype = 5, inherit.aes = FALSE, colour="gray") + geom_point(data=NOGpf2801df, aes(x=meanAge,y=Abundance), inherit.aes = F)+ ylab("CLR-transformed COG2801 Abundance")+ xlab("Mean Age (years)") +ggtitle("")+
  theme_tufte(base_size = 15, base_family = "Arial")+ theme(axis.line = element_line(colour = "black", linetype=1))



tiff("Output/cog2801within.tiff", res =400, units = "in", width = 12, height = 9, bg = "transparent")
ggarrange(COG2801lmwithinplot,COG2801lmwithinmeanplot,labels=c("A","B"))
dev.off()


## genemat ----

genemat <- read.csv("InputTables/Matrix.mat.smaller", sep = "\t") #%>% pivot_wider(names_from = TubeNo, values_from = Reads) %>% mutate(OTU = paste0("OTU", 1:nrow(.)))
genemat2 <- genemat %>% mutate(Genes2 = paste0("Gene", Genes))
# taxa table
genecog2801 <- read.csv("InputTables/COG2801genecat_smaller.tsv", sep = "\t")
genetaxa <- genemat2 %>% dplyr::select(Genes,Genes2) %>% merge(.,genecog2801,by="Genes",all.x=T) %>% mutate(Taxa = replace_na(Taxa, "nogcog2801")) %>% column_to_rownames(., var="Genes2")

# otu table
genesotu <- genemat2  %>% dplyr::select(-Genes) %>% column_to_rownames(., var="Genes2") %>% replace(is.na(.), 0)


genematphyseq <- phyloseq(tax_table(as.matrix(genetaxa)), otu_table(as.matrix(genesotu), taxa_are_rows = T), sample_data(NOGpf))

genematphyseqcore_members <- transform(genematphyseq,"compositional") %>% core_members(.,detection = 0.00001, prevalence = 0.5)

genematphyseqcore <- prune_taxa(genematphyseqcore_members,genematphyseq) 
genematphyseqclr <- genematphyseqcore %>% transform(.,"clr")

#cog2801 only 
genematphyseqcog2801 <- genematphyseq %>% filter_tax_table(Taxa %in% "COG2801")
genematphyseqcog2801core_members <- transform(genematphyseqcog2801,"compositional") %>% core_members(.,detection = 0.001, prevalence = 0.2)
genematphyseqcog2801core <- prune_taxa(genematphyseqcog2801core_members,genematphyseqcog2801)
#genematphyseqcog2801core2 <- aggregate_taxa(genematphyseqcog2801core,Taxa)

genematphyseqancomid <- ancombc2(genematphyseqcog2801core,fix_formula = "capage",rand_formula = "(1|BirdID)" , p_adj_method = "BH")
genematphyseqancomidres <- genematphyseqancomid$res

genematphyseqancomid2 <- ancombc2(genematphyseqcog2801core,fix_formula = "capage+ TerminalYear + season + SampleYear + SexEstimate+ Timeinfridge_scaled + CatchTime_scaled + TQcorrected_scaled",rand_formula = "(1|BirdID)" , p_adj_method = "BH")
genematphyseqancomidres2 <- genematphyseqancomid2$res

genematphyseqancomid2age <- genematphyseqancomidres2 %>% dplyr::select(taxon, ends_with("capage")) %>% mutate(direct = case_when(lfc_capage< 0 & p_capage < 0.05 ~ "Negative_LFC", lfc_capage > 0 & p_capage < 0.05 ~ "p < 0.05", TRUE ~ "p > 0.05"))

genematphyseqancomid2ageplot <- ggplot(genematphyseqancomid2age, aes(x=lfc_capage,y= reorder(taxon, +lfc_capage),color=direct)) +
  geom_point() +
  geom_errorbar(aes(xmin = lfc_capage - 1.96*se_capage, xmax = lfc_capage + 1.96*se_capage)) +
  geom_vline(xintercept = 0, linetype = 3) +
  ggtitle("A. ANCOMBC2")+
  ylab("") +
  xlab("Log fold change with age") +
  scale_color_manual(values = c("p < 0.05"="olivedrab4","p < 0.05"="orangered3","p > 0.05"="seashell3"), name="Significance") +
  theme_tufte(base_size = 10, base_family = "Arial")+ theme(axis.line = element_line(colour = "black", linetype=1))

# genemat gllvm
genematphyseqcog2801coreclr <- transform(genematphyseq,"clr") %>% prune_taxa(genematphyseqcog2801core_members,.)
genematmat <- vegan_otu(genematphyseqcog2801coreclr)
genematphyseqclrst <- as(sample_data(genematphyseqcog2801coreclr),"data.frame") %>% 
  mutate(TerminalYear = as.factor(TerminalYear)) %>% 
  mutate(capage=pmin(SamplingAge,12))%>%     
  dplyr::select(TerminalYear,BirdID,season,SampleYear,capage, SexEstimate,CatchTime_scaled,Timeinfridge_scaled,TQcorrected_scaled) 
genematgllvm <- gllvm(genematmat, genematphyseqclrst, formula = ~ capage+TerminalYear + season + SampleYear + SexEstimate+ Timeinfridge_scaled + CatchTime_scaled + TQcorrected_scaled , family = "gaussian",num.lv = 2,sd.errors = TRUE, row.eff = ~(1|BirdID), seed = 8,starting.val="random")
coefplot.gllvm(genematgllvm,which.Xcoef = 1)
summary(genematgllvm)

genematsummarygllvm <- summary(genematgllvm)
genematgllvmcoef <- as.data.frame(genematsummarygllvm[["Coef.tableX"]]) %>% mutate(terms= row.names(.)) %>% separate_wider_regex(cols = terms, c(terms = ".*", ":", OTU = ".*"))

genematgllvmcoefage <- genematgllvmcoef %>% filter(terms %in% "capage") %>% rename(pval = 'Pr(>|z|)', se = 'Std. Error') %>% mutate(direct = case_when(Estimate < 0 & pval < 0.05 ~ "p < 0.05", Estimate > 0 & pval < 0.05 ~ "p < 0.05", TRUE ~ "p > 0.05"))

genematgllvmcoefageplot <- ggplot(genematgllvmcoefage, aes(x=Estimate,y= reorder(OTU, +Estimate),color=direct)) +
  geom_point() +
  geom_errorbar(aes(xmin = Estimate - 1.96*se, xmax = Estimate + 1.96*se)) +
  geom_vline(xintercept = 0, linetype = 3) +
  ggtitle("B. GLLVM") + 
  ylab("") +
  xlab("Log fold change with age") +
  scale_color_manual(values = c("p < 0.05"="olivedrab4","p < 0.05"="orangered3","p > 0.05"="seashell3"), name="Significance") +
  theme_tufte(base_size = 10, base_family = "Arial")+ theme(axis.line = element_line(colour = "black", linetype=1))

tiff("Output/genecat_Ancomc2GLLVM.tiff", res =400, units = "in", width = 10, height = 15, bg = "transparent")
ggarrange(genematphyseqancomid2ageplot,genematgllvmcoefageplot,common.legend = T,legend="right")
dev.off()

# clr before or after filter?
# clr before prune
cog2801clrbefore <- genematphyseq %>% transform(.,"clr") %>% filter_tax_table(Taxa %in% "COG2801") 

ggplot(genematphyseqcog2801coreclrpsmelt, aes(x=Abundance)) + geom_histogram()

genematphyseqcog2801coreclrpsmelt <- psmelt(cog2801clrbefore) %>% group_by(Sample) %>% mutate(sumabundance = sum(abs(Abundance), na.rm = T)) %>% distinct(Sample,.keep_all = T)

ggplot(genematphyseqcog2801coreclrpsmelt, aes(SamplingAge, sumabundance)) + geom_point()
genematphyseqcog2801coreclrpsmeltlm <- lmer(sumabundance ~ capage+TerminalYear + season + SampleYear + SexEstimate+ Timeinfridge_scaled + CatchTime_scaled + TQcorrected_scaled + (1|BirdID), data= genematphyseqcog2801coreclrpsmelt)
summary(genematphyseqcog2801coreclrpsmeltlm)

genematphyseqcog2801coreclr2 <- genematphyseq  %>% filter_tax_table(Taxa %in% "COG2801") %>% transform(.,"clr")
genematphyseqcog2801coreclrpsmelt2 <- psmelt(genematphyseqcog2801coreclr2) %>% group_by(Sample) %>% mutate(sumabundance = sum(abs(Abundance), na.rm = T)) %>% distinct(Sample,.keep_all = T)
genematphyseqcog2801coreclrpsmeltlm2 <- lmer(sumabundance ~ capage+TerminalYear + season + SampleYear + SexEstimate+ Timeinfridge_scaled + CatchTime_scaled + TQcorrected_scaled + (1|BirdID), data= genematphyseqcog2801coreclrpsmelt2)
summary(genematphyseqcog2801coreclrpsmeltlm2)
