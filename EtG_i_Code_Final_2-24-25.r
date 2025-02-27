#++++++++++++
#+
#+ Analyses and figures for: The proof is in the plumage: dietary ethanol exposure in birds revealed by ethyl glucuronide in feathers
#+ #++++++++++++
library(ape)
library(car)
library(ggplot2)
library(phyr)
library(gridExtra)
library(ggthemes) # Load to get plots with a more simple design

#This includes black axis text and black axis lines
new_theme <- theme_few()  + theme(text = element_text(size = 12))+ theme(legend.position = "none")+ theme(plot.title = element_text(hjust = 0))+ theme(axis.text.x = element_text(color="black"), axis.text.y = element_text(color="black"), panel.border = element_rect(color = "black", fill=NA), axis.ticks = element_line(color="black"))



#Only three hummingbirds out of 366 species in all of Avonet coded as omnivores: Archilochus colubris, Lophornis magnificus, Lophornis ornatus.  Seems odd not to have them all as herbivores as all hummingbirds have insects in their diet. Recoded A. colubris as a herbivore.
#Read in the data.
df <- read.table("EtG_Dataset_A_colubris_as_herbivore.txt", header=TRUE)
df <- data.frame(df)



# Filtering the data
df.postiveFEtGsubset <- subset(df, FiltF_Avg_AbsIncPres %in% c("Present")) # dataset without "Absent" and "Inconclusive" readings - only conclusive positives
df.postiveFEtGsubset$logFEtGconc <- log10(df.postiveFEtGsubset$FiltF_Avg_pgmg) #  log10 transformed EtG, only conclusive positive
df$logFEtGconc_with0s <- log10(df$FiltF_Avg_pgmg+1) # log10 transformed EtG with 0s included
df.noInconclusives <- subset(df, FiltF_Avg_AbsIncPres %in% c("Present", "Absent")) # filter out inconclusives


# Read in tree
filename <- "Tree_EtG_in_Bird_Feathers.nwk"
tree <- ape::read.tree(filename); phytools::read.newick(filename)
#check if the phylogeny and dataframe match in species
setdiff(tree$tip.label,unique(df$species.in.tree))

#++++++++++++

#Check the distribution of the data!
#Assessing the distribution of the prediction data
no.inconclusives.qq.plot <- ggplot(df.noInconclusives, aes(sample = logFEtGconc_with0s)) + stat_qq() + stat_qq_line()
no.inconclusives.histogram <- ggplot(df.noInconclusives, aes(x = logFEtGconc_with0s)) + geom_histogram()

df.postiveFEtGsubset.qq.plot <- ggplot(df.postiveFEtGsubset, aes(sample = logFEtGconc)) + stat_qq() + stat_qq_line()
df.postiveFEtGsubset.histogram <- ggplot(df.postiveFEtGsubset, aes(x = logFEtGconc)) + geom_histogram(bins=10)

grid.arrange(no.inconclusives.qq.plot, no.inconclusives.histogram, df.postiveFEtGsubset.qq.plot, df.postiveFEtGsubset.histogram, nrow=2)

#Any relationship between the EtG measured and the sample weight?  No clear relationship
ggplot(df.noInconclusives, aes(x = FiltF_Avg_Wtmg, y = FiltF_Avg_pgmg)) + geom_point()

#-------------------------

#Making the phylogenies for the tests.
#Have to drop species not included in any particular data set.

#Data with zeros tree
species.to.drop <- setdiff(tree$tip.label,unique(df.noInconclusives$species.in.tree))
tree.noInconclusives <- drop.tip(tree, species.to.drop)
#Data with no zeros tree
species.to.drop2 <- setdiff(tree$tip.label,unique(df.postiveFEtGsubset$species.in.tree))
tree.postiveFEtGsubset <- drop.tip(tree, species.to.drop2)

#-------------------
#Nonphylogenetic tests to assess whether there was variation among species in levels of EtG

# 1. ANOVA test for differences between species with whole dataset (absents as 0s, inconclusives removed)
#Reported result
lm.model_noInconc_with0s <- lm(logFEtGconc_with0s ~ CommonName, data = df.noInconclusives)
summary(lm.model_noInconc_with0s)
Anova(lm.model_noInconc_with0s)

# 2. ANOVA test for differences between species with only positive EtG samples
#Reported result
lm.model_positives <- lm(logFEtGconc ~ CommonName, data = df.postiveFEtGsubset)
summary(lm.model_positives)
Anova(lm.model_positives)


#-------------------
#Phylogenetically corrected tests
#Specified an overall phylogenetic effect using species.in.tree__, which also automatically includes a nonphylogenetic i.i.d. effect of species.
#https://daijiang.github.io/phyr/articles/phyr_example_empirical.html
#I think that we need to go with the REML version of the tests          



#Test for differences in mean FEtG between nectarivorous hummingbirds and all other species.
#Look at distribution of the data
with(df.noInconclusives, boxplot(logFEtGconc_with0s ~ Nect.NonNect))
#Reported result
test1 = phyr::pglmm(logFEtGconc_with0s ~ Nect.NonNect + FiltF_Avg_Wtmg +(1|species.in.tree__) , data = df.noInconclusives, family = "gaussian", REML = TRUE, cov_ranef = list(species.in.tree = tree.noInconclusives))
summary(test1)
fixef(test1)
ranef(test1)   

#Test for differences in mean FEtG between nectarivorous hummingbirds and all other species.
#With only positive EtG samples
test2 = phyr::pglmm(logFEtGconc ~ Nect.NonNect + FiltF_Avg_Wtmg + (1|species.in.tree__) , data = df.postiveFEtGsubset, family = "gaussian", REML = TRUE, cov_ranef = list(species.in.tree = tree.postiveFEtGsubset))
summary(test2) 

#Test for nectarivores + frugivores versus other birds
#Reported result
test1 = phyr::pglmm(logFEtGconc_with0s ~ NectFrug.Other + FiltF_Avg_Wtmg +(1|species.in.tree__) , data = df.noInconclusives, family = "gaussian", REML = TRUE, cov_ranef = list(species.in.tree = tree.noInconclusives))
summary(test1)
   
#Test for nectarivores + frugivores versus other birds
#With only positive EtG samples
test2 = phyr::pglmm(logFEtGconc ~ NectFrug.Other + FiltF_Avg_Wtmg + (1|species.in.tree__) , data = df.postiveFEtGsubset, family = "gaussian", REML = TRUE, cov_ranef = list(species.in.tree = tree.postiveFEtGsubset))
summary(test2) 


#Test for differences in mean FEtG between nectarivores + frugivores, granivores, and all other species.

with(df.noInconclusives, boxplot(logFEtGconc_with0s ~ NectFrug.Gran.Other))

#Trophic level test with all data including zeros (absents as 0s, inconclusives removed)
#Reported result
test1 = phyr::pglmm(logFEtGconc_with0s ~ NectFrug.Gran.Other + FiltF_Avg_Wtmg + (1|species.in.tree__) , data = df.noInconclusives, family = "gaussian", REML = TRUE, cov_ranef = list(species.in.tree = tree.noInconclusives))
summary(test1)   

#Test for differences in mean FEtG between nectarivores + frugivores, granivores, and all other species.
#Trophic level test with only positive EtG samples
test2 = phyr::pglmm(logFEtGconc ~ NectFrug.Gran.Other + FiltF_Avg_Wtmg + (1|species.in.tree__) , data = df.postiveFEtGsubset, family = "gaussian", REML = TRUE, cov_ranef = list(species.in.tree = tree.postiveFEtGsubset))
summary(test2) 



#Test for differences in mean FEtG between trophic levels 

with(df.noInconclusives, boxplot(logFEtGconc_with0s ~ Trophic.Level))


#Trophic level test with all data with zeros (absents as 0s, inconclusives removed)
#Reported result
test1 = phyr::pglmm(logFEtGconc_with0s ~ Trophic.Level + FiltF_Avg_Wtmg + (1|species.in.tree__) , data = df.noInconclusives, family = "gaussian", REML = TRUE, cov_ranef = list(species.in.tree = tree.noInconclusives))
summary(test1)   

#Trophic level test with only positive EtG samples
test2 = phyr::pglmm(logFEtGconc ~ Trophic.Level + FiltF_Avg_Wtmg + (1|species.in.tree__) , data = df.postiveFEtGsubset, family = "gaussian", REML = TRUE, cov_ranef = list(species.in.tree = tree.postiveFEtGsubset))
summary(test2)  
                   
#-------------
#Binomial model tests
#Try the tests with samples scored with presence/absence of EtG 
test3 = phyr::pglmm(FiltF_Avg_EtG_AbsIncPresBinom ~ Nect.NonNect + FiltF_Avg_Wtmg +  (1|species.in.tree__) , data = df.noInconclusives, family = "binomial", REML = TRUE, cov_ranef = list(species.in.tree = tree.noInconclusives))
summary(test3)  

test3 = phyr::pglmm(FiltF_Avg_EtG_AbsIncPresBinom ~ NectFrug.Other + FiltF_Avg_Wtmg +  (1|species.in.tree__) , data = df.noInconclusives, family = "binomial", REML = TRUE, cov_ranef = list(species.in.tree = tree.noInconclusives))
summary(test3)  

test3 = phyr::pglmm(FiltF_Avg_EtG_AbsIncPresBinom ~ NectFrug.Gran.Other + FiltF_Avg_Wtmg +  (1|species.in.tree__) , data = df.noInconclusives, family = "binomial", REML = TRUE, cov_ranef = list(species.in.tree = tree.noInconclusives))
summary(test3) 

#The binomial test won't run because of lack of variation in the dependent variable.
test3 = phyr::pglmm(FiltF_Avg_EtG_AbsIncPresBinom ~ Trophic.Level + FiltF_Avg_Wtmg+ (1|species.in.tree__) , data = df.noInconclusives, family = "binomial", REML = TRUE, cov_ranef = list(species.in.tree = tree.noInconclusives))
summary(test3) 

#--------------------------------------
#Figure 2
#Order the factor levels
df$CommonName2 = factor(df$CommonName, levels = c("Yellow-rumped_Warbler", "Dark-eyed_Junco", "Golden-crowned_Sparrow" , "American_Robin" , "Hermit_Thrush", "House_Wren" , "Bewicks_Wren" ,"Cedar_Waxwing",  "Barn_Swallow", "Acorn_Woodpecker" , "White-tailed_Kite", "Costas_Hummingbird","Annas_Hummingbird","Black-chinned_Hummingbird" ,"Ruby-throated_Hummingbird",  "Common_Murre",   "Band-tailed_Pigeon"))


ggplot(data = df, aes(y = factor(CommonName2), x = logFEtGconc_with0s)) + geom_boxplot(coef=NULL) + geom_jitter(width=0, height=0.2, size = 0.8, color="black") + labs(y="Species", x="log10 Feather EtG Concentration (pg/mg)")+ geom_vline(xintercept=log10(30), linetype="dashed", color = "red") 

#Preliminary plot
ggplot(data = df, aes(x = factor(CommonName2), y = logFEtGconc_with0s)) +
  geom_boxplot(coef=NULL)+ #default is 1.5 (1.5x interquartile range), null = whiskers go to min/max values
  geom_jitter(width=0.2, height=0, size = 0.8) + 
  theme(plot.title = element_text(hjust=0.5),
        
        axis.text.x = element_text(size = 8),
        axis.title = element_text(size = 8),
        legend.position = "right") +
  EnvStats::stat_n_text(size = 2)+
  EnvStats::stat_mean_sd_text(size = 2, hjust = 1)+
  labs(x="Species", y="log10 Feather EtG Concentration (pg/mg)")+
  geom_hline(yintercept=log10(30), linetype="dashed", color = "red") +
  geom_text(aes(0,log10(30),label = "log10(30)", vjust = -0.5, hjust = 1), size = 8*0.36) + #why multiple geom_text by 0.36: https://stackoverflow.com/questions/25061822/ggplot-geom-text-font-size-control#comment121955273_25062509
  coord_flip() + 
  scale_x_discrete(labels = c("Yellow-rumped Warbler",  "Dark-eyed Junco","Golden-crowned Sparrow" , "American Robin" , "Hermit Thrush", "House Wren" , "Bewick's Wren" ,"Cedar Waxwing",  "Barn Swallow", "Acorn Woodpecker" , "White-tailed Kite", "Costa's Hummingbird","Anna's Hummingbird","Black-chinned Hummingbird" ,"Ruby-throated Hummingbird",  "Common Murre",   "Band-tailed Pigeon"))
  
  #Final plot
 final.plot <- ggplot(data = df, aes(x = factor(CommonName2), y = logFEtGconc_with0s)) +
  geom_boxplot(coef=NULL, fill="light gray")+ #default is 1.5 (1.5x interquartile range), null = whiskers go to min/max values
  geom_jitter(width=0.2, height=0, size = 0.8) + 
  theme(plot.title = element_text(hjust=0.5),
        
        axis.text.x = element_text(size = 8),
        axis.title = element_text(size = 8),
        legend.position = "right") +
  EnvStats::stat_n_text(size = 2)+
  EnvStats::stat_mean_sd_text(size = 2, hjust = 1)+
  labs(x="Species", y="log10 Feather EtG Concentration (pg/mg)")+
  geom_hline(yintercept=log10(30), linetype="dashed", color = "red") +
  geom_text(aes(0,log10(30),label = "log10(30)", vjust = -0.5, hjust = 1), size = 8*0.36) + #why multiple geom_text by 0.36: https://stackoverflow.com/questions/25061822/ggplot-geom-text-font-size-control#comment121955273_25062509
  coord_flip() +new_theme + 
   scale_x_discrete(labels = c("Yellow-rumped Warbler",  "Dark-eyed Junco","Golden-crowned Sparrow" , "American Robin" , "Hermit Thrush", "House Wren" , "Bewick's Wren" ,"Cedar Waxwing",  "Barn Swallow", "Acorn Woodpecker" , "White-tailed Kite", "Costa's Hummingbird","Anna's Hummingbird","Black-chinned Hummingbird" ,"Ruby-throated Hummingbird",  "Common Murre",   "Band-tailed Pigeon"))
  
  
dev.new(width=5.5, height=7)
pdf("Fig. 2 Avian species and feather EtG concentration.pdf", width=5.5, height=7)
final.plot
dev.off()

