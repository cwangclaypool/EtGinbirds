#++++++++++++
#+
#+ Analyses and figures for: Wang-Claypool et al.: The proof is in the plumage: dietary ethanol exposure in birds revealed by ethyl glucuronide in feathers
#+ Submitted to ProcB September 2023
#++++++++++++
#+Last updated: 25 July 2023
#+
#+Table of Contents:
#+
#+ Preparing data
#+ ANOVA
#+ Binomial logistic regression
#+ pglmm
#+ Figure 2
#+ Figure S1
#+ Figure S2
#+ Figure S3
#+ Figure S4

#++++++++++++
#+ 
#+ Preparing data
#+ 
#++++++++++++

# Load packages
library(readxl)
library(plyr)
library(tibble)
library(ggplot2)
library(cowplot)
library(ggtree)
library(ggtext)
library(ggpubr)
library(phytools)
library(EnvStats)
library(scales)
library(stats)
library(stringr)
library(tidyr)
library(phyr)
library(dplyr)
library(treeio)
library(ape)
library(aod)


xl_data <- "Dataset_FeatherEtG.xlsx"
excel_sheets(path = xl_data)
tab_names <- excel_sheets(path = xl_data)
list_all <- lapply(tab_names, function(x) read_excel(path = xl_data, sheet = x))

df <- list_all[1]
df <- data.frame(df)

# Filtering the data
df.postiveFEtGsubset <- subset(df, FiltF_Avg_AbsIncPres %in% c("Present")) # dataset without "Absent" and "Inconclusive" readings - only conclusive positives
df.postiveFEtGsubset$logFEtGconc <- log10(df.postiveFEtGsubset$FiltF_Avg_pgmg) #  log10 transformed EtG, only conclusive positive
df$logFEtGconc_with0s <- log10(df$FiltF_Avg_pgmg+1) # log10 transformed EtG with 0s included
df.noInconclusives <- subset(df, FiltF_Avg_AbsIncPres %in% c("Present", "Absent")) # filter out inconclusives

# Read in the data frame for phylogenetic test
df.phylotest <- list_all[2]
df.phylotest <- data.frame(df.phylotest)
row.names(df.phylotest)<-df.phylotest$Newick.label;

df.phylotest$logFEtGconc_with0s <- log10(df.phylotest$Avg.FiltF.Avg.pgmg+1) #+1 to all FEtG concentrations in order to include 0s in log10 transformation
df.phylotest.PostiveVibesOnly <- subset(df.phylotest, EtG.AbsIncPres %in% c("Present"))
df.phylotest.PostiveVibesOnly$logAvgFEtGconc <- log10(df.phylotest.PostiveVibesOnly$Avg.FiltF.Avg.pgmg) #log10 of FEtG concentrations without +1
df.phylotest.NoInconclusives <- subset(df.phylotest, EtG.AbsIncPres %in% c("Present", "Absent"))

# Read in tree
filename <- "Tree_FeatherEtG.nwk"
tree <- ape::read.tree(filename); phytools::read.newick(filename)


#++++++++++++
#+ 
###### ANOVA ######
#+ 
#++++++++++++

# 1. ANOVA test for differences between species with whole dataset (absents as 0s, inconclusives removed)
lm.model_noInconc_with0s <- lm(logFEtGconc_with0s ~ CommonName, data = df.noInconclusives)
summary(lm.model_noInconc_with0s)

# 2A. ANOVA test for differences in mean FEtG between trophic levels (absents as 0s, inconclusives removed)
lm.model.noInconc.with0s.TrophicLevel <- lm(logFEtGconc_with0s ~ Trophic.Level, data = df.noInconclusives)
summary(lm.model.noInconc.with0s.TrophicLevel)
aov.logEtG.noInconc <- ggplot(df.noInconclusives, aes(x = Trophic.Level, y = logFEtGconc_with0s))+
  geom_boxplot() +
  geom_jitter(width=0.2, height=0, size = 0.8) +
  scale_y_continuous(limits = c(-.25, 2)) +
  EnvStats::stat_n_text(size = 5) +
  ggtitle("EtG absent and present results\n (inconclusives dropped)")

# 2B. ANOVA test for differences in mean FEtG between trophic levels, but of only positive FEtG results
lm.model_without0s.TrophicLevel <- lm(logFEtGconc ~ Trophic.Level, data = df.postiveFEtGsubset)
summary(lm.model_without0s.TrophicLevel)
aov.logEtG.postitives <- ggplot(df.postiveFEtGsubset, aes(x = Trophic.Level, y = logFEtGconc))+
  geom_boxplot() +
  geom_jitter(width=0.2, height=0, size = 0.8) +
  scale_y_continuous(limits = c(-.25, 2)) +
  EnvStats::stat_n_text(size = 5) +
  ggtitle("EtG present results\n (absent and inconclusives dropped)")

plot_grid(aov.logEtG.noInconc, aov.logEtG.postitives, labels = c('A', 'B'), label_size = 30)


### Other tests
#ANOVA test for differences between nectarivores vs. non-nectarivores with whole dataset (Abs as 0s, inconclusives removed)
lm.model_with0s.NectNonNect <- lm(logFEtGconc_with0s ~ Nect.NonNect, data = df.noInconclusives)
summary(lm.model_with0s.NectNonNect)

# Tests where the results are driven by the FEtG absent results
lm.model_with0s.Habitat <- lm(logFEtGconc_with0s ~ Habitat, data = df.noInconclusives)
summary(lm.model_with0s.Habitat) #driven by marine which has 0s and inconclusives

lm.model_with0s.TrophicNiche <- lm(logFEtGconc_with0s ~ Trophic.Niche, data = df.noInconclusives)
summary(lm.model_with0s.TrophicNiche) #driven by frugivore which is 0

lm.model_with0s.PrimaryLifestyle<- lm(logFEtGconc_with0s ~ Primary.Lifestyle, data = df.noInconclusives)
summary(lm.model_with0s.PrimaryLifestyle) #driven by aquatic which has 0s

#++++++++++++
#+ 
###### Binomial logistic regression ######
#+ 
#++++++++++++

binomlogreg <- glm(FiltF_Avg_EtG_AbsIncPresBinom ~ Trophic.Level, family = "binomial", data = df.noInconclusives)
summary(binomlogreg)
exp(confint.default(binomlogreg, level = 0.95))
exp(cbind(OR = coef(binomlogreg), confint(binomlogreg)))


#++++++++++++
#+ 
###### pglmm ######
#+ 
#++++++++++++

###Trophic level not significant when analyzing all feather EtG results (absent and inconclusives as 0s), but significant when 0s removed and when log-transformed
bird.phylo.model.TrophicLevel1.AvgFEtG <- pglmm_compare(Avg.FiltF.Avg.pgmg~Trophic.Level, family = "gaussian", phy = tree, data = df.phylotest)
summary(bird.phylo.model.TrophicLevel1.AvgFEtG)
ggplot(df.phylotest, aes(x=`Trophic.Level`, y=`Avg.FiltF.Avg.pgmg`)) +
  geom_boxplot() + 
  EnvStats::stat_n_text() +
  geom_jitter(width=0.2, height=0, size = 0.8)

bird.phylo.model.TrophicLevel1.PositiveAvgFEtG <- pglmm_compare(Avg.FiltF.Avg.pgmg~Trophic.Level, family = "gaussian", phy = tree, data = df.phylotest.PostiveVibesOnly)
summary(bird.phylo.model.TrophicLevel1.PositiveAvgFEtG)
ggplot(df.phylotest.PostiveVibesOnly, aes(x=`Trophic.Level`, y=`Avg.FiltF.Avg.pgmg`)) +
  geom_boxplot() + 
  EnvStats::stat_n_text() +
  geom_jitter(width=0.2, height=0, size = 0.8)

bird.phylo.model.TrophicLevel1.PositiveFEtG <- pglmm_compare(logAvgFEtGconc~Trophic.Level, family = "gaussian", phy = tree, data = df.phylotest.PostiveVibesOnly)
summary(bird.phylo.model.TrophicLevel1.PositiveFEtG) #
ggplot(df.phylotest.PostiveVibesOnly, aes(x=`Trophic.Level`, y=`logAvgFEtGconc`)) +
  geom_boxplot() + 
  EnvStats::stat_n_text() +
  geom_jitter(width=0.2, height=0, size = 0.8)


####Below are other models with significant results, but with low confidence

#Trophic niche, average FEtG vs. log-transformed
bird.phylo.model.TrophicNiche.AvgFEtG  <- pglmm_compare(Avg.FiltF.Avg.pgmg ~ Trophic.Niche, family = "gaussian", phy = tree, data = df.phylotest)
summary(bird.phylo.model.TrophicNiche.AvgFEtG) 
bird.phylo.model.TrophicNiche.AvgFEtG.log  <- pglmm_compare(logFEtGconc_with0s ~ Trophic.Niche, family = "gaussian", phy = tree, data = df.phylotest)
summary(bird.phylo.model.TrophicNiche.AvgFEtG.log) #Frugivore and omnivores significant, but only n=1 for each niche
ggplot(df.phylotest, aes(x=`Trophic.Niche`, y=`logFEtGconc_with0s`)) +
  geom_boxplot() + 
  EnvStats::stat_n_text() + geom_jitter(width=0.2, height=0, size = 0.8)

#Trophic niche, log-transformed absence/presence FEtG only (inconclusives removed)
bird.phylo.model.NoIncLog10.3<- pglmm_compare(logFEtGconc_with0s~Trophic.Niche, family = "gaussian", phy = tree, data = df.phylotest.NoInconclusives) 
summary(bird.phylo.model.NoIncLog10.3) #granivore, nectarivore significant compared to frugivore, but only 1 frugivore
ggplot(df.phylotest.NoInconclusives, aes(x=`Trophic.Niche`, y=`logFEtGconc_with0s`)) +
  geom_boxplot()+ 
  EnvStats::stat_n_text() +
  geom_jitter(width=0.2, height=0, size = 0.8)

#[Nectarivores+frugivores] vs. granivore vs. all other niches, average FEtG vs. log-transformed
bird.phylo.model.NectFrugGranOther.AvgFEtG <- pglmm_compare(Avg.FiltF.Avg.pgmg ~ NectFrug.Gran.Other, family = "gaussian", phy = tree, data = df.phylotest)
summary(bird.phylo.model.NectFrugGranOther.AvgFEtG)
bird.phylo.model.NectFrugGranOther.AvgFEtG.log <- pglmm_compare(logFEtGconc_with0s ~ NectFrug.Gran.Other, family = "gaussian", phy = tree, data = df.phylotest)
summary(bird.phylo.model.NectFrugGranOther.AvgFEtG.log) #significant differences, but likely driven by 0s in [Nectarivores+Frugivores] and Other
ggplot(df.phylotest, aes(x=`NectFrug.Gran.Other`, y=`logFEtGconc_with0s`)) +
  geom_boxplot()+ 
  EnvStats::stat_n_text() +
  geom_jitter(width=0.2, height=0, size = 0.8)


#++++++++++++
#+ 
###### Figure 2 ######
#+ 
#++++++++++++
cols <- c("Aquatic predator" = "#0072B2",
          "Frugivore" = "#660099", 
          "Granivore" = "#FFC425", 
          "Invertivore" = "#333333", 
          "Nectarivore" = "#FF00FF", 
          "Omnivore" = "#663300", 
          "Vertivore" = "#CC0000")
cols <- data.frame(Trophic.Niche = names(cols), color = cols)
df <- merge(df, cols, by = "Trophic.Niche", all.x = TRUE)

df$CommonName2 <- paste0("<span style=\"color: ", df$color, "\">", df$CommonName, "</span>")

species_list2 <- c("<span style=\"color: #333333\">Yellow-rumped Warbler</span>",
                  "<span style=\"color: #663300\">Golden-crowned Sparrow</span>",
                  "<span style=\"color: #FFC425\">Dark-eyed Junco</span>",
                  "<span style=\"color: #333333\">American Robin</span>",
                  "<span style=\"color: #333333\">Hermit Thrush</span>",
                  "<span style=\"color: #333333\">House Wren</span>",
                  "<span style=\"color: #333333\">Bewick's Wren</span>",
                  "<span style=\"color: #660099\">Cedar Waxwing</span>",
                  "<span style=\"color: #333333\">Barn Swallow</span>" ,
                  "<span style=\"color: #FFC425\">Acorn Woodpecker</span>",
                  "<span style=\"color: #CC0000\">White-tailed Kite</span>",
                  "<span style=\"color: #FF00FF\">Costa's Hummingbird</span>",
                  "<span style=\"color: #FF00FF\">Anna's Hummingbird</span>",
                  "<span style=\"color: #FF00FF\">Black-chinned Hummingbird</span>",
                  "<span style=\"color: #FF00FF\">Ruby-throated Hummingbird</span>",
                  "<span style=\"color: #0072B2\">Common Murre</span>",
                  "<span style=\"color: #FFC425\">Band-tailed Pigeon</span>")

ggplot(data = df, aes(x = factor(CommonName2, level= species_list2), y = logFEtGconc_with0s)) +
  geom_boxplot(coef=NULL)+ #default is 1.5 (1.5x interquartile range), null = whiskers go to min/max values
  geom_jitter(width=0.2, height=0, size = 0.8) + 
  theme(plot.title = element_text(hjust=0.5),
        axis.text.y = ggtext::element_markdown(size = 8, face = "bold"), #https://stackoverflow.com/questions/72749018/color-axis-text-to-match-grouping-variable-in-ggplot
        axis.text.x = element_text(size = 8),
        axis.title = element_text(size = 8),
        legend.position = "right") +
  EnvStats::stat_n_text(size = 2)+
  EnvStats::stat_mean_sd_text(size = 2, hjust = 1)+
  labs(x="Species", y="log10 Feather EtG Concentration (pg/mg)")+
  geom_hline(yintercept=log10(30), linetype="dashed", color = "red") +
  geom_text(aes(0,log10(30),label = "log10(30)", vjust = -0.5, hjust = 1), size = 8*0.36) + #why multiple geom_text by 0.36: https://stackoverflow.com/questions/25061822/ggplot-geom-text-font-size-control#comment121955273_25062509
  coord_flip()

#++++++++++++
#+ 
###### Figure S1 ######
#+ 
#++++++++++++

df.flightbodysubset <- subset(df, PreparatorID %in% c("CYWC113","CYWC114","CYWC121","CYWC123")) #just the samples where flight and body feathers were separated
df.flightbodysubset2 <- df.flightbodysubset[c("PreparatorID","FEtGTest1_conc_AbsPres","FEtGTest2_conc_AbsPres")]
colnames(df.flightbodysubset2) <- c("PreparatorID", "Body Feathers","Flight Feathers") # change column names of all the columns in the dataframe print(df)
df.flightbodysubsetreformat <- gather(df.flightbodysubset2, feathertype, concentration, "Body Feathers":"Flight Feathers") #Create long format
df.flightbodysubsetreformat$concentration <- as.numeric(df.flightbodysubsetreformat$concentration)


plotflightbodysubset <- ggplot(df.flightbodysubsetreformat, aes(PreparatorID, concentration, fill=feathertype)) +
  geom_bar(stat = "identity", position = 'dodge',  colour="black") +
  theme(plot.title=element_text(hjust=0.5)) +
  labs(x = "Preparator Number", y = "Concentration (pg/mg)") +
  scale_fill_manual(name = "Feather Type", labels = c("Body","Flight"), values=c("#00BFC4", "#F8766D")) + 
  theme(axis.text.x = element_text(vjust = 1, hjust=1, size = 12),         
        strip.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        axis.title = element_text(size = 18), 
        legend.text=element_text(size=18),
        legend.title=element_text(size=18),
        legend.position=c(.2,.9))+ # Inset axis position
  annotate(geom="text", x=3.25, y=2, label="*",
           color="black", size = 18) + 
  annotate(geom="text", x=.77, y=2, label="*",
           color="black", size = 18)

#Rachis+calamus versus vane+body feathers
df.rachisvanesubset <- subset(df, PreparatorID %in% c("CYWC097","CYWC098")) #just the samples where flight and body feathers were separated
df.rachisvanesubset2 <- df.rachisvanesubset[c("PreparatorID","FEtGTest1_conc_AbsPres","FEtGTest2_conc_AbsPres")]
colnames(df.rachisvanesubset2) <- c("PreparatorID", "Rachis+Calamus","Vane+Body") # change column names of all the columns in the dataframe print(df)
df.rachisvanesubset2reformat <- gather(df.rachisvanesubset2, feathertype, concentration, "Rachis+Calamus":"Vane+Body") #Create long format
df.rachisvanesubset2reformat$concentration <- as.numeric(df.rachisvanesubset2reformat$concentration)

plotrachisvanesubset <- ggplot(df.rachisvanesubset2reformat, aes(PreparatorID, concentration, fill=feathertype)) +
  geom_bar(stat = "identity", position = 'dodge',  colour="black") +
  labs(x = "Preparator Number", y = "Concentration (pg/mg)") +
  scale_fill_manual(name = "Feather Type", labels = c("Rachis+Calamus","Vane+Body"), values=c("#00BFC4", "#F8766D")) + 
  theme(axis.text.x = element_text(vjust = 0, hjust=0, size = 12),         
        strip.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        axis.title = element_text(size = 18), 
        legend.text=element_text(size=18),
        legend.title=element_text(size=18),
        legend.position=c(.25,.9))+ # Inset axis position
  ylim(0, 12) + 
  annotate(geom="text", x=1.22, y=.25, label="*",
           color="black", size = 18)

plot_grid(plotflightbodysubset, plotrachisvanesubset, labels = c('A', 'B'), label_size = 20, ncol = 1)

#++++++++++++
#+ 
###### Figure S2 ######
#+ 
#++++++++++++
###liver samples excluding the problematic ones?
df.allliversubset <- subset(df, L_Dataset %in% c("one_reading")) #just the samples where both readings worked
#view(df.allliversubset)
df.allliversubsetreformat <- gather(df.allliversubset, testtype, concentration, FiltF_Avg_pgmg:LEtG_ngg_plotting) #Create long format
#view(df.allliversubsetreformat)
df.allliversubsetreformat <- subset(df.allliversubsetreformat, testtype %in% c("LEtG_ngg_plotting","FiltF_Avg_pgmg"))
df.allliversubsetreformat$concentration <- as.numeric(df.allliversubsetreformat$concentration)

#str(df.allliversubsetreformat$concentration)
plotAllLivers <- ggplot(df.allliversubsetreformat, aes(PreparatorID, concentration, fill=testtype)) +
  geom_bar(stat = "identity", position = "dodge", colour="black") + #color = black to get black outline around bars: https://stackoverflow.com/questions/23969153/how-to-display-0-value-in-a-bar-chart-using-ggplot2
  scale_y_continuous(  #dual Y axes: https://r-graph-gallery.com/line-chart-dual-Y-axis-ggplot2.html
    name = "Feather EtG Concentration\n (pg/mg)", # Features of the first axis
    sec.axis = sec_axis( trans=~.*1, name="Liver EtG Concentration\n (ng/g)")) + # Add a second axis and specify its features
  theme(plot.title=element_text(hjust=0.5)) +
  labs(x = "Preparator Number by Species") +
  scale_fill_manual(name = "Tissue", labels = c("Feather", "Liver"), values=c("#00BFC4", "#F8766D")) + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = 12),
        strip.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        axis.title = element_text(size = 18)) + 
  facet_grid(~CommonName,  #multi-level labels in R: https://dmitrijskass.netlify.app/2019/06/30/multi-level-labels-with-ggplot2/
             scales = "free_x", # Let the x axis vary across facets.
             space = "free_x",  # Let the width of facets vary and force all bars to have the same width.
             switch = "x",     # Move the facet labels to the bottom.
             labeller = labeller(CommonName = label_wrap_gen(width = 10))) + #wrap facet label in r: https://datavizpyr.com/wrap-really-long-facet-labels-ggplot2/
  theme(strip.placement = "outside",                      # Place facet labels outside x axis labels.
        strip.background = element_rect(fill = "white"),  # Make facet label background white.
        legend.position=c(.9,.9), # Inset axis position
        legend.text=element_text(size=18),
        legend.title=element_text(size=18))
plotAllLivers

#++++++++++++
#+ 
###### Figure S3 ######
#+ 
#++++++++++++
#df2 <- list_all[1]
#df2 <- data.frame(df2)

se <- function(x) sqrt(var(x)/length(x))

assigningmymean <- df %>% group_by(CommonName) %>% dplyr::summarize(my_mean = mean(FiltF_Avg_pgmg))
assigningmyse <- df %>% group_by(CommonName) %>% dplyr::summarize(my_se = se(FiltF_Avg_pgmg))
assigningn <- df %>% group_by(CommonName) %>% dplyr::summarize(n= n())
merged.df <- left_join(df,assigningmymean, by = "CommonName")
merged.df2 <- left_join(merged.df,assigningmyse, by = "CommonName")
merged.df3 <- left_join(merged.df2,assigningn, by = "CommonName")

species_list <- c("Yellow-rumped Warbler",
                  "Golden-crowned Sparrow",
                  "Dark-eyed Junco",
                  "American Robin",
                  "Hermit Thrush",
                  "House Wren",
                  "Bewick's Wren",
                  "Cedar Waxwing",
                  "Barn Swallow" ,
                  "Acorn Woodpecker",
                  "White-tailed Kite",
                  "Costa's Hummingbird",
                  "Anna's Hummingbird",
                  "Black-chinned Hummingbird",
                  "Ruby-throated Hummingbird",
                  "Common Murre",
                  "Band-tailed Pigeon")


#merged barplot and scatter plot
EtG.barscatterplot <- ggplot(merged.df3, aes(x = factor(CommonName, level= species_list), y=FiltF_Avg_pgmg)) + 
  geom_bar(data = merged.df3,
           aes(y = my_mean, x = factor(CommonName, level= species_list),
               ymin = my_mean - my_se,
               ymax = my_mean + my_se), 
           stat="summary", width=0.75,
           fill="darkgrey") + 
  geom_errorbar(data = merged.df3,
                aes(y = my_mean, x = factor(CommonName, level= species_list),
                    ymin = my_mean - my_se,
                    ymax = my_mean + my_se), width=0.75) + 
  geom_point(data = merged.df3, aes(y = FiltF_Avg_pgmg, x = factor(CommonName, level= species_list))) +
  theme(text = element_text(size = 20))+
  scale_x_discrete(labels = function(x) str_wrap(x, width = 15))+
  theme(plot.title=element_text(hjust=0.5)) +
  #EnvStats::stat_n_text(size = 5) +
  geom_text(aes(label = paste0("n=", n), y = -1.5), check_overlap = TRUE, hjust = .75, size = 14*0.36) +
  #ggtitle("Feather EtG concentrations across\n species (with 0s)") + 
  labs(x="Species", y="Feather EtG Concentration (pg/mg)")+
  scale_y_continuous(breaks = scales::pretty_breaks(n = 5)) +
  geom_hline(yintercept=30, linetype="dashed", color = "red") +
  geom_text(aes(0,30,label = "30 pg/mg", vjust = -0.5, hjust = -.1), size = 18*0.36) +
  coord_flip()
EtG.barscatterplot

#++++++++++++
#+ 
###### Figure S4 ######
#+ 
#++++++++++++
monthorder <- c('Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec')
df$daymonth <- factor(format(df$`DateCollected`,"%b"), levels = monthorder)

pd <- position_dodge(.5)
monthspp <- ggplot(df, aes(daymonth, FiltF_Avg_pgmg)) +
  geom_point(aes(color = CommonName), position = pd, size = 2)+
  theme(plot.title=element_text(hjust=0.5), legend.justification = "left") +
  #ggtitle("Feather EtG (pg/mg) during a calendar year") + 
  EnvStats::stat_n_text(size = 3.4)+
  scale_color_manual(name="Species", labels = function(x) str_wrap(x, width = 15), values = c('#3cb44b', '#000000','#e6194b', '#4363d8', '#f58231', '#911eb4', '#46f0f0', '#f032e6', '#bcf60c', '#fabebe', '#008080', '#e6beff', '#9a6324', '#ffe119', '#800000', '#aaffc3', 'deeppink')) +
  labs(x="Calendar Month", y="Feather EtG Concentration (pg/mg)") +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 9)) +
  geom_text(aes(0,30,label = "30 pg/mg", hjust = -12.5, vjust = -.5), color = "black") +
  geom_hline(yintercept=30, linetype="dashed", color = "red")

ANHUdf <- subset(df, CommonName %in% c("Anna's Hummingbird"))
ANHUmonth <- ggplot(ANHUdf, aes(daymonth, FiltF_Avg_pgmg)) +
  geom_point(aes(color = CommonName), position = pd, size = 2)+
  theme(plot.title=element_text(hjust=0.5), legend.justification = "left") +
  #ggtitle("ANHU Feather EtG (pg/mg) during a calendar year") + 
  EnvStats::stat_n_text(size = 3.4)+
  scale_color_discrete(name="Species", labels = function(x) str_wrap(x, width = 15)) +
  labs(x="Calendar Month", y="Feather EtG Concentration (pg/mg)") +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 9)) +
  geom_text(aes(0,30,label = "30 pg/mg", hjust = -12.5, vjust = -.5), color = "black") +
  geom_hline(yintercept=30, linetype="dashed", color = "red")

plot_grid(monthspp, ANHUmonth, ncol = 1, align = "v", labels = c('A', 'B'), label_size = 20) 
