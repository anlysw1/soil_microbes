# Author: Nicole Pietrasiak and Angie Swanson
# Title: SOIL476L data analysis

library(ggplot2)
library(Rmisc)
library(ggpubr)
library(reshape)
library(vegan)
library(RColorBrewer)
library(kableExtra)


# Set the working directory by hitting Ctrl + Shift + H or
setwd("~/Documents/#grad school/##thesis/R data analysis")

# Load the data and check that it loads correctly. 
AS_data <- read.csv("FINAL_DATA_ASthesis.csv")
head(AS_data,2)
dim(AS_data)
str(AS_data)

# Rename the first variable name, adds in a strange character for some reason.
names(AS_data)[1] <- "Site"
names(AS_data)
AS_data$Site <- factor(AS_data$Site)
AS_data$Depth <- factor(AS_data$Depth)
AS_data$Depth_string <- factor(AS_data$Depth_string)
AS_data$Site_depth <- factor(AS_data$Site_depth)
head(AS_data,2)

# Creating a new data set to
# 1) Reorder the sites to be Desert, Turf, Farm
# 2) Add the Adjusted Gram Positive to the final data set for averages

P_AdGramPos <- AS_data$P_GramPos - AS_data$P_Actin

AS_data2 <- cbind(AS_data[,1:6], P_AdGramPos, AS_data[,c(7,9:32)])
names(AS_data2)[7] <- "P_AdGramPos"

names(AS_data2)
str(AS_data2)

# Rename the factors (and their levels) to be in the order I want them to be
AS_data2$Site <- c(rep(1,6),rep(3,6),rep(2,6))
AS_data2$Site <- factor(AS_data2$Site)
levels(AS_data2$Site) <- c("Desert", "Turf Path", "Farm")

AS_data2$Depth_string <- rep(c(rep(1,3),rep(2,3)),3)
AS_data2$Depth_string <- factor(AS_data2$Depth_string)
levels(AS_data2$Depth_string) <- c("surf", "sub")

levels(AS_data2$Depth) <- c("0-1cm", "1-5cm")

AS_data2$Site_depth <- c(1,1,1,2,2,2,5,5,5,6,6,6,3,3,3,4,4,4)
AS_data2$Site_depth <- factor(AS_data2$Site_depth)
levels(AS_data2$Site_depth) <- c("Desert_surf","Desert_sub",
                                 "Turf_surf", "Turf_sub",
                                 "Farm_surf", "Farm_sub")

mydata <- AS_data2[order(AS_data2$Site),]
str(mydata)
head(mydata, 10)


#^#^#^#^#^#^#^#^#^#^#^#^#^#^#^#^#^#^#^#^#^#^#^
# Making confidence intervals for my data
# ?confint()
# The confint method uses a linear model object to generate CIs
# I'm going to need to use the formula for an unknown pop. std dev
# from Dr. Trainor's Lecture 8: Inferences about the Mean
# y_bar +/- t(alpha/2, n-1) * s/sqrt(n)
# Where y_bar is the average for the observations within each variable
# and n is the number of observations per group (site*depth), n=3
# After googling, it looks like I'll need to make my own function for this...
#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#

# Confidence Intervals - Environmental Properties
names(mydata)
head(mydata)
str(mydata)
envm <- mydata[,c(1,3,13,18:32)]
# str(envm)
# dim(envm)
summary(mydata)

env.means <- aggregate(envm[,sapply(envm, is.numeric)], 
                       with(envm, list(Site, Depth)), 
                       mean)
names(env.means)[1] <- "Site"
names(env.means)[2] <- "Depth"

ord<-order(env.means$Site)
env.means <- env.means[ord,]
row.names(env.means) <- NULL

env.stdv <- aggregate(envm[,sapply(envm, is.numeric)], 
                      with(envm, list(Site, Depth)), 
                      sd)
names(env.stdv)[1] <- "Site"
names(env.stdv)[2] <- "Depth"
env.stdv <- env.stdv[ord,]
row.names(env.stdv) <- NULL

t <- qt(p=0.95,df=2)
upper.env <- env.means[,3:18] + (t * env.stdv[,3:18]/sqrt(3))
upper.env <- cbind(env.means[,1:2], upper.env)

lower.env <- env.means[,3:18] - (t * env.stdv[,3:18]/sqrt(3))
lower.env <- cbind(env.means[,1:2], lower.env)



env.CI <- rbind(lower.env, upper.env)
ord1 <- order(env.CI$Depth)
env.CI <- env.CI[ord1,]
ord2 <- order(env.CI$Site)
env.CI <- env.CI[ord2,]
row.names(env.CI) <- NULL
env.CI

# write.csv(env.CI, file = "env_CI.csv", row.names = FALSE)



# Using Andrew Dominguez's ggplot2 code to make some plots for the env data
library(ggplot2)
library(Rmisc)
library(ggpubr)

# Set the graphing parameters that will be carried over to all plots in 
# this group.
point <- 20
c.siz <- 0.3
b.siz <- 10.5

# Colors correspond to Desert, Farm, and Turf
col_sites <- c("blue3", "green3", "sienna2")
p_color <- "Spectral"

# *** Physical soil properties ***

# Dry aggregates MWD
MWDp <- ggplot(mydata, aes(x = Site, y = MWD, fill = Depth)) + 
  geom_boxplot() + theme_light(base_size = b.siz) + 
  scale_fill_brewer(palette = p_color) +
  theme(axis.title.x = element_blank()) +
  labs(y="Mean Weight Diameter (mm)") +
  scale_y_continuous(breaks = round(seq(0,max(mydata$MWD), by = 0.5),1))
MWDp

# Dry aggregates percent fine diameter aggregates
FDAp <- ggplot(mydata, aes(x = Site, y = P_FDA, fill = Depth)) +
  geom_boxplot() + theme_light(base_size = b.siz) +
  scale_fill_brewer(palette = p_color) +
  theme(axis.title.x = element_blank()) +
  labs(y= "% <0.25 mm Diameter") +
  scale_y_continuous(breaks = round(seq(0,100, by = 20),1))
FDAp

# Wet aggregates - percent stable aggregation from 2mm aggregates
WAGp <- ggplot(mydata, aes(x = Site, y = P_wagg, fill = Depth)) +
  geom_boxplot() + theme_light(base_size = b.siz) +
  scale_fill_brewer(palette = p_color) +
  theme(axis.title.x = element_blank()) +
  labs(y= "% Water Stable") +
  scale_y_continuous(breaks = round(seq(-20, 80, by = 20),0))
WAGp

# Does the % wagg differ by depth for Turf?
# t.test(x = AS_data$P_wagg[13:15], y = AS_data$P_wagg[16:18])$p.value
# [1] 2.724549e-05

# Wet aggregates - percent stone of 2 mm aggregates
Stonep <- ggplot(mydata, aes(x = Site, y = P_stone, fill = Depth)) +
  geom_boxplot() + theme_light(base_size = b.siz) +
  scale_fill_brewer(palette = p_color) +
  theme(axis.title.x = element_blank()) +
  labs(y= "% Stone") +
  scale_y_continuous(breaks = round(seq(-10,100, by = 20),0))
Stonep

clayp <- ggplot(mydata, aes(x = Site, y = P_clay, fill = Depth)) +
  geom_boxplot() + theme_light(base_size = b.siz) +
  scale_fill_brewer(palette = p_color) +
  theme(axis.title.x = element_blank()) +
  labs(y= "% Clay") +
  scale_y_continuous(breaks = round(seq(0,30, by = 5),1))
clayp

# Percent moisture
PMp <- ggplot(mydata, aes(x = Site, y = P_moist, fill= Depth)) +
  geom_boxplot() + theme_light(base_size = b.siz) +
  scale_fill_brewer(palette = p_color) + 
  theme(axis.title.x = element_blank()) +
  labs(y="Soil % Moisture")
PMp

phys_figure <- ggarrange(MWDp, FDAp, WAGp, PMp, 
                         ncol = 2, nrow = 2, legend = "top",
                         common.legend = TRUE, label.x = NULL)
phys_figure <- annotate_figure(phys_figure, 
                               top = "Soil Physical Properties by Site and Depth")
phys_figure
# ggsave("phys_figure4.png")


summary(lm(P_moist ~ Site * Depth, data = mydata))
summary(lm(P_OM ~ Site * Depth, data = mydata))
summary(lm(pH ~ Site * Depth, data = mydata))
summary(lm(P_TC ~ Site * Depth, data = mydata))
summary(lm(P_TN ~ Site * Depth, data = mydata))

# *** Chemical soil properties ***

# Electrical Conductivity
ECp <- ggplot(mydata, aes(x = Site, y = EC, fill = Depth))+ 
  geom_boxplot() + theme_light(base_size = b.siz) +
  scale_fill_brewer(palette = p_color) + 
  theme(axis.title.x = element_blank()) +
  labs(y="Electrical Conductivity")
ECp

# t-test to see if the turf surface and subsurface are different
# t.test(x = AS_data$EC[13:15], y = AS_data$EC[16:18])$p.value
# [1] 0.0256

# pH by Site
pHp <- ggplot(mydata, aes(x = Site, y = pH, fill = Depth))+ 
  geom_boxplot() + theme_light(base_size = b.siz) +
  scale_fill_brewer(palette = p_color) + 
  theme(axis.title.x = element_blank()) +
  labs(y="Soil pH")
pHp

# t-test to see if the turf surface and subsurface are different
# t.test(x = AS_data$pH[13:15], y = AS_data$pH[16:18])$p.value
# [1] 0.252

# Total Kjeldahl Nitrogen 
TNp <- ggplot(mydata, aes(x = Site, y = P_TN, fill = Depth)) + 
  geom_boxplot() + theme_light(base_size = b.siz) +
  scale_fill_brewer(palette = p_color) + 
  theme(axis.title.x = element_blank()) +
  labs(y="Total Nitrogen %")
TNp

# Test for difference between surface and subsurface
# t.test(x = AS_data$P_TN[13:15], y = AS_data$P_TN[16:18])$p.value
# [1] 0.0213

# Olsen-P Available Phosphorus mg/L
APp <- ggplot(mydata, aes(x = Site, y = APmgL, fill = Depth)) + 
  geom_boxplot() + theme_light(base_size = b.siz) +
  scale_fill_brewer(palette = p_color) + 
  theme(axis.title.x = element_blank()) +
  labs(y="Available P mg/L")
APp

# % Organic Matter
OMp  <- ggplot(mydata, aes(x = Site, y = P_OM, fill = Depth)) + 
  geom_boxplot() + theme_light(base_size = b.siz) +
  scale_fill_brewer(palette = p_color) + 
  theme(axis.title.x = element_blank()) +
  labs(y="Organic Matter %")
OMp

# Test for difference between surface and subsurface
# t.test(x = AS_data$P_OM[13:15], y = AS_data$P_OM[16:18])$p.value
# [1] 0.0156

# % Organic Carbon
OCp <- ggplot(mydata, aes(x = Site, y = P_OC, fill = Depth)) + 
  geom_boxplot() + theme_light(base_size = b.siz) +
  scale_fill_brewer(palette = p_color) + 
  theme(axis.title.x = element_blank()) +
  labs(y="Organic Carbon %")
OCp

# % Active Carbon
ACp <- ggplot(mydata, aes(x = Site, y = P_AC, fill = Depth)) + 
  geom_boxplot() + theme_light(base_size = b.siz) +
  scale_fill_brewer(palette = p_color) + 
  theme(axis.title.x = element_blank()) +
  labs(y="Active Carbon %")
ACp

# Total % Carbon (Combustion Reactor at 1038C)
TCp <- ggplot(mydata, aes(x = Site, y = P_TC, fill = Depth)) + 
  geom_boxplot() + theme_light(base_size = b.siz) +
  scale_fill_brewer(palette = p_color) + 
  theme(axis.title.x = element_blank()) +
  labs(y="Total % Carbon")
TCp

# Test for difference between surface and subsurface
# t.test(x = AS_data$P_TC[13:15], y = AS_data$P_TC[16:18])$p.value
# [1] 0.03476728

chem_figure <- ggarrange(ECp, pHp, APp, TNp, 
                         OMp, ACp, 
                         ncol = 2, nrow = 3, 
                         common.legend = TRUE, label.x = NULL)
chem_figure <- annotate_figure(chem_figure,
                               top = "Chemical Properties by Site and Depth")
chem_figure
# ggsave("chem_figure6.png")


carbons <- AS_data[,c(5,27,28,30)]

# Load reshape, using the melt() function to move all the variables
# into one column. Probably could do that with some rbind() combo but
# this is an easy cheat
library(reshape)

carbs <- melt(carbons, id.vars = c("Site_depth"), 
              measure.vars = c("P_TC","P_OM","P_OC"))

carbon_figure <- ggplot(carbs) + 
  geom_boxplot(aes(x=Site_depth, y=value, fill=variable))

carbon_figure

# *** Biological soil properties ***

# Actinomycetes Bacteria
Actinop <- ggplot(mydata, aes(x = Site, y = P_Actin, fill=Depth)) +
  geom_boxplot() + theme_light(base_size = b.siz) + 
  scale_fill_brewer(palette = p_color) +
  theme(axis.title.x = element_blank()) +
  labs(y="Actinomycetes %")
Actinop

# Gram Negative Bacteria
Gnegp <- ggplot(mydata, aes(x = Site, y = P_GramNeg, fill=Depth)) +
  geom_boxplot() + theme_light(base_size = b.siz) + 
  scale_fill_brewer(palette = p_color) +
  theme(axis.title.x = element_blank()) +
  labs(y="Gram (-) %")
Gnegp

# Gram Positive Bacteria
Gposp <- ggplot(mydata, aes(x = Site, y = P_AdGramPos, fill=Depth)) +
  geom_boxplot() + theme_light(base_size = b.siz) + 
  scale_fill_brewer(palette = p_color) +
  theme(axis.title.x = element_blank()) +
  labs(y="Other Gram (+) %")
Gposp

# Saprophytic Fungi
SFp <- ggplot(mydata, aes(x = Site, y = P_SapFun, fill=Depth)) +
  geom_boxplot() + theme_light(base_size = b.siz) + 
  scale_fill_brewer(palette = p_color) +
  theme(axis.title.x = element_blank()) +
  labs(y="Saprophyte %")
SFp

# Arbuscular Mycorrhizal Fungi
AMFp <- ggplot(mydata, aes(x = Site, y = P_ArbMyc, fill=Depth)) +
  geom_boxplot() + theme_light(base_size = b.siz) + 
  scale_fill_brewer(palette = p_color) +
  theme(axis.title.x = element_blank()) +
  labs(y="AMF %")
AMFp

# Undifferentiated 
Undiffp <- ggplot(mydata, aes(x = Site, y= P_Undiff, fill=Depth)) +
  geom_boxplot() + theme_light(base_size = b.siz) + 
  scale_fill_brewer(palette = p_color) +
  theme(axis.title.x = element_blank()) +
  labs(y="Undifferentiated %")
Undiffp

# Total PLFA Biomass in ng/g
PLFAp <- ggplot(mydata, aes(x = Site, y = PLFA, fill=Depth)) +
  geom_boxplot() + theme_light(base_size = b.siz) + 
  scale_fill_brewer(palette = p_color) +
  theme(axis.title.x = element_blank()) 
PLFAp + labs(title = "Total PLFA Microbial Biomass")
PLFAp
# ggsave("PLFA_boxplot.png")

# Compile the figure with ggarrange()
species_figure <- ggarrange(Undiffp, Gnegp,
                            Actinop, Gposp,
                            SFp, AMFp,
                            PLFAp,
                            ncol = 2, nrow = 4, 
                            common.legend = TRUE, label.x = NULL)
species_figure <- annotate_figure(species_figure, 
                              top = "Taxonomic Group Composition by Site")
species_figure
# Save the figure (without parameters it will print what is show on screen)
# ggsave("species_figure.png")



###############################################################################

# looking at the data for the ANOVA normality, log transform looks better
# hist(sqrt(AS_data$PLFA), breaks=100)
# 
# library(car)
# Using Type III Sum of Squares (SS) because 
# - Type I SS is a sequential sum of squares, where the order of the factors 
#      influences the significance of the model, giving different results 
#      depending on the sequence of the factors
# - Type III SS tests for the significance of the interaction first, then the
#      main effect, which is better if the interaction between the two factors
#      in a two-way ANOVA is significant, as the sequence of the factors doesn't
#      affect the result of the test. 
# HOWEVER, I'm not running ANOVA because I'm using PERMANOVA instead, which
# is more appropriate for non-parametric data.
# Anova(lm(sqrt(PLFA) ~ Site + Depth_string, data = AS_data), type = "III")
# Anova(lm(PLFA ~ Site + Depth_string, data = AS_data), type = "III")



# Per the suggestions in my manuscript, I'll run a 
# PERMANOVA *~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*
# This is a better choice than the ANOVA probably, because it is a 
# non-parametric test like the NMDS, and uses a distance matrix 
# to look at differences between variables over 999 iterations
library(vegan)

# ?adonis2()

# First find the dissimilarity indices using vegdist()
# The data input object for adonis2() must be a dissimilarity matrix
# defaults:
# vegdist(x, method="bray", binary=FALSE, diag=FALSE, upper=FALSE,
#         na.rm = FALSE, ...) 
names(mydata)
# Dissimiliarity matrix for percent composition
sp.dist <- vegdist(mydata[,6:12])
# Dissimiliarity matrix for total PLFA biomass
PLFA.dist <- vegdist(mydata[,13])

# Use adonis2 to run the PERMANOVA
# defaults: 
# adonis2(formula, data, permutations = 999, method = "bray",
#         sqrt.dist = FALSE, add = FALSE, by = "terms",
#         parallel = getOption("mc.cores"), ...)
# for the formula, using an asterisk * gives the interaction between 
# the factor variables site and depth
adonis2(formula = sp.dist ~ AS_data$Site*AS_data$Depth_string)
adonis2(formula = PLFA.dist ~ AS_data$Site*AS_data$Depth_string)

# No difference really when the order is changed (adonis2 > adonis for this)
# adonis2(formula = sp.dist ~ AS_data$Depth_string*AS_data$Site)
# adonis2(formula = PLFA.dist ~ AS_data$Depth_string*AS_data$Site)

# Based on the dune example
# > data("dune")  <- this is a 20 x 30 dataframe
# > data("dune.env") <- this is a 20 x 5 dataframe
# > adonis(dune ~ Management*A1, data = dune.env)

comp.data <- mydata[,6:12]
env.data <- mydata[,18:32]
adonis(comp.data ~ Site * Depth, data = mydata)


# Composition figure - Stacked geom_plot


# Just the big three- Bac/Fun/Und
names(mydata)
spec.com1 <- mydata[,c(5,14:16)]
str(spec.com1)

comp.means <- by(spec.com1, spec.com1$Site_depth, FUN=function(x){
  means <- colMeans(x[,2:4])
})
# Look at the list of means generated with by()
comp.means

# Transform the list with t() and sapply() to create a data frame with the sites
# on the left, the average of each group variable in the columns
comp.means3 <- t(sapply(comp.means, I))
# Transform again with t() to get the data frame split up by site with the
# variable group names on the left
t(comp.means3)

library(reshape)

spec.com3 <- melt(comp.means3)
names(spec.com3)[1] <- "Site_depth"
names(spec.com3)[2] <- "Group"
names(spec.com3)[3] <- "Percent"

as.numeric(spec.com3$Site_depth)
spec.com3$Site_depth <- c(rep(c(1,2,3,4,5,6),3))
spec.com3$Site_depth <- factor(spec.com3$Site_depth)
levels(spec.com3$Site_depth) <- c("Desert, 0-1cm", "Desert, 1-5cm",
                                  "Turf, 0-1cm", "Turf, 1-5cm",
                                  "Farm, 0-1cm", "Farm, 1-5cm")

spec.com3$Group <- factor(spec.com3$Group)
levels(spec.com3$Group) <- c("Bacteria", "Fungi", "Undifferentiated")

str(spec.com3)

# Phylum composition stacked barplot
library(RColorBrewer)
composition <- ggplot(spec.com3, aes(x = Site_depth, y = Percent)) + 
  geom_col(aes(fill = Group)) + theme_classic() +
  labs(x="Site and Depth", title="Microbial Phylum Composition by Site and Depth") + 
  scale_fill_manual(values=brewer.pal(n=3, name="Dark2"))
composition
# ggsave("composition.png")


# With more detail...
# Adj Gram Positive = Total Gram Pos - Actinomycetes
# Bacteria = Adj Gram Positive + Actino + Gram Negative
# Fungi = Saprophytic Fungi + Arbuscular Mycorrhyzal Fungi

names(mydata)
spec.com2 <- mydata[,c(5:10,12,16)]
str(spec.com2)

comp.means2 <- by(spec.com2, spec.com2$Site_depth, FUN=function(x){
  means <- colMeans(x[,2:8])
})
# Look at the list of means generated with by()
comp.means2

# Transform the list with t() and sapply() to create a data frame with the sites
# on the left, the average of each group variable in the columns
comp.means7 <- t(sapply(comp.means2, I))
# Transform again with t() to get the data frame split up by site with the
# variable group names on the left
t(comp.means7)

spec.com7 <- melt(comp.means7)
names(spec.com7)[1] <- "Site_depth"
names(spec.com7)[2] <- "Group"
names(spec.com7)[3] <- "Percent"
as.numeric(spec.com7$Site_depth)
spec.com7$Site_depth <- c(rep(c(1,2,3,4,5,6),7))
spec.com7$Site_depth <- factor(spec.com7$Site_depth)
levels(spec.com7$Site_depth) <- c("Desert, 0-1cm", "Desert, 1-5cm",
                                  "Turf, 0-1cm", "Turf, 1-5cm",
                                  "Farm, 0-1cm", "Farm, 1-5cm")

spec.com7$Group <- c(rep(1,6),rep(2,6),rep(3,6),rep(4,6),rep(5,6),rep(6,6),
                     rep(7,6))
spec.com7$Group <- factor(spec.com7$Group)
levels(spec.com7$Group) <- c("Actinomycetes", "Other Gram (+)",
                             "Gram (-)", "Saprophytic Fungi",
                             "Arbuscular Mycorrhizal Fungi",
                             "Protozoa", "Undifferentiated")
str(spec.com7)


# Palette for the colorblind, stolen from http://www.cookbook-r.com/Graphs/Colors_(ggplot2)/
cbPalette <- c("#CC79A7", "#E69F00", "#56B4E9", "#009E73", "#F0E442", #"#0072B2", 
               "#D55E00", "#999999")

composition_detail <- ggplot(spec.com7, aes(x = Site_depth, y = Percent)) + 
  geom_col(aes(fill = Group)) + theme_classic() +
  labs(x="Site and Depth", title="Microbial Group Composition by Site and Depth") + 
  scale_fill_manual(values=cbPalette)
composition_detail

# ggsave("composition_detail.png")





#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#


# *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
# Table of environmental properties


# Biological confidence intervals
names(mydata)
biom <- mydata[,c(1,3,6:16)]
# str(envm)
# dim(envm)


bio.means <- aggregate(biom[,sapply(biom, is.numeric)], 
                       with(biom, list(Site, Depth)), 
                       mean)
names(bio.means) <- c("Site", "Depth","Actino", "AdGramPos", "GramNeg", 
                      "SapFun","AMF", "Rhiz", "Proto", "PLFA", "Bac", 
                      "Fun", "Undiff")

ord<-order(bio.means$Site)
bio.means <- bio.means[ord,]
row.names(bio.means) <- NULL

bio.stdv <- aggregate(biom[,sapply(biom, is.numeric)], 
                      with(biom, list(Site, Depth)), 
                      sd)
names(bio.stdv) <- c("Site", "Depth","Actino", "AdGramPos", "GramNeg", 
                     "SapFun","AMF", "Rhiz", "Proto", "PLFA", "Bac", 
                     "Fun", "Undiff")
bio.stdv <- bio.stdv[ord,]

t <- qt(p=0.95,df=2)
upper.bio <- bio.means[3:13] + (t * bio.stdv[3:13]/sqrt(3))
lower.bio <- bio.means[3:13] - (t * bio.stdv[3:13]/sqrt(3))

bio.CI <- cbind(bio.means, bio.stdv, lower.bio, upper.bio)
write.csv(bio.CI, file = "bio_CI.csv", row.names = FALSE)





# need to transform the matrix to get the confidence intervals to look right
# with the table like:
#                   EC             pH
#              lower  upper   lower  upper
# Desert 0-1cm  2.1   1.2       7.2   7.3
#        1-5cm
# Turf path 0-1
#           1-5
# May just make in excel to save some grief (maybe..?)
# Still don't know how to export tables in R.
# install.packages("kableExtra")

library(kableExtra)
kbl(env.means)

# Change the default number of digits to 3 
options(digits = 3)

# Basic example of the table
env.means %>%
  kbl() %>%
  kable_styling(bootstrap_options = "striped", full_width = F, position = "left")

# env.means %>%
#   kbl() %>%
#   kable_paper("hover", full_width = FALSE)
# 
# env.means %>%
#   kbl(caption = "Table 1. Environmental Properties by Site and Depth") %>%
#   kable_classic(full_width = F, html_font = "Calibri") %>%
#   kable_styling(position = "center", font_size = 12)

soil.type <- c("loamy fine sand", " ",
               "sandy loam", " ",
               "clay loam", " ")

phys.means <- cbind(env.means[,c(1:8)], soil.type)
names(phys.means)[9] <- "Soil Class"

# Pretty looking table
phys.means %>%
  kbl(caption = "Table 1. Soil Physical Properties by Site and Depth") %>%
  kable_classic(full_width = F, html_font = "Calibri") %>%
  kable_styling(position = "center", font_size = 12)

# InorgC <- AS_data$P_TC - AS_data$P_OC
# plot(InorgC)

chem.means <- env.means[,c(1:2, 9:11, 16:17, 12,14,13,15)]
chem.means %>%
  kbl(caption = "Table 2. Soil Chemical Properties by Site and Depth") %>%
  kable_classic(full_width = F, html_font = "Calibri") %>%
  kable_styling(position = "center", font_size = 12)

# Could add richness (rich) and Shannon H (H) to the pile too..?
bio.means <- mydata[,c(1,4,13,14,8,6,7,15,9,10,16)]
bio.means <- aggregate(bio.means[,sapply(bio.means, is.numeric)], 
                        with(bio.means, list(Site, Depth_string)), 
                        mean)
names(bio.means)[1] <- "Site"
names(bio.means)[2] <- "Depth"

ord <- order(bio.means$Depth, decreasing = TRUE)
ord2 <- order(bio.means$Site)
bio.means <- bio.means[ord,]
bio.means <- bio.means[ord2,]
rownames(bio.means) <- NULL

nematodes <- c(41,33,5,5,8,12)  # taken from the SOIL 476L data
cyanos <- c(15,0,2,2,10,0)      # taken from the SOIL 476L data
bio.means <- cbind(bio.means, cyanos, nematodes)

names(bio.means) <- c("Site", "Depth", "Total Biomass (ng/g)", "% Bac",
                      "% Gram (-)", "% Actino", "% Gram (+)",
                      "% Fungi", "% Saprophytes", "% AMF",
                      "% Undiff", "Cyano Count", "Nematode Count")
bio.means %>%
  kbl(caption = "Table 3. Soil Biological Properties by Site and Depth") %>%
  kable_classic(full_width = F, html_font = "Calibri") %>%
  kable_styling(position = "center", font_size = 12)

####  Multivariate Analyses  ##################################################
# install.packages("vegan")
# load the vegan package
library(vegan)

# Use the data set that has the most diversity considered but no overlapping
# values (adGramPos = GramPos - Actino)

P_AdGramPos <- AS_data$P_GramPos - AS_data$P_Actin
species <- cbind(AS_data[,c(6,7,9,10,11,12,16)], P_AdGramPos)
head(species)

names(species)
H <- diversity(species)
rich <- specnumber(species)
anova(lm(H ~ AS_data$Site_depth))

# Subset the environmental data by selecting the following variables:
# Column 13 - Total PLFA biomass in ng/g soil
# 18 - MWD, 19 - Percent_fine-diameter-aggregates (<0.25 mm)
# 20 - P_water-stable-2mm-aggregates, 21 - P_stones-in-2mm-aggregates
# 22 - P_sand, 23 - P_clay, 24 - P_moisture, 25 - EC, 26 - pH
# 27 - P_organic-matter, 28 - P_organic-carbon, 29 - P_active-carbon
# 30 - P_total-carbon, 31 - P_total-nitrogen, 
# 32 - total available phosphorus mg/L, 33 - P_available-phosphorus

# Finding the C:N ratio
AS_data$CNrat <- AS_data$P_TC/AS_data$P_TN

names(AS_data)

# Removed OC and C:N because they are directly on top of OM
env <- AS_data[,c(18:27, 29:32)]
# Looking at the variables, deciding what to include
# env <- AS_data[,c(19:27,29,31,32)]
# env <- AS_data[,c(33,26,31,23,20,25,30)]
# env <- AS_data[,c(20,25,26,27,22)]
head(env)

# First, run the unconstrained PCA. 
# env.un.pca <- rda(env)
# env.un.pca

# Plot the unconstrained PCA to look for patterns.
# TKN pops out as a strong vector pulling to the right, strongly associated
# with the x-axis. Sand strongly pulls just off the y-axis to the left. 
# *HOWEVER* - because this is unconstrained (scale = FALSE) the reason the 
# TKN and sand variables have strong vectors is magnitude (these values are
# larger than the rest of the variables).
# biplot(env.un.pca, scaling="symmetric")

# Run the constrained PCA.
# The scaling = TRUE forces the variance for the variables to be equal, or
# rather it uses "unit variance". 
# env.c.pca <- rda(env, scale = TRUE)
# env.c.pca

env.c.pca <- rda(env, scale = TRUE)
summary(env.c.pca)
# The variables that influence PC1
sort(summary(env.c.pca)$species[,1], decreasing = TRUE)
# and PC2
sort(summary(env.c.pca)$species[,2], decreasing = TRUE)


# The summary gives the proportion explained by the given eigenvalue.
# We see the same pattern as the unconstrained PCA with direction of the
# directional vectors but the environmental factors near the farm points
# are all about the same length with standardization. 
summary(eigenvals(env.c.pca))
summary(env.c.pca)
plot(env.c.pca)


# Plot the constrained PCA, make it pretty with points instead of sit#
# Here notice we see the points separating by location, with the desert
# and turf sites on the left and the farm on the right. 
# The turf plot is in the (-x,-y) corner of the ordination, with the main
# drivers of the -x,-y being sand and the x,-y being pH and EC (more 
# correlated with the -y axis). 
# biplot(env.c.pca, scaling="symmetric")
# biplot(env.c.pca, scaling="symmetric", type = c("text", "points"))
biplot(env.c.pca, scaling="symmetric", type = c("text", "points"))


# Using the example from Gavin Simpson's webinar video
# https://youtu.be/tVnnG7mFeqA

# His code for species data:
# pca <- rda(decostand(varespec, method = "hellinger"), scale = TRUE)

# We adapt the standardization decostand() code for our data but
# plotting the env.c.pca2 gives the same plot. This is because the 
# decostand(method = "standardize") does the same thing as scale = TRUE.
# From the vegan help: "standardize: scale x to zero mean and unit variance
# (default MARGIN = 2)"

# env.c.pca2 <- rda(decostand(env, method = "standardize"), scale = TRUE)
# biplot(env.c.pca2, scaling="symmetric", type = c("text", "points"))


# Plotting the location onto the constrained PCA
# Subset metadata:
# Column 1 - Site, 2 - Transect (omitted), 3 - Depth (numeric, omitted),
# 4 - Depth stirng (surface, subsurface)
meta <- AS_data[ , c(1,4)]
head(meta)
str(meta)
dim(meta)

# Set the meta variables as factors
meta$Site<-as.factor(meta$Site)
meta$Depth_string<-as.factor(meta$Depth_string)

## SET GRAPHIC PARAMETERS 
# Adapted from Gavin's code
# Set the graph parameters 
point <- rep(c(16, 17),3)   # for pch =      (point shape)
scl <- "symmetric"          # for scaling =  
sites <- "sites"            # for display = 
size <- 1                   # for cex =      (point size)


# Sites, the color of the points on the PCA
# Desert = burlywood, Farm = tan4, Turf = dark olive green
col_sites <- c("blue3", "sienna2", "green3")
cols <- with(meta, col_sites[Site])
lvl <- with(meta, levels(Site)) 
d.point <- with(meta, point[Depth_string])

# png(filename="sp.all.c.pca.png", width = 600, height = 600, pointsize = 14)


## FROM HERE DOWN PLOTS THE ENVIRONMENTAL VARIABLES PCA #######################
# Set up the ordinal plot for PCA using scaling = "symmetric" and 
# display = "sites"
# Running this line of code will clear the graph!
# plot(env.c.pca, type = "n", scaling = scl, display = sites)

plot(env.c.pca, type = "n", scaling = scl, display = sites)

# Plot the PCA
# biplot() adds to the existing plot. Add the PCA to the plot generated above.
# biplot(env.c.pca, scaling="symmetric", type = c("text", "points"), col="black")

biplot(env.c.pca, scaling="symmetric", type=c("text", "points"), col="black")

# Alternative plot for PCA, with axis range defined and black arrows
# biplot(env.c.pca, col="black",
#        xlim = c(-2.5, 2.5), ylim = c(-2, 1.5),
#        cex.lab = 1, cex.axis = 1,
#        type = c("text", "points"))
# The code stopped working when I tried to add a title, giving the error message:
# Error in x[, 2] : subscript out of bounds
# The biplot() code above worked fine until I added the title... 
# It seems to return a non annoying error like this.
biplot(title(main = "PCA of Environmental Parameters by Site"))

# Add the site points generated by the PCA
# points(env.c.pca, display = sites, scaling = scl, pch = d.point, col = cols,
#        cex = size)
# Here I've assigned the pch (point character) to be d.point, or the depth point,
# showing the difference between depths with two different shapes.
# The cols is going to 
scores(env.c.pca)
points(env.c.pca, display = sites, scaling = scl, pch = d.point, col = cols,
       cex = size)

# Put the legend for the color by sites on the graph
legend("bottomleft", legend = lvl, bty = "n", col = col_sites, pch = 19)
legend("bottomright", legend=c("1-5 cm Depth", "0-1 cm Depth"), bty="n", pch = point)
# legend("bottomleft", legend = lvl, bty = "n", col = col_sites, pch = 19)
# legend("bottomright", legend=c("1-5 cm Depth", "0-1 cm Depth"), bty="n", pch = point)

# This finalizes the save function started above.
# dev.off()



### NMDS  *********************************************************************
# Subset the species
head(AS_data)
# This is all microbial compositional data (% by group)
# 6 - Actinomycetes, 7 - Gram Neg, 8 - Gram pos, 
# 9 - Saprophytes, 10 - AMF, 11 - Rhizobia, 12 - Protozoa
# 13 - Total PLFA biomass
# 14 - % Bac, 15 - % Fun, 16 - % Undiff
# 17 - Fungi to Bacteria ratio

species <- cbind(AS_data[c(6)], P_AdGramPos, AS_data[c(7,9,10,11,12,16)])
head(species)
species

meta <- AS_data[ , c(1,4)]
head(meta)
str(meta)
dim(meta)
# This is for the NMDS ordinal plot to show the Turf isn't actually turf
meta$Site <- c("Desert", "Desert","Desert","Desert","Desert","Desert",
               "Farm","Farm","Farm","Farm","Farm","Farm",
               "Turf Path","Turf Path","Turf Path","Turf Path","Turf Path","Turf Path")

# Set the meta variables as factors
meta$Site<-as.factor(meta$Site)
meta$Depth_string<-as.factor(meta$Depth_string)

levels(meta$Site)

# Compute Nonmetric Multidimensional Scaling (NMDS)
# k is the number of dimensions, k=2 is default
NMDS.sp <- metaMDS(species, k=2)
# View the stressplot
stressplot(NMDS.sp)
names(NMDS.sp)

# Plot the farm OM and PLFA to look for correlations
lm1 <- lm(AS_data[c(7:12),27] ~ AS_data[c(7:12),13], data = AS_data)
plot(AS_data[c(7:12),13], AS_data[c(7:12),27]) + abline(lm1)
summary(lm1)
names(AS_data)

# Draw NMDS - with ellipses and hulls (BELOW IS THE PLOT)
# dev.off()
# png("NMDS_hulls-Sites.png", width = 600, height = 600, pointsize = 18)
# Plotting empty ordination space to choose what is to be displayed
plot(NMDS.sp, disp="sites", type="n")
# Draw the convex "hulls" connecting the vertices of the points associated 
# by each community
ordihull(NMDS.sp, meta$Site, col=col_sites, lwd=2)
# Draw an "ellipse" and using function "ehull" for ellipsoid hull to enclose all 
# points in the group.
ordiellipse(NMDS.sp, meta$Site, col=col_sites, kind = "ehull", lwd=2)
# Draw an ellipse based off the standard deviation of point scores
ordiellipse(NMDS.sp, meta$Site, col=col_sites, draw="polygon", kind = c("sd"))
ordispider(NMDS.sp, meta$Site, col=col_sites, label = TRUE)
biplot(title(main = "NMDS Ordinal Separation of \n Microbial Groups by Site"))
# points(NMDS.sp, disp="sites", pch=21, col="black", bg="white", cex=1)
# dev.off()


# By depth
col_depth <- c("black","red")
# png("NMDS_hulls-Depth.png", width = 600, height = 600, pointsize = 18)
plot(NMDS.sp, disp="sites", type="n")
ordihull(NMDS.sp, meta$Depth_string, col=col_depth, lwd=2)
ordiellipse(NMDS.sp, meta$Depth_string, col=col_depth, kind = "ehull", lwd=2)
ordiellipse(NMDS.sp, meta$Depth_string, col=col_depth, draw="polygon", kind = c("sd"))
ordispider(NMDS.sp, meta$Depth_string, col=col_depth, label = TRUE)
biplot(title(main = "NMDS Ordinal Separation of \n Microbial Groups by Depth"))
# dev.off()


# By Site and by Depth
meta2 <- AS_data[,c(5)]
col_sites2 <- c("#0559be","#11a1e6","#663333","#f0c683","#36831e","#9efe03")

# png("NMDS_hulls-Site_depth.png", width = 600, height = 600, pointsize = 18)
plot(NMDS.sp, disp="sites", type="n")
ordiellipse(NMDS.sp, meta2, col=col_sites2, kind = "ehull", lwd=2)
ordiellipse(NMDS.sp, meta2, col=col_sites2, draw="polygon", kind = c("sd"))
ordispider(NMDS.sp, meta2, col=col_sites2, label = TRUE)
biplot(title(main = "NMDS Ordinal Separation of \n Microbial Group by Site and Depth"))
# dev.off()

#***** Plot from the final script of SOIL 598 for the vector NMDS **************
# Run the environmental fit test; show relationships of env variables to species 
# distribution fits vectors of continuous variables (nutrient content, texture, 
# etc.) and centroids of levels of class variables (textural class, 
# surface/subsurf, etc.) arrow shows the direction of the increasing gradient
# length of the arrow is directly proportional to the correlation between 
# variable and ordination

# Subset the environmental data by selecting the following variables:
# Column 13 - Total PLFA biomass in ng/g soil
# 18 - MWD, 19 - Percent_fine-diameter-aggregates (<0.25 mm)
# 20 - P_water-stable-2mm-aggregates, 21 - P_stones-in-2mm-aggregates
# 22 - P_sand, 23 - P_clay, 24 - P_moisture, 25 - EC, 26 - pH
# 27 - P_organic-matter, 28 - P_organic-carbon, 29 - P_active-carbon
# 30 - P_total-carbon, 31 - P_total-nitrogen, 
# 32 - total available phosphorus mg/L, 33 - P_available-phosphorus

# What we used in the PCA
env <- AS_data[,c(18:27, 29:32)]

# Remove P_clay, P_OM
env <- AS_data[,c(18:22, 24:26, 29:32)]
head(env, 2)      

# NMDS for ENV ~ SPECIES
fit.env <- envfit(ord=NMDS.sp,env=env, perm = 999, na.rm = TRUE) 

scores(fit.env, "vectors")
summary(fit.env)
# Variables influencing PC1
sort(scores(fit.env, "vectors")[,1], decreasing = TRUE)
# Variables influencing PC2
sort(scores(fit.env, "vectors")[,2], decreasing = TRUE)

# Start here for a blank plot *********************************
ordiplot(NMDS.sp,type="n")

# Code to save the plot
# png("NMDS.fit.env.all.png", height = 600, width = 600, pointsize = 18)

orditorp(NMDS.sp,display="species",
         # Trying to remove visible values from the plot but not from the 
         # analysis, like Gavin had mentioned in the webinar. 
         # removed Protozoa, skewed off to the right
         labels = c("Actino", "Gram(-)", "Sapro", "AMF", "Rhizobia","", 
                    "Undiff", "OtherGram(+)"),
         col="#03165a", air=0.01, pch = NA, cex = 0.8)
plot(fit.env, cex = 0.6, col = "black")
plot(fit.env, p.max = 0.05, col = "blue", cex = 0.6)
biplot(title(main = "NMDS of Env and Soil Health Parameters \n by Community Group"))

# dev.off()


# *****************************************************************************


# Attempting to make a table for the biomarkers
# This is way too long and I can't figure out how to break over several columns
# Work on this later....
biomarkers <- read.table(file="biomarkers.csv", sep=",", header = TRUE)
names(biomarkers)[1] <- "Group ID"
names(biomarkers)[2] <- "Biomarker(s)"
names(biomarkers)[3] <- "Taxonomic Group Name"

library(kableExtra)
kbl(biomarkers[,-1])
biomarkers[,-1] %>%
  kbl(caption = "Table X. PLFA Biomarker Taxonomic Group Assignments") %>%
  kable_classic(full_width = F, html_font = "Helvetica") %>%
  kable_styling(position = "center", font_size = 10)
