#To run script automatrically:
  #Must call with Source instead of Run
  #This can be done either by clicking "Source" in the top right of the source panel or "Command+Shift+Enter"
  #After that, can revisit any specific lines of code and "Run" them as usual

#Alternatively, can setwd manually to location of "Data" folder and "Run" in chunks

#All data files must be downladed along with the script for it to run 
  #They must be in the "Data" folder parallel to (in the same location as) this script
this.dir <- dirname(parent.frame(2)$ofile) 
setwd(this.dir) #Set working directory

###########################
#Basic R-prep 
rm(list=ls())  #clear variables
set.seed(42) #Set the starting point for randomized values

library(ggplot2)
library(plyr)
library(dplyr)
library(tidyr)
library(reshape2)
library(gridExtra)
library(grid)
library(agricolae)  
library(foreign)
library(Rmisc)
library(drc)
library(gtools)
library(msm)
library(deSolve)


# The colorblind palette with grey to use for any figures:
cbPalette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7","#999999")
#Reset the ggplot default theme to black and white with text size 25 (user preference)
theme_set(theme_bw(base_size = 25))


#####################################################################################
# Field Diet Survey  ------------------------------------------
#####################################################################################

###########################
#Data
uniqDIET<-read.csv("Data/DataS1.csv")



#Running the Binomial GLM on whether consumption itself was controlled
summary(zeroProb <- glm(Empty~Species*Open*Date, data = uniqDIET, family = binomial(link = "logit")))
#Percentage of diets which were empty
100*sum(uniqDIET$Empty==1)/nrow(uniqDIET) 

#For those exceeding the consumption threshold, what factors controlled their magnitude of consumption
summary(nonZero <- glm(MassTotal~Species*Open*Date, data = subset(uniqDIET, Empty==0), family = Gamma(link="log")))


###########################
uniqEggless<-read.csv("Data/DataS2.csv")


#Empty diets analyzed with a Binomial GLM
summary(zeroProbEggless <- glm(Empty~Species*Open*Date, data = uniqEggless, family = binomial(link = logit)))
#Percentage of diets which were empty
100*sum(uniqEggless$Empty==1)/nrow(uniqEggless) 

#The magnitude of non-egg diets are analyzed with a Gamma GLM
summary(nonZeroEggless <- glm(MassTotalEggless~Species*Open*Date, data = subset(uniqEggless, Empty==0), family = Gamma(link="log")))


###########################
#Plot

#Total Diet first
#First insert Tukey groups into a data frame for printing onto the figure
test<-with(uniqDIET,interaction(Open,Date,Species))
testAOV<-aov(MassTotal~test, data=uniqDIET)
results<-HSD.test(testAOV,"test",group=TRUE)
out<-results$group

#And calculate SEM for each group to insert error bars
stats <- summarySE(uniqDIET, measurevar="MassTotal", groupvars=c("Open","Date","Species"),na.rm=TRUE)
stats$all<-paste(stats$Open,stats$Date)

#Combine mean, SEM, and Tukey groups into one data frame that can be used for plotting
mergeDIET<-merge(stats,out,by.x="MassTotal")

#ggplot graphic creation for whole diets 
#(Figure 2a)
ggplot(mergeDIET, aes(x=all, y=MassTotal,fill=Species)) + 
  geom_col(position=position_dodge(.7),colour="black", width = 0.7) + 
  geom_errorbar(aes(ymin=MassTotal, ymax=MassTotal+se),size=.5, width=.2,position=position_dodge(.7)) + 
  scale_x_discrete(expand=c(0.1,0.1),labels=c("Before","During","Before","During"))+
  scale_fill_manual(values=cbPalette)+
  theme(legend.justification=c(1,1),legend.position=c(0.53,0.888),
        legend.box.background=element_rect(colour = "black",size=0.75),axis.text=element_text(size=25),
        legend.text=element_text(size=20),legend.title=element_text(size=25),legend.key.size=unit(1.5,"cm"))+
  geom_text(aes(x=c(1,1,2,2,3,3,4,4), y=c(1.2,1.2,0.9,0.9,1.3,1.3,5.9,5.9),label=c("a","a","a","a","a","a","b","b")), size=8)+ #Manually adding Tukey groups at variable heights
  geom_text(aes(x=all, y=MassTotal+se+0.4,label=N), position=position_dodge(width=0.7), size=8)+
  ylab("Diet Mass (g)")+ xlab(" No Salmon                                                Salmon    ") #Imperfect means of axis labeling the salmon factor




#Trying with a more b&w friendly palette
ggplot(mergeDIET, aes(x=all, y=MassTotal,fill=Species)) + 
  geom_col(position=position_dodge(.7),colour="black", width = 0.7) + 
  geom_errorbar(aes(ymin=MassTotal, ymax=MassTotal+se),size=.5, width=.2,position=position_dodge(.7)) + 
  scale_x_discrete(expand=c(0.1,0.1),labels=c("Before","During","Before","During"))+
  scale_fill_viridis_d(begin=0.3,end=1)+
  scale_y_continuous(breaks=c(0,2,4,6),labels=c("0","2.0","4.0","6.0"))+
  theme(legend.justification=c(1,1),legend.position=c(0.53,0.888),
        legend.box.background=element_rect(colour = "black",size=0.75),axis.text=element_text(size=25),
        legend.text=element_text(size=20),legend.title=element_text(size=25),legend.key.size=unit(1.5,"cm"))+
  geom_text(aes(x=c(1,1,2,2,3,3,4,4), y=c(1.2,1.2,0.9,0.9,1.3,1.3,5.9,5.9),label=c("a","a","a","a","a","a","b","b")), size=8)+ #Manually adding Tukey groups at variable heights
  geom_text(aes(x=all, y=MassTotal+se+0.4,label=N), position=position_dodge(width=0.7), size=8)+
  geom_label(aes(x=0.75,y=5.57,label="a"),size=15,fill="white",vjust=0.25,hjust=0.5)+
  ylab("Diet Mass (g)")+ xlab(" No Salmon                                                Salmon    ") #Imperfect means of axis labeling the salmon factor




###########################
#Plot

#Creating the same figure for the non-egg diets
#First generating Tukey groupings
eggtest<-with(uniqEggless,interaction(Open,Date,Species))
eggtestAOV<-aov(MassTotalEggless~eggtest, data=uniqEggless)
eggresults<-HSD.test(eggtestAOV,"eggtest",group=TRUE)
eggout<-eggresults$group

#Calculating mean and SEM
eggstats <- summarySE(uniqEggless, measurevar="MassTotalEggless", groupvars=c("Open","Date","Species"),na.rm=TRUE)
eggstats$all<-paste(eggstats$Open,eggstats$Date)

#Combining all together
mergeEggless<-merge(eggout,eggstats,by.x="MassTotalEggless")

#Printing ggplot figure for the non-egg diets 
#(Figure 2b)
ggplot(mergeEggless, aes(x=all, y=MassTotalEggless,fill=Species)) + 
  geom_col(position=position_dodge(.7),colour="black", width = 0.7) + 
  geom_errorbar(aes(ymin=MassTotalEggless, ymax=MassTotalEggless+se),size=.5, width=.2,position=position_dodge(.7)) + 
  scale_x_discrete(expand=c(0.1,0.1),labels=c("Before","During","Before","During"))+
  scale_fill_manual(values=cbPalette)+
  theme(legend.justification=c(1,1),legend.position=c(0.53,0.888),
        legend.box.background = element_rect(colour = "black",size=0.9),axis.text=element_text(size=25),
        legend.text=element_text(size=20),legend.title=element_text(size=25),legend.key.size=unit(1.5,"cm"))+
  geom_text(aes(x=all, y=MassTotalEggless+se+0.07,label=N), position=position_dodge(width=0.7), size=8)+
  geom_text(aes(x=all, y=MassTotalEggless+se+0.14,label=groups), position=position_dodge(width=0.7), size=8)+
  ylab("Diet Mass (g)")+ xlab(" No Salmon                                                Salmon    ") #Imperfect means of axis labeling the salmon factor


#(Figure 2b) B&W friendly
ggplot(mergeEggless, aes(x=all, y=MassTotalEggless,fill=Species)) + 
  geom_col(position=position_dodge(.7),colour="black", width = 0.7) + 
  geom_errorbar(aes(ymin=MassTotalEggless, ymax=MassTotalEggless+se),size=.5, width=.2,position=position_dodge(.7)) + 
  scale_x_discrete(expand=c(0.1,0.1),labels=c("Before","During","Before","During"))+
  scale_fill_viridis_d(begin=0.3)+
  theme(legend.justification=c(1,1),legend.position=c(0.53,0.888),
        legend.box.background = element_rect(colour = "black",size=0.9),axis.text=element_text(size=25),
        legend.text=element_text(size=20),legend.title=element_text(size=25),legend.key.size=unit(1.5,"cm"))+
  geom_text(aes(x=all, y=MassTotalEggless+se+0.07,label=N), position=position_dodge(width=0.7), size=8)+
  geom_text(aes(x=all, y=MassTotalEggless+se+0.14,label=groups), position=position_dodge(width=0.7), size=8)+
  geom_label(aes(x=0.75,y=1.08,label="b"),size=15,fill="white",vjust=0.5,hjust=0.5)+
  ylab("Diet Mass (g)")+ xlab(" No Salmon                                                Salmon    ") #Imperfect means of axis labeling the salmon factor




###########################
#Proportion of stomach contents in each of egg/non-egg categories

Proportions<-read.csv("Data/DataS3.csv")

#Calculating the mean of each diet item type in each species, date, and location cross
PropMeans<-Proportions%>%
  subset(MassTotal>0)%>%
  group_by(OpenDate,Species,Diet_Item)%>%
  summarise(MEAN=mean(Prop)*100,
            SD=sd(Prop)*100,
            SEM=SD/length(Prop)*100)

#Producing a ggplot figure with the proportions of diets due to egg and non-egg items 
#(Figure A3)
ggplot(PropMeans,aes(x=OpenDate,y=MEAN,fill=Diet_Item))+
  geom_col(color="black")+
  scale_fill_manual(values=cbPalette,name="Diet Item")+
  facet_grid(rows=vars(Species))+
  theme(legend.justification=c(0.5,1),legend.position="top",
        legend.box.background=element_rect(colour = "black",size=1.5),axis.text=element_text(size=25))+
  scale_x_discrete(expand=c(0.13,0.13),labels=c("Before","During","Before","During"))+
  ylab("Diet Composition (%)")+ xlab(" No Salmon                                      Salmon    ")


#(Figure A3) B&W friendly
ggplot(PropMeans,aes(x=OpenDate,y=MEAN,fill=Diet_Item))+
  geom_col(color="black")+
  scale_fill_viridis_d(begin=0.3,name="Diet Item")+
  facet_grid(rows=vars(Species))+
  theme(legend.justification=c(0.5,1),legend.position="top",
        legend.box.background=element_rect(colour = "black",size=1.5),axis.text=element_text(size=25))+
  scale_x_discrete(expand=c(0.13,0.13),labels=c("Before","During","Before","During"))+
  ylab("Percent Diet Mass")+ xlab(" No Salmon                                      Salmon    ")




#




#




#




#




#




#

#####################################################################################
# Feeding Experiment ------------------------------------------------------
#####################################################################################

###########################
#Growth Response

#Feeding and Measurements throughout experiment
Growth<-read.csv("Data/DataS4.csv")


#Growth Data during the experiment
Growth<-Growth%>%
  mutate(rateMass=dMass/Days_Alive)

summary(ANOVAgrowth<-aov(rateMass~Species*Treatment_Names+Size_Class, data=Growth))

#Checking ANOVA assumptions
  #listed as an annotation to allow the source call of the script to continue uninterrupted
#residuals<-residuals(object=ANOVAgrowth)
#shapiro.test(residuals)
#plot(ANOVAgrowth)


###########################
#Plot

#Creating a graph for the growth rate in the experiment
#Needs to same Tukey groups as for the previous plots
growthtest<-with(Growth,interaction(Treatment_Names,Species))
growthAOV<-aov(rateMass~growthtest, data=Growth)
growthresults<-HSD.test(growthAOV,"growthtest",group=TRUE)
growthout<-growthresults$group

#And mean, SEM for plotting
growthstats <- summarySE(Growth, measurevar="rateMass", groupvars=c("Treatment_Names","Species"),na.rm=TRUE)
growthstats$Treatment_Names<-factor(growthstats$Treatment_Names,levels = c("Chironomids","Steady","Intermediate","Gorge"))
growthstats$trt<-paste(growthstats$Treatment_Names,growthstats$Species)

#Putting them together
growthmerge<-merge(growthstats,growthout,by.x="rateMass")

#Plotting the growth rate for both species in the experiment 
#(Figure 3a)
ggplot(growthmerge, aes(x=Species, y=rateMass,fill=Treatment_Names)) + 
  geom_col(position=position_dodge(0.7),colour="black", width = 0.7) + 
  geom_errorbar(aes(ymin=rateMass, ymax=rateMass+se),size=0.5, width=0.2,position=position_dodge(0.7)) + 
  scale_x_discrete(expand=c(0.2,0.2))+
  scale_fill_manual(values=cbPalette,name="Treatment")+
  theme(legend.justification=c(1,1),legend.position=c(0.98,0.97),
        axis.text=element_text(size=25),legend.box.background = element_rect(colour = "black",size=0.9),
        legend.text=element_text(size=20),legend.title=element_text(size=25,hjust=0.5),legend.key.size=unit(1.5,"cm"))+
  geom_text(aes(x=Species, y=rateMass+se+0.05,label=N), position=position_dodge(width=0.7), size=8)+
  geom_text(aes(x=Species, y=rateMass+se+0.1,label=groups), position=position_dodge(width=0.7), size=8)+
  geom_label(aes(x=0.67,y=0.54,label="a"),size=15,fill="white",hjust=0.5,vjust=0.25)+
  ylab("Growth Rate (g/day)")+ xlab("Species")+
  guides(fill=guide_legend(nrow=2))


#Trying with a more B&W friendly palette
ggplot(growthmerge, aes(x=Species, y=rateMass,fill=Treatment_Names)) + 
  geom_col(position=position_dodge(0.7),colour="black", width = 0.7) + 
  geom_errorbar(aes(ymin=rateMass, ymax=rateMass+se),size=0.5, width=0.2,position=position_dodge(0.7)) + 
  scale_x_discrete(expand=c(0.2,0.2))+
  scale_fill_viridis_d(begin=0.3,end=1,name="Treatment")+
  theme(legend.justification=c(1,1),legend.position=c(0.98,0.97),
        axis.text=element_text(size=25),legend.box.background = element_rect(colour = "black",size=0.9),
        legend.text=element_text(size=20),legend.title=element_text(size=25,hjust=0.5),legend.key.size=unit(1.5,"cm"))+
  geom_text(aes(x=Species, y=rateMass+se+0.05,label=N), position=position_dodge(width=0.7), size=8)+
  geom_text(aes(x=Species, y=rateMass+se+0.1,label=groups), position=position_dodge(width=0.7), size=8)+
  geom_label(aes(x=0.67,y=0.54,label="a"),size=15,fill="white",hjust=0.5,vjust=0.25)+
  ylab("Growth Rate (g/day)")+ xlab("Species")+
  guides(fill=guide_legend(nrow=2))


###########################
#Conversion Efficiency

Efficiency<-read.csv("Data/DataS5.csv")

#Conversion Efficiency of consumed energy through Experiment
summary(ANOVAEfficiency<-aov(CE~Species*Treatment_Names+Size_Class, data = Efficiency))

#ANOVA assumptions check
  #listed as an annotation to allow the source call of the script to continue uninterrupted
#residuals<-residuals(object=ANOVAEfficiency)
#shapiro.test(residuals)
#plot(ANOVAEfficiency)

###########################
#Plot

#Preparing the Tukey groups
efficiencytest<-with(Efficiency,interaction(Treatment_Names,Species))
efficiencyAOV<-aov(CE~efficiencytest, data=Efficiency)
efficiencyresults<-HSD.test(efficiencyAOV,"efficiencytest",group=TRUE)
efficiencyout<-efficiencyresults$group

#Calculating the mean and SEM for plotting
efficiencystats <- summarySE(Efficiency, measurevar="CE", groupvars=c("Treatment_Names","Species"),na.rm=TRUE)
efficiencystats$Treatment_Names<-factor(efficiencystats$Treatment_Names,levels = c("Chironomids","Steady","Intermediate","Gorge"))
efficiencystats$trt<-paste(efficiencystats$Treatment_Names,efficiencystats$Species)

#Putting them together
efficiencymerge<-merge(efficiencystats,efficiencyout,by.x="CE")
efficiencymerge$trt<-gsub(" B",".B",efficiencymerge$trt)


#Plotting conversion efficiency of trout in experiment 
#(Figure 3b)
ggplot(efficiencymerge, aes(x=Species, y=CE,fill=Treatment_Names)) + 
  geom_col(position=position_dodge(0.7),colour="black", width = 0.7) + ylim(0,0.6)+
  geom_errorbar(aes(ymin=CE, ymax=CE+se),size=.5, width=.2,position=position_dodge(0.7))+
  scale_fill_manual(values=cbPalette,name="Treatment")+scale_x_discrete(expand=c(0.2,0.2))+
  theme(legend.justification=c(1,1),legend.position=c(0.98,0.97),
        axis.text=element_text(size=25),legend.box.background = element_rect(colour = "black",size=0.9),
        legend.text=element_text(size=20),legend.title=element_text(size=25,hjust=0.5),legend.key.size=unit(1.5,"cm"))+
  geom_text(aes(x=Species, y=CE+se+0.033,label=N), position=position_dodge(width=0.7), size=8)+
  geom_text(aes(x=Species, y=CE+se+0.07,label=groups), position=position_dodge(width=0.7), size=8)+
  ylab("Conversion Efficiency (g/g)")+ xlab("Species")+
  guides(fill=guide_legend(nrow=2))


#(Figure 3b) more B&W friendly
ggplot(efficiencymerge, aes(x=Species, y=CE,fill=Treatment_Names)) + 
  geom_col(position=position_dodge(0.7),colour="black", width = 0.7) + ylim(0,0.6)+
  geom_errorbar(aes(ymin=CE, ymax=CE+se),size=.5, width=.2,position=position_dodge(0.7))+
  scale_fill_viridis_d(begin=0.3,name="Treatment")+
  scale_x_discrete(expand=c(0.2,0.2))+
  theme(legend.justification=c(1,1),legend.position=c(0.98,0.97),
        axis.text=element_text(size=25),legend.box.background = element_rect(colour = "black",size=0.9),
        legend.text=element_text(size=20),legend.title=element_text(size=25,hjust=0.5),legend.key.size=unit(1.5,"cm"))+
  geom_text(aes(x=Species, y=CE+se+0.033,label=N), position=position_dodge(width=0.7), size=8)+
  geom_text(aes(x=Species, y=CE+se+0.07,label=groups), position=position_dodge(width=0.7), size=8)+
  geom_label(aes(x=0.67,y=0.582,label="b"),size=15,fill="white",hjust=0.5,vjust=0.5)+
  ylab("Conversion Efficiency (g/g)")+ xlab("Species")+
  guides(fill=guide_legend(nrow=2))


###########################
#Stomach volume of trout over the experiment
DissectionsComplete<-read.csv("Data/DataS6.csv")


###########################
#Comparing an isometric stomach volume to an allometric
summary(stomachISO<-lm(Stomach_V~Mass,data=DissectionsComplete))
summary(stomachALLO<-lm(log(Stomach_V)~log(Mass),data=DissectionsComplete))
#The isometric model is slightly better, however the allometric is still significant

#Generate predictions from the allometric model, then calculate the residuals from those as new volume measures
DissectionsComplete$residStomach<-exp(log(DissectionsComplete$Stomach_V)-predict.lm(stomachALLO,DissectionsComplete))

#Assessing whether the log-log residuals are unbiased with body mass, unlike the isometric
summary(normalV<-lm(residStomach~Mass,data=DissectionsComplete))  
#They are unbiased
#Therefore, this is the best stomach volume response variable to use


#Perform the ANOVA using the residual stomach volume from the allometric, log-log model
summary(ANOVAresidStomach<-aov(residStomach~Species*TOD*Treatment_Names+Size_Class,data=DissectionsComplete))

###########################
#Plot

#Generate a figure for this test, with the Time of Death and Species as the important variables to include
#Preparing the Tukey groups
volumetest<-with(DissectionsComplete,interaction(TOD,Species))
volumeAOV<-aov(residStomach~volumetest, data=DissectionsComplete)
volumeresults<-HSD.test(volumeAOV,"volumetest",group=TRUE)
volumeout<-volumeresults$group

#Calculating the mean and SEM for plotting
volumestats <- summarySE(DissectionsComplete, measurevar="residStomach",groupvars=c("TOD","Species"),na.rm=TRUE)
volumestats$TOD<-factor(volumestats$TOD,levels = c("Before Experiment","After Experiment"))
volumestats$trt<-paste(volumestats$TOD,volumestats$Species)

#Putting them together
volumemerge<-merge(volumestats,volumeout,by.x="residStomach")


#Plotting conversion volume of trout in experiment 
#(Figure A4)
ggplot(volumemerge, aes(x=Species, y=residStomach,fill=TOD)) + 
  geom_col(position=position_dodge(0.7),colour="black", width = 0.7)+
  geom_errorbar(aes(ymin=residStomach, ymax=residStomach+se),size=.5, width=.2,position=position_dodge(0.7)) + 
  scale_fill_manual(values=cbPalette,name="Time of Death")+scale_x_discrete(expand=c(0.2,0.2))+
  theme(legend.justification=c(0.5,1),legend.position="top",
        axis.text=element_text(size=25),legend.box.background = element_rect(colour = "black",size=1.5),
        legend.text=element_text(size=20),legend.title=element_text(size=25,hjust=0.5),legend.key.size=unit(1.5,"cm"))+
  geom_text(aes(x=Species, y=residStomach+se+0.15,label=N), position=position_dodge(width=0.7), size=8)+
  geom_text(aes(x=Species, y=residStomach+se+0.3,label=groups), position=position_dodge(width=0.7), size=8)+
  ylab("Residual Stomach Volume")+ xlab("Species")+
  guides(fill=guide_legend(nrow=1))



#(Figure A4) B&W friendly
ggplot(volumemerge, aes(x=Species, y=residStomach,fill=TOD)) + 
  geom_col(position=position_dodge(0.7),colour="black", width = 0.7)+
  geom_errorbar(aes(ymin=residStomach, ymax=residStomach+se),size=.5, width=.2,position=position_dodge(0.7)) + 
  scale_fill_viridis_d(begin=0.3,name="Time of Death")+
  scale_x_discrete(expand=c(0.2,0.2))+
  theme(legend.justification=c(0.5,1),legend.position="top",
        axis.text=element_text(size=25),legend.box.background = element_rect(colour = "black",size=1.5),
        legend.text=element_text(size=20),legend.title=element_text(size=25,hjust=0.5),legend.key.size=unit(1.5,"cm"))+
  geom_text(aes(x=Species, y=residStomach+se+0.15,label=N), position=position_dodge(width=0.7), size=8)+
  geom_text(aes(x=Species, y=residStomach+se+0.3,label=groups), position=position_dodge(width=0.7), size=8)+
  ylab("Residual Stomach Volume")+ xlab("Species")+
  guides(fill=guide_legend(nrow=1))


#




#




#




#




#




#




#


#####################################################################################
# Bioenergetics Model -----------------------------------------------------
#####################################################################################

set.seed(42)
#Bring in the temperature and discharge data averaged across 10 years at USGS gauge stations in Sherman and Hoxeyville, MI
Temp<-read.csv("Data/DataS7.csv")
Discharge<-read.csv("Data/DataS8.csv")




###########################
#Initialize an empty data frame to fill predictions with
outMass<-data.frame(Species=numeric(length=10000*12),
                    Temp=numeric(length=120000),
                    Time=numeric(length=120000),
                    initMass=numeric(length=120000),
                    finalMass=numeric(length=120000),
                    meanBinge=numeric(length=120000),
                    maxBinge=numeric(length=120000),
                    Replicate=numeric(length=120000),
                    sGrowth=numeric(length=120000))

fish=1

Consump<-matrix(ncol=366,nrow=120000)
CMAX<-matrix(ncol=366,nrow=120000)

###########################
#Model simulations conducted here

for (scenario in 1:10000) { #The whole thing will be simulated 10,000 times
  
  #Create a new randomized vector each replicate simulation series--for salmon and drift
  salmonStored<-rnorm(366,mean=3.79,sd=2.75) #The mass of salmon consumption is reflective of field observations
  salmonStored<-ifelse(salmonStored<0,0,salmonStored) #Don't allow any negatives
  
  #Day 1 is March 23rd
  #Drift rates throughout the year as g/m3 (From Wipfli and Gregovich 2002)
  drift<-c(runif(9,0.75,2.25),runif(30,2.1,9),runif(31,4.2,12),runif(30,3.5,6.2),runif(31,2,8),
           runif(31,1,3.5),runif(30,1.7,3.3),runif(31,0.5,4),runif(30,0.75,14),
           runif(31,0.75,3),runif(31,1.5,12.5),runif(29,1,8.8),runif(22,0.75,2.25))/1000 
  #Adjust relative to discharge pattern from Michigan streams
  foodStored<-Discharge[,2]*drift*0.1 #Proportion a trout can reasonably consume-10%

  for (species in 1:2) { #Species 1 is Brook Trout and Species 2 is Brown Trout
    for (temp in 1:2) { #Temp 1 is variable and temp 2 is stable
      for (time in 1:3) { #time 1 is early, 2 is mid, and 3 is late
        
        food<-foodStored #Restore the food vector for each scenario individual simulation
        ifsalmon<-salmonStored #Restore the salmon vector for each individual simulation
        
        if (species==1) { #Brook trout
          
          #From Hartman and Cox, 2008--all the parameter constants for the bioenergetics model
          #Consumption specific parameters
          CA=0.3013
          CB=-0.3055
          CQ=7.274
          CTO=20.9
          CTM=21
          CTL=24.05
          CK1=0.5
          CK4=0.203
          #Respiration Parameters
          RA=0.0132
          RB=-0.265
          RQ=4.5
          RTO=20.2
          RTL=25 #RTM in the Hartman and Cox and Hanson et al. definition, but just for concision
          RK1=1
          RK4=0.13
          ACT=2.89
          BACT=0.0405
          SDA=0.172
          #Egestion and Excretion Parameters
          FA=0.212
          FB=-0.222
          FG=0.631
          UA=0.0314
          UB=0.58
          UG=-0.299
          
          Jo2=4.63
          ED=rnorm(1,mean=5897,sd=245) #From my lab samples
          
        } else { #Brown Trout
          
          #From Dieterman 2004--same parameters but applicable to Brown Trout
          #Consumption Parameters
          CA=0.2161
          CB=-0.233
          CQ=3.8
          CTO=17.5
          CTM=17.5
          CTL=20.8
          CK1=0.23
          CK4=0.1
          #Respiration Parameters
          RA=0.0013
          RB=-0.269
          RQ=0.0938
          RTO=0.0234
          RTL=25
          RK1=1
          RK4=0.13
          ACT=9.7
          BACT=0.0405
          SDA=0.172
          #Egestion and Excretion Parameters
          FA=0.212
          FB=-0.222
          FG=0.631
          UA=0.0314
          UB=0.58
          UG=-0.299
          
          Jo2=4.63
          ED=rnorm(1,mean=5861,sd=403) #From lab samples
        }
        
        if (temp==1) { #Coolwater
          dailyTemp<-Temp[,2] 
        } else { #Coldwater
          dailyTemp<-Temp[,3] 
        }
        
        Mass<-numeric(length=367) #Need an empty vector with the right length to fill with masses over the simulation
        Mass[1]<-rnorm(1,mean=70,sd=10) #Generate initial mass for trout
        
        sMass<-numeric(length=31)
        binge<-data.frame(ccmax=numeric(length=31),binge=numeric(length=31)) #To keep P-values for binge
        bingeDay<-1
        
        inStomach=0 #Everything starts off as empty for the start of the simulation
        Cmax=1 #start with a dummy Cmax to make the first comparison to for feeding
         
        for (day in 1:366) { #This is where all the feeding happens, 1 day at a time
            
          EDprey=4534.5 #for chironomids
          P=0.52 #For all non-salmon times, they use this proportion of Cmax
          
          if (time==1) { #If time is 1 (Early)
            if (day>=150 & day<=180) { #Then the salmon come 8/19 - 9/18
              salmon=1 
            } else {
              salmon=0
            }
          } else if (time==2) { #If time is 2 (Middle)
            if (day>=180 & day<=210) { #Then salmon come 9/18 - 10/18
              salmon=1
            } else {
              salmon=0
            }
          } else { #If time is neither, ie 3 (Late)
            if (day>=210 & day<=240) { #Then the salmon come 10/18 - 11/17
              salmon=1
            } else {
              salmon=0
            }
          }
            
          if (inStomach<=0.5*Cmax*Mass[day]) { #Using the mass equal to that which these fish could feed daily in the experiment

            #Basic equations just using those parameters, but specific to the species
            G1=(1/(CTO-CQ))*log((0.98*(1-CK1))/(CK1*0.02))
            L1=exp(G1*(dailyTemp[day]-CQ))
            Ka=(CK1*L1)/(1+CK1*(L1-1))
            G2=(1/(CTL-CTM))*log((0.98*(1-CK4))/(CK4*0.02))
            L2=exp(G2*(CTL-dailyTemp[day]))
            Kb=(CK4*L2)/(1+CK4*(L2-1))
            #Respiration Equations--equation 1 for BKT (species 1) and equation 2 for BNT
            if (species==1) {
              V=(RTL-dailyTemp[day])/(RTL-RTO)      
              Z=log(RQ)*(RTL-RTO)
              Y=log(RQ)*(RTL-RTO+2)
              X=(Z^2*(1+(1+40/Y)^0.5)^2)/400
              f=(V^X)*(exp(X*(1-V)))
            } else {
              if (dailyTemp[day]>RTL) {
                VEL=RK1*(Mass[day]^RK4)
              } else {
                VEL=ACT*(Mass[day]^RK4)*exp(BACT*dailyTemp[day])
              }
              ACT=exp(RTO*VEL)
              f=exp(RQ*dailyTemp[day])
            }
            
            #Normal Bioenergetics
            Cmax=CA*(Mass[day]^CB)*Ka*Kb
                
            if (salmon==1) { #Dates of the salmon run and 2=salmon scenario
              EDprey=6533.8 #Salmon eggs
              C=ifsalmon[day]/Mass[day] #Consumption becomes reflective of field sampling
              P=C/Cmax #Calculate our own estimate of P for this time (binge-feed quantity)
              binge[bingeDay,]=c(P,salmon)
              sMass[bingeDay]=Mass[day]
              bingeDay=bingeDay+1
            } else if ((P*Cmax*Mass[day])<food[day]) { #If their bioenergetically allowed consumption is less than what's available...
              C=P*Cmax #Feed to the bioenergetically allowed amount
            } else {
              C=food[day]/Mass[day] #And if it's something less, then eat whatever is there
            }
              
            #Bioenergetic equations based on that consumption...
            R=RA*(Mass[day]^RB)*f*ACT #Respiration...
            Eg=FA*(dailyTemp[day]^FB)*(exp(FG*(P/Mass[day])))*C #Egestion...
            S=SDA*(inStomach+C-Eg) #Specific dynamic action...
            Ex=UA*(dailyTemp[day]^UB)*(exp(UG*(P/Mass[day])))*(C-Eg) #Excretion...
            
            #Growing Equation--putting it all together
            dMass=(((C-Eg-Ex)*EDprey-(R+S)*Jo2)/ED)
            Mass[day+1]=Mass[day]+dMass*Mass[day]
              
            inStomach=C*Mass[day] #The amount that was consumed is now filling the stomach
            
            Consump[fish,day]<-C
            CMAX[fish,day]<-Cmax
               
          } else {
              
            #All temperature dependent for BKT in Sweka 2004--pick the right rate based on current temp
            if (dailyTemp[day]<=6.7) { 
              b=-0.0058
            } else if (dailyTemp[day]>6.7 & dailyTemp[day]<=10.6) {
              b=-0.0135
            } else if (dailyTemp[day]>10.6 & dailyTemp[day]<=13.85) {
              b=-0.0149
            } else if (dailyTemp[day]>13.85 & dailyTemp[day]<=16.3) {
              b=-0.0150
            } else {
              b=-0.0153
            }
              
            inStomach=(1+24*b)*inStomach #Then remove some food in the stomach based on that time
              
            #Then go through all the same bioenergetics growth
            #Basic equations just using those parameters, but specific to the species
            G1=(1/(CTO-CQ))*log((0.98*(1-CK1))/(CK1*0.02))
            L1=exp(G1*(dailyTemp[day]-CQ))
            Ka=(CK1*L1)/(1+CK1*(L1-1))
            G2=(1/(CTL-CTM))*log((0.98*(1-CK4))/(CK4*0.02))
            L2=exp(G2*(CTL-dailyTemp[day]))
            Kb=(CK4*L2)/(1+CK4*(L2-1))
            #Respiration Equations--equation 1 for BNT (species 1) and equation 2 for BKT
            if (species==1) {
              V=(RTL-dailyTemp[day])/(RTL-RTO)      
              Z=log(RQ)*(RTL-RTO)
              Y=log(RQ)*(RTL-RTO+2)
              X=(Z^2*(1+(1+40/Y)^0.5)^2)/400
              f=(V^X)*(exp(X*(1-V)))
            } else {
              if (dailyTemp[day]>RTL) {
                VEL=RK1*(Mass[day]^RK4)
              } else {
                VEL=ACT*(Mass[day]^RK4)*exp(BACT*dailyTemp[day])
              }
              ACT=exp(RTO*VEL)
              f=exp(RQ*dailyTemp[day])
            }
            
            #Consumption and P have to be zero during this time 
            C=0
            Cmax=CA*(Mass[day]^CB)*Ka*Kb
            P=C/Cmax
            R=RA*(Mass[day]^RB)*f*(ACT)
            Eg=FA*(dailyTemp[day]^FB)*(exp(FG*(P/Mass[day])))*C
            S=SDA*(inStomach+C-Eg)
            Ex=UA*(dailyTemp[day]^UB)*(exp(UG*(P/Mass[day])))*(C-Eg)
              
            #Growing Equation
            dMass=(((C-Eg-Ex)*EDprey-(R+S)*Jo2)/ED)
            Mass[day+1]=Mass[day]+dMass*Mass[day]
            
            if (salmon==1) {
              binge[bingeDay,]=c(P,salmon)
              sMass[bingeDay]=Mass[day]
              bingeDay=bingeDay+1
            }
            
            Consump[fish,day]<-C
            CMAX[fish,day]<-Cmax
            
          }
        } 
        
        outMass[fish,]<-c(species,temp,time,Mass[1],Mass[367],
                          mean(subset(binge,ccmax>0)$ccmax),max(binge$ccmax),scenario,sMass[31]-sMass[1])
        fish=fish+1
      }    
    }
  }
  #Progress indicator, every 1,000 simulations an indication of progress is printed until 10,000
  if ((scenario/1000)%%1==0) {
    print(scenario)
  }
  
}

###########################
#Model Output
#Revert the numbers to words for ease of analysis/understanding
#Species 1 is Brook Trout, Species 2 is Brown Trout
#Temp 1 is Coolwater, 2 is Coldwater
#Time 1 is Early, 2 is Average, and 3 is Late

outMass$Species<-ifelse(outMass$Species==1, "Brook Trout", "Brown Trout") 
outMass$Temp<-ifelse(outMass$Temp==1, "Cool-water", "Cold-water")
outMass$Time<-ifelse(outMass$Time==1, "Early Spawn", ifelse(outMass$Time==2,"Average Spawn","Late Spawn"))
outMass$Time<-factor(outMass$Time,levels=c("Early Spawn","Average Spawn","Late Spawn"))
outMass$Env<-paste(outMass$Temp,outMass$Time,sep="\n")
outMass$Env<-factor(outMass$Env,levels=c("Cool-water\nEarly Spawn","Cool-water\nAverage Spawn",
                                         "Cool-water\nLate Spawn","Cold-water\nEarly Spawn",
                                         "Cold-water\nAverage Spawn","Cold-water\nLate Spawn"))


#Calculate growth over the simulation for each individual
outMass$dMass<-outMass$finalMass-outMass$initMass

#And convert to a growth rate as g/day
outMass$rate<-outMass$dMass/366 

#Proportion of growth during the salmon spawning window
outMass$sGrowth_prop<-outMass$sGrowth/outMass$dMass

#Percent of fish that gorged throughout... 
nrow(filter(outMass,meanBinge>1))/nrow(outMass)*100
#...and at least once during spawning
nrow(filter(outMass,maxBinge>1))/nrow(outMass)*100


###########################
#Plotting Output
###########################

#Growth Rates
#Summaries of the different variables for plotting
all<-summarySE(outMass,measurevar="rate",groupvars=c("Time","Temp","Species"))
all$Env<-paste(all$Temp,all$Time,sep="\n")


#Growth rate output--mean and 95% range whiskers
#(Figure 4a)
ggplot()+
  geom_jitter(data=outMass,aes(Env,rate,fill=Species),shape=21,size=0.5,
              position=position_jitterdodge(0.33),alpha=0.2,show.legend=F)+
  geom_linerange(data=all,aes(x=Env,ymin=rate-1.96*sd,ymax=rate+1.96*sd,color=Species),lwd=1.5,
                 position=position_dodge(0.75),show.legend=F)+
  geom_point(data=all,aes(Env,rate,fill=Species),shape=21,size=7,stroke=1.25,
             position=position_dodge(0.75))+
  scale_fill_manual(values=cbPalette)+
  scale_color_manual(values=cbPalette)+
  ylab("Growth Rate (g/day)")+xlab("Scenario")+
  scale_y_continuous(expand=c(0,0),limits=c(0.09,0.411))+
  scale_x_discrete(expand=c(0,0.2))+
  theme(legend.position=c(0.3,0.92),legend.background=element_rect(size=0.75,color="black"),
        legend.title=element_blank(),legend.margin=margin(-5,5,5,5),
        axis.text.y=element_text(size=25),axis.text.x=element_text(size=20),
        legend.text=element_text(size=20),legend.key.size=unit(1.5,"cm"))+
  guides(fill=guide_legend(nrow=1))



#(Figure 4a) B&W friendly
ggplot()+
  geom_jitter(data=outMass,aes(Env,rate,fill=Species),shape=21,size=0.5,
              position=position_jitterdodge(0.33),alpha=0.2,show.legend=F)+
  geom_linerange(data=all,aes(x=Env,ymin=rate-1.96*sd,ymax=rate+1.96*sd,color=Species),lwd=1.5,
                 position=position_dodge(0.75),show.legend=F)+
  geom_point(data=all,aes(Env,rate,fill=Species),shape=21,size=7,stroke=1.25,
             position=position_dodge(0.75))+
  scale_fill_viridis_d(begin=0.3)+
  scale_color_viridis_d(begin=0.3)+
  ylab("Growth Rate (g/day)")+xlab("Scenario")+
  scale_y_continuous(expand=c(0,0),limits=c(0.09,0.411),breaks=c(0.1,0.2,0.3,0.4),
                     labels=c("0.10","0.20","0.30","0.40"))+
  scale_x_discrete(expand=c(0,0.2))+
  theme(legend.position=c(0.3,0.92),legend.background=element_rect(size=0.75,color="black"),
        legend.title=element_blank(),legend.margin=margin(-5,5,5,5),
        axis.text.y=element_text(size=25),axis.text.x=element_text(size=20),
        legend.text=element_text(size=20),legend.key.size=unit(1.5,"cm"))+
  geom_label(aes(x=0.75,y=0.375,label="a"),size=15,fill="white",hjust=0.5,vjust=0.25)+
  guides(fill=guide_legend(nrow=1))



###########################
#Binge scores--quantifying with a C/Cmax in the loop
allB<-summarySE(outMass,measurevar="meanBinge",groupvars=c("Time","Temp","Species"))


#Create the maximum values of each day in the simulation lines for both C and Cmax, then bind together
CMAX2<-apply(CMAX,2,max) #Keep the maximum individual on each day (margin 2)
Consump2<-apply(Consump,2,max) #Keep the maximum individual on each day (margin 2)
C_Cmax<-data.frame(Time=seq(1,366),C=Consump2,Cmax=CMAX2) #Put them together

#C and Cmax comparison
#(Figure 4b)
ggplot(C_Cmax)+
  geom_line(aes(Time,C),lty=2,lwd=1.25)+
  geom_line(aes(Time,Cmax),lwd=1.25)+
  geom_ribbon(data=filter(C_Cmax,Time<150),aes(Time,ymax=Cmax,ymin=C),fill="firebrick3",alpha=0.4)+ #Shade the pre-salmon area in red
  geom_ribbon(data=filter(C_Cmax,Time>=150 & Time<=240),aes(Time,ymax=Cmax,ymin=C),fill="seagreen3",alpha=0.4)+ #Shade salmon area in green
  geom_ribbon(data=filter(C_Cmax,Time>240),aes(Time,ymax=Cmax,ymin=C),fill="firebrick3",alpha=0.4)+ #Shade post-salmon in red
  scale_x_continuous(breaks=c(10,40,71,101,132,163,193,224,254,285,316,345),
                     labels=c("A","M","J","J","A","S","O","N","D","J","F","M"),
                     expand=c(0,5))+ #Label with the month abbreviations
  scale_y_continuous(limits=c(0,0.25),expand=c(0,0.0025))+
  theme(panel.grid.minor.x = element_blank(),
        axis.text=element_text(size=25))+
  ylab(expression("C and C"[max]))+xlab("Month")


#(Figure 4b) new
ggplot(C_Cmax)+
  geom_line(aes(Time,C),lty=2,lwd=1.25)+
  geom_line(aes(Time,Cmax),lwd=1.25)+
  geom_ribbon(data=filter(C_Cmax,Time<150),aes(Time,ymax=Cmax,ymin=C),fill="firebrick3",alpha=0.4)+ #Shade the pre-salmon area in red
  geom_ribbon(data=filter(C_Cmax,Time>=150 & Time<=240),aes(Time,ymax=Cmax,ymin=C),fill="seagreen3",alpha=0.4)+ #Shade salmon area in green
  geom_ribbon(data=filter(C_Cmax,Time>240),aes(Time,ymax=Cmax,ymin=C),fill="firebrick3",alpha=0.4)+ #Shade post-salmon in red
  scale_x_continuous(breaks=c(10,40,71,101,132,163,193,224,254,285,316,345),
                     labels=c("A","M","J","J","A","S","O","N","D","J","F","M"),
                     expand=c(0,5))+ #Label with the month abbreviations
  scale_y_continuous(limits=c(0,0.25),expand=c(0,0.0025))+
  theme(panel.grid.minor.x = element_blank(),
        axis.text=element_text(size=25))+
  geom_label(aes(x=12,y=0.2308,label="b"),size=15,fill="white",hjust=0.5,vjust=0.5)+
  ylab(expression("C and C"[max]))+xlab("Month")

###########################
#Temperature and Spawning Windows Scenarios
#(Figure A2)
ggplot(data=Temp,aes(x=as.numeric(rownames(Temp))))+
  geom_rect(aes(ymin=-Inf,ymax=Inf,xmin=150,xmax=180),fill="#35608DFF")+
  geom_rect(aes(ymin=-Inf,ymax=Inf,xmin=180,xmax=210),fill="#2FB47CFF")+
  geom_rect(aes(ymin=-Inf,ymax=Inf,xmin=210,xmax=240),fill="#FDE725FF")+
  geom_point(aes(y=Mean_Temp_Coldwater),fill="white",shape=21)+ #Coldwater
  geom_point(aes(y=Mean_Temp_Coolwater),fill="black",shape=21)+ #Coolwater
  ylab(expression("Mean Daily Temperature ("^"o"*"C)"))+
  xlab("Month")+
  scale_x_continuous(breaks=c(10,40,71,101,132,163,193,224,254,285,316,345),
                     labels=c("A","M","J","J","A","S","O","N","D","J","F","M"),
                     expand=c(0,5))+
  theme(panel.grid.minor.x = element_blank(),axis.text=element_text(size=25))+
  geom_hline(aes(yintercept=mean(Mean_Temp_Coldwater)))+
  geom_hline(aes(yintercept=mean(Mean_Temp_Coolwater)))



#Now that the whole thing has run--you may return to any particular lines for closer investigation
#Using "Run" should work now
  #Saving the script to a unique directory can enable this from now on, but ensure the data files come with



















