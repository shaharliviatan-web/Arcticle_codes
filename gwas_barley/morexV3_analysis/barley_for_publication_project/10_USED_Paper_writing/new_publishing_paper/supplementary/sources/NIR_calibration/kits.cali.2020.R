##This script is to compare the 3  analyses-kits,calibration, and scans of 2020 in order to
##understand if we can use the data or not. 

##I used reshape2 package and melt function to build the right table- 3 repetitions of each site, each rep. represents the
##type of analyses.

##I need to do posthoc to see if the pvalue in between GI is not muvhak.



library(ggplot2)
library(reshape2)
library(dplyr)
library(GGally)
Protein<-read.csv("/Users/noampintel/Documents/Noam/RESEARCH/Analyses/NIR/NIR2020.vs.NIRcali/Nutrients/Protein.csv")
Starch<-read.csv("/Users/noampintel/Documents/Noam/RESEARCH/Analyses/NIR/NIR2020.vs.NIRcali/Nutrients/Starch.csv")
Betaglucan<-read.csv("/Users/noampintel/Documents/Noam/RESEARCH/Analyses/NIR/NIR2020.vs.NIRcali/Nutrients/Betaglucan.csv")


###Correlations between analyses.

cor.test(Starch$Starch.kits,Starch$Starch.cali) #cor= 0.7848839 
cor.test(Starch$Starch.kits,Starch$Starch.2020) #cor= 0.4659546
cor.test(Starch$Starch.2020,Starch$Starch.cali) #cor= 0.5479516 
s<-ggpairs(data = Starch,columns = 2:4, title = "")

cor.test(Protein$Protein.kits, Protein$Protein.cali) #cor= 0.9149877 
cor.test(Protein$Protein.kits, Protein$Protein.2020) #cor= 0.5563798
cor.test(Protein$Protein.cali, Protein$Protein.2020) #cor= 0.6492421
p<-ggpairs(data = Protein,columns = 2:4, title = "")

cor.test(Betaglucan$Betaglucan.kits, Betaglucan$Betaglucan.cali)  #cor= 1
cor.test(Betaglucan$Betaglucan.kits,Betaglucan$Betaglucan.2020)   #cor= 0.5135945 
cor.test(Betaglucan$Betaglucan.cali, Betaglucan$Betaglucan.2020) #cor= 0.5135945
b<-ggpairs(data = Betaglucan,columns = 2:4, title = "")


###Bocplots to see anova###

######################PROTEIN

Protein.long<-melt(Protein, id.vars = "GI")
ggplot(Protein.long, aes(fill=variable, y=value, x=as.factor(GI))) + 
  geom_bar(position="dodge", stat="identity") 
  

GI.aov.protein<-aov(GI~value, data = Protein.long)
summary(GI.aov)

box1<-ggplot(data= Protein.long,aes(x=as.factor(GI), y=value)) +
  geom_point()+
  geom_boxplot()+
  stat_compare_means(method = "anova")+
  theme(axis.text.x = element_text(angle = -45))+
  ggtitle("Protein")

##colored points to see if one was highest always. if so, maybe do a correction
box11<-ggplot(data= Protein.long,aes(x=as.factor(GI), y=value, color=variable)) +
  geom_point()+
  theme(axis.text.x = element_text(angle = -45))+
  ggtitle("Protein")

######################STARCH

Starch.long<-melt(Starch, id.vars = "GI")
ggplot(Starch.long, aes(fill=variable, y=value, x=as.factor(GI))) + 
  geom_bar(position="dodge", stat="identity")


starch.aov<-aov(value~as.factor(GI), data = Starch.long)
summary(starch.aov)
starch.aov$coefficients


###extract pvalue from anova list
sum_test = unlist(summary(starch.aov))
names(sum_test)
sum_test["Pr(>F)1"]
####

##boxplot to see differences
box2<-ggplot(data= Starch.long,aes(x=as.factor(GI), y=value)) +
  geom_point()+
  geom_boxplot()+
  stat_compare_means(method = "anova")+
  theme(axis.text.x = element_text(angle = -45))+
  ggtitle("Starch")

##colored points to see if one was highest always.
box22<-ggplot(data= Starch.long,aes(x=as.factor(GI), y=value, color=variable)) +
  geom_point()+
  theme(axis.text.x = element_text(angle = -45))+
  ggtitle("Starch")
library(gridExtra)
grid.arrange(p1,p2, nrow=2)



######################BETAGLUCAN

Betaglucan.long<-melt(Betaglucan, id.vars = "GI")
ggplot(Betaglucan.long, aes(fill=variable, y=value, x=as.factor(GI))) + 
  geom_bar(position="dodge", stat="identity")

beta.aov<-aov(value~as.factor(GI), data = Betaglucan.long)
summary(beta.aov)
beta.aov$coefficients

box3<-ggplot(data= Betaglucan.long,aes(x=as.factor(GI), y=value)) +
  geom_point()+
  geom_boxplot()+
  stat_compare_means(method = "anova")+
  theme(axis.text.x = element_text(angle = -45))+
  ggtitle("Betaglucan")

##colored points to see if one was highest always.
box33<-ggplot(data= Betaglucan.long,aes(x=as.factor(GI), y=value, color=variable)) +
  geom_point()+
  theme(axis.text.x = element_text(angle = -45))+
  ggtitle("Betaglucan")



library(gridExtra)
grid.arrange(box1,box2,box3,box11,box22,box33, nrow=2)





