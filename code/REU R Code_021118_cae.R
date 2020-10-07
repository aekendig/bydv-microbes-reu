## R file with reduced code

## Old code (011718 and previous) has models that incorporate molecular work and random effects, but these weren't useful or didn't converge


## Load required packages
library(ggplot2)
library(plyr)
library(Hmisc) #for binconf
library(cowplot)


## Read in data
#dat <- read.csv("~/Desktop/REU Project/Data/REU Infection Data Final CAE 9-24-16.csv")
setwd("~/Google Drive/Lab/REU")
dat=read.csv("REU Infection Data Final CAE 9-24-16.csv",header=T)
str(dat)


## Edit data

# Make new fields that round pav and rpv values of .5 up and down
unique(dat$pav)
unique(dat$rpv)
dat$pav_down=with(dat,ifelse(pav==0.5,0,as.numeric(pav)))
dat$rpv_down=with(dat,ifelse(rpv==0.5,0,as.numeric(rpv)))
dat$pav_up=with(dat,ifelse(pav==0.5,1,as.numeric(pav)))
dat$rpv_up=with(dat,ifelse(rpv==0.5,1,as.numeric(rpv)))

# Make presence a factor (0.5's rounded up, all count)
dat$rpv_pres=as.factor(dat$rpv_up)
dat$pav_pres=as.factor(dat$pav_up)

# Make infection a factor (0.5 rounded down, only 1's count)
dat$rpv_inf=as.factor(dat$rpv_down)
dat$pav_inf=as.factor(dat$pav_down)

# Make a column for interaction of inoculation type
dat$inoc_RPV=as.factor(with(dat,ifelse(disease=="RPV"|disease=="Co",1,0)))
dat$inoc_PAV=as.factor(with(dat,ifelse(disease=="PAV"|disease=="Co",1,0)))

# Identify infections that had contamination. 
# Only count really verified infections (down)
dat$pav_cont_down=with(dat,ifelse((disease=="RPV"|disease=="Healthy")&pav_down==1,1,0))
sum(dat$pav_cont_down,na.rm=T)	#0
dat$rpv_cont_down=with(dat,ifelse((disease=="PAV"|disease=="Healthy")&rpv_down==1,1,0))
sum(dat$rpv_cont_down,na.rm=T)	#7

# Count all detected infections (up)
dat$pav_cont_up=with(dat,ifelse((disease=="RPV"|disease=="Healthy")&pav_up==1,1,0))
sum(dat$pav_cont_up,na.rm=T)	#5
dat$rpv_cont_up=with(dat,ifelse((disease=="PAV"|disease=="Healthy")&rpv_up==1,1,0))
sum(dat$rpv_cont_up,na.rm=T)	#26

# Make new field for average chlorophyll measurements
dat$chlorophyll_avg <- (dat$chlorophyll_1 + dat$chlorophyll_2 + dat$chlorophyll_3)/3

# Reset Soil Contrasts to be More meaningful
#                        	    A   D  H  Sterile
contrasts(dat$soil) <- cbind(c( 1,  1,  1,  -3), #Sterile vs others
                             c(-2,  1,  1,	0),  # D&H vs A   
                             c(-1,  0,  1,	0))  # H vs A
                             
# Subset data to remove untested samples
dat2=subset(dat,!is.na(dat$pav)&!is.na(dat$rpv))
nrow(dat)-nrow(dat2) # lost 6 samples

# Remove contaminated samples based on "up"
dat3=subset(dat2,pav_cont_up==0&rpv_cont_up==0)
nrow(dat2)-nrow(dat3) # lost 31 samples

# Remove contaminated samples based on "down"
dat4=subset(dat2,pav_cont_down==0&rpv_cont_down==0)
nrow(dat2)-nrow(dat4) # lost 7 samples

# Separate PAV and RPV datasets
pDat3<-subset(dat3,disease=="PAV"|disease=="Co")
rDat3<-subset(dat3,disease=="RPV"|disease=="Co")
pDat4<-subset(dat4,disease=="PAV"|disease=="Co")
rDat4<-subset(dat4,disease=="RPV"|disease=="Co")

# Change factor levels
pDat3$disease<-factor(pDat3$disease,levels=c("PAV","Co"))
rDat3$disease<-factor(rDat3$disease,levels=c("RPV","Co"))
pDat4$disease<-factor(pDat4$disease,levels=c("PAV","Co"))
rDat4$disease<-factor(rDat4$disease,levels=c("RPV","Co"))


## Analysis 1: Biomass and chloropyll

# Sample size
ddply(dat3,.(nutrient,soil,disease),summarise,n=length(biomass)) #5-10

# 1.a: Only inoculation type matters (use dat2)
bioModA=lm(biomass~nutrient*soil*inoc_RPV*inoc_PAV,data=dat2)
summary(bioModA)
# N addition, RPV, RPV:PAV
chlorModA=lm(chlorophyll_avg~nutrient*soil*inoc_RPV*inoc_PAV,data=dat2)
summary(chlorModA)
# N addition

# 1.b: Infected samples must be detected (0.5-1, up, pres, dat3)
datB=subset(dat3,(disease=="PAV"&pav_up==1)|(disease=="RPV"&rpv_up==1)|(disease=="Co"&pav_up==1&rpv_up==1))
nrow(dat2)-nrow(datB) # 197 samples excluded (31 from cont.)
bioModB=lm(biomass~nutrient*soil*rpv_pres*pav_pres,data=datB)
summary(bioModB)
# Nothing sig, many comparisons couldn't be made
chlorModB=lm(chlorophyll_avg~nutrient*soil*rpv_pres*pav_pres,data=datB)
summary(chlorModB)
# Nothing sig, many comparisons couldn't be made

# 1.c: Infected samples must be a full band (1, down, inf, dat4)
datC=subset(dat4,(disease=="PAV"&pav_down==1)|(disease=="RPV"&rpv_down==1)|(disease=="Co"&pav_down==1&rpv_down==1))
nrow(dat2)-nrow(datC) # 226 samples excluded (7 from cont.)
bioModC=lm(biomass~nutrient*soil*rpv_inf*pav_inf,data=datC)
summary(bioModC)
# Nothing sig, many comparisons couldn't be made
chlorModC=lm(chlorophyll_avg~nutrient*soil*rpv_inf*pav_inf,data=datC)
summary(chlorModC)
# Sig effect of soil 1, soil 2, interaction between soil 2 and PAV, and soil3 and PAV, many comparisons couldn't be made

# 1.d: Undetected samples may be considered healthy or singly infected (dat3 - still exclude accidental infections because we don't know when they occurred)
bioModD=lm(biomass~nutrient*soil*rpv_pres*pav_pres,data=dat3)
summary(bioModD)
# N addition, RPV
chlorModD=lm(chlorophyll_avg~nutrient*soil*rpv_pres*pav_pres,data=dat3)
summary(chlorModD)
# N addition

# Save models 1d
write.csv(summary(bioModD)$coefficients,"BiomassMod_WeakAndStrong_011818.csv",row.names=T)
write.csv(summary(chlorModD)$coefficients,"ChlorophyllMod_WeakAndStrong_011818.csv",row.names=T)

# 1.e: Weak band samples may be considered healthy or singly infected (dat4 - still exclude accidental infections because we don't know when they occurred)
bioModE=lm(biomass~nutrient*soil*rpv_pres*pav_pres,data=dat4)
summary(bioModE)
# N addition, N:RPV
chlorModE=lm(chlorophyll_avg~nutrient*soil*rpv_pres*pav_pres,data=dat4)
summary(chlorModE)
# N addition

# Save model 1e
write.csv(summary(bioModE)$coefficients,"BiomassMod_JustStrong_011818.csv",row.names=T)


## Analysis 2: Infection success

# Sample size
ddply(pDat3,.(nutrient,soil,disease),summarise,n=length(pav_up)) # 7-10
ddply(rDat3,.(nutrient,soil,disease),summarise,n=length(rpv_up)) #9-10

# 2.a: Infection is anything detected (up, pres, dat3)
pModA=glm(pav_up~nutrient*soil*rpv_pres,data=pDat3,family=binomial)
summary(pModA)
# Soil 3 is sig
rModA=glm(rpv_up~nutrient*soil*pav_pres,data=rDat3,family=binomial)
summary(rModA)
# N:soil1

# Save models 2a
write.csv(summary(rModA)$coefficients,"RPVInfMod_WeakAndStrong_011818.csv",row.names=T)
write.csv(summary(pModA)$coefficients,"PAVInfMod_WeakAndStrong_011818.csv",row.names=T)

# 2.b: Infection is a full band (down, inf, dat4)
pModB=glm(pav_down~nutrient*soil*rpv_inf,data=pDat4,family=binomial)
summary(pModB)
# Nothing sig
rModB=glm(rpv_down~nutrient*soil*pav_inf,data=rDat4,family=binomial)
summary(rModB)
# Nothing sig

# 2.c: Infection is anything detected, but inoculation matters (up, dis, dat3)
pModC=glm(pav_up~nutrient*soil*disease,data=pDat3,family=binomial)
summary(pModC)
# Nothing is sig
rModC=glm(rpv_up~nutrient*soil*disease,data=rDat3,family=binomial)
summary(rModC)
# Coinoculation

# 2.d: Infection is a full band, but inoculation matters (down, dis, dat4)
pModD=glm(pav_down~nutrient*soil*disease,data=pDat4,family=binomial)
summary(pModD)
# Nothing sig
rModD=glm(rpv_down~nutrient*soil*disease,data=rDat4,family=binomial)
summary(rModD)
# Nothing sig


## Figures

# Figure 1: Nutrient, microbes present/absent on RPV infection

rDat3$presence_microbes=ifelse(rDat3$soil=="Sterile","absent","present")

rDatBin<-ddply(rDat3,.(nutrient, presence_microbes),summarise,meanRate=binconf(sum(rpv_up),length(rpv_up))[1],lowerCI=binconf(sum(rpv_up),length(rpv_up))[2],upperCI=binconf(sum(rpv_up),length(rpv_up))[3],nSamps=length(rpv_up))

rDatBin$nutrient=revalue(rDatBin$nutrient,c("Ctrl"="low","N"="high"))

pdf("RPVInf_MicrobesN_011818.pdf",width=3.5,height=3.5)
ggplot(rDatBin, aes(x=presence_microbes,y=meanRate))+geom_point(size=4,aes(shape=nutrient),position=position_dodge(0.25))+geom_errorbar(aes(ymin=lowerCI,ymax=upperCI,group=nutrient),width=0.1,position=position_dodge(0.25))+labs(x = "Soil microbes", y = "RPV infection rate")+scale_shape_manual(values=c(19,17),name="N addition")+theme_bw()+theme(axis.text=element_text(size=10,colour="black"),axis.title=element_text(size=12),strip.background=element_blank(),panel.background=element_blank(),panel.grid.major=element_blank(),panel.grid.minor=element_blank(),axis.line=element_line(colour="black"),panel.border=element_rect(color="black",fill=NA),strip.text=element_text(size=12),legend.text=element_text(size=10),legend.title=element_text(size=12),legend.position=c(0.75,0.85))+ylim(0,1)
dev.off()

# Figure 2: Soil microbial community and PAV

pDatBin<-ddply(pDat3,.(soil),summarise,meanRate=binconf(sum(pav_up),length(pav_up))[1],lowerCI=binconf(sum(pav_up),length(pav_up))[2],upperCI=binconf(sum(pav_up),length(pav_up))[3],nSamps=length(pav_up))

pDatBin$soil=revalue(pDatBin$soil,c("A"="0","D"="34","H"="272","Sterile"="sterile"))
pDatBin$soil=factor(pDatBin$soil,levels=c("sterile","0","34","272"))

pdf("PAVInf_Microbes_011818.pdf",width=3.5,height=3.5)
ggplot(pDatBin, aes(x=soil,y=meanRate))+geom_point(size=4,position=position_dodge(0.25))+geom_errorbar(aes(ymin=lowerCI,ymax=upperCI),width=0.1,position=position_dodge(0.25))+xlab(expression(paste("Soil (kg N ",ha^{-1},yr^{-1},")",sep="")))+ylab("PAV infection rate")+theme_bw()+theme(axis.text=element_text(size=10,colour="black"),axis.title=element_text(size=12),strip.background=element_blank(),panel.background=element_blank(),panel.grid.major=element_blank(),panel.grid.minor=element_blank(),axis.line=element_line(colour="black"),panel.border=element_rect(color="black",fill=NA),strip.text=element_text(size=12),legend.text=element_text(size=10),legend.title=element_text(size=12),legend.position=c(0.75,0.85))+ylim(0,1)
dev.off()

# Proportions for text
mean(subset(pDat3,soil=="H")$pav_up)
mean(subset(pDat3,soil!="H")$pav_up)

# Figure 3: N addition and RPV infection on biomass and chlorophyll

bDatSum<-ddply(dat3,.(nutrient,rpv_pres),summarise,meanBio=mean(biomass),seBio=sd(biomass)/sqrt(length(biomass)),nBio=length(biomass))
cDatSum<-ddply(dat3,.(nutrient),summarise,meanChlor=mean(chlorophyll_avg),seChlor=sd(chlorophyll_avg)/sqrt(length(chlorophyll_avg)),nChlor=length(chlorophyll_avg))

bDatSum$nutrient=revalue(bDatSum$nutrient,c("Ctrl"="low","N"="high"))
cDatSum$nutrient=revalue(cDatSum$nutrient,c("Ctrl"="low","N"="high"))
bDatSum$rpv_pres=revalue(bDatSum$rpv_pres,c("0"="uninfected","1"="infected"))

bPlot=ggplot(bDatSum, aes(x=nutrient,y=meanBio,group=rpv_pres))+geom_point(size=4,position=position_dodge(0.25),aes(shape=rpv_pres))+geom_errorbar(aes(ymin=meanBio-seBio,ymax=meanBio+seBio),width=0.1,position=position_dodge(0.25))+scale_shape_manual(values=c(1,16),name="RPV infection")+xlab("N addition")+ylab("Biomass (g)")+theme_bw()+theme(axis.text=element_text(size=10,colour="black"),axis.title=element_text(size=12),strip.background=element_blank(),panel.background=element_blank(),panel.grid.major=element_blank(),panel.grid.minor=element_blank(),axis.line=element_line(colour="black"),panel.border=element_rect(color="black",fill=NA),strip.text=element_text(size=12),legend.text=element_text(size=10),legend.title=element_text(size=12),legend.position=c(0.27,0.75))
bPlot
cPlot=ggplot(cDatSum, aes(x=nutrient,y=meanChlor))+geom_point(size=4,position=position_dodge(0.25),aes(shape=rpv_pres),shape=15)+geom_errorbar(aes(ymin=meanChlor-seChlor,ymax=meanChlor+seChlor),width=0.1,position=position_dodge(0.25))+xlab("N addition")+ylab("Chlorophyll content (SPAD)")+theme_bw()+theme(axis.text=element_text(size=10,colour="black"),axis.title=element_text(size=12),strip.background=element_blank(),panel.background=element_blank(),panel.grid.major=element_blank(),panel.grid.minor=element_blank(),axis.line=element_line(colour="black"),panel.border=element_rect(color="black",fill=NA),strip.text=element_text(size=12),legend.text=element_text(size=10),legend.title=element_text(size=12),legend.position=c(0.3,0.7))

pdf("PlantTraits_011818.pdf",width=6,height=3)
plot_grid(bPlot,cPlot,labels=c('A','B'),hjust=c(-5,-4.2),vjust=2)
dev.off()


## Aditional figure added CAE 2-11-18

# recalculate rDatBin and pDatBin

rDatBin2 <- ddply(rDat3,.(nutrient, soil),summarise,meanRate=binconf(sum(rpv_up),length(rpv_up))[1],lowerCI=binconf(sum(rpv_up),length(rpv_up))[2],upperCI=binconf(sum(rpv_up),length(rpv_up))[3],nSamps=length(rpv_up))

rDatBin2$nutrient=revalue(rDatBin2$nutrient,c("Ctrl"="low","N"="high"))
rDatBin2$Virus <- "RPV"
soilLUT <- c("Sterile" = "sterile", "A" = "0", "D" = "34", "H" = "272")
rDatBin2$soil <- as.factor(soilLUT[as.character(rDatBin2$soil)]) 
rDatBin2$soil=factor(rDatBin2$soil,levels=c("sterile","0","34","272"))

pDatBin2<-ddply(pDat3,.(nutrient, soil),summarise,meanRate=binconf(sum(pav_up),length(pav_up))[1],lowerCI=binconf(sum(pav_up),length(pav_up))[2],upperCI=binconf(sum(pav_up),length(pav_up))[3],nSamps=length(pav_up))

pDatBin2$soil <- as.factor(soilLUT[as.character(pDatBin2$soil)])
pDatBin2$soil=factor(pDatBin2$soil,levels=c("sterile","0","34","272"))
pDatBin2$nutrient=revalue(pDatBin2$nutrient,c("Ctrl"="low","N"="high"))
pDatBin2$Virus <- "PAV"

DatBin <- rbind(rDatBin2, pDatBin2)

ggplot(DatBin, aes(x=soil,y=meanRate)) + 
  geom_point(size=4,aes(shape=nutrient, color = Virus, group = interaction(nutrient, Virus)), position = position_dodge(width = .5)) + 
  scale_color_grey() + 
  geom_errorbar(data = DatBin, aes(ymin=lowerCI,ymax=upperCI, group = interaction(nutrient, Virus)), width = .01, position = position_dodge(width = .5)) +
  labs(x = "Soil", y = "Infection rate")+scale_shape_manual(values=c(19,17),name="N addition")+theme_bw()+theme(axis.text=element_text(size=10,colour="black"),axis.title=element_text(size=12),strip.background=element_blank(),panel.background=element_blank(),panel.grid.major=element_blank(),panel.grid.minor=element_blank(),axis.line=element_line(colour="black"),panel.border=element_rect(color="black",fill=NA),strip.text=element_text(size=12),legend.text=element_text(size=10),legend.title=element_text(size=12),legend.position=c(.8,.15), legend.key = element_blank(), legend.box = "horizontal") + scale_y_continuous(lim = c(0,1), breaks = c(0, .25, .5, .75, 1))

# chlorophyll by soil

# Figure 3: N addition and RPV infection on biomass and chlorophyll

bDatSum2<-ddply(dat3,.(nutrient,rpv_pres, soil),summarise,meanBio=mean(biomass),seBio=sd(biomass)/sqrt(length(biomass)),nBio=length(biomass))

bDatSum2$nutrient=revalue(bDatSum2$nutrient,c("Ctrl"="low","N"="high"))
bDatSum2$rpv_pres=revalue(bDatSum2$rpv_pres,c("0"="uninfected","1"="infected"))
bDatSum2$soil <- bDatSum2$soil <- as.factor(soilLUT[as.character(bDatSum2$soil)])
bDatSum2$soil=factor(bDatSum2$soil,levels=c("sterile","0","34","272"))

ggplot(bDatSum2, aes(x=soil,y=meanBio, color = nutrient, group = interaction(nutrient, rpv_pres))) + 
  geom_point(size=4,position=position_dodge(0.25),aes(shape=rpv_pres)) +
  scale_color_grey() + geom_errorbar(aes(ymin=meanBio-seBio,ymax=meanBio+seBio, group = interaction(nutrient, rpv_pres)),width=0.1,position=position_dodge(0.25))

#+scale_shape_manual(values=c(1,16),name="RPV infection")+xlab("N addition")+ylab("Biomass (g)")+theme_bw()+theme(axis.text=element_text(size=10,colour="black"),axis.title=element_text(size=12),strip.background=element_blank(),panel.background=element_blank(),panel.grid.major=element_blank(),panel.grid.minor=element_blank(),axis.line=element_line(colour="black"),panel.border=element_rect(color="black",fill=NA),strip.text=element_text(size=12),legend.text=element_text(size=10),legend.title=element_text(size=12),legend.position=c(0.27,0.75))


cDatSum2<-ddply(dat3,.(nutrient, soil),summarise,meanChlor=mean(chlorophyll_avg),seChlor=sd(chlorophyll_avg)/sqrt(length(chlorophyll_avg)),nChlor=length(chlorophyll_avg))

cDatSum2$nutrient=revalue(cDatSum2$nutrient,c("Ctrl"="low","N"="high"))

cDatSum2$soil <- cDatSum2$soil <- as.factor(soilLUT[as.character(cDatSum2$soil)])
cDatSum2$soil=factor(cDatSum2$soil,levels=c("sterile","0","34","272"))

ggplot(cDatSum2, aes(x=soil,y=meanChlor, color = nutrient, group = nutrient)) +
  geom_point(size=4,position=position_dodge(0.25),aes(shape=rpv_pres),shape=15) + 
  scale_color_grey() +
  geom_errorbar(aes(ymin=meanChlor-seChlor,ymax=meanChlor+seChlor, group = nutrient),width=0.1,position=position_dodge(0.25))

#+xlab("N addition")+ylab("Chlorophyll content (SPAD)")+theme_bw()+theme(axis.text=element_text(size=10,colour="black"),axis.title=element_text(size=12),strip.background=element_blank(),panel.background=element_blank(),panel.grid.major=element_blank(),panel.grid.minor=element_blank(),axis.line=element_line(colour="black"),panel.border=element_rect(color="black",fill=NA),strip.text=element_text(size=12),legend.text=element_text(size=10),legend.title=element_text(size=12),legend.position=c(0.3,0.7))

