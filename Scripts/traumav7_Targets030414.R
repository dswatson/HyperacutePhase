## trauma analysis ###

##################################################################################################################################
### 	Version 7: Exported without Samples 226_0hr and 545_0hr & Patient 303 modified to NO_MODS   ##############################
library(limma)
library(sva)
library(VennDiagram)
setwd("C:/Enserio/Barts/Projects/Trauma/limma/rawdata/")

targets <- read.table("Targets_030414.txt",sep="\t",header=T,stringsAsFactors=F,check.names=F)
rawdata <- read.ilmn(files="Sample_probe_Profile_5_12_13.txt",ctrlfiles="Control_probe_Profile_5_12_13.txt")

## arrange the columns of targets as expression data 
x <- targets$Array
y <- colnames(rawdata$E)
z <- as.numeric()
for(i in 1:length(x)){
	z <- c(z, which(y[i] == x))
}
temp <- targets[z,]
targets<-temp
###########
## batch effect correction combat 
edata <- rawdata$E
batch <- targets$Batch
mod <- model.matrix(~as.factor(Covariate3), data=targets)
mod0 <- model.matrix(~1,data=targets)
num.sv(edata,mod,method="leek")

png("PLOTS_030414.png")
combat_edata = ComBat(dat=edata, batch=batch,mod=mod,numCovs=NULL, par.prior=TRUE,prior.plots=T)
dev.off()

setwd("C:/Enserio/Barts/Projects/Trauma/limma/BatchNbckCorData_030414/BatchPlots/")
## Batchs
b1<-which(targets[,3]==1)
b2<-which(targets[,3]==2)
x <- which(colnames(rawdata$E) %in% targets[b1,2])
y <- which(colnames(rawdata$E) %in% targets[b2,2])
cols <- rep("red",118)
cols[x] <- "blue"

png("BxPlot_IntensitiesRaw_030414.png",width=1000,height=700)
boxplot(log2(rawdata$E),col=cols,range=0,ylab="log2 intensity",ylim=c(6,16))
legend(1,16, legend=c("GC", "UCL"), fill=c("blue","red"), cex=0.8)
dev.off()
png("BxPlot_IntensitiesRawbyBatch_030414.png",width=1000,height=700)
boxplot(log2(rawdata$E[,x]),col="blue",at=1:50,range=0,ylab="log2 intensity",ylim=c(6,16),xlim=c(1,118))
boxplot(log2(rawdata$E[,y]),col="red",at=51:118,range=0,ylab="log2 intensity",add=TRUE)
legend(1,16, legend=c("GC", "UCL"), fill=c("blue","red"), cex=0.8)
dev.off()

png("BxPlot_IntensitiesCorBatch_030414.png",width=1000,height=700)
boxplot(log2(combat_edata[,x]),col="blue",at=1:50,range=0,ylab="log2 intensity",ylim=c(6,16),xlim=c(1,118))
boxplot(log2(combat_edata[,y]),col="red",at=51:118,range=0,ylab="log2 intensity",add=TRUE)
legend(1,16, legend=c("GC", "UCL"), fill=c("blue","red"), cex=0.8)
dev.off()

########## 
# place combat data in the set
rawdata$E <- combat_edata
## background and normalization 
bkcor_data<- neqc(rawdata)

png("BxPlot_NormCorIntensitiesBatch_030414.png",width=1000,height=700)
boxplot(log2(bkcor_data$E[,x]),col="blue",at=1:50,range=0,ylab="log2 intensity",ylim=c(2,5),xlim=c(1,118))
boxplot(log2(bkcor_data$E[,y]),col="red",at=51:118,range=0,ylab="log2 intensity",add=TRUE)
legend(1,5, legend=c("GC", "UCL"), fill=c("blue","red"), cex=0.8)
dev.off()
####
rm(b1,b2,batch,cols,combat_edata,edata,i,mod,mod0,temp,x,y,z)
setwd("C:/Enserio/Barts/Projects/Trauma/limma/BatchNbckCorData_030414")
########################################################
pe <- propexpr(rawdata)
# Filter out probes that are not expressed. We keep probes that are expressed in at least
# three arrays according to a detection p-values of 5%:
expressed <- rowSums(bkcor_data$other$Detection < 0.05) >= 3
data <- bkcor_data[expressed,]
> dim(bkcor_data)
[1] 45248   118
> dim(data)
[1] 29385   118

############################
# # > names(targets)
> names(targets)
 [1] "Subject"                        "Array"                          "Batch"                          "Covariate1"                    
 [5] "Covariate2"                     "Covariate3"                     "Severity"                       "adverse_outcome"               
 [9] "MODS"                           "Infection"                      "Gender"                         "Age"                           
[13] "Time_of_admission_post_injury"  "MOI"                            "ISS"                            "ISS_grp"                       
[17] "ISS_grp_h"                      "CSLprior_ACIT_blood_draw"       "Lactate"                        "Shock3Cat"                     
[21] "Shock3Cat_h"                    "ShockYN"                        "ShockYN_h"                      "BD"                            
[25] "BD_grp"                         "BD_grp_h"                       "0_Died_1_Home_2_Transferred"    "Total_LOS"                     
[29] "28_days_outcome_Dead_0_Alive_1" "MOF_Yes_1_No_0"                 "Infection_Yes_1_No_0"           "2h_Hb"                         
[33] "2h_Hct"                         "2h_Neut"                        "2h_Mono"                        "2h_Lymp"                       
[37] "2h_Eosin"                       "24h_Hb"                         "24h_Hct"                        "24h_Neut"                      
[41] "24h_Mono"                       "24h_Lymp"                       "24h_Eosin"                      "48h_Hb"                        
[45] "48h_Hct"                        "48h_Neut"                       "48h_Mono"                       "48h_Lymp"                      
[49] "48h_Eosin"                      "72h_Hb"                         "72h_Hct"                        "72h_Neut"                      
[53] "72h_Mono"                       "72h_Lymp"                       "72h_Eosin"                     


#########################################################################################################
################################## Differential expression ##############################################

############### Controls VS Critical #############################################
## Number of Samples
# > summary(as.factor(targets$Severity))
# Control_0h  Control_24h  Control_72h  Critical_0h Critical_24h Critical_72h  Moderate_0h Moderate_24h Moderate_72h     RNA_ctrl 
           # 5            7            6           27           26           31            5            5            5            1 

# # > targets[which(targets[,1]=="303"),]
    # # Subject   Array Batch  Covariate1 Covariate2      Covariate3     Severity adverse_outcome    MODS Infection Gender Age Time_of_admission_post_injury
# # 7       303  303 0h     1 No_MODS_Inf        T0h  No_MODS_Inf_0h  Critical_0h             Yes NO_MODS       Yes      1  35                           120
# # 28      303 303 24h     1 No_MODS_Inf       T24h No_MODS_Inf_24h Critical_24h             Yes NO_MODS       Yes      1  35                           120
# # 102     303 303 72h     2 No_MODS_Inf       T72h No_MODS_Inf_72h Critical_72h             Yes NO_MODS       Yes      1  35                           120

		   
### cluster factors
setwd("../BatchNbckCorData_030414/ClusterPlots/")
png("MDS_Severity_030414.png",width=800,height=800)
plotMDS(data,labels=targets$Severity)
dev.off()

## Control 0hr Vs Critical 0h	
setwd("C:/Enserio/Barts/Projects/Trauma/limma/BatchNbckCorData_030414")	   
control_critical <- factor(targets$Severity)
# > control_critical
  # [1] Moderate_0h  Control_0h   Control_0h   Critical_0h  Moderate_0h  Critical_0h  ...........
# Levels: Control_0h Control_24h Control_72h Critical_0h Critical_24h Critical_72h Moderate_0h Moderate_24h Moderate_72h RNA_ctrl
design <- model.matrix(~0+control_critical)
colnames(design) <- levels(control_critical)
dupcor<- duplicateCorrelation(data,design,block=targets$Subject,ndups=2)
fit <- lmFit(data,design,block=targets$Subject,correlation=dupcor$consensus.correlation)

#### Critical vs Control General
contrasts <- makeContrasts(Diff24=(Critical_24h-Critical_0h)-(Control_24h-Control_0h),Diff72=(Critical_72h-Critical_24h)-(Control_72h-Control_24h),levels=design)
fit2 <- contrasts.fit(fit, contrasts)
fit2 <- eBayes(fit2, trend=TRUE)
summary(decideTests(fit2, method="separate"))
   # # Diff24 Diff72
# # -1     39      0
# # 0   29339  29385
# # 1       7      0
Critical_Control_tb <- topTable(fit2, number=dim(fit2)[1], adjust="BH")
write.table(Critical_Control_tb,file="Critical_Control_DF.txt",row.names=F,col.names=T,sep="\t",quote=F)

# > Critical_Control_tb[1:10,]
        # # > Critical_Control_tb[1:10,]
        # # ProbeID TargetID     Diff24      Diff72  AveExpr        F      P.Value    adj.P.Val
# # 4180286 4180286   CRISP3 -1.7381395  0.05193343 4.956846 31.05551 1.473141e-11 4.328825e-07
# # 630014   630014   MAN1A1  1.0066108  0.21027028 6.422275 19.43532 5.081743e-08 7.466351e-04
# # 7510132 7510132     OLR1 -2.4900791  0.60041252 5.621202 17.79798 1.763187e-07 1.295838e-03
# # 2260563 2260563     GJB6 -1.5977640  0.04503530 4.822719 17.61945 2.022590e-07 1.295838e-03

Critical_Control_ovr_24_0_tb <- topTable(fit2, number=dim(fit2)[1], adjust="BH",coef=1)
write.table(Critical_Control_ovr_24_0_tb,file="Critical_Control_ovr_24_0.txt",row.names=F,col.names=T,sep="\t",quote=F)
Critical_Control_ovr_72_24_tb <- topTable(fit2, number=dim(fit2)[1], adjust="BH",coef=2)
write.table(Critical_Control_ovr_72_24_tb,file="Critical_Control_ovr_72_24.txt",row.names=F,col.names=T,sep="\t",quote=F)
length(which(Critical_Control_tb$adj.P.Val < 0.05))
# [1] 32
length(which(Critical_Control_ovr_24_0_tb$adj.P.Val <= 0.05))
# [1] 46

########
## Critical Control at time points ####
contrasts <- makeContrasts(Critical_0h-Control_0h, Critical_24h-Control_24h,Critical_72h-Control_72h,levels=design)
fit2 <- contrasts.fit(fit, contrasts)
fit2 <- eBayes(fit2,trend=TRUE)
summary(decideTests(fit2, method="separate"))
> summary(decideTests(fit2, method="separate"))
   # Critical_0h - Control_0h Critical_24h - Control_24h Critical_72h - Control_72h
# -1                      573                       1540                        328
# 0                     28146                      26643                      28690
# 1                       666                       1202                        367

Critical_Control_TP_tb <- topTable(fit2, number=dim(fit2)[1], adjust="BH")
length(which(Critical_Control_TP_tb$adj.P.Val < 0.05))
write.table(Critical_Control_TP_tb,file="Critical_Control_TP_tb.txt",row.names=F,col.names=T,sep="\t",quote=F)
#[1] 3872
### Critical Controls 0hr
Critical_Control_0h_tb <- topTable(fit2, number=dim(fit2)[1], adjust="BH",coef=1)
write.table(Critical_Control_0h_tb,file="Critical_Control_0h_DF.txt",row.names=F,col.names=T,sep="\t",quote=F)
length(which(Critical_Control_0h_tb$adj.P.Val < 0.05))
## [1] 1239

## Critical Controls 24 hr
Critical_Control_24h_tb <- topTable(fit2, number=dim(fit2)[1], adjust="BH",coef=2)
write.table(Critical_Control_24h_tb,file="Critical_Control_24h_DF.txt",row.names=F,col.names=T,sep="\t",quote=F)
length(which(Critical_Control_24h_tb$adj.P.Val < 0.05))
# [1] 2742

## Critical Controls 72hr
Critical_Control_72h_tb <- topTable(fit2, number=dim(fit2)[1], adjust="BH",coef=3)
write.table(Critical_Control_72h_tb,file="Critical_Control_72h_DF.txt",row.names=F,col.names=T,sep="\t",quote=F)
length(which(Critical_Control_72h_tb$adj.P.Val < 0.05))
# [1] 695
################################
## TEST 
#### Critical vs Control General
# # contrasts <- makeContrasts(Diff24=(Critical_24h-Critical_0h)-(Control_24h-Control_0h),Diff72=(Critical_72h-Critical_24h)-(Control_72h-Control_24h),levels=design)
# # fit2 <- contrasts.fit(fit, contrasts)
# # fit2 <- eBayes(fit2, trend=TRUE)
# # summary(decideTests(fit2, method="global"))
# # # > summary(decideTests(fit2, method="global"))
   # # # Diff24 Diff72
# # # -1     18      0
# # # 0   29364  29384
# # # 1       3      1
# # Critical_Control_tb <- topTable(fit2, number=dim(fit2)[1], adjust="BH")
################################



### Criticals Time points  #######
contrasts <- makeContrasts(Critical_24h-Critical_0h,Critical_72h-Critical_24h,Critical_72h-Critical_0h,levels=design)
fit2 <- contrasts.fit(fit, contrasts)
fit2 <- eBayes(fit2, trend=TRUE)
Critical_TP_tb <- topTable(fit2, number=dim(fit2)[1], adjust="BH")

summary(decideTests(fit2, method="separate"))
   # Critical_24h - Critical_0h Critical_72h - Critical_24h Critical_72h - Critical_0h
# -1                       3382                         685                       3373
# 0                       23091                       28211                      23208
# 1                        2912                         489                       2804

length(which(Critical_TP_tb$adj.P.Val < 0.05))
## 7545

Critical_24h_0h_tb <- topTable(fit2, number=dim(fit2)[1], adjust="BH",coef=1)
write.table(Critical_24h_0h_tb,file="Critical_24h_0h_DF.txt",row.names=F,col.names=T,sep="\t",quote=F)
length(which(Critical_24h_0h_tb$adj.P.Val < 0.05))
## [1]  6294

Critical_72h_24h_tb <- topTable(fit2, number=dim(fit2)[1], adjust="BH",coef=2)
write.table(Critical_72h_24h_tb,file="Critical_72h_24h_DF.txt",row.names=F,col.names=T,sep="\t",quote=F)
length(which(Critical_72h_24h_tb$adj.P.Val < 0.05))
# [1] 1174

Critical_72h_0h_tb <- topTable(fit2, number=dim(fit2)[1], adjust="BH",coef=3)
write.table(Critical_72h_0h_tb,file="Critical_72h_0h_DF.txt",row.names=F,col.names=T,sep="\t",quote=F)
length(which(Critical_72h_0h_tb$adj.P.Val < 0.05))
# [1] 6177

################ ADVERSE vs NO ADVERSE OUTCOME ############################
### no adverse only 9 Samples (NO_MODS_NO_INF) ############################
setwd("C:/Enserio/Barts/Projects/Trauma/limma/BatchNbckCorData_030414/ClusterPlots/")
png("MDS_MODS_030414.png",width=800,height=800)
plotMDS(data,labels=targets$MODS)
dev.off()

setwd("C:/Enserio/Barts/Projects/Trauma/limma/BatchNbckCorData_030414")
## SAMPLES #########
targets$adverse_outcome[which(targets$Covariate3=="MODS_Inf_0h")]<-"Yes_0hr"
targets$adverse_outcome[which(targets$Covariate3=="MODS_Inf_24h")]<-"Yes_24hr"
targets$adverse_outcome[which(targets$Covariate3=="MODS_Inf_72h")]<-"Yes_72hr"
targets$adverse_outcome[which(targets$Covariate3=="MODS_No_Inf_0h")]<-"Yes_0hr"
targets$adverse_outcome[which(targets$Covariate3=="MODS_No_Inf_24h")]<-"Yes_24hr"
targets$adverse_outcome[which(targets$Covariate3=="MODS_No_Inf_72h")]<-"Yes_72hr"
targets$adverse_outcome[which(targets$Covariate3=="No_MODS_Inf_0h")]<-"Yes_0hr"
targets$adverse_outcome[which(targets$Covariate3=="No_MODS_Inf_24h")]<-"Yes_24hr"
targets$adverse_outcome[which(targets$Covariate3=="No_MODS_Inf_72h")]<-"Yes_72hr"
targets$adverse_outcome[which(targets$Covariate3=="No_MODS_No_Inf_0h")]<-"No_0hr"
targets$adverse_outcome[which(targets$Covariate3=="No_MODS_No_Inf_24h")]<-"No_24hr"
targets$adverse_outcome[which(targets$Covariate3=="No_MODS_No_Inf_72h")]<-"No_72hr"
targets$adverse_outcome[which(targets$Covariate1=="Control")]<-"No"
targets$adverse_outcome[which(targets$Covariate1=="Moderate")]<-"No"
> targets$adverse_outcome
  [1] "No"       "No"       "No"       "No_0hr"   "No"       "No_0hr"   "No_0hr"   "Yes_0hr"  "No"       "No_0hr"   "No"       "No"       "Yes_0hr"  "Yes_0hr" 
 [15] "Yes_0hr"  "Yes_0hr"  "Yes_0hr"  "Yes_0hr"  "Yes_0hr"  "No"       "Yes_0hr"  "Yes_0hr"  "Yes_0hr"  "Yes_0hr"  "Yes_0hr"  "Yes_0hr"  "No_0hr"   "No_0hr"  
 [29] "No"       "No_0hr"   "No_0hr"   "Yes_0hr"  "Yes_0hr"  "Yes_0hr"  "Yes_0hr"  "Yes_0hr"  "No"       "No"       "No"       "No"       "No_24hr"  "No"      
 [43] "No_24hr"  "Yes_24hr" "No_24hr"  "Yes_24hr" "No"       "No_24hr"  "No_24hr"  "No"       "No"       "Yes_24hr" "Yes_24hr" "Yes_24hr" "Yes_24hr" "Yes_24hr"
 [57] "Yes_24hr" "Yes_24hr" "No"       "Yes_24hr" "Yes_24hr" "Yes_24hr" "Yes_24hr" "Yes_24hr" "No_24hr"  "No_24hr"  "No"       "No"       "No_24hr"  "Yes_24hr"
 [71] "Yes_24hr" "Yes_24hr" "Yes_24hr" "No"       "No"       "No"       "No"       "No"       "No_72hr"  "No"       "No_72hr"  "Yes_72hr" "No_72hr"  "No"      
 [85] "No_72hr"  "No_72hr"  "No"       "No"       "Yes_72hr" "Yes_72hr" "Yes_72hr" "Yes_72hr" "Yes_72hr" "Yes_72hr" "Yes_72hr" "Yes_72hr" "No"       "Yes_72hr"
 [99] "Yes_72hr" "Yes_72hr" "Yes_72hr" "Yes_72hr" "Yes_72hr" "Yes_72hr" "No_72hr"  "No_72hr"  "No"       "No_72hr"  "Yes_72hr" "Yes_72hr" "Yes_72hr" "Yes_72hr"
[113] "Yes_72hr" "No"       "Yes_72hr" "Yes_72hr" "No"       "RNA_ctrl"

adverse_tags <- factor(targets$adverse_outcome)
levels(adverse_tags)
# [1] "No"       "No_0hr"   "No_24hr"  "No_72hr"  "RNA_ctrl" "Yes_0hr"  "Yes_24hr" "Yes_72hr"
> table(targets$adverse_outcome)

      No   No_0hr  No_24hr  No_72hr RNA_ctrl  Yes_0hr Yes_24hr Yes_72hr 
      33        8        8        8        1       19       18       23 

design <- model.matrix(~0+adverse_tags)
colnames(design) <- levels(adverse_tags)
dupcor<- duplicateCorrelation(data,design,block=targets$Subject,ndups=2)
fit <- lmFit(data,design,block=targets$Subject,correlation=dupcor$consensus.correlation)

#### Adverse General
contrasts <- makeContrasts(Adverse0hr=(Yes_0hr - No_0hr), Adverse24hr=(Yes_24hr - No_24hr), Adverse72hr=(Yes_72hr - No_72hr),levels=design)
fit2 <- contrasts.fit(fit, contrasts)
fit2 <- eBayes(fit2, trend=TRUE)
summary(decideTests(fit2, method="separate"))
# > summary(decideTests(fit2, method="separate"))
   # Adverse0hr Adverse24hr Adverse72hr
# -1          1           6          40
# 0       29381       29374       29308
# 1           3           5          37
Adverse_TP_tb <- topTable(fit2, number=dim(fit2)[1], adjust="BH")
write.table(Adverse_TP_tb,file="Adverse_DF.txt",row.names=F,col.names=T,sep="\t",quote=F)
length(which(Adverse_TP_tb$adj.P.Val < 0.05))
# 95
#### Adverse 0hr 
Adverse_0h_tb <- topTable(fit2, number=dim(fit2)[1], adjust="BH",coef=1)
write.table(Adverse_0h_tb,file="Adverse_0h_DF.txt",row.names=F,col.names=T,sep="\t",quote=F)
length(which(Adverse_0h_tb$adj.P.Val < 0.05))
# [1] 4
####### Adverse 24hr 
Adverse_24h_tb <- topTable(fit2, number=dim(fit2)[1], adjust="BH",coef=2)
write.table(Adverse_24h_tb,file="Adverse_24h_DF.txt",row.names=F,col.names=T,sep="\t",quote=F)
length(which(Adverse_24h_tb$adj.P.Val < 0.05))
# [1] 11
####### Adverse 72hr 
Adverse_72h_tb <- topTable(fit2, number=dim(fit2)[1], adjust="BH",coef=3)
write.table(Adverse_72h_tb,file="Adverse_72h_DF.txt",row.names=F,col.names=T,sep="\t",quote=F)
length(which(Adverse_72h_tb$adj.P.Val < 0.05))
# [1] 77
####

#### Adverse overtime
contrasts <- makeContrasts(Diff24=(Yes_24hr-Yes_0hr)-(No_24hr-No_0hr),Diff72=(Yes_72hr-Yes_24hr)-(No_72hr-No_24hr),levels=design)
fit2 <- contrasts.fit(fit, contrasts)
fit2 <- eBayes(fit2, trend=TRUE)
summary(decideTests(fit2, method="separate"))
> summary(decideTests(fit2, method="separate"))
   Diff24 Diff72
-1      0      0
0   29385  29385
1       0      0
################# MODS vs NO_MODS
targets$MODS[which(targets$Covariate3=="MODS_Inf_0h")]<-"MODS_0hr"
targets$MODS[which(targets$Covariate3=="MODS_Inf_24h")]<-"MODS_24hr"
targets$MODS[which(targets$Covariate3=="MODS_Inf_72h")]<-"MODS_72hr"
targets$MODS[which(targets$Covariate3=="MODS_No_Inf_0h")]<-"MODS_0hr"
targets$MODS[which(targets$Covariate3=="MODS_No_Inf_24h")]<-"MODS_24hr"
targets$MODS[which(targets$Covariate3=="MODS_No_Inf_72h")]<-"MODS_72hr"
targets$MODS[which(targets$Covariate3=="No_MODS_Inf_0h")]<-"NO_MODS_0hr"
targets$MODS[which(targets$Covariate3=="No_MODS_Inf_24h")]<-"NO_MODS_24hr"
targets$MODS[which(targets$Covariate3=="No_MODS_Inf_72h")]<-"NO_MODS_72hr"
targets$MODS[which(targets$Covariate3=="No_MODS_No_Inf_0h")]<-"NO_MODS_0hr"
targets$MODS[which(targets$Covariate3=="No_MODS_No_Inf_24h")]<-"NO_MODS_24hr"
targets$MODS[which(targets$Covariate3=="No_MODS_No_Inf_72h")]<-"NO_MODS_72hr"
targets$MODS[which(targets$Covariate1=="Control")]<-"No"
targets$MODS[which(targets$Covariate1=="Moderate")]<-"No"

mods_tags <- factor(targets$MODS)
# > table(mods_tags)
# mods_tags
    # MODS_0hr    MODS_24hr    MODS_72hr           No  NO_MODS_0hr NO_MODS_24hr NO_MODS_72hr     RNA_ctrl 
          # 11           12           14           33           16           14           17            1 

levels(mods_tags)
# [1] "MODS_0hr"     "MODS_24hr"    "MODS_72hr"    "No"           "NO_MODS_0hr"  "NO_MODS_24hr" "NO_MODS_72hr" "RNA_ctrl"    
design <- model.matrix(~0+mods_tags)
colnames(design) <- levels(mods_tags)
dupcor<- duplicateCorrelation(data,design,block=targets$Subject,ndups=2)
fit <- lmFit(data,design,block=targets$Subject,correlation=dupcor$consensus.correlation)

#### MOD  vs no-MODS
contrasts <- makeContrasts(Mods0hr=(MODS_0hr - NO_MODS_0hr), MODS24hr=(MODS_24hr - NO_MODS_24hr), MODS72hr=(MODS_72hr - NO_MODS_72hr),levels=design)
fit2 <- contrasts.fit(fit, contrasts)
fit2 <- eBayes(fit2, trend=TRUE)
summary(decideTests(fit2, method="global"))
> summary(decideTests(fit2, method="separate"))
   Mods0hr MODS24hr MODS72hr
-1     214        5        5
0    29022    29352    29357
1      149       28       23

MODs_TP_tb <- topTable(fit2, number=dim(fit2)[1], adjust="BH")
length(which(MODs_TP_tb$adj.P.Val < 0.05))
write.table(MODs_TP_tb,file="Mods_DF.txt",row.names=F,col.names=T,sep="\t",quote=F)
# 360

#### Mods 0hr 
Mods_0h_tb <- topTable(fit2, number=dim(fit2)[1], adjust="BH",coef=1)
write.table(Mods_0h_tb,file="Mods_0h_DF.txt",row.names=F,col.names=T,sep="\t",quote=F)
length(which(Mods_0h_tb$adj.P.Val < 0.05))
# [1] 363
#### MODS 24 hr
Mods_24h_tb <- topTable(fit2, number=dim(fit2)[1], adjust="BH",coef=2)
write.table(Mods_24h_tb,file="Mods_24h_DF.txt",row.names=F,col.names=T,sep="\t",quote=F)
length(which(Mods_24h_tb$adj.P.Val < 0.05))
# [1] 33
##### MODS 72hr #########
Mods_72h_tb <- topTable(fit2, number=dim(fit2)[1], adjust="BH",coef=3)
write.table(Mods_72h_tb,file="Mods_72h_DF.txt",row.names=F,col.names=T,sep="\t",quote=F)
length(which(Mods_72h_tb$adj.P.Val < 0.05))
# 28

### MODS over time
contrasts <- makeContrasts(Diff24=(MODS_24hr-MODS_0hr)-(NO_MODS_24hr-NO_MODS_0hr),Diff72=(MODS_72hr-MODS_24hr)-(NO_MODS_72hr-NO_MODS_24hr),levels=design)
fit2 <- contrasts.fit(fit, contrasts)
fit2 <- eBayes(fit2, trend=TRUE)
summary(decideTests(fit2, method="global"))
> summary(decideTests(fit2, method="separate"))
   Diff24 Diff72
-1     44      0
0   29315  29385
1      26      0
MODs_overtime_tb <- topTable(fit2, number=dim(fit2)[1], adjust="BH")
# > length(which(MODs_overtime_tb$adj.P.Val < 0.05))
# [1] 63
write.table(MODs_overtime_tb,file="Mods_overtime.txt",row.names=F,col.names=T,sep="\t",quote=F)

#### Mods overtime diff at 24-0
Mods_ovr_24_0_tb <- topTable(fit2, number=dim(fit2)[1], adjust="BH",coef=1)
write.table(Mods_ovr_24_0_tb ,file="Mods_ovr_24_0_DF.txt",row.names=F,col.names=T,sep="\t",quote=F)
length(which(Mods_ovr_24_0_tb$adj.P.Val < 0.05))
# [1] 70

#### Mods overtime diff at 72-24
Mods_ovr_72_24_tb <- topTable(fit2, number=dim(fit2)[1], adjust="BH",coef=2)
write.table(Mods_ovr_72_24_tb ,file="Mods_ovr_72_24_DF.txt",row.names=F,col.names=T,sep="\t",quote=F)
> length(which(Mods_ovr_72_24_tb$adj.P.Val < 0.05))
[1] 0
#################### INFECTION VS NO INFECTION ###################
targets$Infection[which(targets$Covariate3=="MODS_Inf_0h")]<-"Inf_0hr"
targets$Infection[which(targets$Covariate3=="MODS_Inf_24h")]<-"Inf_24hr"
targets$Infection[which(targets$Covariate3=="MODS_Inf_72h")]<-"Inf_72hr"
targets$Infection[which(targets$Covariate3=="MODS_No_Inf_0h")]<-"NoInf_0hr"
targets$Infection[which(targets$Covariate3=="MODS_No_Inf_24h")]<-"NoInf_24hr"
targets$Infection[which(targets$Covariate3=="MODS_No_Inf_72h")]<-"NoInf_72hr"
targets$Infection[which(targets$Covariate3=="No_MODS_Inf_0h")]<-"Inf_0hr"
targets$Infection[which(targets$Covariate3=="No_MODS_Inf_24h")]<-"Inf_24hr"
targets$Infection[which(targets$Covariate3=="No_MODS_Inf_72h")]<-"Inf_72hr"
targets$Infection[which(targets$Covariate3=="No_MODS_No_Inf_0h")]<-"NoInf_0hr"
targets$Infection[which(targets$Covariate3=="No_MODS_No_Inf_24h")]<-"NoInf_24hr"
targets$Infection[which(targets$Covariate3=="No_MODS_No_Inf_72h")]<-"NoInf_72hr"
targets$Infection[which(targets$Covariate1=="Control")]<-"No"
targets$Infection[which(targets$Covariate1=="Moderate")]<-"No"
inf_tags <- factor(targets$Infection)
# > table(targets$Infection)

   # Inf_0hr   Inf_24hr   Inf_72hr         No  NoInf_0hr NoInf_24hr NoInf_72hr   RNA_ctrl 
        # 17         16         21         33         10         10         10          1 
levels(inf_tags)
# [1] "Inf_0hr"    "Inf_24hr"   "Inf_72hr"   "No"         "NoInf_0hr"  "NoInf_24hr" "NoInf_72hr" "RNA_ctrl"  

design <- model.matrix(~0+inf_tags)
colnames(design) <- levels(inf_tags)
dupcor<- duplicateCorrelation(data,design,block=targets$Subject,ndups=2)
fit <- lmFit(data,design,block=targets$Subject,correlation=dupcor$consensus.correlation)

#### INFECTION VS NO INFECTIOn
contrasts <- makeContrasts(Inf0hr=(Inf_0hr - NoInf_0hr), Inf24hr=(Inf_24hr - NoInf_24hr), Inf72hr=(Inf_72hr - NoInf_72hr),levels=design)
fit2 <- contrasts.fit(fit, contrasts)
fit2 <- eBayes(fit2, trend=TRUE)
summary(decideTests(fit2, method="global"))
> summary(decideTests(fit2, method="separate"))
   Inf0hr Inf24hr Inf72hr
-1      0      10      38
0   29383   29363   29309
1       2      12      38
Inf_TP_tb <- topTable(fit2, number=dim(fit2)[1], adjust="BH")
write.table(Inf_TP_tb,file="Inf_TP_DF.txt",row.names=F,col.names=T,sep="\t",quote=F)
length(which(Inf_TP_tb$adj.P.Val < 0.05))
# 110

#### Inf 0hr 
Inf_0h_tb <- topTable(fit2, number=dim(fit2)[1], adjust="BH",coef=1)
write.table(Inf_0h_tb,file="Inf_0h_DF.txt",row.names=F,col.names=T,sep="\t",quote=F)
length(which(Inf_0h_tb$adj.P.Val < 0.05))
# [1] 2

#### INF 24 hr
Inf_24h_tb <- topTable(fit2, number=dim(fit2)[1], adjust="BH",coef=2)
write.table(Inf_24h_tb,file="Inf_24h_DF.txt",row.names=F,col.names=T,sep="\t",quote=F)
length(which(Inf_24h_tb$adj.P.Val < 0.05))
# [1] 22
##### INF 72hr #########
Inf_72h_tb <- topTable(fit2, number=dim(fit2)[1], adjust="BH",coef=3)
write.table(Inf_72h_tb,file="Inf_72h_DF.txt",row.names=F,col.names=T,sep="\t",quote=F)
length(which(Inf_72h_tb$adj.P.Val < 0.05))
# 76
###############################
contrasts <- makeContrasts(Diff24=(Inf_24hr-Inf_0hr)-(NoInf_24hr-NoInf_0hr),Diff72=(Inf_72hr-Inf_24hr)-(NoInf_72hr-NoInf_24hr),levels=design)
fit2 <- contrasts.fit(fit, contrasts)
fit2 <- eBayes(fit2, trend=TRUE)
# # > summary(decideTests(fit2, method="separate"))
   # # Diff24 Diff72
# # -1      0      0
# # 0   29385  29385
# # 1       0      0
Inf_overtime_tb <- topTable(fit2, number=dim(fit2)[1], adjust="BH")
# > length(which(Inf_overtime_tb$adj.P.Val < 0.05))
# [1] 1
write.table(Inf_overtime_tb,file="Inf_overtime.txt",row.names=F,col.names=T,sep="\t",quote=F)

Inf_ovr_24_0_tb <- topTable(fit2, number=dim(fit2)[1], adjust="BH",coef=1)
length(which(Inf_ovr_24_0_tb$adj.P.Val < 0.05))
# 0
Inf_ovr_72_24_tb <- topTable(fit2, number=dim(fit2)[1], adjust="BH",coef=2)
length(which(Inf_ovr_72_24_tb$adj.P.Val < 0.05))
# 0

##############################################################################################################################

#### Survivor vs NOSurvivor
# names(targets)
# [29] "28_days_outcome_Dead_0_Alive_1"
targets[which(targets[,29]==0),29] <- "D"
targets[which(targets[,29]==1),29] <- "S"
x <- paste(targets[,29],targets$Covariate2,sep="_")
targets[,29] <- x
suv_tags <- factor(targets[,29])
levels(suv_tags)
# [1] "D_T0h"       "D_T24h"      "D_T72h"      "NA_RNA_ctrl" "S_T0h"       "S_T24h"      "S_T72h"     
# # > table(targets[,29])

      # # D_T0h      D_T24h      D_T72h NA_RNA_ctrl       S_T0h      S_T24h      S_T72h 
          # # 4           4           6           1          33          34          36 

design <- model.matrix(~0+suv_tags)
colnames(design) <- levels(suv_tags)
dupcor<- duplicateCorrelation(data,design,block=targets$Subject)
fit <- lmFit(data,design,block=targets$Subject,correlation=dupcor$consensus.correlation)

contrasts <- makeContrasts(D_T0h-S_T0h,D_T24h-S_T24h,D_T72h-S_T72h,levels=design)
fit2 <- contrasts.fit(fit, contrasts)
fit2 <- eBayes(fit2, trend=TRUE)
# > summary(decideTests(fit2, method="separate"))
   D_T0h-S_T0h D_T24h-S_T24h D_T72h-S_T72h
-1            12              14             146
0          29349           29309           28843
1             24              62             396

Suv_TP_tb <- topTable(fit2, number=dim(fit2)[1], adjust="BH")
length(which(Suv_TP_tb$adj.P.Val < 0.05))
# [1] 656
write.table(Suv_TP_tb,file="Suv_DF.txt",row.names=F,col.names=T,sep="\t",quote=F)

Suv_0h_tb <- topTable(fit2, number=dim(fit2)[1], adjust="BH",coef=1)
write.table(Suv_0h_tb,file="Suv_0h_DF.txt",row.names=F,col.names=T,sep="\t",quote=F)
length(which(Suv_0h_tb$adj.P.Val < 0.05))
#[1] 36
Suv_24h_tb <- topTable(fit2, number=dim(fit2)[1], adjust="BH",coef=2)
write.table(Suv_24h_tb,file="Suv_24h_DF.txt",row.names=F,col.names=T,sep="\t",quote=F)
length(which(Suv_24h_tb$adj.P.Val < 0.05))
#[1] 76
Suv_72h_tb <- topTable(fit2, number=dim(fit2)[1], adjust="BH",coef=3)
write.table(Suv_72h_tb,file="Suv_72h_DF.txt",row.names=F,col.names=T,sep="\t",quote=F)
length(which(Suv_72h_tb$adj.P.Val < 0.05))
#[1] 542

#################################### Prolonged mods ############################################################
##define MODS
targets$MODS[which(targets$Covariate3=="MODS_Inf_0h")]<-"MODS_0hr"
targets$MODS[which(targets$Covariate3=="MODS_Inf_24h")]<-"MODS_24hr"
targets$MODS[which(targets$Covariate3=="MODS_Inf_72h")]<-"MODS_72hr"
targets$MODS[which(targets$Covariate3=="MODS_No_Inf_0h")]<-"MODS_0hr"
targets$MODS[which(targets$Covariate3=="MODS_No_Inf_24h")]<-"MODS_24hr"
targets$MODS[which(targets$Covariate3=="MODS_No_Inf_72h")]<-"MODS_72hr"
targets$MODS[which(targets$Covariate3=="No_MODS_Inf_0h")]<-"NO_MODS_0hr"
targets$MODS[which(targets$Covariate3=="No_MODS_Inf_24h")]<-"NO_MODS_24hr"
targets$MODS[which(targets$Covariate3=="No_MODS_Inf_72h")]<-"NO_MODS_72hr"
targets$MODS[which(targets$Covariate3=="No_MODS_No_Inf_0h")]<-"NO_MODS_0hr"
targets$MODS[which(targets$Covariate3=="No_MODS_No_Inf_24h")]<-"NO_MODS_24hr"
targets$MODS[which(targets$Covariate3=="No_MODS_No_Inf_72h")]<-"NO_MODS_72hr"
targets$MODS[which(targets$Covariate1=="Control")]<-"No"
targets$MODS[which(targets$Covariate1=="Moderate")]<-"No"
################### excluding 303 ############################################################
prmods <- targets$MODS
> prmods[which(targets$Subject==314)]
[1] "MODS_0hr"  "MODS_24hr" "MODS_72hr"
> prmods[which(targets$Subject==327)]
[1] "MODS_0hr"  "MODS_24hr" "MODS_72hr"
> prmods[which(targets$Subject==429)]
[1] "MODS_0hr"  "MODS_24hr" "MODS_72hr"
> prmods[which(targets$Subject==443)]
[1] "MODS_0hr"  "MODS_24hr" "MODS_72hr" "MODS_72hr"
prmods[which(targets$Subject==314)]<- c("PR_MODS_0hr","PR_MODS_24hr","PR_MODS_72hr")
prmods[which(targets$Subject==327)]<- c("PR_MODS_0hr","PR_MODS_24hr","PR_MODS_72hr")
prmods[which(targets$Subject==429)]<- c("PR_MODS_0hr","PR_MODS_24hr","PR_MODS_72hr")
prmods[which(targets$Subject==443)]<- c("PR_MODS_0hr","PR_MODS_24hr","PR_MODS_72hr","PR_MODS_72hr")

prmods_tags <- factor(prmods)
levels(prmods_tags)
 # [1] "MODS_0hr"     "MODS_24hr"    "MODS_72hr"    "No"           "NO_MODS_0hr"  "NO_MODS_24hr" "NO_MODS_72hr" "PR_MODS_0hr"  "PR_MODS_24hr" "PR_MODS_72hr"
# [11] "RNA_ctrl"    
# # > table(prmods)
# # prmods
    # # MODS_0hr    MODS_24hr    MODS_72hr           No  NO_MODS_0hr NO_MODS_24hr NO_MODS_72hr  PR_MODS_0hr PR_MODS_24hr PR_MODS_72hr     RNA_ctrl 
           # # 7            8            9           33           16           14           17            4            4            5            1 
design <- model.matrix(~0+prmods_tags)
colnames(design) <- levels(prmods_tags)
dupcor<- duplicateCorrelation(data,design,block=targets$Subject,ndups=2)
fit <- lmFit(data,design,block=targets$Subject,correlation=dupcor$consensus.correlation)
#### MOD  vs PR-MODS
contrasts <- makeContrasts(PRMods0hr=(PR_MODS_0hr - MODS_0hr), PRMODS24hr=(PR_MODS_24hr - MODS_24hr), PRMODS72hr=(PR_MODS_72hr - MODS_72hr),levels=design)
fit2 <- contrasts.fit(fit, contrasts)
fit2 <- eBayes(fit2, trend=TRUE)
summary(decideTests(fit2, method="separate"))
   PRMods0hr PRMODS24hr PRMODS72hr
-1        31          1        128
0      29312      29378      29162
1         42          6         95


PR_MODs_TPno303_tb <- topTable(fit2, number=dim(fit2)[1], adjust="BH")
length(which(PR_MODs_TPno303_tb$adj.P.Val < 0.05))
#154
write.table(PR_MODs_TPno303_tb,file="PR_Mods_DF_no303.txt",row.names=F,col.names=T,sep="\t",quote=F)

#### PR_Mods 0hr 
PR_Modsno303_0h_tb <- topTable(fit2, number=dim(fit2)[1], adjust="BH",coef=1)
write.table(PR_Modsno303_0h_tb,file="PR_Mods_0h_DF_no303.txt",row.names=F,col.names=T,sep="\t",quote=F)
length(which(PR_Modsno303_0h_tb$adj.P.Val < 0.05))
# [1] 73

####  PR MODS 24 hr
PR_Modsno303_24h_tb <- topTable(fit2, number=dim(fit2)[1], adjust="BH",coef=2)
write.table(PR_Modsno303_24h_tb,file="PR_Mods_24h_DF_no303.txt",row.names=F,col.names=T,sep="\t",quote=F)
length(which(PR_Modsno303_24h_tb$adj.P.Val < 0.05))
[1] 7
##### PR MODS 72hr #########
PR_Modsno303_72h_tb <- topTable(fit2, number=dim(fit2)[1], adjust="BH",coef=3)
write.table(PR_Modsno303_72h_tb,file="PR_Mods_72h_DF_no303.txt",row.names=F,col.names=T,sep="\t",quote=F)
length(which(PR_Modsno303_72h_tb$adj.P.Val < 0.05))
#223
!! SAVE

######################################################################################################
### Compare genes in PR_MODS vs MODS against genes in MODS vs NO_MODS 
prmods_tags <- factor(prmods)
levels(prmods_tags)
 # [1] "MODS_0hr"     "MODS_24hr"    "MODS_72hr"    "No"           "NO_MODS_0hr"  "NO_MODS_24hr" "NO_MODS_72hr" "PR_MODS_0hr"  "PR_MODS_24hr" "PR_MODS_72hr"
# [11] "RNA_ctrl"  
design <- model.matrix(~0+prmods_tags)
colnames(design) <- levels(prmods_tags)
dupcor<- duplicateCorrelation(data,design,block=targets$Subject,ndups=2)
fit <- lmFit(data,design,block=targets$Subject,correlation=dupcor$consensus.correlation)

contrasts <- makeContrasts(Mods0hr=(MODS_0hr - NO_MODS_0hr), MODS24hr=(MODS_24hr - NO_MODS_24hr), MODS72hr=(MODS_72hr - NO_MODS_72hr),levels=design)
fit2 <- contrasts.fit(fit, contrasts)
fit2 <- eBayes(fit2, trend=TRUE)
# # > summary(decideTests(fit2, method="separate"))
# # Mods0hr MODS24hr MODS72hr
# # -1      16        9        0
# # 0    29333    29335    29375
# # 1       36       41       10
> summary(decideTests(fit2, method="separate"))
   Mods0hr MODS24hr MODS72hr
-1      11        0        0
0    29349    29385    29385
1       25        0        0

### ! different results from the previous comparison of MODs and no MODs because some samples now are categorized as Prolonged MODs
MODs_NOMods_withPRMODstb <- topTable(fit2, number=dim(fit2)[1], adjust="BH")
Mods_NOMods_withPRMODs_0h_tb <- topTable(fit2, number=dim(fit2)[1], adjust="BH",coef=1)
Mods_NOMods_withPRMODs_24h_tb <- topTable(fit2, number=dim(fit2)[1], adjust="BH",coef=2)
Mods_NOMods_withPRMODs_72h_tb <- topTable(fit2, number=dim(fit2)[1], adjust="BH",coef=3)

x1 <- MODs_NOMods_withPRMODstb[which(MODs_NOMods_withPRMODstb$adj.P.Val < 0.05),2]
x2 <- PR_MODs_TP_tb[which(PR_MODs_TP_tb$adj.P.Val < 0.05),2]

library(VennDiagram)
## prmods(with303) vs mods (old mods)
venn.diagram(x=list(MODs=x1,PR_MODs=x2),lwd=3,filename="PRMODSno303_MODSold.tiff",col="black",fill=c("cornflowerblue", "green"),cex = 2.5)
## prmods(with 303 vs Mods (new mods)
venn.diagram(x=list(MODs=x1,PR_MODs=x2),lwd=3,filename="PRMODSno303_MODS.tiff",col="black",fill=c("cornflowerblue", "green"),cex = 2.5)
# prmods(no 303) vs mods(old mods)
venn.diagram(x=list(MODs=x1,PR_MODs=x2),lwd=3,filename="PRMODS_MODS.tiff",col="black",fill=c("cornflowerblue", "green"),cex = 2.5)
# prmods(no 303) vs mods(new mods)
venn.diagram(x=list(MODs=x1,PR_MODs=x2),lwd=3,filename="PRMODS_MODS.tiff",col="black",fill=c("cornflowerblue", "green"),cex = 2.5)



####### compare with MODS vs NO_MODs (when only two categories) ##################























######################### TABLE for GENEGO ############################################################

> ls(pattern="_tb")

MODS: "Mods_ovr_24_0_tb","Mods_ovr_72_24_tb","MODs_TP_tb","Mods_0h_tb","Mods_24h_tb","Mods_72h_tb","MODs_overtime_tb"

ADVERSE: "Adverse_0h_tb","Adverse_24h_tb","Adverse_72h_tb","Adverse_TP_tb",
CRITICAL:,"Critical_24h_0h_tb","Critical_72h_0h_tb","Critical_72h_24h_tb","Critical_Control_0h_tb","Critical_Control_24h_tb","Critical_Control_72h_tb", "Critical_Control_ovr_24_0_tb","Critical_Control_ovr_72_24_tb","Critical_Control_tb","Critical_Control_TP_tb","Critical_TP_tb"
INFECTION:,"Inf_0h_tb","Inf_24h_tb","Inf_72h_tb","Inf_overtime_tb","Inf_ovr_24_0_tb","Inf_ovr_72_24_tb","Inf_TP_tb",
SURVIVOR:,"Suv_0h_tb","Suv_24h_tb","Suv_72h_tb","Suv_overtime_tb","Suv_ovr_24_0_tb","Suv_ovr_72_24_tb","Suv_TP_tb",

#### ORDER ALL DATAFRAMES
x<-ls(pattern="_tb")
for(i in x){
	temp <- get(i)
	y <- paste("t",i,sep="_")
	z <- paste("Adj.PVal",i,sep="_")
	colnames(temp)[5] <- y
	colnames(temp)[7] <- z
	assign(i,temp[order(temp[,1]),])
}

trauma_exp <- data.frame(MODs_overtime_tb[,c(1,2,5,7)],Mods_ovr_24_0_tb[,c(5,7)],Mods_0h_tb[,c(5,7)],Mods_24h_tb[,c(5,7)],Mods_72h_tb[,c(5,7)],Adverse_0h_tb[,c(5,7)],Adverse_24h_tb[,c(5,7)],Inf_ovr_24_0_tb[,c(5,7)],Inf_0h_tb[,c(5,7)],Inf_24h_tb[,c(5,7)],Inf_72h_tb[,c(5,7)],check.names=TRUE)
write.table(trauma_exp,file="Exp.txt",sep="\t",col.names=T,row.names=F,quote=F)

trauma_vals <- data.frame(MODs_overtime_tb[,c(1,2,7)],Mods_ovr_24_0_tb[,7],Mods_0h_tb[,7],Mods_24h_tb[,7],Mods_72h_tb[,7],Adverse_0h_tb[,7],Adverse_24h_tb[,7],Inf_ovr_24_0_tb[,7],Inf_0h_tb[,7],Inf_24h_tb[,7],Inf_72h_tb[,7],check.names=TRUE)
write.table(trauma_vals,file="Vals.txt",sep="\t",col.names=T,row.names=F,quote=F)



x <- order(trauma_exp$Adj.PVal_MODs_overtime_tb,trauma_exp$Adj.PVal_Mods_0h_tb)
trauma_exp<-trauma_exp[x,]
 




 
 ### CELLMIX
 ######
inj_tab <- unique(targets[,c(1,15,32:55)])
inj_tab[1:3,]
write.table(inj_tab, file="targets_tmpoints.txt",col.names=T,row.names=F,sep="\t")
 
 
 
 
 
 
 
 
 
 
 
 
 