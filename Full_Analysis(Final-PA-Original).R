setwd("C://Users/Angus/Downloads/MasterThesis/EPI Revision/ANALYSIS/Original/")
source("C:/Users/Angus/Downloads/MasterThesis/EPI Revision/ANALYSIS/DataFunctions(Final).R")
source("C:/Users/Angus/Downloads/MasterThesis/EPI Revision/ANALYSIS/EstFunctions(Final).R")


library('ltmle')
library('SuperLearner')
library("tidyverse")
library("xtable")
library("patchwork")
library("data.table")

##########################################21 Days(aka. PRIMARY ANALSYSIS)#########################################
filepath <- "C://Users/Angus/Downloads/MasterThesis/EPI Revision/DATA/CaseAsOutcome/PrimaryAnalysis/COVID-19_Dataset_primaryAnalysis7.csv"
endDate <- "2020-12-01" 
this.date <- endDate
FILE <- "Race"
mobility <- T
pre_policy <- T
relativeGrowth <- F
DC <- F
NY <- T
screen=T
SL.library <- list( "SL.mean", "SL.gam", "SL.earth" , "SL.rpart", "SL.xgboost")
A = "A_lv3_targetDateplus0"
Y = "relativeChange_7days"


filepath <- "C://Users/Angus/Downloads/MasterThesis/EPI Revision/DATA/CaseAsOutcome/PrimaryAnalysis/COVID-19_Dataset_primaryAnalysis21.csv"
Y = "relativeChange_21days"

primary_analysis_21days <- RUN(filepath = filepath, this.date=this.date, mobility=mobility, 
                               pre_policy=pre_policy,DC=DC, NY=NY, screen=screen, relativeGrowth=relativeGrowth,
                               SL.library = SL.library, A=A, Y=Y)

##########################################30 Days(aka. PRIMARY ANALSYSIS)#########################################
filepath <- "C://Users/Angus/Downloads/MasterThesis/EPI Revision/DATA/CaseAsOutcome/PrimaryAnalysis/COVID-19_Dataset_primaryAnalysis30.csv"
Y = "relativeChange_30days"

primary_analysis_30days <- RUN(filepath = filepath, this.date=this.date, mobility=mobility, 
                               pre_policy=pre_policy,DC=DC, NY=NY, screen=screen, relativeGrowth=relativeGrowth,
                               SL.library = SL.library, A=A, Y=Y)

##########################################45 Days(aka. PRIMARY ANALSYSIS)#########################################
filepath <- "C://Users/Angus/Downloads/MasterThesis/EPI Revision/DATA/CaseAsOutcome/PrimaryAnalysis/COVID-19_Dataset_primaryAnalysis45.csv"
Y = "relativeChange_45days"

primary_analysis_45days <- RUN(filepath = filepath, this.date=this.date, mobility=mobility, 
                               pre_policy=pre_policy,DC=DC, NY=NY, screen=screen, relativeGrowth=relativeGrowth,
                               SL.library = SL.library, A=A, Y=Y)

##########################################60 Days(aka. PRIMARY ANALSYSIS)#########################################
filepath <- "C://Users/Angus/Downloads/MasterThesis/EPI Revision/DATA/CaseAsOutcome/PrimaryAnalysis/COVID-19_Dataset_primaryAnalysis60.csv"
Y = "relativeChange_60days"

primary_analysis_60days <- RUN(filepath = filepath, this.date=this.date, mobility=mobility, 
                               pre_policy=pre_policy,DC=DC, NY=NY, screen=screen, relativeGrowth=relativeGrowth,
                               SL.library = SL.library, A=A, Y=Y)


############################################Compile Results into table############################################
unadjusted_PA_result <- as.data.frame(rbind(primary_analysis_21days$Result_Unadj_21days,
                                                  primary_analysis_30days$Result_Unadj_30days,
                                                  primary_analysis_45days$Result_Unadj_45days,
                                                  primary_analysis_60days$Result_Unadj_60days))%>%
  mutate(Adjustment = "Unadjusted Analysis") %>%
  mutate("Outcome at" = c(21,30,45,60))

primary_analysis_result <- as.data.frame(rbind(primary_analysis_21days$Result_TMLE_21days,
                                               primary_analysis_30days$Result_TMLE_30days,
                                               primary_analysis_45days$Result_TMLE_45days,
                                               primary_analysis_60days$Result_TMLE_60days))%>%
  mutate(Adjustment = "Primary Analysis") %>%
  mutate("Outcome at" = c(21,30,45,60))

result_PA_all <- as.data.frame(rbind(primary_analysis_result,unadjusted_PA_result)) %>% 
  mutate("Txt(95% CI)" = paste0(round(txt.pt,2)," (", round(txt.CI.lo,2),", ", round(txt.CI.hi,2), ")")) %>%
  mutate("Con(95% CI)" = paste0(round(con.pt,2)," (", round(con.CI.lo,2),", ", round(con.CI.hi,2), ")")) %>%
  mutate("RD(95% CI)" = paste0(round(RD.pt,2)," (", round(RD.CI.lo,2),", ", round(RD.CI.hi,2), ")")) %>%
  mutate("RD p-value" = format(signif(RD.pval, digits=2)))%>%
  mutate("RR(95% CI)" = paste0(round(RR.pt,2)," (", round(RR.CI.lo,2),", ", round(RR.CI.hi,2), ")")) %>%
  mutate("RR p-value" = format(signif(RR.pval, digits=2)))

##RR Table
result_PA_RR <- result_PA_all %>%
  select("Adjustment","Outcome at", "Txt(95% CI)", "Con(95% CI)", "RR(95% CI)", "RR p-value") %>%
  mutate("Outcome at"= paste(as.character(`Outcome at`),"days"))
#write.csv(result_PA_RR, paste0("(FINAL)result-RR(PA))",'.csv'),row.names = F)

##RD Table
result_PA_RD <- result_PA_all %>%
  select("Adjustment","Outcome at", "Txt(95% CI)", "Con(95% CI)", "RD(95% CI)", "RD p-value") %>%
  mutate("Outcome at"= paste(as.character(`Outcome at`),"days"))
#write.csv(result_PA_RD, paste0("(FINAL)result-RD(PA))",'.csv'),row.names = F)

## Write result to csv for reference and create tables/figures
write.csv(result_PA_all, paste0("result-PA(Original)",'.csv'),row.names = F)


##############################################Supplementary Material##############################################
##Supplmentary Table 1
supTable1 <- result_PA_all %>% 
  select("Adjustment","Outcome at", "Txt(95% CI)","Con(95% CI)","RR(95% CI)","RD(95% CI)")
write.csv(supTable1, paste0("Supplementary Table 1 (PA-Original)",'.csv'),row.names = F)


##Supplmentary Figure 1 (PScore Histogram)
PScore <- as.data.frame(cbind(primary_analysis_21days$Propensity_Score$pscore_TMLE,
                              primary_analysis_30days$Propensity_Score$pscore_TMLE,
                              primary_analysis_45days$Propensity_Score$pscore_TMLE,
                              primary_analysis_60days$Propensity_Score$pscore_TMLE))
colnames(PScore) <- c("21 days","30 days","45 days","60 days")

PScore <- PScore %>% gather(Days, "PropensityScore")
PScore$Days <- factor(PScore$Days, levels=c("21 days","30 days","45 days","60 days"))


supFig1 <- ggplot(data=PScore, aes(x=PropensityScore)) +
  geom_histogram(bins=50)+
  facet_grid(Days ~ .)+
  scale_x_continuous(breaks=seq(0,1,0.1))+
  scale_y_continuous(breaks=seq(0,12,2))

tiff("PScore(PA-Original).tiff", units="in", width=7.5, height= 6, res=1200)
supFig1
dev.off()

setDT(PScore)
PScore_table <- PScore[, as.list(summary(PropensityScore)), by = Days]
write.csv(PScore_table, paste0("PropensityScoreSummary(PA-Original)",'.csv'),row.names = F)
##################################################################################################################



