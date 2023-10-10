library("survival")
library("ggplot2")
library("survminer")
library("ggforce")
library("ggfortify")
library("ggpubr")
library(readxl)


df = read_excel("F:/ferroptosis/TCGA/survival_data/Survival_Table_RKIP_BACH1.xlsx",sheet = 2)
#Group the expression data of RKIP and BACH1 such that if the expression value is greater than the median, its R+/B+, else R-/B- respectively.
x1=df %>%group_by(`cancer type abbreviation`) %>%mutate(RKIP = ifelse(RKIP>median(RKIP),"R+",ifelse(RKIP<median(RKIP),"R-",NA)))
x1=x1 %>%group_by(`cancer type abbreviation`) %>%mutate(BACH1 = ifelse(BACH1>median(BACH1),"B+",ifelse(BACH1<median(BACH1),"B-",NA)))
x1<-as.data.frame(x1)
x1$Type<-paste(x1[,"RKIP"], x1[,"BACH1"])  #Group the RKIP and BACH1 columns to make different combinations.
cancer<-c("ACC","BLCA","BRCA","CESC","CHOL","COAD","DLBC","ESCA","GBM","HNSC","KICH","KIRC","KIRP","LAML","LGG","LIHC","LUAD","LUSC","MESO","OV","PAAD","PCPG","PRAD","READ","SARC","SKCM","STAD","TGCT","THCA","THYM","UCEC","UCS","UVM")

for(i in 1:length(cancer))
{ 
  df<-x1[x1[,c("cancer type abbreviation")]==cancer[i],]
  y = df[,c("sample", "DFI", "DFI.time", "Type")]  #change OS to DSS,PFI or DFI if needed
  y$Type = factor(y$Type, levels = c( "R- B+", "R+ B+", "R+ B-", "R- B-"))
# y<-na.omit(y)
  fit.coxph <- coxph(Surv(DFI.time, DFI) ~ Type, data = y)
  s<-summary(fit.coxph)
  p.value<-signif(s$coefficients[, "Pr(>|z|)"],2)
  HR <-signif(s$coefficients[,"exp(coef)"], 2)
  log2HR<-log2(HR)
  HR.confint.lower <- signif(s$conf.int[,"lower .95"], 2)
  HR.confint.upper <- signif(s$conf.int[,"upper .95"],2)
  Log_Rank_pvalue = signif(s$logtest[['pvalue']],2)
  table<-data.frame(HR,HR.confint.lower,HR.confint.upper,p.value)
  labels<-substr(row.names(table),nchar(row.names(table))-4,nchar(row.names(table)))
  table$labels<-labels
  
  ggplot(data=table, aes(labels, HR)) +
    geom_pointrange(aes(ymin=HR.confint.lower, ymax=HR.confint.upper,colour=ifelse(HR>1, "red", "blue")), position = position_dodge(0.5), shape=20, size= 1.2) + 
    geom_errorbar(aes(ymin=HR.confint.lower, ymax=HR.confint.upper,colour=ifelse(HR>1, "red", "blue")), size = 2, width=0.2) +
    theme_classic()+geom_text(aes(label = ifelse(p.value<0.05 & p.value>=0.01, "*", ifelse(p.value<0.01 & p.value>0.001, "**", ifelse(p.value<0.001, "***", ""))), y= ifelse(HR.confint.lower<=0 & HR.confint.upper<=0, HR.confint.lower-0.5, HR.confint.upper+0.5)), 
                              position = position_dodge(width = 1.3), vjust= 0.7, size = 25 / .pt)+xlab(" ")+ylab("Hazard Ratio")+
    geom_hline(yintercept=1, lty=3, size = 1, colour= 'black') + annotate(geom="text", x=3.4, y=3.5, label=paste0("Reference : R- B+ \n Log-Rank p: ",Log_Rank_pvalue), size=8,fontface="bold")+
    theme(legend.position="none",axis.ticks = element_line(colour = "black", size = 0.7,), axis.ticks.length=unit(.07, "cm"), axis.text.y =element_text(size=20, face = "bold"), axis.text.x =element_text(size=20,colour = "black", face = "bold"),axis.title.x = element_text(size = 22, face="bold"))+geom_hline(yintercept=1, lty=3, size = 0.88, colour= 'black') + coord_flip()
  

  ggsave(file = paste0("path/",cancer[i],".pdf"))
}

dev.off()