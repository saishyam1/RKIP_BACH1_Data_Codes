library("survival")
library("ggplot2")
library("survminer")
library("ggforce")
library("ggfortify")
library("ggpubr")



df = read.table(path_to_file_with_survival_data)
#Group the expression data of RKIP and BACH1 such that if the expression value is greater than the median, its R+/B+, else R-/B- respectively.
x1=x %>%group_by(`cancer type abbreviation`) %>%mutate(RKIP = ifelse(RKIP>median(RKIP),"R+",ifelse(RKIP<median(RKIP),"R-",NA)))
x1=x1 %>%group_by(`cancer type abbreviation`) %>%mutate(BACH1 = ifelse(BACH1>median(BACH1),"B+",ifelse(BACH1<median(BACH1),"B-",NA)))
x1<-as.data.frame(x1)
x1$Type<-paste(x1[,"RKIP"], x1[,"BACH1"])  #Group the RKIP and BACH1 columns to make different combinations.
cancer<-c("ACC","BLCA","BRCA","CESC","CHOL","COAD","DLBC","ESCA","GBM","HNSC","KICH","KIRC","KIRP","LAML","LGG","LIHC","LUAD","LUSC","MESO","OV","PAAD","PCPG","PRAD","READ","SARC","SKCM","STAD","TGCT","THCA","THYM","UCEC","UCS","UVM")

for(i in 1:length(cancer))
{ 
  df<-x1[x1[,c("cancer type abbreviation")]==cancer[i],]
  y = df[,c("sample", "OS", "OS.time", "Type")]  #change OS to DSS,PFI or RFI if needed
  y$Type = factor(y$Type, levels = c( "R+ B-", "R+ B+", "R- B+", "R- B-"))
  df<-na.omit(df)
  fit.coxph <- coxph(Surv(OS.time, OS) ~ Type, data = df)
  ggforest(fit.coxph, data = df2,fontsize=1.2)
  ggsave(file = paste0("path/",cancer[i],".pdf"))
}
}
dev.off()