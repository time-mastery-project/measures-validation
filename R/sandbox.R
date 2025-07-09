
################################

rm(list=ls())
library(readxl)
library(lavaan)
library(semTools)
library(psych)
library(ggplot2)
library(tidyr)
library(dplyr)
library(corrplot)
source("R/custom-functions.R")
fi = c("rmsea","cfi","nnfi")

df = data.frame(read_excel("data/Database_Time_16.05.25_PD&MI.xlsx")); nrow(df)
df = df[!is.na(df$Class),]; nrow(df)

dd = data.frame(read_excel("data/Data_Dictionary.xlsx"))
dd$Var_name = gsub("[ °-]",".",dd$Var_name)

################################

# list of constructs and ratios of complete observations

tb = as.data.frame(table(dd$ScaleName,dd$Construct))
tb = tb[tb$Freq!=0,1:2]
names(tb) = c("Scale","Construct")
tb$completeObs = NA
for(i in 1:nrow(tb)){
  x = dd$Var_name[dd$ScaleName==tb$Scale[i] & dd$Construct==tb$Construct[i]]; x = x[!is.na(x)]
  tb$completeObs[i] = round(mean(rowMeans(!is.na(df[,x]))),3)
}; tb
# write.csv(tb,"tb.csv",row.names=F)

################################

# basic descriptives 

length(unique(df$School)); length(unique(df$Class))
table(df$Grade)

################################

# check missing values and compute preliminary correlations

mean(apply(is.na(df[,paste0("OTm_",1:16)]),1,sum)==0)
mean(apply(is.na(df[,paste0("OTp_",1:10)]),1,sum)==0)
mean(apply(is.na(df[,paste0("QSTm_",1:10,"_t")]),1,sum)==0)
mean(apply(is.na(df[,paste0("QSTp_",1:10,"_p")]),1,sum)==0)

df$OTm = apply(df[,paste0("OTm_",1:16)],1,sum)
df$QSTmt = apply(df[,paste0("QSTm_",1:10,"_t")],1,sum)
df$QSTmp = apply(df[,paste0("QSTm_",1:10,"_p")],1,sum)

mean(is.na(df$TEdd_Mean))
mean(is.na(df$AverageDevAbs_TR))
mean(is.na(df$RatioTD))

hist(df$OTm,breaks=30)
hist(df$RatioTD,breaks=20)
hist(log(df$RatioTD),breaks=20)

mean(is.na(df$SDAI_dis))

names(df)

cor(df$RatioTD,df$SDAI_dis,use="pairwise.complete",method="spearman")

cor(df$AverageDevAbs_TR,df$SDAI_dis,use="pairwise.complete",method="spearman")
cor(df$AverageDevAbs_TR,df$SDAI_dis,use="pairwise.complete",method="spearman")

cor(df$OTm,df$SDAI_dis,use="pairwise.complete",method="spearman")
cor(df$QSTmt,df$SDAI_dis,use="pairwise.complete",method="spearman")
cor(df$QSTmp,df$SDAI_dis,use="pairwise.complete",method="spearman")
cor(df$QSTmp,df$QSTmt,use="pairwise.complete",method="spearman")

################################

#### MEASUREMENT MODELS 

unique(dd$Construct)
unique(dd$ScaleName)
vars = c("Time Knowledge","Time Orientation",
         "QSTpadua_parent",
         "QSTmilan_teacher")

var = "OTm"
items = dd$Var_name[dd$Level=="item"&dd$ScaleName==var]
items = items[!is.na(items)]
model = paste0(gsub(" ","",var),"_lat=~",paste0(items,collapse="+"))
print(model)
fit = cfa(model=model, data=df, ordered=T)
fit@Data
fitMeasures(fit,fit.measures=fi)
modificationIndices(fit,sort.=T)[1:10,]

################################

# TIME REPRODUCTION 

items = dd$Var_name[dd$ScaleName=="TR"&dd$Level=="item"]
items = items[!is.na(items) & grepl("PercDevAbs_TR_",items)]
items
itemsNorm = paste0("TR",2:12)
df[,itemsNorm] = NA
for(i in 1:length(items)) df[,itemsNorm[i]] = normalize(abs(df[,items[i]]))
model = paste("TR_lat =~",paste(itemsNorm,collapse="+"))
model
fit = cfa(model=model, data=df)
fit@Data
fitMeasures(fit,fit.measures=fi)
modificationIndices(fit,sort.=T)[1:10,]
psych::alpha(df[,itemsNorm])
corrplot(cor(df[,itemsNorm],use="pairwise.complete",method="pearson"),
         method = "color",  
         type = "full",     
         addCoef.col = "black", 
         tl.col = "black",  
         tl.cex = 1, 
         number.cex = 1)

dfLong = df %>% pivot_longer(cols=all_of(items), names_to="item", values_to="value")
ggplot(dfLong, aes(x = "", y = value)) +
  geom_violin() +
  facet_wrap(~ item, scales = "free") +
  theme_minimal() +
  labs(x = NULL, y = "Value", title = "Violin Plots for Selected Items")

################################

