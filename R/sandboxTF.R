
# Data and Libraries ####

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
fi = c("rmsea","srmr","cfi","nnfi")

df = data.frame(read_excel("data/Database_Time_16.05.25_PD&MI.xlsx")); nrow(df)
df = df[!is.na(df$Class),]; nrow(df)

dd = data.frame(read_excel("data/Data_Dictionary.xlsx"))
dd$Var_name = gsub("[ °-]",".",dd$Var_name)

# Constructs and NA ####

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

# basic descriptives ####

length(unique(df$School)); length(unique(df$Class))
table(df$Grade)

# check missing values and compute preliminary correlations ####

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

# ------------------------------------------------------------------------------ #
# Questionnaire parents ####
var = "QSTpadua_parent"
items = dd$Var_name[dd$Level=="item"&dd$ScaleName==var]
itemsQstParent = items[!is.na(items)]
modelQstParent = paste0(gsub(" ","",var),"_lat=~",paste0(itemsQstParent,collapse="+"))
print(modelQstParent)
fitQstParent = cfa(model=modelQstParent, data=df, ordered=T)
fitQstParent@Data
fitMeasures(fitQstParent,fit.measures=fi)
modificationIndices(fitQstParent,sort.=T)[1:10,]

mi_line <- "
QSTp_4_p ~~ QSTp_5_p
QSTp_2_p ~~ QSTp_6_p
"
modelQstParentRev <- paste(c(modelQstParent, mi_line), collapse = "\n")
fitQstParentRev = cfa(model=modelQstParentRev, data=df, ordered=T)
fitQstParentRev@Data
fitMeasures(fitQstParentRev,fit.measures=fi)
modificationIndices(fitQstParentRev,sort.=T)[1:10,]

# ------------------------------------------------------------------------------ #
# Questionnaire teachers ####
var = "QSTmilan_teacher"
items = dd$Var_name[dd$Level=="item"&dd$ScaleName==var]
itemsQstTeacher = items[!is.na(items)]
modelQstTeacher = paste0(gsub(" ","",var),"_lat=~",paste0(itemsQstTeacher,collapse="+"))
print(modelQstTeacher)
fitQstTeacher = cfa(model=modelQstTeacher, data=df, ordered=T)
fitQstTeacher@Data
fitMeasures(fitQstTeacher,fit.measures=fi)
modificationIndices(fitQstTeacher,sort.=T)[1:10,]

# ------------------------------------------------------------------------------ #
# Time Reproduction ####
dd$Var_name[dd$Construct == "Time Reproduction"]

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
psych::alpha(df[,itemsNorm])$total$raw_alpha
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

normTot <- rowMeans(df[,itemsNorm], na.rm=T)
baseTot <- rowMeans(abs(df[,items]), na.rm=T)
normTotNorm <- normalize(normTot)
cor(normTot,baseTot, use = "complete.obs")
cor(normTotNorm,baseTot, use = "complete.obs")
cor(normTot,normTotNorm, use = "complete.obs")

# ------------------------------------------------------------------------------ #
# Time Orientation #### 
var = "OTm"
items = dd$Var_name[dd$Level=="item"&dd$ScaleName==var]
itemsTimeOt = items[!is.na(items)]
modelTimeOt = paste0(gsub(" ","",var),"_lat=~",paste0(itemsTimeOt,collapse="+"))
print(modelTimeOt)
fitTimeOt = cfa(model=modelTimeOt, data=df, ordered=T)
fitTimeOt@Data
fitMeasures(fitTimeOt,fit.measures=fi)
modificationIndices(fitTimeOt,sort.=T)[1:10,] # Residui correlati item 12 e 14 - unici con risposta sì-no

# ------------------------------------------------------------------------------ #
# Time Discrimination ####
dd$Var_name[dd$Construct == "Time Discrimination"]
df$RatioTD

# ------------------------------------------------------------------------------ #
# Time Estimation ####
dd$Var_name[dd$Construct == "Time Estimation"]
df$TEdd_Mean
cor(df$TE_Barca,df$TE_Ladro, use = "complete.obs")
cor(df$TE_Barca_dd,df$TE_Ladro_dd, use = "complete.obs")


################################

# Correlations ####
dcor <- subset(df, select = c(QSTp_Total_parent, QSTm_Total_teacher,
                              OTm, RatioTD,
                              TEdd_Mean, AverageDevAbs_TR,
                              Grade, Gender))
dcor$AverageDevAbs_TR <- normalize(dcor$AverageDevAbs_TR)
(corTab<-round(cor(dcor, use = "pairwise.complete"),2))
corrplot::corrplot(corTab, addCoef.col = "black")

# GAMSS ####
library(gamlss)
library(ggplot2)

# install.packages(c("gamlss", "gamlss.dist", "gamlss.add"))  # if needed
library(gamlss)
library(gamlss.dist)  # LOGNO, qLOGNO, etc.
library(gamlss.add)   # smoothers like pb()

# --- EXAMPLE DATA (remove this block and plug in your own 'dat') ---
set.seed(1)
n   <- 600
dat <- data.frame(
  age = runif(n, 0, 18)
)
# generate a log-normal outcome with age-varying mean & sd on log-scale
mu_true    <- 1 + 0.12*sin(dat$age/3)          # mean of log(Y)
sigma_true <- 0.3 + 0.01*(dat$age - 9)         # sd of log(Y)
y          <- rLOGNO(n, mu = mu_true, sigma = pmax(0.15, sigma_true))
dat$y <- y
# -------------------------------------------------------------------

# 1) FIT: log-normal GAMLSS with smooths in mu and sigma
# - mu link is identity on log-scale mean
# - sigma link is log by default (keeps sigma > 0)
fit <- gamlss(
  y ~ pb(age),               # smooth effect for mu
  sigma.fo = ~ pb(age),      # allow spread to change with age
  family   = LOGNO,          # log-normal
  data     = dat,
  trace    = FALSE
)

# 2) PERCENTILE CURVES over an age grid (e.g., 3rd..97th)
ages_grid <- seq(min(dat$age), max(dat$age), by = 0.25)
centiles_out <- centiles.pred(
  fit,
  xname  = "age",
  xvalues = ages_grid,
  cent   = c(3, 10, 25, 50, 75, 90, 97),
  plot   = TRUE,             # set to FALSE if you only want the values
  show   = TRUE,
  ylab   = "Outcome (y)",
  xlab   = "Age"
)

# centiles_out is a data.frame with columns: x, C3, C10, ..., C97
head(centiles_out)

# 3) EXACT PERCENTILES at specific ages you care about
ages_of_interest <- c(2, 5, 10, 15)
newdat <- data.frame(age = ages_of_interest)

# Get fitted distribution parameters at those ages
pa <- predictAll(fit, newdata = newdat)  # returns mu, sigma (for LOGNO)

# Choose percentiles you want
p <- c(.03, .10, .25, .50, .75, .90, .97)

# Build a tidy table of percentile values by age
perc_mat <- sapply(seq_len(nrow(newdat)), function(i) {
  qLOGNO(p, mu = pa$mu[i], sigma = pa$sigma[i])
})
dimnames(perc_mat) <- list(paste0(p*100, "%"), paste0("age=", newdat$age))
percentiles_by_age <- t(perc_mat)
percentiles_by_age




df1 <- df[is.na(df$QSTp_Total_parent)==FALSE,c("Grade","QSTp_Total_parent","Age")]
df1$QSTp_scaled <- df1$QSTp_Total_parent / 30

model <- gamlss(
  QSTp_scaled ~ pb(Age),             # Smooth function for mean
  sigma.fo = ~ pb(Age),               # Smooth function for SD
  data = df1,
  family = BEINF                            # Box-Cox Power Exponential, flexible for skewed data
)

# New grade values to predict over
newdata <- data.frame(Age = seq(min(df$Age), max(df$Age), by = 0.1))

# Predict model parameters
predData <- predictAll(model, newdata, type = "response")

predData <- cbind.data.frame(newdata,predData)

ggplot(predData, aes(x = Age, y = mu)) +
  geom_point()

## https://rdrr.io/cran/gamlss/man/centiles.pred.html
centiles(model,xvar = df1$Age)

centiles.pred(model,
              xname = "Age",
              xvalues = seq(min(df1$Age), max(df1$Age), by = 0.1),
              cent = c(2, 10, 25, 50, 75, 90, 98))


newx<-seq(5,12,1)
mat <- centiles.pred(model, xname="Age", xvalues=newx )
mat

## now plot the centile curves  
mat <- centiles.pred(model, xname="Age",xvalues=newx, plot=TRUE )


## bring the data and fit the model
data(abdom)
a<-gamlss(y~pb(x),sigma.fo=~pb(x), data=abdom, family=BCT)
## plot the centiles
centiles(a,xvar=abdom$x)
##-----------------------------------------------------------------------------
## the first use of the function centiles.pred()
## to calculate the centiles at new x values
##-----------------------------------------------------------------------------
newx<-seq(12,40,2)
mat <- centiles.pred(a, xname="x", xvalues=newx )
mat
## now plot the centile curves  
mat <- centiles.pred(a, xname="x",xvalues=newx, plot=TRUE )
##-----------------------------------------------------------------------------
## the second use of the function centiles.pred()
## to calculate (nornalised) standard-centiles for new x
## values using the fitted model
##-----------------------------------------------------------------------------
newx <- seq(12,40,2)
mat <- centiles.pred(model, xname="x",xvalues=newx, type="standard-centiles" )
mat
## now plot the standard centiles  
mat <- centiles.pred(a, xname="x",xvalues=newx, type="standard-centiles",
                     plot = TRUE )
##-----------------------------------------------------------------------------
## the third use of the function centiles.pred()
##  if we have new x and y values what are their z-scores?
##-----------------------------------------------------------------------------
# create new y and x values and plot them in the previous plot
newx <- c(20,21.2,23,20.9,24.2,24.1,25)
newy <- c(130,121,123,125,140,145,150)
for(i in 1:7) points(newx[i],newy[i],col="blue")
## now calculate their z-scores
znewx <- centiles.pred(a, xname="x",xvalues=newx,yval=newy, type="z-scores" )
znewx
## Not run: 
##-----------------------------------------------------------------------------
## What we do if the x variables is transformed?
##----------------------------------------------------------------------------
##  case 1 : transformed x-variable within the formula
##----------------------------------------------------------------------------
## fit model
aa <- gamlss(y~pb(x^0.5),sigma.fo=~pb(x^0.5), data=abdom, family=BCT)
## centiles is working in this case
centiles(aa, xvar=abdom$x, legend = FALSE)
## get predict for values of x at 12, 14, ..., 40
mat <- centiles.pred(aa, xname="x", xvalues=seq(12,40,2), plot=TRUE )
mat
# plot all prediction points
xx <- rep(mat[,1],9)
yy <- unlist(mat[,2:10])
points(xx,yy,col="red")
##----------------------------------------------------------------------------
##  case 2 : the x-variable is previously transformed 
##----------------------------------------------------------------------------
nx <- abdom$x^0.5
aa <- gamlss(y~pb(nx),sigma.fo=~pb(nx), data=abdom, family=BCT)
centiles(aa, xvar=abdom$x)
# equivalent to fitting
newd<-data.frame( abdom, nx=abdom$x^0.5)
aa1 <- gamlss(y~pb(nx),sigma.fo=~pb(nx), family=BCT, data=newd)
centiles(aa1, xvar=abdom$x)
# getting the centiles at x equal to 12, 14, ...40
mat <-  centiles.pred(aa, xname="nx", xvalues=seq(12,40,2), power=0.5, 
                      data=newd, plot=TRUE)
# plot all prediction points         
xxx <- rep(mat[,1],9)
yyy <- unlist(mat[,2:10])
points(xxx,yyy,col="red")
# the idea is that if the transformed x-variable is used in the fit
# the power argument has to used in centiles.pred()

## End(Not run)