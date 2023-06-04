
######################################################### statistics for nirS project #########################################################

install.packages('nortest')

library(ggplot2)
library(readxl)
library(car)
library(nortest)


####################### read files #######################
nirSproject.data<-read_excel("./datasum_nirS project.xlsx")

nirSproject.data$sym.density<-as.numeric(nirSproject.data$sym.density)
nirSproject.data$Fv.Fm<-as.numeric(nirSproject.data$Fv.Fm)
nirSproject.data$nirS.fold.change<-as.numeric(nirSproject.data$nirS.fold.change)
nirSproject.data$CN<-as.numeric(nirSproject.data$CN)
nirSproject.data$prop.A<-as.numeric(nirSproject.data$prop.A)
nirSproject.data$prop.B<-as.numeric(nirSproject.data$prop.B)

nirSproject.data$Host<-as.factor(nirSproject.data$Host)
nirSproject.data$Sym<-as.factor(nirSproject.data$Sym)


####################### coral physiologt #######################

####################### symbiont density #######################
symb_anova <- aov(log10 (sym.density) ~ Host * Sym, data=nirSproject.data)
plot(fitted(symb_anova), residuals(symb_anova))
# data normality
ad.test(residuals(symb_anova)) # P=0,6895, should >0,05
cvm.test(residuals(symb_anova)) # P=0,7364, should >0,05
shapiro.test(residuals(symb_anova)) # P=0,6063, should > 0,05
# statistical analysis output
summary(symb_anova)
TukeyHSD(symb_anova)


####################### Fv/Fm #######################
Fv.Fm_anova <- aov(log10 (Fv.Fm) ~ Host * Sym, data=nirSproject.data)
plot(fitted(Fv.Fm_anova), residuals(Fv.Fm_anova))

ad.test(residuals(Fv.Fm_anova)) # P=0,6728, should >0,05
cvm.test(residuals(Fv.Fm_anova)) # P=0,5413, should >0,05
shapiro.test(residuals(Fv.Fm_anova)) # P=0,9404, should > 0,05

summary(Fv.Fm_anova)
TukeyHSD(Fv.Fm_anova)


####################### C/N ratio #######################
CN_anova <- aov(log10 (CN) ~ Host * Sym, data=nirSproject.data)
plot(fitted(CN_anova), residuals(CN_anova))

ad.test(residuals(CN_anova)) # P=0,06147, should >0,05
cvm.test(residuals(CN_anova)) # P=0,05929, should >0,05
shapiro.test(residuals(CN_anova)) # P=0,158, should > 0,05

summary(CN_anova)
TukeyHSD(CN_anova)


####################### relative nirS abundance #######################
nirS.fold.change_anova <- aov(log10 (GeneCopy) ~ Host * Sym, data=nirSproject.data)
plot(fitted(nirS.fold.change_anova), residuals(nirS.fold.change_anova))

ad.test(residuals(nirS.fold.change_anova)) # P=0,1973, should >0,05
cvm.test(residuals(nirS.fold.change_anova)) # P=0,2258, should >0,05
shapiro.test(residuals(nirS.fold.change_anova)) # P=0,2684, should > 0,05

summary(nirS.fold.change_anova)
TukeyHSD(nirS.fold.change_anova)


####################### qPCR delta Ct #######################
ab_nirS=read.table("./qPCR_nirS_deltaCalApril22.txt", header = TRUE, sep ='\t')

ab_nirS$ID<-paste(ab_nirS$Host, ab_nirS$Symbiont, sep = "_")
ab_nirS$deltaCt<-as.numeric(ab_nirS$deltaCt)
ab_nirS$Host <-as.factor(ab_nirS$Host)
ab_nirS$Symbiont <-as.factor(ab_nirS$Symbiont)
ab_nirS=ab_nirS %>% drop_na()

sum.qPCR <- ddply(ab_nirS, c("ID"), summarise,
                  N = length (deltaCt),
                  mean = mean(deltaCt),
                  sd   = sd(deltaCt),
                  se   = sd / sqrt(N))

sum.qPCR = separate(sum.qPCR, ID, into = c("Host", "Symbiont"), sep = "_", remove = F)
sum.qPCR$Host=as.factor(sum.qPCR$Host)
sum.qPCR$Symbiont=as.factor(sum.qPCR$Symbiont)



####################### Non-normal distributed data #######################
#Kruskal–Wallis Test
#if(!require(psych)){install.packages("psych")}
if(!require(FSA)){install.packages("FSA")}
if(!require(lattice)){install.packages("lattice")}

if(!require(coin)){install.packages("coin")}
if(!require(multcompView)){install.packages("multcompView")}
if(!require(rcompanion)){install.packages("rcompanion")}
if(!require(PMCMRplus)){install.packages("PMCMRplus")}

str(ab_nirS)
kruskal.test(deltaCt ~ Host,
             data = ab_nirS)
kruskal.test(deltaCt ~ Symbiont,
             data = ab_nirS)
