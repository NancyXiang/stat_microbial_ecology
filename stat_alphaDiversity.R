
######################################################### Stats for Beta diversity #########################################################


####################### two-way ANOVA for Alpha Diversity #######################

# Chao1
chao_AN <- aov( S.chao1 ~ host * symbiont, data=alpha_six)
plot(fitted(chao_AN), residuals(chao_AN))

library(nortest)
ad.test(residuals(chao_AN)) # P=0.5345, should >0,05, Anderson-Darling test for normality
cvm.test(residuals(chao_AN)) # P=0,6758, should >0,05, Cramer-Von Mises Test for normal distribution
shapiro.test(residuals(chao_AN)) # P=0,3243, should > 0,05, Normality Test

summary(chao_AN)
TukeyHSD(chao_AN)


# Shannon
Shannon_AN <- aov( Shannon ~ host * symbiont, data=alpha_six)

library(nortest)
ad.test(residuals(Shannon_AN)) # P=0.08
cvm.test(residuals(Shannon_AN)) # P=0,09
shapiro.test(residuals(Shannon_AN)) # P=0,08

summary(Shannon_AN)
TukeyHSD(Shannon_AN)


# Simpson
simpson_AN <- aov(log10(simpson) ~ host * symbiont, data=alpha_six) #log transformation to meet normality

ad.test(residuals(simpson_AN)) # P=0.2332
cvm.test(residuals(simpson_AN)) # P=0,227
shapiro.test(residuals(simpson_AN)) # P=0,2131

summary(simpson_AN)
TukeyHSD(simpson_AN)


####################### t-test for Alpha Diversity #######################
# Chao1
t.test(S.chao1~ symbiont, alpha_Sym)# p-value = 0.296

# Shannon
t.test(Shannon~ symbiont, alpha_Sym)#p-value = 0.8011

# Simpson
t.test(simpson~ symbiont, alpha_Sym)#p-value = 0.7664
