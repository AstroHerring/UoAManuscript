#–––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––-
# Unit of analysis issues in laboratory based research: 
# a review of concepts and guidance on study design and reporting
#–––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––-
# Nick R. Parsons, M. Dawn Teare and Alice J. Sitch
# 12th December 2017
#–––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––-
# R code to reproduce analyses in Appendix 3
#–––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––-
# Code was developed in R 3.3.1, but should run in any version
# provided pacakges are updated appropriately
# Users need to install the packages nlme and lme4:
install.packages("nlme")
library(nlme)
install.packages("lme4")
library(lme4)
#–––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––-
# Example 1
# Adjuvant radiotherapy and lymph node size in colorectal cancer
#–––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––-
#
# Data
LNsize <- c(1.71, 1.72, 1.98, 1.98, 1.88, 1.85, 2.51, 2.98,
            2.55, 3.20, 2.65, 2.80, 1.69, 1.82, 1.72, 1.97,
            1.80, 1.73, 1.72, 2.50, 1.78, 2.65, 2.04, 2.77,
            3.32, 3.11, 3.27, 3.03, 3.07, 3.11, 2.33, 2.86,
            2.48, 2.87, 2.53, 2.52, 2.37, 2.36, 2.36, 2.62,
            2.20, 2.60, 1.33, 1.90, 1.35, 1.87, 1.15, 1.85,
            1.70, 2.07, 1.78, 1.76, 1.78, 1.85, 2.23, 2.50,
            2.14, 2.33, 2.21, 2.16, 2.10, 2.11, 1.89, 2.16,
            1.75, 2.12, 2.58, 2.77, 2.54, 2.65, 2.59, 2.60)
Subject <- factor(rep(1:12, each = 6), levels = 1:12)
Sample <- factor(rep(1:2, times = 36), levels = 1:2)
Slice <- factor(rep(rep(1:3, each = 2), times = 12), levels = 1:3)
RadioTherapy <- factor(rep(1:2, each = 36), levels = 1:2,
labels = c("None", "RTShort"))
LymphNode <- data.frame(Subject, Sample, Slice, RadioTherapy, LNsize)
#
#–––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––-
#
# Strip plot
par(pty = "s", mar = c(4.5,4.5,2,2))
stripchart(LNsize ~ Subject, data = LymphNode, vertical = TRUE, las = 1,
     ylab = "Lymph Node Size (mm)", xlab = "Subject", ylim = c(1, 3.5),
     pch = 21, col = 1, bg = 1, cex.axis = 1.2, subset = Sample == 1)
stripchart(LNsize ~ Subject, data = LymphNode, vertical = TRUE, 
     pch = 21, col = 1, bg = 0, subset = Sample == 2, add = TRUE)
legend(x = 1, y = 1.25, legend = c("Sample 1", "Sample 2"), pch = 21, 
      pt.bg = c(1, 0), col =c(1, 1), bty = "n")
lines(x = c(6.5, 6.5), y =c(0, 4), lty = 2, col = 1)
text(x = c(3.5, 9.5), y = c(3.5, 3.5), c("None", "Short RT"), cex = 1.2)
#
#–––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––-
#
# Naive analysis
t.test(LNsize ~ RadioTherapy, var.equal = TRUE, data = LymphNode)
summary(mod <- lm(LNsize ~ RadioTherapy, data = LymphNode))
cbind(coef(mod), confint(mod))
#
#–––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––-
#
# Analysis by subject
LNsize.means <- tapply(LymphNode$LNsize, list(LymphNode$Subject),
                            mean, na.rm = TRUE)
RT.means <- factor(rep(1:2, each = 6), levels = 1:2,
labels = c("None", "RTShort"))
t.test(LNsize.means ~ RT.means, var.equal = TRUE)
summary(mod.lm <- lm(LNsize.means ~ RT.means))
anova(mod.lm)
cbind(coef(mod.lm), confint(mod.lm))
#
#–––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––-
#
# Mixed-effect model
mod.lme <- lme(LNsize ~ RadioTherapy, random = ~ 1 | Subject / Sample,
data = LymphNode)
summary(mod.lme)
anova(mod.lme)
# Estimates
intervals(mod.lme, which = "fixed")
intervals(mod.lme, which = "var-cov")
# Compare models
mod0.lme <- update(mod.lme, random = ~ 1 | Subject)
anova(mod0.lme, mod.lme)
mod0.lme <- update(mod.lme, random = ~ 1 | Sample)
anova(mod0.lme, mod.lme)
# Diagnostic plots
# Standardized residuals versus fitted values by subject
plot(mod.lme, resid(., type = "response") ~ fitted(.) | Subject,
                                    abline = 0)
# Observed versus fitted values by subject
plot(mod.lme, LNsize ~ fitted(.) | Subject, abline = c(0,1))
# Box-plots of residuals by Subject
plot(mod.lme, Subject ~ resid(., type = "response"), aspect = "fill", 
 scales = list(cex = 1.1), xlab = list(label = "Residuals", cex = 1.2), 
 ylab = list(label = "Subject", cex = 1.2),
 main = "Boxplots of residuals")
# Quantile-quantile plot
qqnorm(resid(mod.lme, type = "response"), pch = 19, col = 1,
                main = "Q-Q plot of residuals", las = 1)
qqline(resid(mod.lme, type = "response"), lty = 2)
#
#–––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––-
#
# Unbalanced data analysis
# Randomly remove some data and re-fit model
set.seed(8845391)
remove.cells <- sample(1:72, 36, replace=FALSE)
Unbalanced.LymphNode <- LymphNode[setdiff(1:72, remove.cells),]
umod.lme <- lme(LNsize ~ RadioTherapy, random = ~ 1 | Subject / Sample,
data = Unbalanced.LymphNode)
summary(umod.lme)
anova(umod.lme)
# Estimates
intervals(umod.lme, which = "fixed")
intervals(umod.lme, which = "var-cov")
# Analysis by subject
tapply(Unbalanced.LymphNode$LNsize, list(Unbalanced.LymphNode$Subject),
mean, na.rm = TRUE)
UB.LNsize.means <- tapply(Unbalanced.LymphNode$LNsize,
list(Unbalanced.LymphNode$Subject), mean, na.rm = TRUE)
t.test(UB.LNsize.means ~ RT.means, var.equal = TRUE)
summary(umod.lm <- lm(UB.LNsize.means ~ RT.means))
confint(umod.lm)
#
#–––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––-
# Example 2
# Grouped binary lymph node count data
#–––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––-
#
# Data
LN.ind <- c(4, 3, 2, 2, 3, 2, 1, 1, 1, 2, 4, 3,
            4, 4, 3, 4, 4, 5, 0, 2, 0, 1, 2, 4,
            NA, 5, 3, 1, 4, 5, 0, NA, 1, 4, 4, 3,
            NA, 2, 2, 2, 3, 3, 0, NA, 0, 0, 3, NA,
            NA, NA, NA, 1, 5, 3, 0, NA, 2, 2, 3, NA)
Subject <- factor(rep(1:12, times = 5), levels = 1:12)
Sample <- factor(rep(1:5, each = 12), levels = 1:5)
nRoI <- rep(5, 60)
RadioTherapy <- factor(rep(rep(1:2, each = 6),times=5), levels = 1:2,
labels = c("None", "RTShort"))
gLymphNodeInd <- data.frame(Subject, Sample, RadioTherapy,
gLNind=as.numeric(LN.ind)/5,nRoI)
gLymphNodeInd <- subset(gLymphNodeInd, complete.cases(gLymphNodeInd))
#
#–––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––-
#
# Logistic regression
tapply(gLymphNodeInd$gLNind, list(gLymphNodeInd$RadioTherapy),
function(x){5 * sum(x, na.rm = TRUE)})
log.reg <- glm(gLNind ~ RadioTherapy, data = gLymphNodeInd,
family=binomial("logit"),weight=nRoI)
summary(log.reg)
anova(log.reg)
exp(coef(log.reg))
exp(confint(log.reg))
#
#–––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––-
#
# ME logistic regression
gmod.lme4 <- glmer(gLNind ~ RadioTherapy + (1 | Subject),
data = gLymphNodeInd, family=binomial("logit"),weight=nRoI)
summary(gmod.lme4)
anova(gmod.lme4)
# Standard errors
se <- sqrt(diag(vcov(gmod.lme4)))
# 95% CIs
par.CI <- confint(gmod.lme4, method = "Wald")
cbind(exp(fixef(gmod.lme4)), exp(par.CI[2:3,]))
# ranef(gmod.lme4)
# Predictions for each RT group
predict(gmod.lme4, newdata = data.frame(RadioTherapy = "RTShort"),
                     type = "response",re.form = ~ 0)
predict(gmod.lme4, newdata = data.frame(RadioTherapy = "None"),
                     type = "response",re.form = ~ 0)
#
#–––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––-
# end 
#–––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––-