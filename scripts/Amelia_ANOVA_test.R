library(ggplot2)
library(ggpubr)
library(tidyverse)
library(broom)
#install.packages("AICcmodavg")
library(AICcmodavg)


# #ANOVA tutorial ---------------------------------------------------------


#tutorial from https://www.scribbr.com/statistics/anova-in-r/
#load in data
crop.data <- read.csv("C:\\Users\\ameli\\Downloads\\crop.data_.anova_\\crop.data.csv", header = TRUE, colClasses = c("factor", "factor", "factor", "numeric"))

#check that the data seems correct
summary(crop.data)

#first we can do a one way ANOVA, model crop yield as a function of the fertilizer type
#can see if one type of fertilizer gives more yield than others
ggplot(crop.data, aes(x = fertilizer, y = yield)) + geom_boxplot() + geom_point() + geom_jitter()

#run the anova
one.way <- aov(yield ~ fertilizer, data = crop.data)

#display the results
summary(one.way)

#we can see there are 2 degrees of freedom since there are three fertilizer types (df = n -1, variables free to vary minus one)
#residuals, this is the variation not explained by the independent variables (ie fertilizer)
#degrees of freedom for residuals = 93 (number of observations minus one and minus the number of fertilizer variables)
#sum sq = sum of squares, total variation between group means and overall mean
#mean sq = mean sum of squares, found by dividing the sum of squares by the degrees of freedom for each parameter
#Fvalue = test statistic from the F test, the mean square of each independent variable divided by the mean square of residiuals.
#Large F value = more likely the variation caused by independent varibale is real and not due to chance
#Pr(>F) = p value of the F statistic. How likely the F value is if the null hypothesis (no difference between any fertilizer group) is true

#so this anova shows a high F value 7.863 with a low p value 7e-04, meaning we can reject the null hypothesis, the group means are statistically different

#in a Two-way anova we can model if the dependent variable is the result of two indepdendent variables 
#for example, is crop yield a funciton of fertilizer type and planting density?

two.way <- aov(yield ~ fertilizer + density, data = crop.data)

summary(two.way)

#we can see, there is a high F value and low p value for both fertilizer and density.
#Meaning they have a statistically significant effect on crop yield.
#this still doesn't explain all of the variance but it reduced the sum sq of variance from 35 to 30

#we can also test to see if two variables have an interaction effect instead of an additive effect
#for example, does planting density affect the ability to uptake fertilizer?

interaction <- aov(yield ~ fertilizer*density, data = crop.data)

summary(interaction)

#the fertilizer demsity interaction has a low f value and a low p value 
#so this interaction doesn't appear to be significant

#a blocking variable can be added if there appears to be a confounding variable to control for
#for example, farmers plant crops in different blocks of fields that have different soil, sunlight etc
#we may want to control for block in our ANOVA

blocking <- aov(yield ~ fertilizer + density + block, data = crop.data)

summary(blocking)

#block has a low f value and high p value so its probably not affecting the model 
#adding block also doesn't change the sum of squares for either independent variable
#meaning it isn't affecting the amount of variation those variables explain. Can probably be excluded

#we've now made four ANOVA models to explain our crop data
#we can chose the best model using AIC scores, which weight the amount of variation explained by the number of parameters
#lowest number explains the most variance aka better fit
#we will use aictab function to calculate the AIC scores for each of the models

model.set <- list(one.way, two.way, interaction, blocking)
model.names <- c("one.way", "two.way", "interaction", "blocking")

aictab(model.set, modnames = model.names)

#in this case the lowest AIC score is for the two.way model with a weight of 0.71
#this model explains 71% of variance in the dependent variable

#check for homoscedasticity (an assumption of ANOVA that is an assumption of equal or similar variances in different groups being compared) in the two-way model
par(mfrow=c(2,2))
plot(two.way)
par(mfrow=c(1,1))
#mean of residuals (red line) should be roughly horizontal and centered on zero, meaning no large outliers
#the Q-Q plot should have a slope of 1, it plots the observed residuals against perfectly homoscedastic residuals

#if your model doesn't fit this assumption you can do the kruskall-wallis test instead

#post-hoc test
#ANOVA tells us if there are differences among group means, not what the differences are
#to find out which groups are statistically different from each other we can do a Tukey's HSD test for pairwise comparisons
#Tukey HSD: Honestly significantly different

turkey.two.way <- TukeyHSD(two.way)
turkey.two.way

#the differences between fertilizer group 3-1 and 3-2 is significant but not between group 2-1
#the two density groups also have a statistically significant difference

#plot the results in a graph
#we want to know which combinations of fertilizer and planting density are statistically different from each other
#we can do an ANOVA-TukeyHSD using interaction of fertilizer and planting density to do this
#we know this model isn't as good as the two way but we can use it to add info to the graph
tukey.plot.aov<-aov(yield ~ fertilizer:density, data=crop.data)
tukey.plot.test<-TukeyHSD(tukey.plot.aov)
#display the results as a graph instead of as a table
plot(tukey.plot.test, las = 1)

#significant groupwise differences are any where the 95% CI doesn't include zero
#3:1-1:1 (read as “fertilizer type three + planting density 1 contrasted with fertilizer type 1 + planting density type 1”)
#and 1:2-1:1, 2:2-1:1, 3:2-1:1, and 3:2-2:1.
#we can divide out group into group A (for 1:1), B(intermediate combinations), C(for 3:2)

#first summarize by density and fertilizer
mean.yield.data <- crop.data %>%
  group_by(fertilizer, density) %>%
  summarise(
    yield = mean(yield)
  )

#add the new group labels
mean.yield.data$group <- c("a","b","b","b","b","c")
mean.yield.data

#plot the raw data 
two.way.plot <- ggplot(crop.data, aes(x = density, y = yield, group=fertilizer)) +
  geom_point(cex = 1.5, pch = 1.0,position = position_jitter(w = 0.1, h = 0))
two.way.plot

#add standard error and means
two.way.plot <- two.way.plot +
  stat_summary(fun.data = 'mean_se', geom = 'errorbar', width = 0.2) +
  stat_summary(fun.data = 'mean_se', geom = 'pointrange') +
  geom_point(data=mean.yield.data, aes(x=density, y=yield))
two.way.plot

#facet wrap by fertilizer so they're not stack on top of each other
two.way.plot <- two.way.plot +
  geom_text(data=mean.yield.data, label=mean.yield.data$group, vjust = -8, size = 5) +
  facet_wrap(~ fertilizer)
two.way.plot

#Example of writeup
#We found a statistically-significant difference in average crop yield by both fertilizer type (F(2)=9.018, p < 0.001) and by planting density (F(1)=15.316, p < 0.001).
#A Tukey post-hoc test revealed that fertilizer mix 3 resulted in a higher yield on average than fertilizer mix 1 (0.59 bushels/acre), and a higher yield on average than fertilizer mix 2 (0.42 bushels/acre). Planting density was also significant, with planting density 2 resulting in an higher yield on average of 0.46 bushels/acre over planting density 1.
#A subsequent groupwise comparison showed the strongest yield gains at planting density 2, fertilizer mix 3, suggesting that this mix of treatments was most advantageous for crop growth under our experimental conditions.


#practice with real data
filename <- "whippomorpha_four_state_max_crep_ER_SYM_ARD_bridge_only_models.rds"

#requires the filename and the number of Mk models (3: ER, SYM, ARD or 4: ER, SYM, ARD, bridge_only)
#function returns a dataframe of the likelihoods for all 1k trees x number of Mk models
df_full <- plot1kLikelihoods(readRDS(here(paste0("finalized_1k_models/", filename))), 4)

#One way ANOVA, does likelihood vary with the Mk model type
one.way <- aov(likelihoods ~ model, data = df_full)
summary(one.way)
#df = 4 (number of variable levels free to vary - 1, four MK models -1)
#high f value, very low p value. So model does have a strong effect on likelihood

#Tukey_HSD
tukey.one.way <- TukeyHSD(one.way)
tukey.one.way

#shows bridge-ARD are not statistically different but all other model comparisons are

#check if variance  of all models is equal
par(mfrow=c(2,2))
plot(one.way)
par(mfrow=c(1,1))

# # Mean comparison on plots tutorial -------------------------------------

#tutorial from https://www.sthda.com/english/articles/24-ggpubr-publication-ready-plots/76-add-p-values-and-significance-levels-to-ggplots/
#install.packages("ggpubr")
library(ggpubr)

#we'll use the built in toothgrowth dataset
data("ToothGrowth")
head(ToothGrowth)

#T test is used to compare the means of two groups (parametric)
#Wilcoxon test is also used to compare two groups (non-parametric)
#ANOVA is used to compare the means of multiple groups (more than two) (parametric)
#Kruskal-Wallis is also used to compare multiple groups (non-parametric)

#the function compare_means() can be used to input your data and use any of the above mean comparisons
#function stat_compare_means() is used add these mean comparisons to ggplot

#example 1: Compare two unpaired samples
#By default compare means will do a wilcoxon test
#we can compare length vs sup

compare_means(len ~ supp, data = ToothGrowth)
#there is a non significant association (p = 0.064)

#we can plot the association between these two variables (len and supp)
p <- ggboxplot(ToothGrowth, x = "supp", y = "len",
               color = "supp", palette = "jco",
               add = "jitter")

#then we use stat_compare_means to add the p value (are the means statistically different)
p + stat_compare_means()
#we can also easily change the method we use to compare means, for instance switch to t test
p + stat_compare_means(method = "t.test")

#by default stat_compare_means will display the method and p value, can also be adjusted to include other stats
#for example can include the significance level and chance the position of the label
p + stat_compare_means(aes(label = ..p.signif..), label.x = 1.5, label.y = 40)

#example 2: compare two paired samples
#we will run the wilcoxon test on len and supp but this time assume the samples are paired (from the same tooth)
compare_means(len ~ supp, data = ToothGrowth, paired = TRUE)
#this time we can see that the difference is significant, p = 0.0043

#we can visualize this with ggpaired
ggpaired(ToothGrowth, x = "supp", y = "len",
         color = "supp", line.color = "gray", line.size = 0.4,
         palette = "jco")+
  stat_compare_means(paired = TRUE)

#example 3: global test
#we can first do a global test with ANOVA, to compare len as a function of dose
compare_means(len ~ dose, data = ToothGrowth, method = "anova")

#we can see globally there is a significant affect of dose on length. p = 9.53e-16
#we can use stat_compare_means to add these comparisons to a boxplot of len vs dose

# Default method = "kruskal.test" for multiple groups
ggboxplot(ToothGrowth, x = "dose", y = "len",
          color = "dose", palette = "jco")+
  stat_compare_means()
# Change method to anova
ggboxplot(ToothGrowth, x = "dose", y = "len",
          color = "dose", palette = "jco")+
  stat_compare_means(method = "anova")

#can see for both methods the p value is significant

#if the grouping variable (ie dose) contains more than two levels (ie 0.5, 1, 2), the function will automatically do a pairwise comparison
#will default t wilcoxon, can specify t test as well
compare_means(len ~ dose,  data = ToothGrowth)

#can see that it first compares 0.5 dose to 1 dose, then 0.5 dose to 2 dose, then 1 dose to 2 dose

#we can visualize this
my_comparisons <- list( c("0.5", "1"), c("1", "2"), c("0.5", "2") )
ggboxplot(ToothGrowth, x = "dose", y = "len",
          color = "dose", palette = "jco")+ 
  stat_compare_means(comparisons = my_comparisons)+ # Add pairwise comparisons p-value
  stat_compare_means(label.y = 50)     # Add global p-value
#can see that globally the p value is significant, and for each comparison it is also significant

#example with my data :o
filename <- "whippomorpha_four_state_max_crep_ER_SYM_ARD_bridge_only_models.rds"
df_full <- plot1kAIC(readRDS(here(paste0("finalized_1k_models/", filename))), 4)

my_comparisons <- list( c("ER", "SYM"), c("ER", "ARD"), c("ER", "bridge_only"), c("SYM", "ARD"), c("SYM", "bridge_only"), c("ARD", "bridge_only"))
ggboxplot(df_full, x = "model", y = "AIC_score",
          color = "model", palette = "jco")+ 
  stat_compare_means(method = "t.test", comparisons = my_comparisons, aes(label = ..p.signif..))+ # Add pairwise comparisons p-value
  stat_compare_means(label.y = 255)+     # Add global p-value
  ggtitle(filename)

#add to original plots
if(length(unique(df_full$model)) == 4){
  my_comparisons <- list( c("ER", "SYM"), c("ER", "ARD"), c("ER", "bridge_only"), c("SYM", "ARD"), c("SYM", "bridge_only"), c("ARD", "bridge_only"))
} 

if(length(unique(df_full$model)) == 3){
  my_comparisons <- list( c("ER", "SYM"), c("ER", "ARD"), c("SYM", "ARD"))
}

ggplot(df_full, aes(x = fct_inorder(model), y = AIC_score)) + geom_jitter(alpha = 0.6, color = "#F8766D") + 
  geom_boxplot(alpha = 0.5, outlier.shape = NA, colour = "black") + ggtitle(filename) +
  stat_compare_means(label = "p.signif", method = "t.test", comparisons = my_comparisons) +
  stat_compare_means(method = "anova", label.y = (max(df_full$AIC_score)+10))

#check if the data has a normal distribution
ggplot(df_full, aes(y = AIC_score)) + geom_histogram(alpha = 0.5, colour = "black")  + ggtitle(filename) + facet_wrap(~model)

#ER is not normal, SYM
df_model <- df_full %>% filter(model == "ARD") 
qqnorm(df_model$AIC_score, main='Normal')


#julia analysis
library(car)
library(tidyverse)

# verify equality of variances
lav_test <- leveneTest(AIC_score ~ model,
                       data = df_full,
)
# The p-value being larger than the significance level of 0.05, 
# we do not reject the null hypothesis, so we cannot reject the hypothesis 
# that variances are equal between groups (p-value = 0.4889499).

# normality (ANOVA is required to check the normality assumption)
res_aov <- aov(AIC_score ~ model,
               data = df_full
)

hist(res_aov$residuals)
qqPlot(res_aov$residuals,
       id = FALSE # id = FALSE to remove point identification
)
shapiro.test(res_aov$residuals)
# normality is met, one-way ANOVA can be used

#perform the one-way ANOVA, gives the global p-value
oneway.test(AIC_score ~ model,
            data = df_full,
            var.equal = TRUE
)

#p-value = < 2.2e-16 so the AIC_score is a function of the Mk model

#now we can see how the means of each group (mean AIC_scores of each Mk model) compare
#which have a statistically significant difference?

# Tukey's post-hoc test to study each pair of treatment
model=lm(df_full$AIC_score ~ df_full$model)
ANOVA=aov(model)

# Tukey test 
TUKEY <- TukeyHSD(x=ANOVA, 'df_full$model', conf.level=0.95)

# Tukey test visual representation :
#any comparison that doesn't have a 95% CI through zero is statistically significant
plot(TUKEY , las=1 , col="brown", cex.axis = 0.8)

# to see the difference in means and their p-values numerically:
TUKEY[["df_full$model"]]


#write this as a function:
anovaAssumptions <- function(results_df = df_full, metric = AIC_scores) {
  # takes the dataframe of all model results, determines if it means the assumptions of ANOVA
  # verify equality of variances
  lav_test <- leveneTest(metric ~ model,
                         data = results_df,
  )
  #if the p value is less than 0.05 the residuals fit the assumption of equal variance
  if(av_test$`Pr(>F)` < 0.05){
    residal_variance <- "Equal-variance"  
  } else{
    #error variance not equal
    residual_variance <- "Not-equal-variance"
    }
  
  #run an ANOVA to check if the residuals are normally distributed
  res_aov <- aov(metric ~ model,
                 data = results_df
  )
  #check for normal distribution with the shapiro test
  norm_test <- shapiro.test(res_aov$residuals)
  if(norm_test$p.value < 0.05){
    normality <- "normal"
  } else{
    normality <- "not-normal"
  }
  
  model=lm(results_df$metric ~ df_full$model )
  ANOVA=aov(model)
  
  # Tukey test 5o compare means between test groups
  TUKEY <- TukeyHSD(x=ANOVA, 'df$model', conf.level=0.95)

  # to see the difference in means and their p-values numerically:
  tukey_test <- TUKEY[["df$model"]]

  results_list <- c(residual_variance, normality, tukey_test)
  return(results_list)
}