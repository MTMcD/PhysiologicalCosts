#################################################################
# Molly McDermott
# created 1/11/21
# edited 4/22/22
# Analysis of BSM field experiment - female physiology
################################################################

#### SET UP ####

#load libraries
library(dplyr)
library(tidyr)
library(ggplot2)
library(lme4)
library(car)
library(predictmeans)
library(lmerTest)
library(chron)
library(MuMIn)
library(skimr)

#read in data file
phys_wide <- read.csv("BSM_FemalePhys_wide.csv")

#set variable structures
phys_wide$Nest2 <- as.factor(phys_wide$Nest2)
phys_wide$Tagged <- as.factor(phys_wide$Tagged)
phys_wide$Site <- as.factor(phys_wide$Site)
phys_wide$SiteSize <- as.factor(phys_wide$SiteSize)
phys_wide$Region <- as.factor(phys_wide$Region)
phys_wide$Year <- as.factor(phys_wide$Year)
phys_wide$Time_1 <- chron(times = phys_wide$Time_1)
phys_wide$Time_2 <- chron(times = phys_wide$Time_2)
phys_wide$CI_1_j <- as.integer(format(as.Date(phys_wide$CI_1, "%m/%d/%Y"), "%j"))
phys_wide$CI_2_j <- as.integer(format(as.Date(phys_wide$CI_2, "%m/%d/%Y"), "%j"))

#reorder levels so model summaries use reduced brood as base
#phys_wide$BroodTrt <- factor(phys_wide$BroodTrt)
phys_wide$BroodTrt <- factor(phys_wide$BroodTrt, levels = c('R', 'E'))

#create treatment variable for plotting and summarizing data
phys_wide$Treatment <- as.factor(paste(phys_wide$Tagged, phys_wide$BroodTrt, sep = "_"))
summary(phys_wide$Treatment)

#exclude bird that was recaptured during second brood attempt for physiology analyses
phys_wide_2 <- phys_wide %>% 
  filter(Band != "2850-57445")

# exclude birds that did not re-lay from egg analysis
egg_df <- phys_wide %>%
  filter(Nest2 == "Y")

#create dataset for legitimate measures of feather growth (excludes birds whose feathers were collected prior to the start of the experiment or not collected at all)
phys_wide_R1 <- filter(phys_wide, R1_coll_1 == "Y")


#### DIFFERENCES AT BASELINE MEASURE ####

#calculate number of days between captures for each treatment group
phys_wide_2 %>% 
  group_by(Treatment) %>%
  summarize(duration = mean(Day_diff, na.rm = TRUE),
            std = sd(Day_diff, na.rm = TRUE),
            max = max(Day_diff),
            min = min(Day_diff))

#model number of days between captures, using same model structure as outcomes of interest
dmod <- lmer(Day_diff ~ Tagged + BroodTrt + Year + (1|Region), data = phys_wide_2)
residplot(dmod)
summary(dmod)

#summarize mass for each treatment group
phys_wide_2 %>% 
  group_by(Treatment) %>%
  summarize(mass_mean = mean(Mass_1, na.rm = TRUE),
            mass_sd = sd(Mass_1, na.rm = TRUE),
            mass_max = max(Mass_1),
            mass_min = min(Mass_1))

#model mass using same model structure as outcomes of interest
mass_start <- lmer(Mass_1 ~ Tagged + BroodTrt + Year + (1|Region), data = phys_wide_2)
residplot(mass_start)
summary(mass_start) 

#summarize blood glucose at first capture for each treatment group
phys_wide_2 %>% 
  group_by(Treatment) %>%
  summarize(BloodGluc_mean = mean(BloodGluc_1, na.rm = TRUE),
            BloodGluc_sd = sd(BloodGluc_1, na.rm = TRUE),
            BloodGluc_max = max(BloodGluc_1),
            BloodGluc_min = min(BloodGluc_1))

#model blood glucose at first capture using same model structure as outcomes of interest
BloodGluc_start <- lmer(BloodGluc_1 ~ Tagged + BroodTrt + Year + (1|Region), data = phys_wide_2)
residplot(BloodGluc_start)
summary(BloodGluc_start)

#summarize H:L ratio at first capture for each treatment group
phys_wide_2 %>% 
  group_by(Treatment) %>%
  summarize(HL_ratio_1_mean = mean(HL_ratio_1, na.rm = TRUE),
            HL_ratio_1_sd = sd(HL_ratio_1, na.rm = TRUE),
            HL_ratio_1_max = max(HL_ratio_1),
            HL_ratio_1_min = min(HL_ratio_1))

#model H:L ratio at first capture using same model structure as outcomes of interest
HL_start <- lmer(HL_ratio_1 ~ Tagged + BroodTrt + Year + (1|Region), data = phys_wide_2)
residplot(HL_start)
summary(HL_start)


#turn dates and times into chron object, sin and cos to linearize time
time2 <- chron(times = phys_wide_2$Time_2)
tsin <- sin(2*pi*as.numeric(time2))
tcos <- cos(2*pi*as.numeric(time2))

#visualize relationship between time of day and blood glucose
Sys.setenv(TZ='GMT')
ggplot(phys_wide, aes(Time_2, BloodGluc_2)) +
  geom_point() +
  theme_bw() +
  scale_x_chron(format = "%H:%M", n = 6) +
  theme(text = element_text(size=16)) +
  xlab("Time of Capture") +
  ylab("Blood Glucose (mg/dL)")


############### summary stats ##############################

# summarize physiological variables by treatment group
a <- phys_wide_2 %>%
  group_by(Treatment) %>%
  summarise(across(c(BroodSize_1, BroodSize_2, Mass_1, Mass_2, BloodGluc_1, BloodGluc_2, 
                    HL_ratio_1, HL_ratio_2, Day_diff),
                   .fns = list(mean = mean, std = sd, min = min, max = max), na.rm = TRUE,
                   .names = "{col}_{fn}"))

#summarize feather growth by treatment group
b <- phys_wide_R1 %>%
  group_by(Treatment) %>%
  summarise(across(R1_mm_day, 
                   .fns = list(mean = mean, std = sd, min = min, max = max), na.rm= TRUE,
                   .names = "R1_{fn}"))

#summarize number of eggs by treatment group
c <- egg_df %>%
  group_by(Treatment) %>%
  summarise(across(Nest2Eggs, 
                   .fns = list(mean = mean, std = sd, min = min, max = max), na.rm= TRUE,
                   .names = "eggs_{fn}"))

#summarize inter-clutch interval by treatment group
d <- ICI_df %>%
  group_by(Treatment) %>%
  summarise(across(ICI, 
                   .fns = list(mean = mean, std = sd, min = min, max = max), na.rm= TRUE,
                   .names = "ICI_{fn}"))

sm <- cbind(a,b,c,d)
#save summary data to reporst in ESM
write.csv(sm, "SummaryStatsBSM.csv")

#### MODEL FITTING ####
# note that the worse-fitting model (usually with an interaction term) is commented out

# MASS
#mass_mod <- lmer(Mass_diff ~ Tagged*BroodTrt + Year + (1|Region), data = phys_wide_2)
#residplot(mass_mod)
mass_mod2 <- lmer(Mass_diff ~ Tagged + BroodTrt + Year + (1|Region), data = phys_wide_2)
#residplot(mass_mod2)
#anova(mass_mod, mass_mod2) #compare model with and without an interaction term
summary(mass_mod2)
r.squaredGLMM(mass_mod2)

# FEATHER GROWTH
#R1_mod <- lmer(R1_mm_day ~Tagged*BroodTrt + Year + (1|Region), data = phys_wide_R1)
#residplot(R1_mod)
R1_mod2 <- lmer(R1_mm_day ~Tagged + BroodTrt + Year + (1|Region), data = phys_wide_R1)
#residplot(R1_mod2)
#anova(R1_mod, R1_mod2) #compare model with and without an interaction term
summary(R1_mod2)
r.squaredGLMM(R1_mod2)

# BLOOD GLUCOSE
#gluc_mod <- lmer(BloodGluc_2 ~ Tagged*BroodTrt + scale(Mass_2) + tsin + tcos + Year + (1|Region), data = phys_wide_2)
#residplot(gluc_mod)
gluc_mod2 <- lmer(BloodGluc_2 ~ Tagged + BroodTrt + scale(Mass_2) + tsin + tcos + Year + (1|Region), data = phys_wide_2)
#residplot(gluc_mod2)
#anova(gluc_mod, gluc_mod2) #compare model with and without an interaction term
summary(gluc_mod2)
r.squaredGLMM(gluc_mod2)

# H:L RATIO
#remove outlier (2640-97500) before running this model
HL_df <- phys_wide_2 %>%
  filter(Band != "2640-97500") %>%
  mutate(HL_diff = HL_ratio_2-HL_ratio_1)

#HL_mod <- lmer(HL_diff ~ Tagged*BroodTrt + Year + (1|Region), data = HL_df)
#residplot(HL_mod)
HL_mod2 <- lmer(HL_diff ~ Tagged + BroodTrt + Year + (1|Region), data = HL_df)
#residplot(HL_mod2)
#anova(HL_mod, HL_mod2) #compare model with and without an interaction term
summary(HL_mod2)
r.squaredGLMM(HL_mod2)

# SECOND NEST (Y/N)
#B2_mod <- glmer(Nest2 ~ Tagged*BroodTrt + Year + CI_1_j + (1|Region), family = "binomial", data = phys_wide)
#residplot(B2_mod)
B2_mod2 <- glmer(Nest2 ~ Tagged + BroodTrt + Year + CI_1_j + (1|Region), family = "binomial", data = phys_wide)
#residplot(B2_mod2)
#anova(B2_mod, B2_mod2) #compare model with and without an interaction term
summary(B2_mod2)
r.squaredGLMM(B2_mod2)

#check for overdispersion
overdisp_fun <- function(model) {
  rdf <- df.residual(model)
  rp <- residuals(model,type="pearson")
  Pearson.chisq <- sum(rp^2)
  prat <- Pearson.chisq/rdf
  pval <- pchisq(Pearson.chisq, df=rdf, lower.tail=FALSE)
  c(chisq=Pearson.chisq,ratio=prat,rdf=rdf,p=pval)
}
#overdisp_fun(B2_mod2)


# NUMBER OF EGGS
#egg_mod <- glmer(Nest2Eggs ~ Tagged*BroodTrt + Year + (1|Region), family = "poisson", data = egg_df)
#residplot(egg_mod)
#overdisp_fun(egg_mod)
egg_mod2 <- glmer(Nest2Eggs ~ Tagged + BroodTrt + Year + (1|Region), family = "poisson", data = egg_df)
#residplot(egg_mod2)
#overdisp_fun(egg_mod2)
#anova(egg_mod, egg_mod2)#compare model with and without an interaction term
summary(egg_mod2)
r.squaredGLMM(egg_mod2)

# INTER-CLUTCH INTERVAL
# exclude outlier
ICI_df <- phys_wide %>%
  filter(Band != "2640-97598")

#ici_mod <- lmer(ICI ~ Tagged*BroodTrt + Year + (1|Region),  data = ICI_df)
#residplot(ici_mod)
ici_mod2 <- lmer(ICI ~ Tagged + BroodTrt + Year + (1|Region), data = ICI_df)
#residplot(ici_mod2)
#anova(ici_mod, ici_mod2) #compare model with and without interaction term
summary(ici_mod2)
r.squaredGLMM(ici_mod2)




#### MARGINAL MEANS ####

library ('emmeans')

#create labels for x-axis
trt_labs <- c("No tag","No tag", "Tag", "Tag")


#estimate marginal means for mass
mass_all_mean <- as.data.frame(emmeans(mass_mod2, specs = c("Tagged","BroodTrt")))
mass_all_mean$Treatment <- as.factor(paste(mass_all_mean$Tagged, mass_all_mean$BroodTrt, sep = "_"))
mass_all_mean


#plot effect of treatment on mass with marginal means and confidence intervals
pdf(file = "mass_emmean_bw.pdf", width = 3, height = 4)
ggplot() + 
  geom_point(data = phys_wide_2, 
             aes(x=Treatment, y=Mass_diff, fill = BroodTrt),
             shape = 21, size = 5, alpha = 0.75) +
  geom_point(data = mass_all_mean, 
             aes(x=Treatment, y=emmean), 
             shape = 23, size = 8) +
  geom_errorbar(data = mass_all_mean, 
                aes(x=Treatment, y=emmean, ymin = lower.CL, ymax = upper.CL), 
                width = 0.1) +
  geom_hline(yintercept = 0) +
  theme_bw() +
  theme(text = element_text(size=12),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  xlab("Treatment") + 
  ylab("Within-individual difference in mass (g)") +
  scale_fill_manual(values=c("#999999", "#111111"),
                      name="", # Legend label, use darker colors
                      breaks=c("E", "R"),
                      labels=c("Brood enlarged", "Brood reduced")) +
  theme(legend.position = c(0.85,0.9))+
  ylim(-2.5, 2)+
  scale_x_discrete(labels = trt_labs) 
dev.off()




# estimate marginal means for feather growth
R1_all_mean <- as.data.frame(emmeans(R1_mod2, specs = c("Tagged","BroodTrt")))
R1_all_mean$Treatment <- as.factor(paste(R1_all_mean$Tagged, R1_all_mean$BroodTrt, sep = "_"))
R1_all_mean

#plot effect of treatment on feather growth with marginal means and confidence intervals
pdf(file = "R1_emmean_bw.pdf", width = 3, height = 4)

ggplot() + 
  geom_point(data = phys_wide_R1, 
             aes(x=Treatment, y=R1_mm_day, fill = BroodTrt),
             shape = 21, size = 5, alpha = 0.75) +
  geom_point(data = R1_all_mean, 
             aes(x=Treatment, y=emmean), 
             size = 8, shape = 23) +
  geom_errorbar(data = R1_all_mean, 
                aes(x=Treatment, y=emmean, ymin = lower.CL, ymax = upper.CL), 
                width = 0.1) +
  geom_hline(yintercept = 0) +
  theme_bw() +
  theme(text = element_text(size=12),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  xlab("Treatment") + 
  ylab("Feather growth per day (mm)") +
  scale_fill_manual(values=c("#999999", "#111111"),
                      name="", # Legend label, use darker colors
                      breaks=c("E", "R"),
                      labels=c("Brood enlarged", "Brood reduced")) +
  theme(legend.position = c(0.85,0.9))+
  ylim(-0.1, 1)+
  scale_x_discrete(labels = trt_labs) 

dev.off()




# estimate marginal means for blood glucose
bg_all_mean <- as.data.frame(emmeans(gluc_mod2, specs = c("Tagged","BroodTrt")))
bg_all_mean$Treatment <- as.factor(paste(bg_all_mean$Tagged, bg_all_mean$BroodTrt, sep = "_"))
bg_all_mean

#plot effect of treatment on feather growth with marginal means and confidence intervals
pdf(file = "bg_emmean_bw.pdf", width = 3, height = 4)

ggplot() + 
  geom_point(data = phys_wide_2, 
             aes(x=Treatment, y=BloodGluc_2, fill = BroodTrt),
             shape = 21, size = 5, alpha = 0.75) +
  geom_point(data = bg_all_mean, 
             aes(x=Treatment, y=emmean), 
             size = 8, shape = 23) +
  geom_errorbar(data = bg_all_mean, 
                aes(x=Treatment, y=emmean, ymin = lower.CL, ymax = upper.CL), 
                width = 0.1) +
  theme_bw() +
  theme(text = element_text(size=12),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  xlab("Treatment") + 
  ylab("Blood glucose (mg/dL)") +
  scale_fill_manual(values=c("#999999", "#111111"),
                      name="", # Legend label, use darker colors
                      breaks=c("E", "R"),
                      labels=c("Brood enlarged", "Brood reduced")) +
  theme(legend.position = c(0.2,0.9))+
  ylim(150, 350)+
  scale_x_discrete(labels = trt_labs) 

dev.off()





# estimate marginal means for H:L ratio
hl_all_mean <- as.data.frame(emmeans(HL_mod2, specs = c("Tagged","BroodTrt")))
hl_all_mean$Treatment <- as.factor(paste(hl_all_mean$Tagged, hl_all_mean$BroodTrt, sep = "_"))
hl_all_mean

#reorder levels to arrange graph from low to high cost treatment
HL_df$Treatment <- factor(HL_df$Treatment, levels = c("N_R", "N_E", "Y_R", "Y_E"))
hl_all_mean$Treatment <- factor(hl_all_mean$Treatment, levels = c("N_R", "N_E", "Y_R", "Y_E"))

#plot effect of treatment on feather growth with marginal means and confidence intervals
#pdf(file = "hl_emmean_bw.pdf", width = 3, height = 4)
pdf(file = "hl_emmean_bw_exittalk.pdf", width = 3, height = 3)

ggplot() + 
  geom_point(data = HL_df, 
             aes(x=Treatment, y=HL_diff, fill = BroodTrt),
             shape = 21, size = 5, alpha = 0.75) +
  geom_point(data = hl_all_mean, 
             aes(x=Treatment, y=emmean), 
             size = 8, shape = 23) +
  geom_errorbar(data = hl_all_mean, 
                aes(x=Treatment, y=emmean, ymin = lower.CL, ymax = upper.CL), 
                width = 0.1) +
  geom_hline(yintercept = 0) +
  theme_bw() +
  theme(text = element_text(size=12),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  xlab("") + 
  ylab("Difference in H:L ratio") +
  scale_fill_manual(values=c("#999999", "#111111"),
                      name="", # Legend label, use darker colors
                      breaks=c("E", "R"),
                      labels=c("Brood enlarged", "Brood reduced")) +
  theme(legend.position = "none")+
  ylim(1.75, -1.75) +
  scale_x_discrete(labels = trt_labs) 

dev.off()


# graph of percent females who laid a second brood


nest_df <- phys_wide %>% 
  filter(!is.na(Nest2)) %>%
  group_by(Treatment, Tagged, BroodTrt) %>% 
  summarize(per = table(Nest2)/length(Nest2)*100,
            val = names(table(Nest2))) %>%
  filter(val == "Y")
  
pdf(file = "nest2_bw.pdf", width = 3, height = 4)

ggplot(nest_df, aes(x = Treatment, y = per, fill= BroodTrt)) +
  geom_col(alpha = 0.75) +
  xlab("Treatment") + 
  ylab("Second brood (%)") +
  theme_bw() +
  theme(text = element_text(size=12),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  scale_fill_manual(values=c("#999999", "#111111"),
                    name="", # Legend label, use darker colors
                    breaks=c("E", "R"),
                    labels=c("Brood enlarged", "Brood reduced")) +
  theme(legend.position = c(0.8,0.9))+
  ylim(0, 100) +
  scale_x_discrete(labels = trt_labs) 

dev.off()




# estimate marginal means for num eggs
egg_all_mean <- as.data.frame(emmeans(egg_mod2, specs = c("Tagged","BroodTrt")))
egg_all_mean$Treatment <- as.factor(paste(egg_all_mean$Tagged, egg_all_mean$BroodTrt, sep = "_"))
egg_all_mean

#reorder factor levels
egg_df$Treatment <- factor(egg_df$Treatment, levels = c("N_R", "N_E", "Y_R", "Y_E"))
egg_all_mean$Treatment <- factor(egg_all_mean$Treatment, levels = c("N_R", "N_E", "Y_R", "Y_E"))


#plot effect of treatment on num eggs with marginal means and confidence intervals
pdf(file = "egg_emmean_bw.pdf", width = 3, height = 4)

ggplot() + 
  geom_jitter(data = egg_df, 
              aes(x=Treatment, y=Nest2Eggs, fill = BroodTrt),
              shape = 21, size = 5, alpha = 0.75, width = 0.1, height = 0.15) +
  geom_point(data = egg_all_mean, 
             aes(x=Treatment, y=exp(emmean)), 
             size = 8, shape = 23) +
  geom_errorbar(data = egg_all_mean, 
                aes(x=Treatment, y=exp(emmean), ymin = exp(asymp.LCL), ymax = exp(asymp.UCL)), 
                width = 0.1) +
  theme_bw() +
  theme(text = element_text(size=12),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  xlab("") + 
  ylab("Eggs in second brood") +
  scale_fill_manual(values=c("#999999", "#111111"),
                      name="", # Legend label, use darker colors
                      breaks=c("E", "R"),
                      labels=c("Brood enlarged", "Brood reduced")) +
  theme(legend.position = "none") +
  ylim(2.5, 7) +
  scale_x_discrete(labels = trt_labs) 

dev.off()




# estimate marginal means for H:L ratio
ici_all_mean <- as.data.frame(emmeans(ici_mod2, specs = c("Tagged","BroodTrt")))
ici_all_mean$Treatment <- as.factor(paste(ici_all_mean$Tagged, ici_all_mean$BroodTrt, sep = "_"))
ici_all_mean

#plot effect of treatment on feather growth with marginal means and confidence intervals
pdf(file = "ICI_emmean_bw.pdf", width = 3, height = 4)

ggplot() + 
  geom_point(data = ICI_df, 
             aes(x=Treatment, y=ICI, fill = BroodTrt),
             shape = 21, size = 5, alpha = 0.75) +
  geom_point(data = ici_all_mean, 
             aes(x=Treatment, y=emmean), 
             size = 8, shape = 23) +
  geom_errorbar(data = ici_all_mean, 
                aes(x=Treatment, y=emmean, ymin = lower.CL, ymax = upper.CL), 
                width = 0.1) +
  theme_bw() +
  theme(text = element_text(size=12),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  xlab("Treatment") + 
  ylab("Inter-clutch interval (days)") +
  scale_fill_manual(values=c("#999999", "#111111"),
                    name="", # Legend label, use darker colors
                    breaks=c("E", "R"),
                    labels=c("Brood enlarged", "Brood reduced")) +
  theme(legend.position = c(0.75,0.9))+
  ylim(42, 60) +
  scale_x_discrete(labels = trt_labs) 

dev.off()


