#################################################################
# Molly McDermott
# created 9/16/21
# processes field data, and analyzed physiological and reproductive outcomes associated with brood size manipulation and range size
####################################################################

#### SET UP ####
#load libraries
library(tidyr)
library(ggplot2)
library(lme4)
library(lmerTest)
library(chron)
library(plyr)
library(dplyr)
library(predictmeans)
library(MuMIn)

# load physiology dataset
phys_wide <- read.csv("BSM_FemalePhys_wide.csv")

#set variable structure
phys_wide$Nest2 <- as.factor(phys_wide$Nest2)
phys_wide$Tagged <- as.factor(phys_wide$Tagged)
phys_wide$Site <- as.factor(phys_wide$Site)
phys_wide$Year <- as.factor(phys_wide$Year)
phys_wide$Time_1 <- chron(times = phys_wide$Time_1)
phys_wide$Time_2 <- chron(times = phys_wide$Time_2)
phys_wide$CI_1_j <- as.integer(format(as.Date(phys_wide$CI_1, "%m/%d/%Y"), "%j"))
phys_wide$CI_2_j <- as.integer(format(as.Date(phys_wide$CI_2, "%m/%d/%Y"), "%j"))


#reorder levels so model summaries use reduced brood as base
#phys_wide$BroodTrt <- factor(phys_wide$BroodTrt)
phys_wide$BroodTrt <- factor(phys_wide$BroodTrt, levels = c('R', 'E'))


#drop bird that was tagged but not part of BSM (2850-57510),  bird that was a male (2640-97570), and for some analyses, bird that was recaptured during second brood (2850-57445)
phys_wide <- phys_wide %>%
  filter(Band != "2850-57510" & Band != "2640-97570") 

#create Treatment variable for plotting
phys_wide$Treatment <- as.factor(paste(phys_wide$Tagged, phys_wide$BroodTrt, sep = "_"))
summary(phys_wide$Treatment)


# load GPS dataset
gps <- read.csv("BARS_GPS_area.csv")

#set variable structure
gps$Band <- as.factor(gps$band)


# is area related to number of fixes?
eq <- gps %>%
  mutate(area_km2 = foraging_area100/1000000) %>%
  filter(time == 'After')

fix <- lm(area_km2 ~ points, data = eq)
summary(fix)
confint(fix)

#reformat dataset so that each bird corresponds to 1 row of data
#for change/outcome models - feather growth, second brood, ICI, #eggs
df_gps_change <- gps %>%
  filter(tag == "Y") %>%
  group_by(Band, time) %>%
  summarize(area = mean(foraging_area100, na.rm = TRUE)) %>%
  pivot_wider(id_cols = Band,
              names_from = time,
              values_from = area) %>%
  mutate(area_diff = After - Before) %>%
  left_join(phys_wide, by = "Band")

#convert area from m2 to km2
df_gps_change$area_km2 <- df_gps_change$After / 1000000
summary(df_gps_change$area_km2)


# outcome models
# calculate number of days between captures
df_gps_change$Day_diff <- df_gps_change$ChickAge_2 - df_gps_change$ChickAge_1

# calculate feather growth rate
df_gps_change$R1_mm_day <- df_gps_change$R1_mm_2/df_gps_change$Day_diff

# filter feather growth dataset to only those with legitimate measurements
df_gps_R1 <- filter(df_gps_change, R1_coll_1 == "Y" & Band != "2850-57445")

#format time of day variables
df_gps_change$Time_2 <- chron(times = df_gps_change$Time_2)

#linearize time by taking sin and cos
df_gps_change$tsin_bg <- sin(2*pi*as.numeric(df_gps_change$Time_2))
df_gps_change$tcos_bg <- cos(2*pi*as.numeric(df_gps_change$Time_2))



######### SUMMARY STATS ###############
a <- df_gps_change %>%
  group_by(Treatment) %>%
  summarise(across(area_km2,
                   .fns = list(mean = mean, std = sd, min = min, max = max), na.rm = TRUE,
                   .names = "{col}_{fn}"))
a



############## MODELS ####################
# note: models with worse fit (usually those with an interaction term) are commented out

#FEATHER GROWTH RATE
#R1_mod <- lmer(R1_mm_day ~ scale(After)*BroodTrt + (1|Region), data = df_gps_R1)
#residplot(R1_mod)
R1_mod2 <- lmer(R1_mm_day ~ scale(After) + BroodTrt + (1|Region), data = df_gps_R1)
#residplot(R1_mod2) #residuals look fine; singular fit
#compare models with and without interaction term
#anova(R1_mod, R1_mod2)
summary(R1_mod2)
r.squaredGLMM(R1_mod2)

#SECOND NEST (Y/N)
#B2_mod <- glmer(Nest2 ~ scale(After)*BroodTrt + CI_1_j + (1|Region), family = "binomial", data = df_gps_change)
#residplot(B2_mod)
B2_mod2 <- glmer(Nest2 ~ scale(After) + BroodTrt + CI_1_j + (1|Region), family = "binomial", data = df_gps_change)
#residplot(B2_mod2)#residuals look fine; singular fit
#compare models with and without interaction term
#anova(B2_mod, B2_mod2)
summary(B2_mod2)
r.squaredGLMM(B2_mod2)

#NUMBER OF EGGS
egg_gps <- df_gps_change %>%
  filter(Nest2 == "Y")
#egg_mod <- glmer(Nest2Eggs ~ scale(After)*BroodTrt + (1|Region), family = "poisson", data = egg_gps)
#residplot(egg_mod)
egg_mod2 <- glmer(Nest2Eggs ~ scale(After) + BroodTrt+ (1|Region), family = "poisson", data = egg_gps)
#residplot(egg_mod2) # residuals are ok, singular fit
#compare models with and without interaction term
#anova(egg_mod, egg_mod2)
summary(egg_mod2) #consider dropping this entirely; so few data
r.squaredGLMM(egg_mod2)

# INTER-CLUTCH INTERVAL
df_gps_change$ICI <- df_gps_change$CI_2_j - df_gps_change$CI_1_j
#hist(df_gps_change$ICI)
#plot(ICI ~ BroodSize_1,  data = df_gps_change)
#plot(ICI ~ BroodSize_2,  data = df_gps_change)

ici_mod <- lmer(ICI ~ scale(After)*BroodTrt + (1|Region),  data = df_gps_change)
#residplot(ici_mod) #residuals look fine
#ici_mod2 <- lmer(ICI ~ scale(After) + BroodTrt + (1|Region), data = df_gps_change)
#residplot(ici_mod2)
#anova(ici_mod, ici_mod2) #compare models with and without interaction term
summary(ici_mod)
r.squaredGLMM(ici_mod)

#MASS
#ma_mod <- lmer(Mass_diff ~ scale(After)*BroodTrt + (1|Region), data = df_gps_change)
#residplot(ma_mod)
ma_mod2 <- lmer(Mass_diff ~ scale(After) + BroodTrt + (1|Region), data = df_gps_change)
#residplot(ma_mod2) #small curve in both fitted and random residuals
#anova(ma_mod, ma_mod2) #compare models with and without interaction term
summary(ma_mod2)
r.squaredGLMM(ma_mod2)

#BLOOD GLUCOSE
#bg_mod <- lmer(BloodGluc_2 ~ scale(After)*BroodTrt + scale(Mass_2) + tsin_bg + tcos_bg + (1|Region), data = df_gps_change)
#residplot(bg_mod)
bg_mod2 <- lmer(BloodGluc_2 ~ scale(After) + BroodTrt + scale(Mass_2) + tsin_bg + tcos_bg + (1|Region), data = df_gps_change)
#residplot(bg_mod2)
#anova(bg_mod, bg_mod2) #compare models with and without interaction term
summary(bg_mod2)
r.squaredGLMM(bg_mod2)


# H:L RATIO
HL_gps_df <- filter(df_gps_change, Band != "2640-97500")

#hl_mod <- lmer(HL_diff ~ scale(After)*BroodTrt + (1|Region), data = HL_gps_df)
#residplot(hl_mod)
hl_mod2 <- lmer(HL_diff ~ scale(After) + BroodTrt + (1|Region), data = HL_gps_df)
#residplot(hl_mod2) #much better without outlier
#anova(hl_mod, hl_mod2) #compare models with and without interaction term
summary(hl_mod2)
r.squaredGLMM(hl_mod2)




#### PLOTS ####
library(ggpubr)

#plot effect of treatment on mass
ma <- ggplot(df_gps_change, aes(x=area_km2, y=Mass_diff, color=BroodTrt)) + 
  geom_point(size = 3, alpha = 0.75) +
  xlab("") +
  ylab("Difference in mass (g)") +
  scale_color_manual(values=c("#999999","#111111"),
                     name="", # Legend label, use darker colors
                     breaks=c("E", "R"),
                     labels=c("Brood enlarged", "Brood reduced")) +
  theme_bw() +
  theme(legend.position = c(0.75,0.85),
        text = element_text(size=12),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

#blood glucose
bg <- ggplot(df_gps_change, aes(x=area_km2, y=BloodGluc_2, color=BroodTrt)) + 
  geom_point(size = 3, alpha = 0.75) +
  xlab("Range size (km2)") +
  ylab("Blood Glucose (mg/dL)") +
  scale_color_manual(values=c("#999999","#111111"),
                     name="", # Legend label, use darker colors
                     breaks=c("E", "R"),
                     labels=c("Brood enlarged", "Brood reduced")) +
  theme_bw() +
  theme(legend.position = c(0.75,0.85),
                    text = element_text(size=12),
                    axis.line = element_blank(),
                    axis.ticks = element_blank(),
                    panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank())

#HL ratio
hl <- ggplot(HL_gps_df, aes(x=area_km2, y=HL_diff, color=BroodTrt)) + 
  geom_point(size = 3, alpha = 0.75) +
  xlab("") +
  ylab("Difference in H:L ratio") +
  scale_color_manual(values=c("#999999","#111111"),
                     name="", # Legend label, use darker colors
                     breaks=c("E", "R"),
                     labels=c("Brood enlarged", "Brood reduced")) +
  theme_bw() +
  theme(legend.position = c(0.75,0.85),
        text = element_text(size=12),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())


# outcome models
#feather growth
#range size was significant, add model prediction line
R1_mod2 <- lm(R1_mm_day ~ scale(area_km2) + BroodTrt, data = df_gps_R1)

df_plot <- df_gps_R1 %>%
  select(After, BroodTrt, R1_mm_day, area_km2) %>%
  filter(!is.na(After), !is.na(R1_mm_day))
R1.predict <- data.frame(predict(R1_mod2, interval = 'prediction'))

plot_gps_R1 <- cbind(df_plot, R1.predict)


r1 <- ggplot(plot_gps_R1, aes(x=area_km2, y=R1_mm_day, color=BroodTrt)) + 
  geom_point(size = 3, alpha = 0.75) +
  geom_line(aes(y = fit)) +
  geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.1) +
  xlab("Range size (km2)") +
  ylab("Feather growth (mm/day)") +
  scale_color_manual(values=c("#999999","#111111"),
                     name="", # Legend label, use darker colors
                     breaks=c("E", "R"),
                     labels=c("Brood enlarged", "Brood reduced")) +
  theme_bw() +
  theme(legend.position = c(0.75,0.85),
        text = element_text(size=12),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())


pdf(file = "range_mass_bw.pdf", width = 2.5, height = 3.5)
ma + xlab("Range size (km2)")
dev.off()


pdf(file = "range_hl_bw.pdf", width = 2.5, height = 3.5)
hl + xlab("Range size (km2)")
dev.off()


pdf(file = "range_bg_bw.pdf", width = 2.5, height = 3.5)
bg
dev.off()


pdf(file = "range_r1_bw.pdf", width = 2.5, height = 3.5)
r1
dev.off()


pdf(file = "range_phys.pdf", width = 6.5, height = 5)
ggarrange(ma, hl, bg, r1,
          ncol = 2, nrow = 2, common.legend = TRUE)
dev.off()






#second brood
plot_nest2 <- filter(df_gps_change, Nest2 != "NA")
n2 <- ggplot(plot_nest2, aes(x=area_km2, y=Nest2, color=BroodTrt)) + 
  geom_point(size = 3, alpha = 0.75) +
  xlab("") +
  ylab("Second brood") +
  scale_color_manual(values=c("#999999","#111111"),
                     name="", # Legend label, use darker colors
                     breaks=c("E", "R"),
                     labels=c("Brood enlarged", "Brood reduced")) +
  theme_bw() +
  theme(legend.position = c(0.75,0.85),
        text = element_text(size=12),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())



#number of eggs
eg <- ggplot(df_gps_change, aes(x=area_km2, y=Nest2Eggs, color=BroodTrt)) + 
  geom_point(size = 3, alpha = 0.75) +
  xlab("") +
  ylab("Number of eggs in second brood") +
  scale_color_manual(values=c("#999999","#111111"),
                     name="", # Legend label, use darker colors
                     breaks=c("E", "R"),
                     labels=c("Brood enlarged", "Brood reduced")) +
  ylim(2.5, 5.5)+
  theme_bw() +
  theme(legend.position = c(0.75,0.85),
        text = element_text(size=12),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())


#ICI
#range size was significant, add model line
ici_lm <- lm(ICI ~ scale(area_km2)*BroodTrt,  data = df_gps_change) #for now, leaving out random effect (main effects are virtually unchanged in lm vs lmer)

df_ici_plot <- df_gps_change %>%
  select(After, BroodTrt, ICI, area_km2, Region) %>%
  filter(!is.na(After), !is.na(ICI))

df_ici_lm<- data.frame(predict(ici_lm, interval = 'prediction'))

plot_gps_ici <- cbind(df_ici_plot, df_ici_lm)

ic <- ggplot(data = plot_gps_ici, 
       aes(x=area_km2, y=ICI, color=BroodTrt)) + 
  geom_point(size = 3, alpha = 0.75) +
  geom_line(aes(y = fit)) +
  geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.1) +
  xlab("Foraging area (km2)") +
  ylab("Inter-clutch interval") +
  scale_color_manual(values=c("#999999","#111111"),
                     name="", # Legend label, use darker colors
                     breaks=c("E", "R"),
                     labels=c("Brood enlarged", "Brood reduced")) +
  theme_bw() +
  theme(legend.position = c(0.75,0.85),
        text = element_text(size=12),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())



pdf(file = "range_rep.pdf", width = 3.25, height = 7.5)
ggarrange(n2, eg, ic, ncol = 1, common.legend = TRUE)
dev.off()

pdf(file = "range_n2_bw.pdf", width = 2.5, height = 3.5)
n2
dev.off()

pdf(file = "range_egg_bw.pdf", width = 2.5, height = 3.5)
eg
dev.off()

pdf(file = "range_ICI_bw.pdf", width = 2.5, height = 3.5)
ic
dev.off()

