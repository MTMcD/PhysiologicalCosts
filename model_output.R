################################################################################
# Molly McDermott
# Created 9/3/21
# script tidies model output (lmer or glmer) and produces a table in excel and a figure with model estimates
################################################################################


#### SET UP ####
# load packages
library(lme4)
library(lmerTest)
library(chron)
library(plyr)
library(dplyr)
library(tidyr)
library(purrr)
library(ggplot2)

# read in main data file
phys_wide <- read.csv("BSM_FemalePhys_wide.csv")

# set variables as factors
phys_wide$Nest2 <- as.factor(phys_wide$Nest2)
phys_wide$Tagged <- as.factor(phys_wide$Tagged)
phys_wide$Site <- as.factor(phys_wide$Site)
phys_wide$SiteSize <- as.factor(phys_wide$SiteSize)
phys_wide$Region <- as.factor(phys_wide$Region)
phys_wide$Year <- as.factor(phys_wide$Year)
#convert time columns using chron package
phys_wide$Time_1 <- chron(times = phys_wide$Time_1)
phys_wide$Time_2 <- chron(times = phys_wide$Time_2)
#linearize time by taking sin and cos
phys_wide$tsin <- sin(2*pi*as.numeric(phys_wide$Time_2))
phys_wide$tcos <- cos(2*pi*as.numeric(phys_wide$Time_2))
#convert clutch initiation dates to julian day
phys_wide$CI_1_j <- as.integer(format(as.Date(phys_wide$CI_1, "%m/%d/%Y"), "%j"))
phys_wide$CI_2_j <- as.integer(format(as.Date(phys_wide$CI_2, "%m/%d/%Y"), "%j"))

#reorder levels so model summaries use reduced brood as base
#phys_wide$BroodTrt <- factor(phys_wide$BroodTrt)
phys_wide$BroodTrt <- factor(phys_wide$BroodTrt, levels = c('R', 'E'))

gps <- read.csv("BARS_GPS_area.csv")
phys <- read.csv("BSM_FemalePhys.csv")
phys$Year <- as.factor(format(as.Date(phys$Date, "%m/%d/%Y"), "%Y"))
phys$Tagged <- as.factor(phys$Tagged)
phys$Site <- as.factor(phys$Site)
phys$BroodTrt <- factor(phys$BroodTrt)
gps$Band <- as.factor(gps$band)
phys$Band <- as.factor(phys$Band)

sites <- read.csv("SiteKey.csv")
phys$SiteSize <- mapvalues(phys$Site, from = sites$Site, to = sites$SizeGroup)
phys$Region <- mapvalues(phys$Site, from = sites$Site, to = sites$Region)





#### CREATE DATA FRAMES ####

# Analyses of all birds
#exclude bird that was recaptured during second brood attempt for physiology analyses
phys_wide_2 <- phys_wide %>% 
  filter(Band != "2850-57445")

#create dataset for legitimate measures of feather growth (excludes birds whose feathers were collected prior to the start of the experiment or not collected at all)
phys_wide_R1 <- filter(phys_wide, R1_coll_1 == "Y")

#remove HL outlier (2640-97500), calculate HL ratio within-ind difference
HL_df <- phys_wide %>%
  filter(Band != "2640-97500") %>%
  mutate(HL_diff = HL_ratio_2-HL_ratio_1)

HL_df %>%
  group_by(Tagged) %>%
  summarise(mean = mean(HL_diff, na.rm = TRUE))

#exclude birds that did not re-lay from egg analysis
egg_df <- phys_wide %>%
  filter(Nest2 == "Y")

#remove outlier from ICI dataset
ICI_df <- phys_wide %>%
  filter(Band != "2640-97598")


# Analyses of tagged birds

#drop bird that was tagged but not part of BSM (2850-57510),  bird that was a male (2640-97570), and for some analyses, bird that was recaptured during second brood (2850-57445)
phys <- phys %>%
  mutate(capture = as.factor(ifelse(ChickAge < 8, "Before", "After"))) %>%
  filter(Band != "2640-97570" & Band != "2850-57510")

phys$BSM_trt <- "C"
phys[phys$capture == "After" & phys$BroodTrt == "E", "BSM_trt"] <- "E"
phys[phys$capture == "After" & phys$BroodTrt == "R", "BSM_trt"] <- "R"
phys$BSM_trt <- factor(phys$BSM_trt)

#for repeated measure models with all data - mass, HL ratio, blood glucose
df_gps <- gps %>%
  group_by(Band, time) %>%
  summarize(area = mean(foraging_area100, na.rm = TRUE)) %>%
  left_join(phys, by = c("Band", "time" = "capture")) %>%
  filter(Tagged == "Y" & Band != "2850-57445") #exclude bird recaptured during second brood attempt
summary(df_gps$BSM_trt)

#convert area from m2 to km2, check distribution for outlier
df_gps$area_km2 <- df_gps$area / 1000000
summary(df_gps$area_km2)
df_gps_noo <- filter(df_gps, area_km2 < 20)
df_gps_noo$Time <- chron(times = df_gps_noo$Time)
df_gps_noo$tsin_bg <- sin(2*pi*as.numeric(df_gps_noo$Time))
df_gps_noo$tcos_bg <- cos(2*pi*as.numeric(df_gps_noo$Time))


#for change/outcome models - feather growth, second brood, ICI, #eggs
df_gps_change <- gps %>%
  filter(tag == "Y") %>%
  group_by(Band, time)%>%
  summarize(area = mean(foraging_area100, na.rm = TRUE)) %>%
  pivot_wider(id_cols = Band,
              names_from = time,
              values_from = area) %>%
  mutate(area = After) %>%
  left_join(phys_wide, by = "Band")

df_gps_change$area_km2 <- df_gps_change$After / 1000000
df_gps_change$Time_2 <- chron(times = df_gps_change$Time_2)
df_gps_change$tsin_bg <- sin(2*pi*as.numeric(df_gps_change$Time_2))
df_gps_change$tcos_bg <- cos(2*pi*as.numeric(df_gps_change$Time_2))


# outcome models
df_gps_change$Day_diff <- df_gps_change$ChickAge_2 - df_gps_change$ChickAge_1
df_gps_change$R1_mm_day <- df_gps_change$R1_mm_2/df_gps_change$Day_diff

df_gps_R1 <- filter(df_gps_change, R1_coll_1 == "Y" & Band != "2850-57445")

HL_gps_df <- filter(df_gps_change, Band != "2640-97500")

egg_gps <- df_gps_change %>%
  filter(Nest2 == "Y")



#### RUN MODELS ####
# bind_cols_fill function, similar to cbind but allows columns of unequal length
bind_cols_fill <- function(df_list) {
  
  max_rows <- map_int(df_list, nrow) %>% max()
  
  map(df_list, function(df) {
    if(nrow(df) == max_rows) return(df)
    first <- names(df)[1] %>% sym()
    df %>% add_row(!!first := rep(NA, max_rows - nrow(df)))
  }) %>% bind_cols()
}

#list of models - main physiology models
modlist <- list(
              mod1 <- lmer(Mass_diff ~ Tagged + BroodTrt + Year + (1|Region), 
                           data = phys_wide_2),
              mod2 <- lmer(R1_mm_day ~Tagged + BroodTrt + Year + (1|Region), 
                           data = phys_wide_R1),
              mod3 <- lmer(BloodGluc_2 ~ Tagged + BroodTrt + scale(Mass_2) + tsin + tcos + Year + (1|Region), 
                           data = phys_wide_2),
              mod4 <- lmer(HL_diff ~ Tagged + BroodTrt + Year + (1|Region), 
                           data = HL_df),
              mod5 <- glmer(Nest2 ~ Tagged + BroodTrt + Year + scale(CI_1_j) + (1|Region), 
                            family = "binomial", 
                            data = phys_wide),
              mod6 <- glmer(Nest2Eggs ~ Tagged + BroodTrt + Year + (1|Region), 
                            family = "poisson", 
                            data = egg_df),
              mod7 <- lmer(ICI ~ Tagged + BroodTrt + Year + (1|Region), 
                           data = ICI_df)
)

# extract and rearrange summary statistics from models
for(i in 1:length(modlist)) {

  a <- modlist[[i]] # set focal model
  coeffs <- as.data.frame(coef(summary(a))) # extract estimates, SE, df, P
  rownames(coeffs) <- NULL # remove rownames from output
  re <- as.data.frame(VarCorr(a)) # extract variance and residual for random effects
  coeffsp <- bind_cols_fill(list(tibble((response = names(a@frame)[1])),
                                 tibble(fixed = rownames(coef(summary(a)))), 
                                 tibble(coeffs),
                                 tibble(re))) # bind together the info extracted from the model
  
  #to create a new file
  #write.csv(coeffsp, "coeffsp.csv") 
  
  #to add to an existing file
  #rename or else information will be added to old file
  write.table(coeffsp, "BSM_model_results_9_28.csv", 
              sep = ",", 
              append = TRUE, 
              row.names = FALSE) # export
}


#list of models - foraging area models
tagmods <- list(
  mod1 <- lmer(Mass_diff ~ scale(After) + BroodTrt + (1|Region), data = df_gps_change),
  mod2 <- lmer(BloodGluc_2 ~ scale(After) + BroodTrt + scale(Mass_2) + tsin_bg + tcos_bg +                      (1|Region), data = df_gps_change),
  mod3 <- lmer(HL_diff ~ scale(After) + BroodTrt + (1|Region), data = HL_gps_df),
  mod4 <- lmer(R1_mm_day ~ scale(After) + BroodTrt + (1|Region), 
               data = df_gps_R1),
  mod5 <- glmer(Nest2 ~ scale(After) + BroodTrt + scale(CI_1_j) + (1|Region), 
                family = "binomial", 
                data = df_gps_change),
  mod6 <- glmer(Nest2Eggs ~ scale(After) + BroodTrt + (1|Region), 
                family = "poisson", 
                data = egg_gps),
  mod7 <- lmer(ICI ~ scale(After)*BroodTrt + (1|Region),  
               data = df_gps_change)
)

for(i in 1:length(tagmods)) {

  a <- tagmods[[i]]
  coeffs <- as.data.frame(coef(summary(a))) # get estimates, etc...
  rownames(coeffs) <- NULL
  re <- as.data.frame(VarCorr(a))
  coeffsp <- bind_cols_fill(list(tibble((response = names(a@frame)[1])),
                               tibble(fixed = rownames(coef(summary(a)))), 
                              tibble(coeffs),
                              tibble(re))) 

#to create a new file
#write.csv(coeffsp, "coeffsp.csv") 
#to add to an existing file
  write.table(coeffsp, "BSM_foraging_results_12_13.csv", sep = ",", append = TRUE, row.names = FALSE) # export
}




################################################################################
#### COEFFICIENT PLOTS ####
## All birds ##
# Using code from Zach Laubach

library(ggpubr)


modlist <- list(
  mod1 <- lmer(Mass_diff ~ Tagged + BroodTrt + Year + (1|Region), 
               data = phys_wide_2),
  mod2 <- lmer(R1_mm_day ~Tagged + BroodTrt + Year + (1|Region), 
               data = phys_wide_R1),
  mod3 <- lmer(BloodGluc_2 ~ Tagged + BroodTrt + Mass_2 + tsin + tcos + Year + (1|Region), 
               data = phys_wide_2),
  mod4 <- lmer(HL_diff ~ Tagged + BroodTrt + Year + (1|Region), 
               data = HL_df),
  mod5 <- glmer(Nest2 ~ Tagged + BroodTrt + Year + CI_1_j + (1|Region), 
                family = "binomial", 
                data = phys_wide),
  mod6 <- glmer(Nest2Eggs ~ Tagged + BroodTrt + Year + (1|Region), 
                family = "poisson", 
                data = egg_df),
  mod7 <- lmer(ICI ~ Tagged + BroodTrt + Year + (1|Region), 
               data = ICI_df)
)


# initialize list of plots
plots <- list ()

#list new labels for plot panels and corresponding response variable names
mylabs <- c("Difference in mass (g)","Feather growth (mm/day)","Blood glucose (mg/dL)","Difference in H:L ratio",
            "Likelihood of second clutch","Eggs in second clutch","Inter-clutch interval (days)")
names(mylabs) <- c("Mass_diff", "R1_mm_day","BloodGluc_2","HL_diff",
                   "Nest2","Nest2Eggs","ICI")


for(i in 1:length(modlist)) {

  a <- modlist[[i]]
  df.mod <- data.frame(summary(a)$coefficients) #extract model coefficients
  df.mod$response <- (response = names(a@frame)[1])#extract name of response variable
  df.mod$term <- rownames(df.mod) #keep predictor names as a variable (not row names)
  rownames(df.mod) <- NULL 
  CI <- data.frame(confint(a, parm = df.mod$term)) #calculate CIs for estimates
  df.mod$low.conf <- CI$X2.5.. #save lower CI bound
  df.mod$hi.conf <- CI$X97.5..#save upper CI bound
  df.mod <- filter(df.mod, term == 'TaggedY' | term == 'BroodTrtE')
  
  p_out <- ggplot(df.mod, aes(y = Estimate, x = term, color = term, shape = term)) +
    geom_hline(yintercept = 0, color = 'black',
               linetype = 2) + # line at null behind coefs
    geom_point(size = 6, position=position_dodge(width = 0.5), alpha =0.8) + 
    
    geom_errorbar(aes(ymin=low.conf, ymax=hi.conf), width=.1,
                  position=position_dodge(.5)) +
    theme_bw() + #use black and white theme
    theme(text = element_text(size=12)) +
    #theme(legend.position = c(0.75, 0.75)) +
    theme(axis.ticks = element_blank()) + # remove axis ticks
    # remove background color
    theme(axis.text.y = element_blank(), #format axis text
          axis.text.x = element_text(color='black', 
                                     size=12, angle=0, 
                                     margin = margin(t = 0, r = 0, 
                                                     b = 0, l = 10))) +
    theme(axis.line = element_blank()) +
    scale_colour_manual(values=c("#000000", "#666666"), #add manual colors and format legend
                        name="", # Legend label
                        breaks=c("BroodTrtE", "TaggedY"),
                        labels=c("Brood Enlarged", "Tagged"))+
    scale_shape_discrete(name="", # Legend label
                        breaks=c("BroodTrtE", "TaggedY"),
                        labels=c("Brood Enlarged", "Tagged")) +
    xlab("") + #no label on x axis for term
    coord_flip() +# flip x and y axes
    facet_wrap(vars(response),
               labeller = labeller(response = mylabs)) # set response variables in individual panels, relabel panels
  
  plots[[i]] <- p_out #store each plot in list plots

}

# calculate exponentiated coeffs and CIs for logistic regression model
i = 5
a <- modlist[[i]]
df.mod <- data.frame(summary(a)$coefficients) #extract model coefficients
df.mod$response <- (response = names(a@frame)[1])#extract name of response variable
df.mod$term <- rownames(df.mod) #keep predictor names as a variable (not row names)
rownames(df.mod) <- NULL 
CI <- data.frame(confint(a, parm = df.mod$term)) #calculate CIs for estimates
df.mod$low.conf <- exp(CI$X2.5..) #save lower CI bound
df.mod$hi.conf <- exp(CI$X97.5..)#save upper CI bound
df.mod$Estimate_exp <- exp(df.mod$Estimate)
df.mod <- filter(df.mod, term == 'TaggedY' | term == 'BroodTrtE')

p_out <- ggplot(df.mod, aes(y = Estimate_exp, x = term, color = term, shape = term)) +
  geom_hline(yintercept = 1, color = 'black',
             linetype = 2) + # line at null behind coefs
  geom_point(size = 6, position=position_dodge(width = 0.5), alpha =0.8) + 
  
  geom_errorbar(aes(ymin=low.conf, ymax=hi.conf), width=.1,
                position=position_dodge(.5)) +
  theme_bw() + #use black and white theme
  theme(text = element_text(size=12)) +
  #theme(legend.position = c(0.75, 0.75)) +
  theme(axis.ticks = element_blank()) + # remove axis ticks
  # remove background color
  theme(axis.text.y = element_blank(), #format axis text
        axis.text.x = element_text(color='black', 
                                   size=12, angle=0, 
                                   margin = margin(t = 0, r = 0, 
                                                   b = 0, l = 10))) +
  theme(axis.line = element_blank()) +
  scale_colour_manual(values=c("#000000", "#666666"), #add manual colors and format legend
                      name="", # Legend label
                      breaks=c("BroodTrtE", "TaggedY"),
                      labels=c("Brood Enlarged", "Tagged"))+
  scale_shape_discrete(name="", # Legend label
                       breaks=c("BroodTrtE", "TaggedY"),
                       labels=c("Brood Enlarged", "Tagged")) +
  xlab("") + #no label on x axis for term
  coord_flip() +# flip x and y axes
  facet_wrap(vars(response),
             labeller = labeller(response = mylabs)) # set response variables in individual panels, relabel panels

plots[[i]] <- p_out #store each plot in list plots



#set up multipanel plots for physiological outcomes
pdf(file = "mod_est_all_exp.pdf", width = 6, height = 6)

ggarrange(plots[[1]] + 
            ylim(-1,1) + 
            theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
            ylab(""), #supress label
          plots[[5]] + 
            ylim(0,20) +
            theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
            ylab(""), #supress label
          plots[[2]] + 
            ylim(-0.5,0.5) +
            theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
            ylab(""),  #supress label
          plots[[6]] + 
            ylim(-0.5,0.5) +
            theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
            ylab(""),#supress label
          plots[[3]] + 
            ylim(45,-45) +
            theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
            ylab(""),  #supress label
          plots[[7]] + 
            ylim(6,-6) +
            theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()),
          plots[[4]] + 
            ylim(1,-1) +
            theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()),
          ncol = 2, nrow = 4, align = "hv", common.legend = TRUE)
          
dev.off()

#set up multipanel plots for reproductive outcomes

pdf(file = "mod_est_rep_all.pdf", width = 3, height = 4.5)
ggarrange(plots[[5]] + ylim(-3,3)+
            theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
            ylab(""), #supress label
          plots[[6]] + ylim(-0.5,0.5)+
            theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
            ylab(""),#supress label
          plots[[7]] + ylim(6,-6) +
          theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()), 
          ncol = 2, common.legend = TRUE)
dev.off()







##################################################################################
## Tagged birds models ##

#list of models - foraging area models
tagmods <- list(
  mod1 <- lmer(Mass_diff ~ scale(After) + BroodTrt + (1|Region), data = df_gps_change),
  mod2 <- lmer(BloodGluc_2 ~ scale(After) + BroodTrt + scale(Mass_2) + tsin_bg + tcos_bg +                      (1|Region), data = df_gps_change),
  mod3 <- lmer(HL_diff ~ scale(After) + BroodTrt + (1|Region), data = HL_gps_df),
  mod4 <- lmer(R1_mm_day ~ scale(After) + BroodTrt + (1|Region), 
               data = df_gps_R1),
  mod5 <- glmer(Nest2Eggs ~ scale(After) + BroodTrt + (1|Region), 
                family = "poisson", 
                data = egg_gps),
  mod6 <- lmer(ICI ~ scale(After)*BroodTrt + (1|Region),  
               data = df_gps_change)
)

#labels for y-axes and corresponding variable names
taglabs <- c("Mass (g)","Feather growth (mm/day)","Blood glucose (mg/dL)","H:L ratio",
            "Likelihood of second clutch","Eggs in second clutch","Inter-clutch interval (days)")
names(taglabs) <- c("Mass_diff", "R1_mm_day","BloodGluc_2","HL_diff",
                   "Nest2","Nest2Eggs","ICI")

# initialize list of plots
tagplots <- list ()

for(i in 1:length(tagmods)) {

  a <- tagmods[[i]]
  df.mod <- data.frame(summary(a)$coefficients) #extract model coefficients
  df.mod$response <- (response = names(a@frame)[1])#extract name of response variable
  df.mod$term <- rownames(df.mod) #keep predictor names as a variable (not row names)
  rownames(df.mod) <- NULL 
  CI <- data.frame(confint(a, parm = df.mod$term)) #calculate CIs for estimates
  df.mod$low.conf <- CI$X2.5.. #save lower CI bound
  df.mod$hi.conf <- CI$X97.5..#save upper CI bound
  df.mod <- filter(df.mod, term == 'scale(After)' | term == 'BroodTrtE')
  
  p_out <- ggplot(df.mod, aes(y = Estimate, x = term, color = term, shape = term)) +
    geom_hline(yintercept = 0, color = 'black',
               linetype = 2) + # line at null behind coefs
    geom_point(size = 6, position=position_dodge(width = 0.5), alpha =0.8) + 
    geom_errorbar(aes(ymin=low.conf, ymax=hi.conf), width=.1,
                  position=position_dodge(.5)) +
    theme_bw() + #use black and white theme
    theme(text = element_text(size=12)) +
    #theme(legend.position = c(0.75, 0.75)) +
    theme(axis.ticks = element_blank()) + # remove axis ticks
    # remove background color
    theme(axis.text.y = element_blank(), #format axis text
          axis.text.x = element_text(color='black', 
                                     size=12, angle=0, 
                                     margin = margin(t = 0, r = 0, 
                                                     b = 0, l = 10))) +
    theme(axis.line = element_blank()) +
    xlab("") + #no label on x axis for term
    coord_flip() +# flip x and y axes
    facet_wrap(vars(response),
               labeller = labeller(response = taglabs)) # set response variables in individual panels, relabel panels
  
  tagplots[[i]] <- p_out #store each plot in list plots
  
}

#####################################################################
#  CIS for tagged birds second clutch model
a <- glm(Nest2 ~ scale(After) + BroodTrt + CI_1_j, 
         family = "binomial", 
         data = df_gps_change)

df.mod <- data.frame(summary(a)$coefficients) #extract model coefficients
df.mod$response <- "Nest2" # name of response variable
df.mod$term <- rownames(df.mod) #keep predictor names as a variable (not row names)
rownames(df.mod) <- NULL 
CI <- data.frame(confint(a, parm = df.mod$term)) #calculate CIs for estimates
df.mod$low.conf <- exp(CI$X2.5..) #save lower CI bound
df.mod$hi.conf <- exp(CI$X97.5..)#save upper CI bound
df.mod <- filter(df.mod, term == 'scale(After)' | term == 'BroodTrtE')
df.mod$Estimate_exp <- exp(df.mod$Estimate)

p_out <- ggplot(df.mod, aes(y = Estimate_exp, x = term, color = term, shape = term)) +
  geom_hline(yintercept = 1, color = 'black',
             linetype = 2) + # line at null behind coefs
  geom_point(size = 6, position=position_dodge(width = 0.5), alpha =0.8) + 
  
  geom_errorbar(aes(ymin=low.conf, ymax=hi.conf), width=.1,
                position=position_dodge(.5)) +
  theme_bw() + #use black and white theme
  theme(text = element_text(size=12)) +
  #theme(legend.position = c(0.75, 0.75)) +
  theme(axis.ticks = element_blank()) + # remove axis ticks
  # remove background color
  theme(axis.text.y = element_blank(), #format axis text
        axis.text.x = element_text(color='black', 
                                   size=12, angle=0, 
                                   margin = margin(t = 0, r = 0, 
                                                   b = 0, l = 10))) +
  theme(axis.line = element_blank()) +
  xlab("") + #no label on x axis for term
  coord_flip() +# flip x and y axes
  facet_wrap(vars(response),
             labeller = labeller(response = taglabs)) # set response variables in individual panels, relabel panels

p_out

tagplots[[7]] <- p_out









#################################################################
#set up multi plot panel using ggpubr package
pdf(file = "mod_est_tag_exp.pdf", width = 6, height = 6)

ggarrange(tagplots[[1]] + 
            ylim(-1,1) +
            theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
            scale_colour_manual(values=c("#111111", "#000000"), #add manual colors and format legend
                                name="", # Legend label
                                breaks=c("scale(After)","BroodTrtE"),
                                labels=c("Range Size", "Brood Enlarged")) +
            scale_shape_manual(values = c(0,19),
                               name="", # Legend label
                                 breaks=c("scale(After)","BroodTrtE"),
                                 labels=c("Range Size", "Brood Enlarged")) +
            ylab(""),
          tagplots[[7]] +
            ylim(0,25) +
            theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
            scale_colour_manual(values=c("#111111", "#000000"), #add manual colors and format legend
                                name="", # Legend label
                                breaks=c("scale(After)", "BroodTrtE"),
                                labels=c("Range Size", "Brood Enlarged"))  +
            scale_shape_manual(values = c(0,19),
                               name="", # Legend label
                               breaks=c("scale(After)","BroodTrtE"),
                               labels=c("Range Size", "Brood Enlarged")) +
            ylab(""), 
          tagplots[[4]] + 
            ylim(-0.5,0.5)+
            theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
            scale_colour_manual(values=c("#111111", "#000000"), #add manual colors and format legend
                                name="", # Legend label
                                breaks=c("scale(After)","BroodTrtE"),
                                labels=c("Range Size", "Brood Enlarged")) +
            scale_shape_manual(values = c(0,19),
                               name="", # Legend label
                               breaks=c("scale(After)","BroodTrtE"),
                               labels=c("Range Size", "Brood Enlarged")) +
            ylab(""),  
          tagplots[[5]] + 
            ylim(-1,1) +
            theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
            scale_colour_manual(values=c("#111111", "#000000"), #add manual colors and format legend
                                name="", # Legend label
                                breaks=c("scale(After)", "BroodTrtE"),
                                labels=c("Range Size", "Brood Enlarged")) +
            scale_shape_manual(values = c(0,19),
                               name="", # Legend label
                               breaks=c("scale(After)","BroodTrtE"),
                               labels=c("Range Size", "Brood Enlarged")) +
            ylab(""),
          tagplots[[2]] + 
            ylim(50,-50)+
            theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
            scale_colour_manual(values=c("#111111", "#000000"), #add manual colors and format legend
                                name="", # Legend label
                                breaks=c("scale(After)","BroodTrtE"),
                                labels=c("Range Size", "Brood Enlarged")) +
            scale_shape_manual(values = c(0,19),
                               name="", # Legend label
                               breaks=c("scale(After)","BroodTrtE"),
                               labels=c("Range Size", "Brood Enlarged")) +
            ylab(""),  
          tagplots[[6]] +
            ylim(6,-6) +
            theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
            scale_colour_manual(values=c("#111111", "#000000"), #add manual colors and format legend
                                name="", # Legend label
                                breaks=c("scale(After)", "BroodTrtE"),
                                labels=c("Range Size", "Brood Enlarged")) +
            scale_shape_manual(values = c(0,19),
                               name="", # Legend label
                               breaks=c("scale(After)","BroodTrtE"),
                               labels=c("Range Size", "Brood Enlarged")),
          tagplots[[3]] + 
            ylim(1.5,-1.5)+
            theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
            scale_colour_manual(values=c("#111111", "#000000"), #add manual colors and format legend
                                name="", # Legend label
                                breaks=c("scale(After)","BroodTrtE"),
                                labels=c("Range Size", "Brood Enlarged")) +
            scale_shape_manual(values = c(0,19),
                               name="", # Legend label
                               breaks=c("scale(After)","BroodTrtE"),
                               labels=c("Range Size", "Brood Enlarged")), 
          ncol = 2, nrow = 4, common.legend = TRUE)

dev.off()
