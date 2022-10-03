################################################################################
# Molly McDermott
# Created 9/3/21
# There are many useful tools for summarizing and visualizing model results in R, however, I often have a need for customization to make publication quality figures and tables. This script tidies model output (lmer or glmer) and produces a table and a figure with coefficient estimates. This script was developed for use with repeated-measures (longitudinal) data.
################################################################################


#### SET UP ####

# load packages
library(lme4)
library(lmerTest)
library(plyr)
library(dplyr)
library(tidyr)
library(purrr)
library(ggplot2)
library(ggpubr)

# create example data
df <- data.frame(
  mass = rnorm(50, mean = 18, sd = 1), #mass of birds in grams
  WBC_count = rpois(50, lambda = 44), #white blood cell count
  tagged = rep(c("Y","N"), times=25, each=1), #experimental treatment: tagged or not?
  eggs = rpois(n = 50, lambda = 4), #reproductive output: number of eggs
  date = rep(c(160,180), times=1, each=25), #julian day
  ID = rep(seq(from = 1, to = 25, by = 1), times = 2)  #unique ID for each individual, two observations of each
)

#set variable types
df$tagged <- as.factor(df$tagged)
df$ID <- as.factor(df$ID)


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

#list of models 
modlist <- list(
  mod1 <- lmer(mass ~ tagged*date + (1|ID), 
               data = df),
  mod2 <- lmer(WBC_count ~ tagged*date + (1|ID), 
               data = df),
  mod3 <- glmer(eggs ~ tagged*date + (1|ID), 
               family = "poisson",
               data = df)
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
  write.table(coeffsp, "model_results.csv", 
              sep = ",", 
              append = TRUE, 
              row.names = FALSE) # export
}



################################################################################
#### COEFFICIENT PLOTS ####

# initialize list of plots
plots <- list ()

#list new labels for plot panels and corresponding response variable names
mylabs <- c("Mass (g)",
            "White Blood Cell count",
            "Eggs in second clutch"
            )

names(mylabs) <- c("mass", 
                   "WBC_count",
                   "eggs"
                   )


for(i in 1:length(modlist)) {

  a <- modlist[[i]]
  df.mod <- data.frame(summary(a)$coefficients) #extract model coefficients
  df.mod$response <- (response = names(a@frame)[1])#extract name of response variable
  df.mod$term <- rownames(df.mod) #keep predictor names as a variable (not row names)
  rownames(df.mod) <- NULL 
  CI <- data.frame(confint(a, parm = df.mod$term)) #calculate CIs for estimates
  df.mod$low.conf <- CI$X2.5.. #save lower CI bound
  df.mod$hi.conf <- CI$X97.5..#save upper CI bound
  df.mod <- filter(df.mod, term != '(Intercept)') #filter out intercept parameters
  
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
               labeller = labeller(response = mylabs)) # set response variables in individual panels, relabel panels
  
  plots[[i]] <- p_out #store each plot in list plots
  
}


#set up multipanel plots for outcomes
pdf(file = "model_estimates.pdf", width = 6, height = 6)

ggarrange(plots[[1]] + 
            ylim(-1,1) + 
            theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
            scale_colour_manual(values=c("#000000", "#666666"), #add manual colors and format legend
                                name="", # Legend label
                                breaks=c("date", "taggedY"),
                                labels=c("Date (Julian Day)", "Tagged"))+
            scale_shape_discrete(name="", # Legend label
                                 breaks=c("date", "taggedY"),
                                 labels=c("Date (Julian Day)", "Tagged")) +
            ylab(""), #supress label
          plots[[2]] + 
            ylim(-0.5,0.5) +
            theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
            ylab(""),  #supress label
          plots[[3]] + 
            ylim(45,-45) +
            theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
            ylab(""),  #supress label
)

dev.off()


