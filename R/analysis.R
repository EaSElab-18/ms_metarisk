## Code for: 

## --- Clear out --- ##
rm(list = ls())

## --- Libraries --- ##
library(ggpubr)
library(patchwork)
library(gt)
library(metafor)
library(multcomp)
library(metaAidR)
library(tidyverse)

## --- Functions --- ##

# MLMR builder with common parameters, for simplicity
mlmr <- function(dat, variable, vcv){
  rma.mv(yi, 
         vcv,
         mods = ~ variable - 1,
         random = list(~ 1 | obs, ~1 | study_id),
         method = "REML", 
         data = dat)
}

r2mod <- function(model) {
  
  # fixed effect variance
  fix <- var(as.numeric(as.vector(model$b) %*% t(as.matrix(model$X))))
  
  # marginal
  R2m <- fix / (fix + sum(model$sigma2))
  
  # conditional
  R2c <- (fix + sum(model$sigma2) - model$sigma2[length(model$sigma2)]) /
    (fix + sum(model$sigma2))
  
  R2s <- c(R2_marginal = R2m, R2_coditional = R2c)
  return(R2s)
}

# Funnel plotter
funnel_plot <- function(dat, ...){
  funnel(residuals(dat),
         dat$vi,
         yaxis = "seinv", 
         steps = 15, 
         digits = 1, 
         back = "white", 
         shade = "white", 
         hlines = "white", 
         pch = 19, 
         col = rgb(0, 100, 200, max = 255), 
         bg = rgb(255, 255, 255, max = 255, alpha = 150), 
         cex = 0.9,
         ...)
  abline(v = 0, lty = 2)
}

## --- Data --- ##

# Eavesdropper-preference data
dat_pref <- read.csv('../data/dat_predpref.csv')

# Eavesdropping-risk data
dat_risk <- read.csv('../data/dat_predrisk.csv')

## --- Processing --- ##

# Calculate inverse s.e for later convenience
dat_pref$precision <- sqrt(1/dat_pref$vi) 
dat_risk$precision <- sqrt(1/dat_risk$vi) 

## --- Meta-analysis models --- ##

# Create a variance-covariance matrix to account for repeated measurements from groups
control_vcv_pref <- make_VCV_matrix(dat_pref, V = 'vi', cluster = 'grp', rho = 0.5)
control_vcv_risk <- make_VCV_matrix(dat_risk, V = 'vi', cluster = 'grp', rho = 0.5)

# Models
models_pref <- list()
models_risk <- list()

# Predator preference - null
models_pref$null <- rma.mv(yi, 
                           control_vcv_pref,
                           random = list(~ 1 | obs, ~1 | study_id),
                           method = "REML",
                           slab = paste(author, year, sep=", "),
                           data = dat_pref)
models_pref$null_predict <- predict(models_pref$null)
models_pref$null_isq <- I2(models_pref$null, dat_pref$vi, obs = 'obs')   

# Predation risk - null
models_risk$null <- rma.mv(yi, 
                           control_vcv_risk,
                           random = list(~ 1 | obs, ~1 | study_id),
                           method = "REML",
                           slab = paste(author, year, sep=", "),
                           data = dat_risk)
models_risk$null_predict <- predict(models_risk$null)
models_risk$null_isq <- I2(models_risk$null, dat_pref$vi, obs = 'obs')   

## Publication bias & related checks

# Funnel plots
png('../figs/fig_funnels.png', width = 19, height = 11, res = 300, units = 'in')
par(mfrow = c(1, 2))
funnel_plot(models_pref$null, 
            xlab = 'Conditional residuals')
title('Eavesdropper preference', adj = 0, line = 1)
funnel_plot(models_risk$null, 
            xlab = 'Conditional residuals',
            xlim = c(-4, 4),
            ylab = '')
title('Eavesdropping risk', adj = 0, line = 1)
dev.off()

# Regression tests
bias_preftest <- regtest(rma(yi = residuals(models_pref$null), sei = 1/sqrt(1/models_pref$null$vi)), model = "lm")
bias_pref <- trimfill(rma(yi = residuals(models_pref$null), sei = 1/sqrt(1/models_pref$null$vi)), estimator = 'R0')

bias_risktest <- regtest(rma(yi = residuals(models_risk$null), sei = 1/sqrt(1/models_risk$null$vi)), model = "lm")
bias_risk <- trimfill(rma(yi = residuals(models_risk$null), sei = 1/sqrt(1/models_risk$null$vi)), estimator = 'R0')

## --- Moderators --- ##

# Modality (visual/auditory/vibratory/olfactory)
models_pref$modality <- mlmr(dat_pref, dat_pref$modality, control_vcv_pref)
models_pref$modality_predict <- 
  predict(models_pref$modality,
          newmods = matrix(c(1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1), 
                           nrow = 4), 
          addx = TRUE)
models_pref$modality_r2 <- r2mod(models_pref$modality)
models_pref$modality_i2 <- I2(models_pref$modality, dat_risk$vi, obs = 'obs')
models_pref$modality_comp <- summary(glht(models_pref$modality, 
                                          linfct = contrMat(table(dat_pref$modality), 
                                                            type = "Tukey")))

models_risk$modality <- mlmr(dat_risk, dat_risk$modality, control_vcv_risk)
models_risk$modality_predict <- 
  predict(models_risk$modality,
          newmods = matrix(c(1, 0, 0, 0, 1, 0, 0, 0, 1), 
                           nrow = 3), 
          addx = TRUE)
models_risk$modality_r2 <- r2mod(models_risk$modality)
models_risk$modality_i2 <- I2(models_risk$modality, dat_risk$vi, obs = 'obs')
models_risk$modality_comp <- 
  summary(glht(models_risk$modality, linfct = contrMat(table(dat_risk$modality),
                                                       type = "Tukey")))

# Receiver (predator/parasite/parasitoid)
models_pref$receiver <- mlmr(dat_pref, dat_pref$pred_type, control_vcv_pref)
models_pref$receiver_predict <- 
  predict(models_pref$receiver,
          newmods = matrix(c(1, 0, 0, 1), 
                           nrow = 2), 
          addx = TRUE)
models_pref$receiver_r2 <- r2mod(models_pref$receiver)
models_pref$receiver_i2 <- I2(models_pref$receiver, dat_risk$vi, obs = 'obs')
models_pref$receiver_comp <- 
  summary(glht(models_pref$receiver, linfct = contrMat(table(dat_pref$pred_type),
                                                       type = "Tukey")))

models_risk$receiver <- mlmr(dat_risk, dat_risk$pred_type, control_vcv_risk)
models_risk$receiver_predict <- 
  predict(models_risk$receiver,
          newmods = matrix(c(1, 0, 0, 0, 1, 0, 0, 0, 1), 
                           nrow = 3), 
          addx = TRUE)
models_risk$receiver_r2 <- r2mod(models_risk$receiver)
models_risk$receiver_i2 <- I2(models_risk$receiver, dat_risk$vi, obs = 'obs')
models_risk$receiver_comp <- 
  summary(glht(models_risk$receiver, linfct = contrMat(table(dat_risk$pred_type),
                                                       type = "Tukey")))

# Magnitude of manipulation (discrete/continuous)
models_pref$manipulation <- mlmr(dat_pref, dat_pref$continuous_or_discrete_manipulation, control_vcv_pref)
models_pref$manipulation_predict <- 
  predict(models_pref$manipulation,
          newmods = matrix(c(1, 0, 0, 1), 
                           nrow = 2), 
          addx = TRUE)
models_pref$manipulation_r2 <- r2mod(models_pref$manipulation)
models_pref$manipulation_i2 <- I2(models_pref$manipulation, dat_risk$vi, obs = 'obs')
models_pref$manipulation_comp <- 
  summary(glht(models_pref$manipulation, linfct = contrMat(table(dat_pref$continuous_or_discrete_manipulation),
                                                           type = "Tukey")))

models_risk$manipulation <- mlmr(dat_risk, dat_risk$continuous_or_discrete_manipulation, control_vcv_risk)
models_risk$manipulation_predict <- 
  predict(models_risk$manipulation,
          newmods = matrix(c(1, 0, 0, 1), 
                           nrow = 2), 
          addx = TRUE)
models_risk$manipulation_r2 <- r2mod(models_risk$manipulation)
models_risk$manipulation_i2 <- I2(models_risk$manipulation, dat_risk$vi, obs = 'obs')
models_risk$manipulation_comp <- 
  summary(glht(models_risk$manipulation, linfct = contrMat(table(dat_risk$continuous_or_discrete_manipulation),
                                                           type = "Tukey")))

## --- Figures --- ##

# Forest plot - null & modality

# Summary data
dat_pref$intercept <- "Intercept-only"
dat_risk$intercept <- "Intercept-only"
int_pref <- data.frame(intercept = "Intercept-only",
                       yi = transf.ilogit(models_pref$null$b),
                       lower = transf.ilogit(models_pref$null$ci.lb),
                       upper = transf.ilogit(models_pref$null$ci.ub),
                       lower.pi = transf.ilogit(models_pref$null_predict$pi.lb),
                       upper.pi = transf.ilogit(models_pref$null_predict$pi.ub))
int_risk <- data.frame(intercept = "Intercept-only",
                       yi = models_risk$null$b,
                       lower = models_risk$null$ci.lb,
                       upper = models_risk$null$ci.ub,
                       lower.pi = models_risk$null_predict$pi.lb,
                       upper.pi = models_risk$null_predict$pi.ub)

mod_est_pref <- data.frame(modality = c('Auditory', 'Olfactory', 'Vibratory', 'Visual'),
                           yi = transf.ilogit(models_pref$modality$b),
                           lower = transf.ilogit(models_pref$modality$ci.lb),
                           upper = transf.ilogit(models_pref$modality$ci.ub),
                           lower.pi = transf.ilogit(models_pref$modality_predict$pi.lb),
                           upper.pi = transf.ilogit(models_pref$modality_predict$pi.ub))
dat_pref$modality <- factor(dat_pref$modality, levels = c('Vibratory', 'Auditory', 'Olfactory', 'Visual'))

mod_est_risk <- data.frame(modality = c('Auditory', 'Olfactory', 'Visual', 'Vibratory'),
                           yi = c(models_risk$modality$b, NA),
                           lower = c(models_risk$modality$ci.lb, NA),
                           upper = c(models_risk$modality$ci.ub, NA),
                           lower.pi = c(models_risk$modality_predict$pi.lb, NA),
                           upper.pi = c(models_risk$modality_predict$pi.ub, NA))
levels(dat_risk$modality) <- c(levels(dat_risk$modality), "Vibratory")
dat_risk$modality <- factor(dat_risk$modality, levels = c('Vibratory', 'Auditory', 'Olfactory', 'Visual'))

# Plots
(forest_pref_null <-
    ggplot(dat_pref, aes(x = intercept, y = transf.ilogit(yi))) + 
    geom_point(aes(size = precision), 
               shape = 21, 
               fill = "white",
               colour = '#00B9E3',
               stroke = 0.8, 
               alpha = 0.7,
               position = position_jitter(width = 0.2, height = 0, seed = 1115)) +
    scale_size(range = c(0.5, 6)) +
    ylab('') +
    xlab('') +
    ggtitle('Eavesdropper preference') +
    coord_flip(ylim = c(0, 1)) +
    theme_classic() +
    labs(tag = "Mean") +
    theme(legend.position = "none",
          axis.text.x = element_blank(),
          axis.line.x = element_blank(),
          axis.ticks.x = element_blank(),
          plot.margin = unit(c(1, 1, 0, 1), "cm"),
          plot.tag.position = c(0.05, 0.9),
          plot.tag = element_text(size = 11)) +
    geom_hline(yintercept = 0.5, lty = 2) +
    geom_point(data = int_pref, mapping = aes(x = intercept, y = yi), size = 2, col = "black") +
    geom_errorbar(data = int_pref, mapping = aes(x = intercept, ymin = lower, ymax = upper),
                  width = 0, size = 1, color = "black") +
    geom_errorbar(data = int_pref, mapping = aes(x = intercept, ymin = lower.pi, ymax = upper.pi),
                  width = 0, size = 0.5, lty = 2, color = "black"))

(forest_risk_null <- 
    ggplot(dat_risk, aes(x = intercept, y = yi)) + 
    geom_point(aes(size = precision), 
               shape = 21, 
               fill = "white",
               colour = '#00B9E3',
               stroke = 0.8, 
               alpha = 0.7,
               position = position_jitter(width = 0.2, height = 0, seed = 1115)) +
    scale_size(range = c(0.5, 6)) +
    scale_x_discrete(drop = FALSE) +
    scale_y_continuous(breaks = seq(-8, 8, by = 2)) +
    xlab("") +  
    ylab("") +
    ggtitle('Eavesdropping risk') +
    coord_flip(ylim = c(-2, 8)) +
    theme_classic() +
    theme(axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.line.x = element_blank(),
          axis.ticks.x = element_blank(),
          plot.margin = unit(c(1, 1, 0, 1), "cm")) +
    geom_hline(yintercept = 0, lty = 2) +
    geom_point(data = int_risk, size = 2, col = "black") +
    geom_errorbar(data = int_risk, mapping = aes(x = intercept, ymin = lower, ymax = upper),
                  width = 0, size = 1, color = "black") +
    geom_errorbar(data = int_risk, mapping = aes(x = intercept, ymin = lower.pi, ymax = upper.pi),
                  width = 0, size = 0.5, lty = 2, color = "black"))

(forest_pref <- 
    ggplot(dat_pref, aes(x = modality, y = transf.ilogit(yi))) + 
    geom_point(aes(size = precision), 
               shape = 21, 
               fill = "white",
               colour = '#00B9E3',
               stroke = 0.8, 
               alpha = 0.7,
               position = position_jitter(width = 0.2, height = 0, seed = 1115)) +
    scale_size(range = c(0.5, 6)) +
    xlab('') +
    ylab("Proportion of choices") +
    coord_flip(ylim = c(0, 1)) +
    theme_classic() +
    labs(tag = "Modality") +
    theme(legend.position = "none",
          plot.margin = unit(c(0, 1, 1, 1), "cm"),
          plot.tag.position = c(0.1, 1),
          plot.tag = element_text(size = 11)) +
    geom_hline(yintercept = 0.5, lty = 2) +
    geom_point(data = mod_est_pref, mapping = aes(x = modality, y = yi), size = 2, col = "black") +
    geom_errorbar(data = mod_est_pref, mapping = aes(x = modality, ymin = lower, ymax = upper),
                  width = 0, size = 1, color = "black") +
    geom_errorbar(data = mod_est_pref, mapping = aes(x = modality, ymin = lower.pi, ymax = upper.pi),
                  width = 0, size = 0.5, lty = 2, color = "black"))

(forest_risk <- 
    ggplot(dat_risk, aes(x = modality, y = yi)) + 
    geom_point(aes(size = precision), 
               shape = 21, 
               fill = "white",
               colour = '#00B9E3',
               stroke = 0.8, 
               alpha = 0.7,
               position = position_jitter(width = 0.2, height = 0, seed = 1115)) +
    scale_size(range = c(0.5, 6)) +
    scale_x_discrete(drop = FALSE) +
    scale_y_continuous(breaks = seq(-8, 8, by = 2)) +
    xlab("") +  
    ylab("SMD") +
    coord_flip(ylim = c(-2, 8)) +
    theme_classic() +
    theme(legend.position = "none",
          axis.text.y = element_blank(),
          plot.margin = unit(c(0, 1, 1, 1), "cm")) +
    geom_hline(yintercept = 0, lty = 2) +
    geom_point(data = mod_est_risk, size = 2, col = "black") +
    geom_errorbar(data = mod_est_risk, mapping = aes(x = modality, ymin = lower, ymax = upper),
                  width = 0, size = 1, color = "black") +
    geom_errorbar(data = mod_est_risk, mapping = aes(x = modality, ymin = lower.pi, ymax = upper.pi),
                  width = 0, size = 0.5, lty = 2, color = "black"))

(forest_pref_null + forest_risk_null) / 
  (forest_pref + forest_risk)

ggsave('../figs/fig_forest_modality.tiff', width = 24, height = 18, units = 'cm')
ggsave('../figs/fig_forest_modality.png', width = 24, height = 18, units = 'cm')

# Forest plot - receiver and manipulation

levels(dat_pref$pred_type) <- c(levels(dat_pref$receiver), "Parasite")
dat_pref$pred_type <- factor(dat_pref$pred_type, levels = c('Parasite', 'Parasitoid', 'Predator'))

rec_est_pref <- data.frame(eavesdropper = c('Parasitoid', 'Predator', 'Parasite'),
                           yi = c(transf.ilogit(models_pref$receiver$b), NA),
                           lower = c(transf.ilogit(models_pref$receiver$ci.lb), NA),
                           upper = c(transf.ilogit(models_pref$receiver$ci.ub), NA),
                           lower.pi = c(transf.ilogit(models_pref$receiver_predict$pi.lb), NA),
                           upper.pi = c(transf.ilogit(models_pref$receiver_predict$pi.ub), NA))

#levels(dat_risk$pred_type) <- c(levels(dat_risk$receiver), "Parasite")
dat_risk$pred_type <- factor(dat_risk$pred_type, levels = c('Parasite', 'Parasitoid', 'Predator'))

rec_est_risk <- data.frame(eavesdropper = c('Parasite', 'Parasitoid', 'Predator'),
                           yi = models_risk$receiver$b,
                           lower = models_risk$receiver$ci.lb,
                           upper = models_risk$receiver$ci.ub,
                           lower.pi = models_risk$receiver_predict$pi.lb,
                           upper.pi = models_risk$receiver_predict$pi.ub)

manip_est_pref <- data.frame(manipulation = c('Continuous', 'Discrete'),
                             yi = transf.ilogit(models_pref$manipulation$b),
                             lower = transf.ilogit(models_pref$manipulation$ci.lb),
                             upper = transf.ilogit(models_pref$manipulation$ci.ub),
                             lower.pi = transf.ilogit(models_pref$manipulation_predict$pi.lb),
                             upper.pi = transf.ilogit(models_pref$manipulation_predict$pi.ub))

manip_est_risk <- data.frame(manipulation = c('Continuous', 'Discrete'),
                             yi = models_risk$manipulation$b,
                             lower = models_risk$manipulation$ci.lb,
                             upper = models_risk$manipulation$ci.ub,
                             lower.pi = models_risk$manipulation_predict$pi.lb,
                             upper.pi = models_risk$manipulation_predict$pi.ub)

# Plots
(forest_rec_pref <- 
    ggplot(dat_pref, aes(x = pred_type, y = transf.ilogit(yi))) + 
    geom_point(aes(colour = pred_type, size = precision), 
               shape = 21, 
               fill = "white", 
               colour = '#00B9E3',
               stroke = 0.8, 
               alpha = 0.7,
               position = position_jitter(width = 0.2, height = 0, seed = 1115)) +
    scale_size(range = c(0.5, 6)) +
    xlab('') +
    ylab("") +
    ggtitle('Eavesdropper preference') +
    coord_flip(ylim = c(0, 1)) +
    theme_classic() +
    labs(tag = "Eavesdropper") +
    scale_x_discrete(drop = FALSE) +
    theme(legend.position = "none",
          axis.line.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.x = element_blank(),
          plot.margin = unit(c(1, 1, 0, 1), "cm"),
          plot.tag.position = c(0.05, 0.95),
          plot.tag = element_text(size = 11)) +
    geom_hline(yintercept = 0.5) +
    geom_point(data = rec_est_pref, mapping = aes(x = eavesdropper, y = yi), size = 2, col = "black") +
    geom_errorbar(data = rec_est_pref, mapping = aes(x = eavesdropper, ymin = lower, ymax = upper),
                  width = 0, size = 1, color = "black") +
    geom_errorbar(data = rec_est_pref, mapping = aes(x = eavesdropper, ymin = lower.pi, ymax = upper.pi),
                  width = 0, size = 0.5, lty = 2, color = "black"))

(forest_manip_pref <- 
    ggplot(dat_pref, aes(x = continuous_or_discrete_manipulation, y = transf.ilogit(yi))) + 
    geom_point(aes(colour = continuous_or_discrete_manipulation, size = precision), 
               shape = 21, 
               fill = "white",
               colour = '#00B9E3',
               stroke = 0.8, 
               alpha = 0.7,
               position = position_jitter(width = 0.2, height = 0, seed = 1115)) +
    scale_size(range = c(0.5, 6)) +
    xlab('') +
    ylab("Proportion of choices") +
    coord_flip(ylim = c(0, 1)) +
    theme_classic() +
    labs(tag = "Manipulation") +
    theme(legend.position = "none",
          plot.margin = unit(c(0, 1, 1, 1), "cm"),
          plot.tag.position = c(0.05, 0.95),
          plot.tag = element_text(size = 11)) +
    geom_hline(yintercept = 0.5) +
    geom_point(data = manip_est_pref, mapping = aes(x = manipulation, y = yi), size = 2, col = "black") +
    geom_errorbar(data = manip_est_pref, mapping = aes(x = manipulation, ymin = lower, ymax = upper),
                  width = 0, size = 1, color = "black") +
    geom_errorbar(data = manip_est_pref, mapping = aes(x = manipulation, ymin = lower.pi, ymax = upper.pi),
                  width = 0, size = 0.5, lty = 2, color = "black"))

(forest_rec_risk <- 
    ggplot(dat_risk, aes(x = pred_type, y = yi)) + 
    geom_point(aes(colour = pred_type, size = precision), 
               shape = 21, 
               fill = "white", 
               colour = '#00B9E3',
               stroke = 0.8, 
               alpha = 0.7,
               position = position_jitter(width = 0.2, height = 0, seed = 1115)) +
    scale_size(range = c(0.5, 6)) +
    scale_y_continuous(breaks = seq(-8, 8, by = 2)) +
    xlab("") +  
    ylab("") +
    ggtitle('Eavesdropping risk') +
    coord_flip(ylim = c(-2, 8)) +
    theme_classic() +
    theme(#legend.position = "none",
      axis.text.y = element_blank(),
      axis.line.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.text.x = element_blank(),
      plot.margin = unit(c(1, 1, 0, 1), "cm")) +
    geom_hline(yintercept = 0) +
    geom_point(data = rec_est_risk, mapping = aes(x = eavesdropper, y = yi), size = 2, col = "black") +
    geom_errorbar(data = rec_est_risk, mapping = aes(x = eavesdropper, ymin = lower, ymax = upper),
                  width = 0, size = 1, color = "black") +
    geom_errorbar(data = rec_est_risk, mapping = aes(x = eavesdropper, ymin = lower.pi, ymax = upper.pi),
                  width = 0, size = 0.5, lty = 2, color = "black"))

(forest_manip_risk <- 
    ggplot(dat_risk, aes(x = continuous_or_discrete_manipulation, y = yi)) + 
    geom_point(aes(colour = continuous_or_discrete_manipulation, size = precision), 
               shape = 21, 
               fill = "white", 
               colour = '#00B9E3',
               stroke = 0.8, 
               alpha = 0.7,
               position = position_jitter(width = 0.2, height = 0, seed = 1115)) +
    scale_size(range = c(0.5, 6)) +
    scale_y_continuous(breaks = seq(-8, 8, by = 2)) +
    xlab("") +  
    ylab("SMD") +
    coord_flip(ylim = c(-2, 8)) +
    theme_classic() +
    theme(legend.position = "none",
          axis.text.y = element_blank(),
          plot.margin = unit(c(0, 1, 1, 1), "cm")) +
    geom_hline(yintercept = 0) +
    geom_point(data = manip_est_risk, mapping = aes(x = manipulation, y = yi), size = 2, col = "black") +
    geom_errorbar(data = manip_est_risk, mapping = aes(x = manipulation, ymin = lower, ymax = upper),
                  width = 0, size = 1, color = "black") +
    geom_errorbar(data = manip_est_risk, mapping = aes(x = manipulation, ymin = lower.pi, ymax = upper.pi),
                  width = 0, size = 0.5, lty = 2, color = "black"))


(forest_rec_pref + forest_rec_risk) / 
  (forest_manip_pref + forest_manip_risk)

ggsave('../figs/fig_forest_rec_manip.tiff', width = 24, height = 24, units = 'cm')
ggsave('../figs/fig_forest_rec_manip.png', width = 24, height = 24, units = 'cm')