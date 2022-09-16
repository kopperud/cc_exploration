## Suppose we set up a reference model, with
## 
## \lambda(t) = sigmoidally increasing
## \mu(t) = 0.28


## And conversely,
## \lambda(t) = 0.28
## \mu(t) = sigmoidally increasing

library(ggplot2)
library(ape)
library(CRABS)
library(magrittr)
library(dplyr)
library(patchwork)

source("scripts/utils.R")


logistic_shift <- function(min, max, steepness, x0){
  res <- function(x) {
    ((max-min) / (1 + steepness^(-1*(x - x0)))) + min
  }
  return(res)
}
####################################
## 
##            Settings
## 
####################################
phy <- read.tree("trees/condamine_etal_2019_Picidae.tre")
height <- max(node.depth.edgelength(phy))
times <- seq(0, height, length.out = 500)
times_fine <- seq(0, height, length.out = 500)
nsamples <- 1000

####################################
## 
##         Reference models
## 
####################################
mu0 <- logistic_shift(0.1, 0.6, 0.1, height/2)
lambda0 <- function(t) 0.35 + 0*t
sigmoidal_mu <- CRABS::create.model(func_spec0 = lambda0, mu0, times = times)

mu1 <- function(t) 0.35 + 0*t
lambda1 <- logistic_shift(0.1, 0.6, 0.1, height/2)
sigmoidal_lambda <- CRABS::create.model(func_spec0 = lambda1, mu1, times = times)

####################################
## 
##       Alternative rates
## 
####################################

p1a <- model2df(sigmoidal_lambda) %>%
  filter(rate %in% c("Speciation", "Extinction")) %>%
  ggplot(aes(x = Time, y = value, color = rate)) +
  geom_line() +
  theme_classic() +
  scale_x_reverse() +
  scale_color_manual(values = c("orange", "black")) +
  theme(legend.position = c(0.2, 0.8),
        legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5)) +
  xlab("Time (Ma)") +
  ylim(c(0.0, 0.6)) +
  ggtitle("Reference model")

foo <- function(reference, proposal, proposed_rate = "lambda", nsamples = 1000){
  ## Sample rates according to proposals
  
  proposals <- list()
  for (i in 1:nsamples){
    proposals[[i]] <- proposal()
  }
  
  ## Set up the congruent models
  if (proposed_rate == "lambda"){
    model_set <- congruent.models(reference, lambdas = proposals, keep_ref = FALSE)
    rate_color <- "orange"
    title <- "Proposed λ'"
    ylab <- "SD[µ']"
  }else if (proposed_rate == "mu"){
    model_set <- congruent.models(reference, mus = proposals, keep_ref = FALSE)
    rate_color <- "black"
    title <- "Proposed µ'"
    ylab <- "SD[λ']"
  }else{
    stop("must be lambda or mu")
  }
  
  flip <- c("lambda" = "Extinction", 
            "mu" = "Speciation")
  rate_name <- flip[proposed_rate]
  
  rate_df <- models2df(model_set) %>%
    filter(rate == rate_name) %>% 
    mutate("n" = as.numeric(stringr::str_extract(name, "[0-9]+"))) %>% 
    filter(n <= 10) ## Don't plot more than 10 rate trajectories
  
  rate_plot <- rate_df %>%
    ggplot(aes(x = Time, y = value, color = name)) +
    geom_line() +
    theme_classic() +
    scale_x_reverse() +
    ylab("Rate") +
    theme(plot.title = element_text(hjust = 0.5),
          axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          legend.position = "none",
          axis.title.y = element_text(angle = 90)) +
    ggtitle(title)
  
  if (proposed_rate == "lambda"){
    col1 <- colorspace::sequential_hcl(palette = "orange", n = length(unique(rate_df$name)))
    rate_plot <- rate_plot +
      ylab("µ'") +
      scale_color_manual(values = col1)
  }else{
    col1 <- head(colorspace::sequential_hcl(palette = "grays", n = length(unique(rate_df$name))*3), n = length(unique(rate_df$name)))
    rate_plot <- rate_plot +
      ylab("λ'") +
      scale_color_manual(values = col1)
    
  }
  
  ## Compute the variance at each time bin
  var_df <- models2df(model_set) %>%
    filter(rate == rate_name) %>%
    group_by(Time) %>%
    summarize("v" = var(value)) %>%
    mutate("sd" = sqrt(v))
  var_df$rate <- rate_name
  ## Generate plots
  
  p <- var_df %>%
    filter(Time > 0.32) %>%
    filter(Time < 77.0) %>%
    ggplot(aes(x = Time, y = sd)) +
    geom_line(color = rate_color) +
    theme_classic() +
    scale_x_reverse() +
    theme(legend.position = c(0.2, 0.7),
          legend.background = element_blank(),
          legend.title = element_blank(),
          plot.title = element_text(hjust = 0.5)) +
    xlab("Time (Ma)") +
    ylab(ylab)
    
  
  ## return plots
  res <- list("rprime" = rate_plot,
              "v" = p,
              "models" = model_set,
              "var_df" = var_df)
  return(res)
}

####################################
## 
##     Generate plots
## 
####################################

ps <- list()
ps[["mu"]] <- list()
ps[["lambda"]] <- list()

####################################
## 
##    Brownian motion
## 
####################################

proposal_bm <- function(){
  dt <- times_fine[[2]] - times_fine[[1]]
  sigma <- 0.02
  
  noise <- c(0.6, rnorm(length(times_fine)-1, sd = sigma * sqrt(dt)))
  y <- cumsum(noise)
  
  f <- approxfun(times_fine, y)
  return(f)
}

ps[["lambda"]][["bm"]] <- foo(sigmoidal_lambda, proposal_bm, proposed_rate = "lambda", nsamples = nsamples)
ps[["mu"]][["bm"]] <- foo(sigmoidal_lambda, proposal_bm, proposed_rate = "mu", nsamples = nsamples)

####################################
## 
##    Ornstein-Uhlenbeck rate
## 
####################################

sim_ou <- function(times, sigma, theta, alpha){
  x <- list()
  x[[1]] <- theta
  dt <- abs(times[[2]] - times[[1]])
  noise <- rnorm(length(times), mean = 0, sd = sigma * sqrt(dt))
  
  for (i in tail(seq_along(times), n = -1)){
    ## dx = - alpha * (x - theta) * dt + sigma * dW
    x[[i]] <- x[[i-1]] - alpha * (x[[i-1]] - theta) * dt + noise[i]
  }
  return(unlist(x))
}

proposal_lambda_ou <- function(){
  alpha <- 0.5
  sigma <- 0.05
  theta <- lambda1(0.0)
  
  y <- sim_ou(times_fine, sigma, theta, alpha)
  f <- approxfun(times_fine, y)
  return(f)
}

proposal_mu_ou <- function(){
  alpha <- 0.5
  sigma <- 0.05
  theta <- lambda1(0.0)
  
  y <- sim_ou(times_fine, sigma, theta, alpha)
  f <- approxfun(times_fine, y)
  return(f)
}

ps[["lambda"]][["ou"]] <- foo(sigmoidal_lambda, proposal_lambda_ou, proposed_rate = "lambda", nsamples = nsamples)
ps[["mu"]][["ou"]] <- foo(sigmoidal_lambda, proposal_mu_ou, proposed_rate = "mu", nsamples = nsamples)


#ps$lambda$ou

####################################
## 
##    IID lognormal
## 
####################################

proposal_mu_lniid <- function(){
  m <- 0.6
  sdlog <- 0.05
  
  y <- c(m, rlnorm(length(times_fine)-1, meanlog = log(m), sdlog = sdlog))
  f <- approxfun(times_fine, y)
  return(f)
}

ps[["lambda"]][["lniid"]] <- foo(sigmoidal_lambda, proposal_mu_lniid, proposed_rate = "lambda", nsamples = nsamples)
ps[["mu"]][["lniid"]] <- foo(sigmoidal_lambda, proposal_mu_lniid, proposed_rate = "mu", nsamples = nsamples)

########################################################################
##
##        Calculate likelihood of the congruent model sets
##
########################################################################
trees <- TESS::tess.sim.taxa.age(5, 25, age = 70, lambda = 1, mu = 0.5)

library(progress)

calc_logLs <- function(tree, model_set){
  logLs <- list()
  pb <- progress_bar$new(
    format = "  calculating logL [:bar] :percent eta: :eta",
    total = length(model_set), clear = FALSE, width= 60)
  pb$tick(0)
  for (i in seq_along(model_set)){
    pb$tick()
    logLs[[i]] <- crabs.loglikelihood(tree, model_set[[i]])
  }
  return(unlist(logLs))
}

logLs_mu <- calc_logLs(trees[[1]], ps$mu$ou$models)
logLs_lambda <- calc_logLs(trees[[1]], ps$lambda$ou$models)

var(logLs_mu)
var(logLs_lambda)

####################################
## 
##    Multi-panel plot
## 
####################################

## Turn off axis labels and adjust things etc
for (rate in c("mu", "lambda")){
  for (proposal in c("ou", "bm")){
    ps[[rate]][[proposal]][["v"]] <- ps[[rate]][[proposal]][["v"]] +
      theme(axis.text.x = element_blank(),
            axis.title.x = element_blank())
  }
  
  ## turn off the title
  for (proposal in c("ou", "lniid")){
    ps[[rate]][[proposal]][["v"]] <- ps[[rate]][[proposal]][["v"]] +
      theme(plot.title = element_blank())
    ps[[rate]][[proposal]][["rprime"]] <- ps[[rate]][[proposal]][["rprime"]] +
      theme(plot.title = element_blank())
  }
}

## BM
ps$lambda$bm$v <- ps$lambda$bm$v +
  ylim(c(0.0, 0.28)) + ylab("SD[µ']")
ps$lambda$bm$rprime <- ps$lambda$bm$rprime +
  ylab("µ'") +
  ylim(c(-0.5, 2.0))
  theme(plot.title = element_text(hjust = 0.5))
ps$mu$bm$rprime <- ps$mu$bm$rprime + 
  ylim(c(-0.5, 2.0)) +
  theme(plot.title = element_text(hjust = 0.5))
ps$mu$bm$v <- ps$mu$bm$v + 
  ylim(c(0.0, 0.28))

## OU
ps$lambda$ou$v <- ps$lambda$ou$v +
  ylim(c(0.0, 0.25)) + 
  ylab("SD[µ']")
ps$lambda$ou$rprime <- ps$lambda$ou$rprime +
  ylab("µ'") +
  ylim(c(-0.5, 2.0))
ps$mu$ou$rprime <- ps$mu$ou$rprime +
  ylim(c(-0.5, 2.0))
ps$mu$ou$v  <- ps$mu$ou$v + 
  ylim(c(0.0, 0.25))

## ln iid
ps$lambda$lniid$v <- ps$lambda$lniid$v +
  ylim(c(0.0, 0.35)) + ylab("SD[µ']")
ps$lambda$lniid$rprime <- ps$lambda$lniid$rprime +
  ylab("µ'") +
  ylim(c(-0.5, 2.0))
ps$mu$lniid$rprime <- ps$mu$lniid$rprime +
  ylim(c(-0.5, 2.0))
ps$mu$lniid$v <-  ps$mu$lniid$v +
  ylim(c(0.0, 0.05)) + ylab("SD[λ']*")

ptbm <- ggplot() + 
  annotate("text", x = 2, y = 4, size=4, label = "Brownian motion", angle = 90) + 
  theme_void()
ptou <- ggplot() + 
  annotate("text", x = 2, y = 4, size=4, label = "Ornstein-Uhlenbeck", angle = 90) + 
  theme_void()
ptln <- ggplot() + 
  annotate("text", x = 2, y = 4, size=4, label = "Lognormal IID", angle = 90) + 
  theme_void()

bm <- ptbm + ((ps$lambda$bm$rprime + ps$mu$bm$rprime) / (ps$lambda$bm$v + ps$mu$bm$v) ) +
  plot_layout(widths = c(0.05, 0.9))
ou <- ptou + ((ps$lambda$ou$rprime + ps$mu$ou$rprime) / (ps$lambda$ou$v + ps$mu$ou$v) ) +
  plot_layout(widths = c(0.05, 0.9))
lniid <- ptln + ((ps$lambda$lniid$rprime + ps$mu$lniid$rprime) / (ps$lambda$lniid$v + ps$mu$lniid$v) ) +
  plot_layout(widths = c(0.05, 0.9))

top_panel <- (plot_spacer() + p1a  + ylab("Rate") + plot_layout(ncol = 2, widths = c(0.05, 0.9)))
variance_plot <- top_panel / bm / ou / lniid

ggsave("figures/suppmat/variance_plot.pdf", variance_plot, width = 180, height = 220, units = "mm", device = cairo_pdf)

####################################
## 
##   Simpler plot for main text
## 
####################################

bm_df <- bind_rows(
  ps$lambda$bm$var_df,
  ps$mu$bm$var_df
)

ou_df <- bind_rows(
  ps$lambda$ou$var_df,
  ps$mu$ou$var_df
)

lniid_df <- bind_rows(
  ps$lambda$lniid$var_df,
  ps$mu$lniid$var_df
)
compact_dfs <- list(bm_df, ou_df, lniid_df)

ps2 <- list()
ylab2 <- c("BM", "OU", "Ln")
for (i in 1:3){
  ps2[[i]] <- compact_dfs[[i]] %>%
    filter(Time > 0.32) %>%
    filter(Time < 77.0) %>%
    ggplot(aes(x = Time, y = sd, color = rate)) +
    geom_line() +
    theme_classic() +
    scale_x_reverse() +
    scale_color_manual(values = c("orange", "black")) +
    theme(legend.position = c(0.15, 0.5),
          legend.title = element_blank(),
          legend.background = element_blank()) +
    xlab("Time (Ma)") +
    ylab(paste0(ylab2[i], "\nSD[rate]"))
  if (i > 1){
    ps2[[i]] <- ps2[[i]] +
      theme(legend.position = "none")
  }
  if(i < 3){
    ps2[[i]] <- ps2[[i]] + 
      theme(axis.text.x = element_blank(),
            axis.title.x = element_blank())
  }
}

varplot_tiny <- p1a + ylab("Rate") + theme(legend.background = element_blank(),
                                           legend.position = "none") +
  ps2[[1]] + ggtitle("Variation in congruent rates") + theme(plot.title = element_text(hjust = 0.5),
                                                           legend.position = "none") +
  ps2[[2]] +
  ps2[[3]] + theme(legend.position = c(0.5, 0.5)) +
  plot_layout(ncol = 1)
#varplot_tiny
ggsave("figures/ms/varplot_tiny.pdf", plot = varplot_tiny, units = "mm", height = 110, width = 85)

Reduce("+", ps2) +
  plot_layout(ncol = 2)


# p2a <- plot(model_set1_frommu)[[2]] +
#   ylim(c(0.0, 0.7)) +
#   theme(plot.title = element_text(hjust = 0.5),
#         axis.title.x = element_blank(),
#         axis.text.x = element_blank()) +
#    ggtitle("Proposed rates")
# 
# var_sp <- models2df(model_set1_frommu) %>%
#   filter(rate == "Speciation") %>%
#   group_by(Time) %>%
#   summarize("v" = var(value))
# var_sp$rate <- "Speciation"
# 
# # n_bins <- 50
# # time_bin_start <- head(seq(-0.1, height+0.01, length.out = n_bins+1), n = -1)
# # time_bin_end <- tail(seq(-0.1, height+0.01, length.out = n_bins+1), n = -1)
# # bin_names <- paste0("bin", 1:n_bins)
# 
# # which_time_bin <- function(x){
# #   bin_names[x > time_bin_start & x < time_bin_end]
# # }
# 
# var_ex <- models2df(model_set1_fromlambda) %>%
#   filter(rate == "Extinction") %>%
#   group_by(Time) %>%
#   summarize("v" = var(value))# %>%
# var_ex$rate <- "Extinction"
# var_df <- bind_rows(var_sp, var_ex)
# 
# p3a <- ggplot(var_df, aes(x = Time, y = v, color = rate)) +
#   geom_line() +
#   theme_classic() +
#   scale_x_reverse() +
#   theme(legend.position = c(0.2, 0.7),
#         legend.background = element_blank(),
#         legend.title = element_blank(),
#         plot.title = element_text(hjust = 0.5),
#         axis.title.y = element_blank()) +
#   scale_color_manual(values = c("orange", "black")) +
#   xlab("Time (Ma)") +
#   ggtitle("Variance in inferred rate")
# 
# ## Right panel
# 
# p1b <- model2df(sigmoidal_mu) %>%
#   dplyr::filter(rate %in% c("Speciation", "Extinction")) %>%
#   ggplot(aes(x = Time, y = value, color = rate)) +
#   geom_line() +
#   theme_classic() +
#   scale_x_reverse() +
#   scale_color_manual(values = c("orange", "black")) +
#   theme(legend.position = c(0.2, 0.8),
#         legend.title = element_blank(),
#         plot.title = element_text(hjust = 0.5),
#         axis.title.x = element_blank(),
#         axis.text.x = element_blank()) +
#   ylab("Rate") +
#   xlab("Time (Ma)") +
#   ylim(c(0.0, 0.6)) +
#   ggtitle("Reference model")
#   
# # p1b <- plot(model_set0_frommu)[[2]] + 
# #   ylim(c(0.0, 0.7)) + 
# #   ylab("Rate") +
# #   theme(plot.title = element_text(hjust = 0.5),
# #         axis.title.x = element_blank(),
# #         axis.text.x = element_blank(),
# #         axis.title.y = element_text(angle = 90)) +
# #   ggtitle("Proposed rates")
# 
# var_sp <- models2df(model_set0_frommu) %>%
#   filter(rate == "Speciation") %>%
#   group_by(Time) %>%
#   summarize("v" = var(value))
# var_sp$rate <- "Speciation"
# 
# var_ex <- models2df(model_set0_fromlambda) %>%
#   filter(rate == "Extinction") %>%
#   group_by(Time) %>%
#   summarize("v" = var(value))
# var_ex$rate <- "Extinction"
# var_df <- bind_rows(var_sp, var_ex)
# 
# p3b <- ggplot(var_df, aes(x = Time, y = v, color = rate)) +
#   geom_line() +
#   theme_classic() +
#   scale_x_reverse() +
#   ylab("Var[rate]") +
#   theme(legend.position = c(0.2, 0.7),
#         legend.background = element_blank(),
#         plot.title = element_text(hjust = 0.5),
#         legend.title = element_blank()) +
#   scale_color_manual(values = c("orange", "black")) +
#   xlab("Time (Ma)") +
#   ggtitle("Variance in inferred rate")
# 
# 
# p0a + p0b + 
#   p1a + p1b + 
#   p3a + p3b +
#   plot_layout(ncol = 2)



########################################################
##
##
##
##
##
########################################################

# source("scripts/hypothetical_models_frommu.R")
# plot(cgs_frommu$up$reference)

# source("scripts/hypothetical_models_fromlambda.R")
# plot(cgs_fromlambda$up$reference)

######
# 
# up_models_frommu <- models2df(cgs_frommu$up)
# up_models_fromlambda <- models2df(cgs_fromlambda$up)
# 
# p2 <- up_models_frommu %>%
#   filter(name == "reference") %>%
#   filter(rate %in% c("Speciation", "Extinction")) %>%
#   ggplot(aes(x = Time, y = value, color = rate)) +
#   geom_line() +
#   theme_classic() +
#   scale_x_reverse() +
#   scale_color_manual(values = c("orange", "black")) +
#   theme(legend.position = c(0.2, 0.8),
#         legend.title = element_blank()) +
#   ylab("Rate") +
#   xlab("Time (Ma)") +
#   ylim(c(0.0, 0.6)) +
#   ggtitle("Proposed extinction rate")
# 
# p3 <- up_models_fromlambda %>%
#   filter(name == "reference") %>%
#   filter(rate %in% c("Speciation", "Extinction")) %>%
#   ggplot(aes(x = Time, y = value, color = rate)) +
#   geom_line() +
#   theme_classic() +
#   scale_x_reverse() +
#   scale_color_manual(values = c("orange", "black")) +
#   theme(legend.position = c(0.2, 0.8),
#         legend.title = element_blank()) +
#   ylab("Rate") +
#   xlab("Time (Ma)") +
#   ylim(c(0.0, 0.6)) +
#   ggtitle("Proposed speciation rate")
# 
# (p2 | p3) / p1
# 
# 
# abs_dev <- function(df, name, rate = "Speciation") {
#   df_ref <- df %>%
#     filter(rate == rate) %>%
#     filter(name == "reference")
#   
#   df2 <- df %>%
#     filter(name == name) %>%
#     filter(rate == rate)
#   
#   abs_deviation <- df2 %>%
#     (function(e) e$value - df_ref$value) %>%
#     abs()
#   
#   rel_deviation <- df2 %>%
#     (function(e) (e$value - df_ref$value)/e$value) %>%
#     abs()
#   
#   res <- tibble(
#     "Time" = df2$Time,
#     "abs_deviation" = abs_deviation,
#     "rel_deviation" = rel_deviation,
#     "name" = name
#   )
#   return(res)
# }
# 
# deviation_sp  <- lapply(unique(up_models_frommu$name), function(name) abs_dev(up_models_frommu, name)) %>%
#   do.call(bind_rows, .)
# 
# deviation_ex  <- lapply(unique(up_models_fromlambda$name), function(name) abs_dev(up_models_fromlambda, name)) %>%
#   do.call(bind_rows, .)
# 
# ggplot(deviation_sp, aes(x = Time, y = name, fill = abs_deviation)) + 
#   geom_raster()
# 


