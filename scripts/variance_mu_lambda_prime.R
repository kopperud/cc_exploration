## Suppose we set up a reference model, with
## 
## \lambda(t) = sigmoidally increasing
## \mu(t) = 0.28



## And conversely,
## \lambda(t) = 0.28
## \mu(t) = sigmoidally increasing




logistic_shift <- function(min, max, steepness, x0){
  res <- function(x) {
    ((max-min) / (1 + steepness^(-1*(x - x0)))) + min
  }
  return(res)
}

## settings
phy <- read.tree("trees/condamine_etal_2019_Picidae.tre")
height <- max(node.depth.edgelength(phy))
times <- seq(0, height, length.out = 500)

mu0 <- logistic_shift(0.1, 0.6, 0.1, height/2)
lambda0 <- function(t) 0.35 + 0*t
sigmoidal_mu <- CRABS::create.model(func_spec0 = lambda0, mu0, times = times)

mu1 <- function(t) 0.35 + 0*t
lambda1 <- logistic_shift(0.1, 0.6, 0.1, height/2)
sigmoidal_lambda <- CRABS::create.model(func_spec0 = lambda1, mu1, times = times)

lambda_proposals0 <- list()
times_fine <- seq(0, height, length.out = 50)
for (i in 1:100){
  lambda_proposals0[[i]] <- sample.basic.models(times = times_fine, 
                                               rate0 = lambda0(0.0), 
                                               "MRF", 
                                               MRF.type = "GMRF", 
                                               fc.mean = 2.0, 
                                               min.rate = 0.05, 
                                               max.rate = 0.65,
                                               mrf.sd.scale = 5.0)
  
}

model_set0_fromlambda <- congruent.models(sigmoidal_mu, lambdas = lambda_proposals0, keep_ref = FALSE)
model_set0_frommu <- congruent.models(sigmoidal_mu, mus = lambda_proposals0, keep_ref = FALSE)

p0a <- model2df(sigmoidal_mu) %>%
  filter(rate %in% c("Speciation", "Extinction")) %>%
  ggplot(aes(x = Time, y = value, color = rate)) +
  geom_line() +
  theme_classic() +
  scale_x_reverse() +
  scale_color_manual(values = c("orange", "black")) +
  theme(legend.position = c(0.2, 0.8),
        legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5),
        axis.title.x = element_blank(),
        axis.text.x = element_blank()) +
  ylab("Rate") +
  xlab("Time (Ma)") +
  ylim(c(0.0, 0.6)) +
  ggtitle("Reference model")
  

p1a <- plot(model_set0_frommu)[[2]] + 
  ylim(c(0.0, 0.7)) + 
  ylab("Rate") +
  theme(plot.title = element_text(hjust = 0.5),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.y = element_text(angle = 90)) +
  ggtitle("Proposed rates")

var_sp <- models2df(model_set0_frommu) %>%
  filter(rate == "Speciation") %>%
  group_by(Time) %>%
  summarize("v" = var(value))
var_sp$rate <- "Speciation"

var_ex <- models2df(model_set0_fromlambda) %>%
  filter(rate == "Extinction") %>%
  group_by(Time) %>%
  summarize("v" = var(value))
var_ex$rate <- "Extinction"
var_df <- bind_rows(var_sp, var_ex)

p3a <- ggplot(var_df, aes(x = Time, y = v, color = rate)) +
  geom_line() +
  theme_classic() +
  scale_x_reverse() +
  #scale_y_log10() +
  ylab("Var[rate]") +
  theme(legend.position = c(0.2, 0.7),
        legend.background = element_blank(),
        plot.title = element_text(hjust = 0.5),
        legend.title = element_blank()) +
  scale_color_manual(values = c("orange", "black")) +
  xlab("Time (Ma)") +
  ggtitle("Variance in inferred rate")


#p0a / p2a / p3a

p0b <- model2df(sigmoidal_lambda) %>%
  filter(rate %in% c("Speciation", "Extinction")) %>%
  ggplot(aes(x = Time, y = value, color = rate)) +
  geom_line() +
  theme_classic() +
  scale_x_reverse() +
  scale_color_manual(values = c("orange", "black")) +
  theme(legend.position = c(0.2, 0.8),
        legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5),
        axis.title.x = element_blank(),
        axis.text.x = element_blank()) +
  xlab("Time (Ma)") +
  ylim(c(0.0, 0.6)) +
  ggtitle("Reference model")

lambda_proposals1 <- list()
times_fine <- seq(0, height, length.out = 50)
for (i in 1:100){
  lambda_proposals1[[i]] <- sample.basic.models(times = times_fine, 
                                                rate0 = lambda1(0.0), 
                                                "MRF", 
                                                MRF.type = "GMRF", 
                                                fc.mean = 2.0, 
                                                min.rate = 0.05, 
                                                max.rate = 0.65,
                                                mrf.sd.scale = 5.0)
  
}


model_set1_fromlambda <- congruent.models(sigmoidal_lambda, lambdas = lambda_proposals1, keep_ref = FALSE)
model_set1_frommu <- congruent.models(sigmoidal_lambda, mus = lambda_proposals1, keep_ref = FALSE)

#p2b <- plot(model_set1_frommu)[[2]] + ylim(c(0.0, 0.7))
p1b <- plot(model_set1_frommu)[[2]] +
  ylim(c(0.0, 0.7)) +
  theme(plot.title = element_text(hjust = 0.5),
        axis.title.x = element_blank(),
        axis.text.x = element_blank()) +
   ggtitle("Proposed rates")

var_sp <- models2df(model_set1_frommu) %>%
  filter(rate == "Speciation") %>%
  group_by(Time) %>%
  summarize("v" = var(value))
var_sp$rate <- "Speciation"

n_bins <- 50
time_bin_start <- head(seq(-0.1, height+0.01, length.out = n_bins+1), n = -1)
time_bin_end <- tail(seq(-0.1, height+0.01, length.out = n_bins+1), n = -1)
bin_names <- paste0("bin", 1:n_bins)

which_time_bin <- function(x){
  bin_names[x > time_bin_start & x < time_bin_end]
}
sapply(x, which_time_bin)

var_ex <- models2df(model_set1_fromlambda) %>%
  filter(rate == "Extinction") %>%
  #mutate("bin" = sapply(Time, which_time_bin)) %>%
  #group_by(bin) %>%
  group_by(Time) %>%
  summarize("v" = var(value))# %>%
  #mutate("nbin" = readr::parse_number(bin))
var_ex$rate <- "Extinction"
var_df <- bind_rows(var_sp, var_ex)

p3b <- ggplot(var_df, aes(x = Time, y = v, color = rate)) +
#p3b <- ggplot(var_df, aes(x = nbin, y = v, color = rate)) +
  geom_line() +
  theme_classic() +
  scale_x_reverse() +
  #scale_y_log10() +
  theme(legend.position = c(0.2, 0.7),
        legend.background = element_blank(),
        legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5),
        axis.title.y = element_blank()) +
  scale_color_manual(values = c("orange", "black")) +
  xlab("Time (Ma)") +
  ggtitle("Variance in inferred rate")

p0a + p0b + 
  p1a + p1b + 
  p3a + p3b +
  plot_layout(ncol = 2)


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

p2 <- up_models_frommu %>%
  filter(name == "reference") %>%
  filter(rate %in% c("Speciation", "Extinction")) %>%
  ggplot(aes(x = Time, y = value, color = rate)) +
  geom_line() +
  theme_classic() +
  scale_x_reverse() +
  scale_color_manual(values = c("orange", "black")) +
  theme(legend.position = c(0.2, 0.8),
        legend.title = element_blank()) +
  ylab("Rate") +
  xlab("Time (Ma)") +
  ylim(c(0.0, 0.6)) +
  ggtitle("Proposed extinction rate")

p3 <- up_models_fromlambda %>%
  filter(name == "reference") %>%
  filter(rate %in% c("Speciation", "Extinction")) %>%
  ggplot(aes(x = Time, y = value, color = rate)) +
  geom_line() +
  theme_classic() +
  scale_x_reverse() +
  scale_color_manual(values = c("orange", "black")) +
  theme(legend.position = c(0.2, 0.8),
        legend.title = element_blank()) +
  ylab("Rate") +
  xlab("Time (Ma)") +
  ylim(c(0.0, 0.6)) +
  ggtitle("Proposed speciation rate")

(p2 | p3) / p1


abs_dev <- function(df, name, rate = "Speciation") {
  df_ref <- df %>%
    filter(rate == rate) %>%
    filter(name == "reference")
  
  df2 <- df %>%
    filter(name == name) %>%
    filter(rate == rate)
  
  abs_deviation <- df2 %>%
    (function(e) e$value - df_ref$value) %>%
    abs()
  
  rel_deviation <- df2 %>%
    (function(e) (e$value - df_ref$value)/e$value) %>%
    abs()
  
  res <- tibble(
    "Time" = df2$Time,
    "abs_deviation" = abs_deviation,
    "rel_deviation" = rel_deviation,
    "name" = name
  )
  return(res)
}

deviation_sp  <- lapply(unique(up_models_frommu$name), function(name) abs_dev(up_models_frommu, name)) %>%
  do.call(bind_rows, .)

deviation_ex  <- lapply(unique(up_models_fromlambda$name), function(name) abs_dev(up_models_fromlambda, name)) %>%
  do.call(bind_rows, .)

ggplot(deviation_sp, aes(x = Time, y = name, fill = abs_deviation)) + 
  geom_raster()



