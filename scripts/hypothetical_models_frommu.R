###################################
#
#    SYNTHETIC DATA
#
###################################

library(ACDC)
library(ape)
library(dplyr)
library(tidyr)
library(pracma)
library(patchwork)

logistic_shift <- function(min, max, steepness, x0){
  res <- function(x) {
    ((max-min) / (1 + steepness^(-1*(x - x0)))) + min
  }
  return(res)
}

lseq <- function (from = 1, to = 1e+05, length.out = 6) {
  exp(seq(log(from), log(to), length.out = length.out))
}

expseq <- function (from = 1, to = 1e+05, length.out = 6) {
  log(seq(exp(from), exp(to), length.out = length.out))
}

concat.factor <- function (...) {
  as.factor(do.call(c, lapply(list(...), as.character)))
}


## settings
phy <- read.tree("trees/condamine_etal_2019_Picidae.tre")
height <- max(node.depth.edgelength(phy))
times <- seq(0, height, length.out = 500)
n <- 5

foo_rate <- list()

## Constants
bs <- seq(0.1, 0.6, length.out = n)
constant <- lapply(bs, function(b) function(t) b)
names(constant) <- paste0("constant", seq_along(constant))
foo_rate[["constant"]] <- constant

## Modal
bs <- seq(0.1, 0.5, length.out = n)
modal <- lapply(bs, function(b) function(t) b * exp((-(t - height*0.5)^2) / 6) + 0.1)
names(modal) <- paste0("modal", seq_along(modal))
foo_rate[["modal"]] <- modal

## Linear-up
bs <- seq(0.001, 0.007, length.out = n)
linear_up <- sapply(bs, function(b) function(t) 0.6 - b * t)
names(linear_up) <- paste0("linear_up", seq_along(linear_up))
foo_rate[["linear_up"]] <- linear_up

## Linear-down
bs <- seq(0.001, 0.009, length.out = n)
linear_down <- sapply(bs, function(b) function(t) 0.1 + b * t)
names(linear_down) <- paste0("linear_down", seq_along(linear_down))
foo_rate[["linear_down"]] <- linear_down

## Up-shift
steepnesses <- lseq(0.1, 0.9, length.out = n)
up <- lapply(steepnesses, function(s) logistic_shift(min=0.1, max=0.6, steepness = s, x0 = height/2.0))
names(up) <- paste0("up", seq_along(up))
foo_rate[["up"]] <- up

## Down-shift
steepnesses <- lseq(1.12, 5, length.out = n)
down <- lapply(steepnesses, function(s) logistic_shift(min=0.1, max=0.6, steepness = s, x0 = height/2.0))
names(down) <- paste0("down", seq_along(down))
foo_rate[["down"]] <- down

## Exponentially increasing
bs <- seq(0.05, 0.2, length.out = n)
exp_up <- sapply(bs, function(b) function(t) 0.1 + 0.3*exp(-b*t))
names(exp_up) <- paste0("exp_up", seq_along(exp_up))
foo_rate[["exp_up"]] <- exp_up

## Exponential decreasing
exp_down <- sapply(bs, function(b) function(t) 0.1 + 0.3*exp(-b*(height-t)))
names(exp_down) <- paste0("exp_down", seq_along(exp_down))
foo_rate[["exp_down"]] <- exp_down

l <- list()
for (i in seq_along(foo_rate)){
  item <- names(foo_rate)[i]
  print(item)
  foos <- foo_rate[[i]]
  
  df <- as_tibble(sapply(foos, function(foo) sapply(times, foo)))
  df$time <- times
  
  df <- gather(df, key = "subitem", value = "rate", -time)
  df$item <- item
  df$subitemidx <- readr::parse_number(df$subitem)
  df$is_three <- df$subitemidx == 3
  
  l[[i]] <- df
}
pdata <- bind_rows(l)

p <- ggplot(pdata, aes(x = time, y = rate, linetype = factor(subitemidx), color = factor(is_three))) +
  facet_wrap(.~ item, nrow = 2) +
  geom_line() +
  theme_classic() +
  scale_x_reverse() +
  scale_color_manual(values = c("gray", "blue")) +
  scale_linetype_manual(values = rep("solid", 6))
plot(p)

#
  

## Reference model
## What is the average extinction rate of the proposed rate functions?
avg_mu <- sapply(unlist(foo_rate), function(f) quadgk(f, 0, height)/height)

median(avg_mu)
mean(avg_mu)
## let's choose 0.28

ref_foos <- list(
  "constant" = foo_rate$constant$constant3,
  "modal" = foo_rate$modal$modal3,
  "linear_up" = foo_rate$linear_up$linear_up3,
  "linear_down" = foo_rate$linear_down$linear_down3,
  "up" = foo_rate$up$up3,
  "down" = foo_rate$down$down3,
  "exp_up" = foo_rate$exp_up$exp_up3,
  "exp_down" = foo_rate$exp_down$exp_down3
)
mu <- function(t) 0.28
#mu <- function(t) 0.12
references <- lapply(ref_foos, function(foo) create.model(func_spec0 = foo, func_ext0 = mu, times = times))


## CONGRUENCE CLASSES
cgs <- lapply(references, function(ref) congruent.models(ref, mus = unlist(foo_rate), keep_ref = TRUE))
group_names <- factor(c(names(foo_rate), "reference"), levels = c(names(foo_rate), "reference"))
ps <- lapply(cgs,
             function(cg) summarize.trends(cg, threshold = 0.02, group_names = group_names)); names(ps) <- names(foo_rate)

## Test exhibit with smaller mu
mu2 <- list(
  function(t) 0.08,
  function(t) 0.18,
  function(t) 0.28,
  function(t) 0.38
)
ref2 <- lapply(mu2, function(mu) create.model(func_spec0 = ref_foos$down, 
                                              func_ext0 = mu, 
                                              times = times))
cgs2 <- lapply(ref2, function(ref) congruent.models(ref, mus = unlist(foo_rate), keep_ref = TRUE))
ps2 <- lapply(cgs2,
              function(cg) summarize.trends(cg, threshold = 0.02, group_names = group_names)); names(ps) <- names(foo_rate)

l <- list()
for (i in seq_along(cgs2)){
  df <- model2df(cgs2[[i]]$reference)
  df[["mu_ref"]] <- cgs2[[i]]$reference$mu(0.0)
  l[[i]] <- df
}; res5 <- bind_rows(l)

res5 %>%
  filter(rate == "Pulled net-diversification") %>%
  ggplot(aes(x = Time, y = value, color = factor(mu_ref))) +
  geom_line() +
  scale_x_reverse() +
  theme_classic()

## Labeller function
lf <- c("reference" = "*",
        "modal" = "M",
        "exp_up" = "E+",
        "exp_down" = "E-",
        "constant" = "C",
        "up" = "S+",
        "down" = "S-",
        "linear_up" = "L+",
        "linear_down" = "L-")

p5s <- list()
for (i in seq_along(ps2)){
  p5 <- ps2[[i]][[2]] + 
    facet_grid(group_name~., scales="free_y", 
               space="free_y", switch = "y", 
               labeller = labeller(group_name = lf)) +
    theme(legend.title = element_blank(),
          plot.title = element_text(hjust = 0.5),
          axis.text.x = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          legend.position = "none") +
    labs(y = "Congruent models",
         x = "Time (Ma)")
         #title = paste0(latex2exp::TeX("$\\mu$"), mu2[[i]](0.0)))
  p5s[[i]] <- p5
}

sigmoid_down_variable_mu <- p5s[[1]] + labs(title = latex2exp::TeX("Reference $\\mu = 0.08$")) + theme(axis.title.y = element_text(angle = 90),
                                                                           legend.position = c(0.2, 0.5)) +
  p5s[[2]] + labs(title = latex2exp::TeX("Reference $\\mu = 0.18$")) +
  p5s[[3]] + labs(title = latex2exp::TeX("Reference $\\mu = 0.28$")) + theme(axis.title.y = element_text(angle = 90),
                                                                             axis.text.x = element_text(),
                                                                             axis.title.x = element_text()) + 
  p5s[[4]] + labs(title = latex2exp::TeX("Reference $\\mu = 0.38$")) + theme(axis.text.x = element_text(),
                                                                             axis.title.x = element_text()) +
  plot_layout(ncol = 2)

ggsave("figures/suppmat/sigmoidal_down_variable_mu.pdf", sigmoid_down_variable_mu, width = 150, height = 150, units = "mm")



ll <- list()
for (i in seq_along(ps)){
  item <- names(ps)[i]
  tmp_plot <- ps[[i]][[2]]
  tmp_plot <- tmp_plot + 
    facet_grid(group_name~., scales="free_y", 
               space="free_y", switch = "y", 
               labeller = labeller(group_name = lf)) +
    #labs(title = paste0("Reference = ", group_names[i])) +
    theme(axis.text.x = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          legend.position = "none")
  
  if (item %in% c("linear_up", "constant")){
    tmp_plot <- tmp_plot + ylab("Models") +
      theme(axis.title.y = element_text(angle = 90))
  }else{
    tmp_plot <- tmp_plot + 
      theme(axis.title.y = element_blank())
  }
  
  if (item %in% c("linear_down", "linear_up", "exp_up", "exp_down")){
    tmp_plot <- tmp_plot +
      theme(axis.title.x = element_text(),
            axis.text.x = element_text()) +
      xlab("Time (Ma)")
  }
  
  ll[[i]] <- tmp_plot
}; names(ll) <- names(ps)

########################
## Rateplots

foobar <- function(pdata, item1, ref_only = FALSE){
  d <- pdata %>%
    dplyr::filter(item == item1)
  
  if (ref_only){
    d <- d %>%
      dplyr::filter(subitem == paste0(item1, "3"))
  }

  p <- 
    ggplot(d, aes(x = time, y = rate, linetype = factor(subitemidx), color = factor(is_three))) +
    geom_line() +
    theme_classic() +
    scale_x_reverse() +
    scale_linetype_manual(values = rep("solid", 6)) +
    labs(title = lf[item1]) +
    theme(legend.position = "none",
          axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          plot.title = element_text(hjust = 0.5)) +
    ylim(c(0.0, 0.8))
  
  if(ref_only){
    p <- p +
      scale_color_manual(values = c("black"))
  }else{
    p <- p +
      scale_color_manual(values = c("gray", "blue"))
  }
  
  if (item1 %in% c("linear_up", "constant")){
    p <- p + 
      ylab("Rate")
  }else{
    p <- p + 
      theme(axis.title.y = element_blank(),
            axis.text.y = element_blank(),
            )
  }
  
  return(p)
}

pdata %>%
  filter(item == "constant") %>%
  (function(e) e$subitem) %>%
  unique()

ps2 <- lapply(names(foo_rate), function(rate) foobar(pdata, rate)); names(ps2) <- names(foo_rate)
ps3 <- lapply(names(foo_rate), function(rate) foobar(pdata, rate, ref_only = TRUE)); names(ps3) <- names(foo_rate)

p_combined <- ps2$constant + ps2$modal + ps2$up + ps2$down +
  ll$constant + ll$modal + ll$up + ll$down + 
  ps2$linear_up + ps2$linear_down + ps2$exp_up + ps2$exp_down +
  ll$linear_up + ll$linear_down + ll$exp_up + 
  theme(legend.position = c(0.45, 0.8),
        legend.title = element_blank(),
        legend.key.size = unit(0.8,"line"),
        legend.background = element_blank()) + ll$exp_down +
  plot_layout(ncol = 4)

ggsave("figures/suppmat/artificial_hypotheses_frommu.pdf", p_combined, width = 170, height = 170, units = "mm")

p_small <- ps3$linear_up + ylab("Speciation rate") +
  ps3$exp_down + ps3$down +
  ll$linear_up + theme(legend.position = c(0.38, 0.67),
                       legend.title = element_blank(),
                       legend.background = element_blank(),
                       legend.key.size = unit(3, "mm"),
                       legend.key = element_rect()) +
  ll$exp_down + 
  ll$down + theme(axis.title.x = element_text(),
                                               axis.ticks.x = element_line(),
                                               axis.text.x = element_text()) + xlab("Time (Ma)") +
  plot_layout(ncol = 3)
p_small
ggsave("figures/ms/artificial_hypotheses_frommu_reduced.pdf", p_small, width = 140, height = 100, units = "mm")




source("scripts/utils.R")

for (cg_name in names(cgs)){
  print(paste("Making spaghetti plot for", cg_name,"."))
  
  p <- spaghetti_plots(cgs[[cg_name]], columns = c("Extinction", "Speciation"))
  ggsave(paste0("figures/suppmat/spaghetti_", cg_name, "_frommu.pdf"),
         p, width = 100, height = 130, units = "mm")
}

