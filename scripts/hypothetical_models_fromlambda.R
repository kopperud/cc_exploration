###################################
#
#    SYNTHETIC DATA
#       proposals from lambda
#
###################################

library(CRABS)
library(ape)
library(dplyr)
library(tidyr)
library(pracma)
library(patchwork)

logistic_shift <- function(min, max, steepness, x0){
  res <- function(x) {
    ((max-min) / (1 + 1.0*steepness^(-1*(x - x0)))) + min
  }
  return(res)
}

lseq <- function (from = 1, to = 1e+05, length.out = 6) {
  exp(seq(log(from), log(to), length.out = length.out))
}

expseq <- function (from = 1, to = 1e+05, length.out = 6) {
  log(seq(exp(from), exp(to), length.out = length.out))
}


## settings
setwd("~/projects/cc_exploration")
phy <- read.tree("trees/condamine_etal_2019_Picidae.tre")
height <- max(node.depth.edgelength(phy))
times <- seq(0, height, length.out = 500)
n <- 5

lambda0 <- 0.28
lambda_reference <- function(t) 0.0*t + lambda0


foo_rate <- list()

## Constants
# bs <- seq(0.1, 0.6, length.out = n)
# constant <- lapply(bs, function(b) function(t) b)
# names(constant) <- paste0("constant", seq_along(constant))
# foo_rate[["constant"]] <- constant

## Modal
bs <- seq(0.1, 0.5, length.out = n)
modal <- lapply(bs, function(b) function(t) b * exp((-(t - height*0.5)^2) / 6) + lambda0)
names(modal) <- paste0("modal", seq_along(modal))
foo_rate[["modal"]] <- modal

## Linear-up
bs <- seq(0.001, 0.003, length.out = n)
linear_up <- sapply(bs, function(b) function(t) lambda0 - b * t)
names(linear_up) <- paste0("linear_up", seq_along(linear_up))
foo_rate[["linear_up"]] <- linear_up

## Linear-down
bs <- seq(0.001, 0.003, length.out = n)
linear_down <- sapply(bs, function(b) function(t) lambda0 + b * t)
names(linear_down) <- paste0("linear_down", seq_along(linear_down))
foo_rate[["linear_down"]] <- linear_down

## Up-shift
steepnesses <- lseq(0.1, 0.9, length.out = n)
up <- lapply(steepnesses, function(s) logistic_shift(min=0.02, max=0.28, steepness = s, x0 = height/2.0))
up <- lapply(up, function(f) function(t) f(t) + (lambda0 - f(0.0)))
names(up) <- paste0("up", seq_along(up))
foo_rate[["up"]] <- up

## Down-shift
steepnesses <- lseq(1.1, 2.5, length.out = n)
down <- lapply(steepnesses, function(s) logistic_shift(min=0.28, max=0.54, steepness = s, x0 = height/2.0))
down <- lapply(down, function(f) function(t) f(t) + (lambda0 - f(0.0)))
names(down) <- paste0("down", seq_along(down))
foo_rate[["down"]] <- down

## Exponentially increasing
bs <- seq(0.05, 0.2, length.out = n)
exp_up <- sapply(bs, function(b) function(t) lambda0*exp(-b*t))
names(exp_up) <- paste0("exp_up", seq_along(exp_up))
foo_rate[["exp_up"]] <- exp_up

## Exponential decreasing
bs <- seq(0.05, 0.2, length.out = n)
exp_down <- sapply(bs, function(b) function(t) exp(-b*(height-t))*lambda0 + lambda0)
exp_down <- lapply(exp_down, function(f) function(t) f(t) + (lambda0 - f(0.0)))
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

# p <- ggplot(pdata, aes(x = time, y = rate, linetype = factor(subitemidx), color = factor(is_three))) +
#   facet_wrap(.~ item, nrow = 2) +
#   geom_line() +
#   theme_classic() +
#   scale_x_reverse() +
#   scale_color_manual(values = c("gray", "red")) +
#   scale_linetype_manual(values = rep("solid", 6))
# plot(p)

sapply(foo_rate, function(models) lapply(models, function(f) f(0.0)))




## Reference model
## What is the average extinction rate of the proposed rate functions?
#avg_mu <- sapply(unlist(foo_rate), function(f) quadgk(f, 0, height)/height)
#median(avg_mu)
#mean(avg_mu)
## let's choose 0.28

ref_foos <- list(
  #"constant" = foo_rate$constant$constant3,
  "modal" = foo_rate$modal$modal3,
  "linear_up" = foo_rate$linear_up$linear_up3,
  "linear_down" = foo_rate$linear_down$linear_down3,
  "up" = foo_rate$up$up3,
  "down" = foo_rate$down$down3,
  "exp_up" = foo_rate$exp_up$exp_up3,
  "exp_down" = foo_rate$exp_down$exp_down3
)

lambda <- function(t) 0.28
references <- lapply(ref_foos, function(foo) create.model(func_spec0 = lambda, func_ext0 = foo, times = times))

## CONGRUENCE CLASSES
cgs_fromlambda <- lapply(references, function(ref) congruent.models(ref, lambdas = unlist(foo_rate), keep_ref = TRUE))
cgs <- cgs_fromlambda
group_names <- factor(c(names(foo_rate), "reference"), levels = c(names(foo_rate), "reference"))
ps <- lapply(cgs,
             function(cg) summarize.trends(cg, threshold = 0.02, 
                                           rate_name = "mu",
                                           group_names = group_names)); names(ps) <- names(foo_rate)

## How many models included negative rates?
n_negative <- function(models){
  time1 <- models[[1]]$times
  res <- 0
  
  for (i in seq_along(models)){
    if (any(models[[i]]$mu(time1) < 0)){
      res <- res + 1
    }
  }
  return(res)
}
sum(sapply(cgs, n_negative)) / sum(sapply(cgs, length)-1)

## Labeller function
lf <- c("reference" = "*",
        "modal" = "M",
        "exp_up" = "E+",
        "exp_down" = "E-",
        #"constant" = "C",
        "up" = "S+",
        "down" = "S-",
        "linear_up" = "L+",
        "linear_down" = "L-")

ll <- list()
for (i in seq_along(ps)){
  item <- names(ps)[i]
  tmp_plot <- ps[[i]][[2]]
  tmp_plot <- tmp_plot + 
    facet_grid(group_name~., scales="free_y", 
               space="free_y", switch = "y", 
               labeller = labeller(group_name = lf)) +
    theme(axis.text.x = element_blank(),
          axis.title.x = element_blank(),
          legend.position = "none")
  
  if (item %in% c("linear_up", "modal")){
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

lf2 <- c(
  "modal" = "Modal", 
  "linear_up" = "Linear down", 
  "linear_down" = "Linear up", 
  "up" = "Sigmoidal up", 
  "down" = "Sigmoidal down", 
  "exp_up" = "Exponential up", 
  "exp_down" = "Exponential down"
  )

foobar <- function(pdata, item1){
  d <- pdata %>%
    dplyr::filter(item == item1)
  
  p <- 
    ggplot(d, aes(x = time, y = rate,
                  linetype = factor(subitemidx), 
                  color = factor(is_three))) +
    geom_line() +
    theme_classic() +
    scale_x_reverse() +
    scale_color_manual(values = c("gray", "red")) +
    scale_linetype_manual(values = rep("solid", 6)) +
    labs(title = lf[item1]) +
    theme(legend.position = "none",
          axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          plot.title = element_text(hjust = 0.5)) +
    ylim(c(0.0, 0.8))
  
  if (item1 %in% c("linear_up", "modal")){
    p <- p + ylab("Rate")
  }else{
    p <- p + 
      theme(axis.title.y = element_blank(),
            axis.text.y = element_blank())
  }
  
  return(p)
}

ps2 <- lapply(names(foo_rate), function(rate) foobar(pdata, rate)); names(ps2) <- names(foo_rate)


p_combined <- ps2$modal + ps2$up + ps2$down +plot_spacer() + 
  #ll$modal + ll$down + ll$down + plot_spacer() +
  ll$modal + ll$up + ll$down + plot_spacer() + 
  ps2$linear_up + ps2$linear_down + ps2$exp_up + ps2$exp_down +
  ll$linear_up + ll$linear_down + ll$exp_up + ll$exp_down +
  plot_layout(ncol = 4)
  

ggsave("figures/suppmat/artificial_fromlambda_hypotheses.pdf", p_combined, width = 170, height = 170, units = "mm")


source("scripts/utils.R")

for (cg_name in names(cgs)){
  print(paste("Making spaghetti plot for", cg_name,"."))
  
  p <- spaghetti_plots(cgs[[cg_name]])
  ggsave(paste0("figures/suppmat/spaghetti_", cg_name, "_fromlambda.pdf"),
         p, width = 100, height = 130, units = "mm")
}



