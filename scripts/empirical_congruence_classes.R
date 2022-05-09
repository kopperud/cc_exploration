## Congruence classes
## On Picidae -- The Woodpeckers
library(readr)
library(ACDC)
library(tibble)
library(dplyr)
library(RColorBrewer)
library(deSolve)
library(pracma)
library(ape)
library(ggplot2)
library(tidyr)
library(ggtree)
library(patchwork)
library(ggtree)

#setwd("~/projects/diversify/")
source("scripts/utils.R")

reduced_datasets <- c("onagraceae_final_map",
                      "condamine_etal_2019_Pipridae",
                      "condamine_etal_2019_Fringillidae",
                      "condamine_etal_2019_Picidae2",
                      "condamine_etal_2019_Viduidae",
                      "condamine_etal_2019_Procellariidae",
                      "condamine_etal_2019_Tyrannidae",
                      "condamine_etal_2019_Muridae",
                      "condamine_etal_2019_Psittacidae",
                      "condamine_etal_2019_Zosteropidae",
                      "condamine_etal_2019_Parulidae",
                      "condamine_etal_2019_Molossidae");


datasets <- lapply(reduced_datasets, read_dataset); names(datasets) <- reduced_datasets

bar <- function(dataset_name, datasets, threshold = 0.02){
  dataset <- datasets[[dataset_name]]
  
  ###############################################################
  #                                                             #
  #          The Diversification rates fitted using RB          #
  #                                                             #
  ###############################################################
  
  scalexformat <- function(x) sprintf("%.0f", abs(round(x, 1)))
  scaleyformat <- function(x) sprintf("%.1f", abs(round(x, 2)))

  ci_plotdata <- bind_rows(dataset$speciation,
                           dataset$extinction)
  ci_plotdata$item <- dplyr::recode(ci_plotdata$item, 
                                    extinction_rate = "Extinction rate", 
                                    speciation_rate = "Speciation rate")
  
  p_revbayesrates <- ci_plotdata %>%
    ggplot(aes(x = time, y = median)) +
    geom_line(aes(color = item)) +
    scale_x_reverse() +
    theme_classic() +
    geom_ribbon(aes(ymin = lower, ymax = upper, fill = item), 
                color = NA, alpha = 0.2) +
    ylab("Rate") +
    xlab("Time (Ma)") +
    theme(legend.position = c(0.3, 0.8),
          legend.title = element_blank(),
          legend.key.size = unit(5, "mm")) +
    #coord_cartesian(xlim = c(max(node.depth.edgelength(dataset$phy)),-0.1)) +
    scale_color_manual(values = c("orange", "black")) +
    scale_fill_manual(values = c("orange", "black"))
  
  ################################################
  #                                              #
  #      Specific hypotheses: rate functions     #
  #                                              #
  ################################################
  
  ylim <- c(0, 0.7)
  
  p_rbfit <- ci_plotdata %>%
    ggplot(aes(x = time, y = median)) +
    #scale_y_continuous(breaks = scales::pretty_breaks(n = 4)) +
    geom_line() +
    scale_x_reverse() +
    facet_grid(cols = vars(item)) +
    theme_classic() +
    geom_ribbon(aes(ymin = lower, ymax = upper), fill = "black", alpha = 0.2) +
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank()) +
    ylab("Reference")
  
  dataset_constants <- c("onagraceae_final_map" = 0.3, 
                         "condamine_etal_2019_Pipridae" = 0, 
                         "condamine_etal_2019_Fringillidae" = 0, 
                         "condamine_etal_2019_Picidae2" = 0, 
                         "condamine_etal_2019_Viduidae" = 0.25, 
                         "condamine_etal_2019_Procellariidae" = 0, 
                         "condamine_etal_2019_Tyrannidae" = 0, 
                         "condamine_etal_2019_Muridae" = 0, 
                         "condamine_etal_2019_Psittacidae" = 0, 
                         "condamine_etal_2019_Zosteropidae" = 0, 
                         "condamine_etal_2019_Parulidae" = 0,
                         "condamine_etal_2019_Molossidae" = 0
  )
  addconstant <- dataset_constants[dataset_name]
  
  m1 <- alternative_lambdas(dataset)
  m2 <- alternative_mus(dataset, addconstant)
  
  times <- m2$constants$reference$times
  

  ps2 <- rateplots(m2, ylim = ylim)
  
  ################################################
  #                                              #
  #      From process: rate functions            #
  #                                              #
  ################################################
  
  HSMRF_mus <- list()
  n_samples <- 20
  set.seed(123)
  
  xtimes <- dataset$reference$times
  for (i in 1:n_samples){
    mu <- ACDC::sample.basic.models(times = xtimes, model = "MRF", MRF.type = "HSMRF", 
                                   min.rate = 0.0, max.rate = 0.3, rate0.median = 0.1 + addconstant, fc.mean = 3.0)
    HSMRF_mus[[i]] <- mu
  }
  
  HSMRF_set <- congruent.models(dataset$reference, mus = HSMRF_mus)
  names(HSMRF_set) <- gsub("model", "HSMRF_", names(HSMRF_set))
  
  plotdata5 <- do.call(bind_rows, mapply(f, HSMRF_set, paste0("HSMRF_", 1:(n_samples+1)), SIMPLIFY = FALSE)) %>%
    pivot_longer(`Speciation rate`:`Extinction rate`, names_to = "rate", values_to = "value")
  
  cols1 <- c("black", head(colorspace::sequential_hcl(palette = "orange", n = length(HSMRF_set)+2), n = -3))
  cols2 <- c("black", head(colorspace::sequential_hcl(palette = "blue", n = length(HSMRF_set)+2), n = -3))

  p3a <- plotdata5 %>%
    dplyr::filter(rate == "Extinction rate") %>%
    dplyr::filter(readr::parse_number(Model) <= 10) %>%
    ggplot(aes(x = times, y = value, color = Model)) +
    geom_line() +
    scale_x_reverse() +
    theme_classic() +
    ylab("HSMRF") +
    xlab("Time (Ma)") +
    scale_color_manual(values = cols1) +
    theme(
      strip.background = element_blank(),
      strip.text.x = element_blank(),
      legend.position =  "none",
      legend.title = element_blank(),
      legend.background = element_blank()
    )
  
  p3b <- plotdata5 %>%
    dplyr::filter(rate == "Speciation rate") %>%
    ggplot(aes(x = times, y = value, color = Model)) +
    geom_line() +
    scale_x_reverse() +
    theme_classic() +
    ylab("HSMRF") +
    xlab("Time (Ma)")+
    scale_color_manual(values = cols2) +
    theme(
      strip.background = element_blank(),
      strip.text.x = element_blank(),
      legend.position =  "none",
      legend.title = element_blank(),
      legend.background = element_blank(),
      axis.title.y = element_blank(),
      axis.text.y = element_blank(),
      axis.line.y = element_blank(),
      axis.ticks.y = element_blank()
    )

  ps2 <- rateplots(m2, 
                   ylabs = c("Constants", "Increasing", "Decreasing", "Modal1", "Modal2", "Modal3"),
                   set_names = c("constants", "exp1", "exp2", "modal1", "modal2", "modal3"),
                   ylim = ylim)
  
  
  f3 <- (p_rbfit) /
    (ps2$constants$`Extinction rate` + ps2$constants$`Speciation rate`) /
    (ps2$exp1$`Extinction rate` + ps2$exp1$`Speciation rate`) /
    (ps2$exp2$`Extinction rate` + ps2$exp2$`Speciation rate`) /
    (ps2$modal1$`Extinction rate` + ps2$modal1$`Speciation rate`) /
    (ps2$modal2$`Extinction rate` + ps2$modal2$`Speciation rate`) /
    (ps2$modal3$`Extinction rate` + 
       theme(axis.text.x = element_blank(),
             axis.title.x = element_blank()) + 
       ps2$modal3$`Speciation rate` + 
       theme(axis.text.x = element_blank(),
             axis.title.x = element_blank())) / 
    (p3a + p3b) & 
    scale_y_continuous(labels = scaleyformat, 
                       breaks = scales::pretty_breaks(n = 4))
  
  ggsave(paste0("figures/suppmat/", dataset_name, "/hyp_from_mu.pdf"), f3, width = 90, height = 170, units = "mm")
  
  
  ################################################
  #                                              #
  #     Proposed models:     lambda -> mu        #
  #                                              #
  ################################################
  
  ci_plotdata$item <- factor(ci_plotdata$item, levels = c("Speciation rate", "Extinction rate"))
  p_rbfit2 <- ci_plotdata %>%
    ggplot(aes(x = time, y = median)) +
    geom_line() +
    scale_x_reverse() +
    facet_grid(cols = vars(item)) +
    theme_classic() +
    geom_ribbon(aes(ymin = lower, ymax = upper), fill = "black", alpha = 0.2) +
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank()) +
    ylab("Reference") +
    coord_cartesian(ylim = ylim) 
  
  ybottom <- -0.5
  
  ps4 <- rateplots(m1, ylim = NA, set_names = c("exp1", "two_inc", "two_dec"),
                   rate_names = c("Speciation rate", "Extinction rate"),
                   ylabs = c("Increasing", "Up-shift", "Down-shift"))

  
  ################################################
  #                                              #
  #    From HSMRF-distribution: lambda -> mu     #
  #                                              #
  ################################################
  
  
  HSMRF_lambdas <- list()
  n_samples <- 20
  set.seed(123)
  
  xtimes <- dataset$reference$times
  for (i in 1:n_samples){
    lambda <- sample.basic.models(times = xtimes, model = "MRF", MRF.type = "HSMRF", min.rate = 0, 
                             rate0 = dataset$reference$lambda(0.0), 
                             max.rate = dataset$reference$lambda(0.0)*4.0, fc.mean = 3.0)
    HSMRF_lambdas[[i]] <- lambda
  }
  
  
  
  HSMRF_set2 <- congruent.models(dataset$reference, lambdas = HSMRF_lambdas)
  names(HSMRF_set2) <- gsub("model", "HSMRF_", names(HSMRF_set2))
  
  plotdata5 <- do.call(bind_rows, mapply(f, HSMRF_set2, paste0("HSMRF_", 1:(n_samples+1)), SIMPLIFY = FALSE)) %>%
    pivot_longer(`Speciation rate`:`Extinction rate`, names_to = "rate", values_to = "value")
  cols <- c("black", head(colorspace::sequential_hcl(palette = "orange", n = length(HSMRF_set2)+2), n = -3))
  
  hsmrf_a <- plotdata5 %>%
    dplyr::filter(rate == "Speciation rate") %>%
    dplyr::filter(readr::parse_number(Model) <= 10) %>%
    ggplot(aes(x = times, y = value, color = Model)) +
    geom_line() +
    scale_x_reverse() +
    theme_classic() +
    ylab("HSMRF") +
    xlab("Time (Ma)") +
    scale_color_manual(values = cols2) +
    theme(
      strip.background = element_blank(),
      strip.text.x = element_blank(),
      legend.position =  "none",
      legend.title = element_blank(),
      legend.background = element_blank()
    )
  
  hsmrf_b <- plotdata5 %>%
    dplyr::filter(rate == "Extinction rate") %>%
    ggplot(aes(x = times, y = value, color = Model)) +
    geom_line() +
    scale_x_reverse() +
    theme_classic() +
    ylab("HSMRF") +
    xlab("Time (Ma)")+
    scale_color_manual(values = cols1) +
    theme(
      strip.background = element_blank(),
      strip.text.x = element_blank(),
      legend.position =  "none",
      legend.title = element_blank(),
      legend.background = element_blank(),
      axis.title.y = element_blank()
    )

  p4a <- (p_rbfit2 + labs(tag = "b)")) & coord_cartesian(ylim = c(0, 0.6))
  p4b <- (ps4$exp1$`Speciation rate` + ps4$exp1$`Extinction rate`) & coord_cartesian(ylim = c(-0.2, 0.6))
  p4c <- (ps4$two_inc$`Speciation rate` + ps4$two_inc$`Extinction rate`) & coord_cartesian(ylim = c(-0.2, 0.5))
  p4d <- (ps4$two_dec$`Speciation rate` + ps4$two_dec$`Extinction rate`) & coord_cartesian(ylim = c(-0.2, 1.5)) &
    theme(axis.title.x = element_blank(), axis.text.x = element_blank())
  p4e <- (hsmrf_a + coord_cartesian(ylim = c(-0.1, 0.5)) +
            hsmrf_b + coord_cartesian(ylim = c(-1.5, 1.5))) 
  
  f4 <-  p4a /
    p4b /
    p4c /
    p4d /
    p4e &
    scale_y_continuous(labels = scaleyformat,
                       breaks = scales::pretty_breaks(n = 4))
  ggsave(paste0("figures/suppmat/", 
                dataset_name, 
                "/hyp_from_lambda.pdf"), 
         f4, width = 90, height = 130, units = "mm")

  ################################################
  #                                              #
  #     Specific hypotheses: summary plots       #
  #                                              #
  ################################################
  
  l <- unlist(m1, recursive = FALSE); class(l) <- c("list", "ACDCset")
  l$reference <- m1$constant$reference
  l$reference <- NULL
  l$linear.reference <- NULL
  l$linear_dec.reference <- NULL
  l$constant.reference <- NULL
  l$two_inc.reference <- NULL
  l$two_dec.reference <- NULL
  l$three.reference <- NULL
  l$exp1.reference <- NULL
  
  l2 <- unlist(m2, recursive = FALSE); class(l2) <- c("list", "ACDCset")
  l2$reference <- NULL
  l2$modal3.reference <- NULL
  l2$constants.reference <- NULL
  l2$exp1.reference <- NULL
  l2$exp2.reference <- NULL
  l2$modal1.reference <- NULL
  l2$modal2.reference <- NULL
  names(l2) <- gsub("\\.model", "", names(l2))
  
  l2 <- c(HSMRF_set, l2); class(l) <- c("list", "ACDCset")
  
  thresholds <- c(0.01, threshold, 0.05)
  
  lf <- c("reference" = "*",
          "modal1" = "M1",
          "modal2" = "M2",
          "modal3" = "M3",
          "exp1" = "I",
          "exp2" = "D",
          "constants" = "C",
          "HSMRF" = "HSMRF")
  
  trend_summaries <- list()
  rate_hypotheses <- list()
  
  for (thresh in thresholds){
    sthresh <- as.character(thresh)
    
    group_names <- factor(c("reference", "constants", "exp1", "exp2", "modal1", "modal2", "modal3", "HSMRF"), 
                          levels = c("reference", "constants", "exp1", "exp2", "modal1", "modal2", "modal3", "HSMRF"))
    p10 <- summarize.trends(l2, threshold = thresh, rm_singleton = TRUE, 
                            group_names = group_names)
    
    p10[[1]] <- p10[[1]] + theme(legend.position = "none",
                                 plot.title = element_text(hjust = 0.5)) + 
      labs(title = "Specific hypotheses", y = latex2exp::TeX("$\\Delta\\lambda$"))
    
    p10[[2]] <- p10[[2]] + 
      facet_grid(group_name~., scales="free_y", space="free_y", switch = "y", labeller = labeller(group_name = lf)) +
      theme(legend.position = "none") + 
      ylab("Models") + 
      xlab("Time (Ma)") +
      plot_layout(ncol = 1)
    
    rate_hypotheses[[sthresh]] <- p10[[2]]
    
    
    p12 <- summarize.trends(HSMRF_set, threshold = thresh, 
                            rm_singleton = FALSE, group_names = c("reference", "HSMRF"))
    
    
    lf2 <- c("reference" = "*",
             "HSMRF" = "HSMRF")
    
    p12[[1]] <- p12[[1]] + theme(legend.position = "none",
                                 plot.title = element_text(hjust = 0.5),
                                 axis.title.y = element_blank()) + 
      labs(title = "Drawn randomly")
    p12[[2]] <- p12[[2]] + 
      facet_grid(group_name~., scales="free_y", space="free_y", switch = "y", labeller = labeller(group_name = lf2)) +
      theme(legend.position = c(0.3, 0.4),
            legend.key.height = unit(3, "mm"),
            legend.background = element_blank(),
            legend.key = element_rect(colour = "black", fill = NA, size = unit(0.1, "mm")),
            legend.title = element_blank(),
            axis.title.y = element_blank()) + 
      xlab("Time (Ma)")

    trend_summary <- 
      p10[[1]] + coord_cartesian(ylim = c(-0.09, 0.1)) + 
      p12[[1]] + coord_cartesian(ylim = c(-0.09, 0.1)) +
      p10[[2]] + 
      p12[[2]] +
      plot_layout(ncol = 2)
    ggsave(paste0("figures/suppmat/", dataset_name, "/trend_summary", thresh, ".pdf"), trend_summary, width = 130, height = 80, units = "mm")
    trend_summaries[[sthresh]] <- trend_summary
  }
  
  print(paste0("finished:  ", dataset_name))
  
  figures <- list(
    "p_revbayesrates" = p_revbayesrates,
    "f3" = f3,
    "f4" = f4,
    "trend_summaries" = trend_summaries,
    "rate_hypotheses" = rate_hypotheses,
    "m1" = m1,
    "m2" = m2
  )
  return(figures)
}


#rd2 <- c("condamine_etal_2019_Psittacidae", "condamine_etal_2019_Tyrannidae", "condamine_etal_2019_Picidae2")
#figures <- lapply(rd2, function(name) bar(name, datasets)); names(figures) <- rd2
figures <- lapply(reduced_datasets, function(name) bar(name, datasets)); names(figures) <- reduced_datasets

## How many trajectories have negative rates?
n_negative <- function(models){
  models1 <- unlist(models, recursive = FALSE)
  time1 <- models1[[1]]$times
  res <- 0
  
  for (i in seq_along(models1)){
    if (any(models1[[i]]$mu(time1) < 0)){
      res <- res + 1
    }
  }
  return(res)
}

c("condamine_etal_2019_Viduidae", "condamine_etal_2019_Parulidae", 
  "condamine_etal_2019_Molossidae", "condamine_etal_2019_Fringillidae", 
  "condamine_etal_2019_Pipridae", "condamine_etal_2019_Procellariidae",
  "condamine_etal_2019_Psittacidae", "condamine_etal_2019_Tyrannidae", 
  "condamine_etal_2019_Picidae2", "condamine_etal_2019_Muridae", 
  "condamine_etal_2019_Zosteropidae", "onagraceae_final_map") %>%
  sapply(function(e) n_negative(figures[[e]]$m1)) %>%
  sum() %>%
  (function(e) e / (12*25))
## Answer is 73%
  


### Multi data summary
plot_multidata <- function(d1, d2, d3,
                          titles = c("Psittacidae", "Tyrannidae", "Picidae"),
                          thresholds = c("0.02", "0.02", "0.02")){
  p1a <- d1$p_revbayesrates +
    labs(y = "Reference*", title = titles[1]) +
    theme(plot.title = element_text(vjust = 0, hjust = 0.5),
          legend.position = c(0.5, 0.8),
          legend.background = element_blank(),
          axis.title.x = element_blank(),
          axis.text.x = element_blank())
  p2a <- d2$p_revbayesrates +
    labs(title = titles[2]) +
    theme(plot.title = element_text(vjust = 0, hjust = 0.5),
          legend.position = "none",
          axis.title.y = element_blank(),
          axis.title.x = element_blank(),
          axis.text.x = element_blank())
  p3a <- d3$p_revbayesrates +
    labs(title = titles[3]) +
    theme(plot.title = element_text(vjust = 0, hjust = 0.5),
          legend.position = "none",
          axis.title.y = element_blank(),
          axis.title.x = element_blank(),
          axis.text.x = element_blank())
  
  p1b <- d1$rate_hypotheses[[thresholds[1]]] + 
    theme(legend.position = "none") +
    labs(y = "Congruent models")
  
  p2b <- d2$rate_hypotheses[[thresholds[2]]] + 
    theme(strip.background = element_blank(),
          strip.text = element_blank(),
          axis.title.y = element_blank(),
          legend.position = "none")
  p3b <- d3$rate_hypotheses[[thresholds[3]]] + 
    theme(strip.background = element_blank(),
          strip.text = element_blank(),
          axis.title.y = element_blank(),
          legend.position = c(0.4, 0.2),
          legend.title = element_blank(),
          legend.background = element_blank(),
          legend.key = element_rect(color="black"))
  
  p <- 
    p1a + p2a + p3a +
    p1b + p2b + p3b +
    plot_layout(ncol = 3,
                heights = c(0.35, 0.65))
  return(p)
}

### Picidae, Tyrannidae, Psittacidae
multidataset <- plot_multidata(figures$condamine_etal_2019_Psittacidae,
                               figures$condamine_etal_2019_Tyrannidae,
                               figures$condamine_etal_2019_Picidae)

## Make them a bit prettier for the publication
multidataset[[1]] <- multidataset[[1]] +
  theme(legend.position = "none")
multidataset[[3]] <- multidataset[[3]] +
  theme(legend.position = c(0.45, 0.9),
        legend.background = element_blank())

multidataset[[4]] <- multidataset[[4]] +
  theme(legend.position = c(0.4, 0.2),
        legend.title = element_blank(),
        legend.key = element_rect(color="black"))
multidataset[[6]] <- multidataset[[6]] +
  theme(legend.position = "none")

ggsave("figures/ms/multidataset.pdf", multidataset, width = 150, height = 120, units = "mm")


rd3 <- c("condamine_etal_2019_Viduidae", "condamine_etal_2019_Parulidae", "condamine_etal_2019_Molossidae")
figures_biggert <- lapply(rd3, function(name) bar(name, datasets, threshold = 0.08)); names(figures_biggert) <- rd3
rd4 <- c("condamine_etal_2019_Fringillidae", "condamine_etal_2019_Pipridae", "condamine_etal_2019_Procellariidae")
figures_smallert <- lapply(rd4, function(name) bar(name, datasets, threshold = 0.01)); names(figures_smallert) <- rd4
figures_smallert$condamine_etal_2019_Procellariidae <- bar("condamine_etal_2019_Procellariidae", datasets, threshold = 0.005)

###############################################################################
#                                                                             #
#       Setting up congruence class exploration for additional datasets       #
#                                                                             #
###############################################################################

md1 <- plot_multidata(figures_smallert$condamine_etal_2019_Fringillidae,
                      figures_smallert$condamine_etal_2019_Pipridae,
                      figures_smallert$condamine_etal_2019_Procellariidae,
                      titles = c("Fringillidae", "Pipridae", "Procellariidae"),
                      thresholds = c("0.01", "0.01", "0.005"))

md2 <- plot_multidata(figures$condamine_etal_2019_Muridae,
                      figures$condamine_etal_2019_Zosteropidae,
                      figures$onagraceae_final_map,
                      titles = c("Muridae", "Zosteropidae", "Onagraceae"))
md2[[1]] <- md2[[1]] + coord_cartesian(y = c(0.0, 1.0))

md3 <- plot_multidata(figures_biggert$condamine_etal_2019_Viduidae,
                      figures_biggert$condamine_etal_2019_Parulidae,
                      figures_biggert$condamine_etal_2019_Molossidae,
                      titles = c("Viduidae", "Parulidae", "Molossidae"),
                      threshold = c("0.08", "0.08", "0.08"))
md3[[4]] <- md3[[4]] + theme(legend.position = c(0.4, 0.3),
                             legend.title = element_blank())
md3[[6]] <- md3[[6]] + theme(legend.position = "none")
md3[[2]] <- md3[[2]] + coord_cartesian(y = c(0.0, 1.2))
md3[[3]] <- md3[[3]] + coord_cartesian(y = c(0.0, 1.2))

ggsave("figures/suppmat/md1.pdf", md1, width = 150, height = 120, units = "mm")
ggsave("figures/suppmat/md2.pdf", md2, width = 150, height = 120, units = "mm")
ggsave("figures/suppmat/md3.pdf", md3, width = 150, height = 120, units = "mm")

###########################################
#                                         #
#       Reading in the posterior          #
#                                         #
###########################################

post_trends <- function(ds, window_size){
  scalexformat <- function(x) sprintf("%.0f", abs(round(x, 1)))
  
  l <- list()
  for (i in 1:4){
    df <- read.table(paste0("output/HSMRF_", ds, "_run_", i, ".log"), header = TRUE)
    
    l[[i]] <- df
  }
  
  df <- bind_rows(l)
  phy <- read.tree(paste0("trees/", ds, ".tre"))
  max_t <- max(node.depth.edgelength(phy))
  
  posterior <- read.RevBayes(df, n_times = 500, n_samples = 110, max_t = max_t)
  
  max_mu0 <- max(sapply(posterior, function(e) e$mu(0.0)))
  print(paste0("max_mu0 (", ds, "): ", max_mu0))
  
  set.seed(1234)
  samples <- sample.congruence.class.posterior(posterior, 10, rate.type = "extinction", 
                                               model = "MRF", MRF.type = "HSMRF", 
                                               mu0.equal = TRUE,
                                               max.rate = max_mu0 + 0.1, min.rate = 0.00001)
  
  trend_data <- summarize.posterior(samples, threshold = 0.02, return_data = TRUE)
  
  bad_samples <- trend_data %>%
    filter(is.na(direction)) %>%
    summarise(unique(paste0(posterior, "_", name))) %>%
    (function(e) e[[1]])

  if (bad_samples[[1]] != "_"){
    bad_samples <- strsplit(bad_samples, "_")

    for (item in bad_samples){
      sample_name <- paste0("posterior", readr::parse_number(item[[1]]))
      model_name <- item[[2]]

      samples[[sample_name]][[model_name]] <- NULL
    }
  }
  
  ## remove nulls
  for (i in seq_along(samples)){
    samples[[i]] <- Filter(Negate(is.null), samples[[i]])
  }
  
  ls <- list()
  for (j in seq_along(window_sizes)){
    window_size <- window_sizes[j]
    p <- summarize.posterior(samples, threshold = 0.02, window_size = window_size)
    
    if (max_t < 30){
      breaks1 <- seq(from = 0, by = 10, length.out = 10)
    }else{
      breaks1 <- seq(from = 0, by = 20, length.out = 10)
    }
    breaks1 <- rev(breaks1)
    
    p <- p +
      xlab("Time (Ma)") +
      scale_x_continuous(labels = scalexformat,
                         breaks = breaks1) +
      coord_cartesian(xlim = c(max_t, 0))
    
    l <- list("p" = p,
              "trend_data" = trend_data,
              "samples" = samples)
    ls[[j]] <- l
  }
  return(ls)
}

window_sizes <- c(10, 25, 55, 80)

posteriors <- list()
for (i in seq_along(rd2)){
  name <- rd2[[i]]
  posteriors[[i]] <- post_trends(name, window_sizes)
}; names(posteriors) <- rd2


ps <- lapply(posteriors, function(e) e[[3]]$p)

titlemap <- c("condamine_etal_2019_Picidae2" = "Picidae", 
              "condamine_etal_2019_Tyrannidae" = "Tyrannidae", 
              "condamine_etal_2019_Psittacidae" = "Psittacidae")

########
posterior_trends <- ps[[1]] + 
  labs(title = titlemap[rd2[[1]]],
       y = "Coverage") +
  theme(legend.position = "none") +
  ps[[2]] + labs(title = titlemap[rd2[[2]]],
                 y = "Coverage") + 
  theme(axis.title.y = element_blank(), 
        legend.position = c(0.5, 0.6),
        axis.text.y = element_blank(),
        legend.key.height = unit(3, "mm"),
        legend.key = element_rect(colour = "black", fill = NA, size = unit(0.5, "mm"))) +
  ps[[3]] + labs(title = titlemap[rd2[[3]]],
                 y = "Coverage") + 
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        legend.position = "none") +
  plot_layout(ncol = 3) &
  theme(plot.title = element_text(hjust = 0.5))

fname <- paste0("figures/ms/posterior-trends.pdf")
ggsave(fname, posterior_trends, width = 120, height = 80, units = "mm")

posteriors2 <- posteriors
for (i in 1:3){
  for (j in 1:4){
    if (j == 1){
      posteriors2[[i]][[j]]$p <- posteriors2[[i]][[j]]$p +
        ggtitle(titlemap[names(posteriors2)[i]]) +
        theme(plot.title = element_text(hjust = 0.5))
    }
    if (i == 1){
      posteriors2[[i]][[j]]$p <- posteriors2[[i]][[j]]$p +
        ylab(paste0("Coverage (k = ", window_sizes[j], ")"))
    }
    if (j < 4){
      posteriors2[[i]][[j]]$p <- posteriors2[[i]][[j]]$p +
        theme(axis.text.x = element_blank(),
              axis.title.x = element_blank())
    }
    if (i > 1){
      posteriors2[[i]][[j]]$p <- posteriors2[[i]][[j]]$p +
        theme(axis.text.y = element_blank(),
              axis.title.y = element_blank())
    }
    if (!(j == 1 && i == 1)){
      posteriors2[[i]][[j]]$p <- posteriors2[[i]][[j]]$p +
        theme(legend.position = "none")
    }else{
      posteriors2[[i]][[j]]$p <- posteriors2[[i]][[j]]$p +
        theme(legend.position = c(0.5, 0.5),
              legend.background = element_blank())
    }
  }
}

plots5 <- list(); k <- 1
for (i in 1:3){
  for (j in 1:4){
    plots5[[k]] <- posteriors2[[i]][[j]]$p
    k <- k + 1
  }
}

posterior_trend_windowsize <- Reduce("+", plots5) +
  plot_layout(ncol = 3, byrow = FALSE)
ggsave("figures/suppmat/posterior_trends_varyingwindowsize.pdf", 
       posterior_trend_windowsize,
       width = 150, height = 150, units = "mm")

#### Plot the posterior spaghetti plots for Picidae

ls <- list()
for (j in seq_along(posteriors)){
  ds_name <- names(posteriors)[j]
  l <- list()
  for (i in seq_along(posteriors[[j]][[3]]$samples)){
    df2 <- model2df(posteriors[[j]][[3]]$samples[[i]]$reference, 
                    compute.pulled.rates = FALSE)
    df2$name <- paste0("posterior",i)
    l[[i]] <- df2
  }
  df_posterior <- bind_rows(l)
  df_posterior$dataset <- titlemap[ds_name]
  ls[[j]] <- df_posterior
}

df_post <- bind_rows(ls)

posterior_spaghetti <- df_post %>%
  filter(rate %in% c("Speciation", "Extinction")) %>%
  filter(parse_number(name) < 20) %>%
  ggplot(aes(x = Time, y = value, color = rate, linetype = name)) +
  facet_grid(rows = vars(rate),
             cols = vars(dataset),
             scales = "free_x") +
  geom_line(alpha = 0.35) +
  scale_x_reverse() +
  theme_classic() +
  theme(legend.position = "none") +
  coord_cartesian(ylim = c(0.0, 0.6)) +
  scale_color_manual(values = c("orange", "black")) +
  scale_linetype_manual(values = rep(1, 110)) +
  xlab("Time (Ma)") +
  ylab("Rate")

ggsave("figures/suppmat/posterior-spaghetti.pdf", 
       posterior_spaghetti, 
       width = 100, height = 100,
       units = "mm") 
