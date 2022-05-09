# convergence_check_Picidae

library(tibble)
library(ape)
library(ggplot2)
library(dplyr)
library(patchwork)
library(tidyr)

make_quantiles <- function(dfs, tree){
  l <- list()
  for ( i in seq_along(dfs) ){
    df <- quantiles_inner(dfs[[i]], tree)
    df$replicate <- paste0("Run ", i)
    l[[i]] <- df
  }
  res <- bind_rows(l)
  return(res)
}

quantiles_inner <- function(samples, tree){

  
  extinction_prefix <- "extinction_rate."
  speciation_prefix <- "speciation_rate."
  
  speciation <- samples[, startsWith(names(samples), speciation_prefix)]
  extinction <- samples[, startsWith(names(samples), extinction_prefix)]
  nsamples <- nrow(speciation)
  
  q_sp <- as_tibble(t(apply(speciation, 2, function(episode) quantile(episode, probs = c(0.025, 0.5, 0.975))))); colnames(q_sp) <- c("lower", "median", "upper")
  q_sp$mean <- apply(speciation, 2, "mean")
  q_sp$time <- seq(0, max(node.depth.edgelength(tree)), length.out = ncol(speciation))
  q_sp$rate_name <- "Speciation rate"
  
  q_ex <- as_tibble(t(apply(extinction, 2, function(episode) quantile(episode, probs = c(0.025, 0.5, 0.975))))); colnames(q_ex) <- c("lower", "median", "upper")
  q_ex$mean <- apply(extinction, 2, "mean")
  q_ex$time <- seq(0, max(node.depth.edgelength(tree)), length.out = ncol(extinction))
  q_ex$rate_name <- "Extinction rate"
  
  q <- bind_rows(q_sp, q_ex)
  
  return(q)
}

make_facetplot <- function(q, limit = c(0, 1)){
  
  p <- q %>%
    ggplot(aes(time, median, fill = rate_name, color = rate_name)) +
    geom_line() +
    geom_ribbon(aes(ymin = lower, ymax = upper), color = NA, alpha = 0.1) +
    coord_cartesian(ylim = limit) +
    scale_x_reverse() +
    scale_fill_manual(values = c("orange", "black")) +
    scale_color_manual(values = c("orange", "black")) +
    facet_wrap(vars(replicate)) +
    theme_classic() +
    theme(legend.position = c(0.15, 0.9),
          legend.title = element_blank()) +
    labs(x = "time before present (Ma)",
         y = "rate")
  
  return(p)
}




make_ksdata <- function(fpaths){
  dfs <- lapply(fpaths, function(path) tibble(read.table(path, header = TRUE, sep = "\t")))
  n_samples <- sapply(dfs, nrow)
  
  rate_names <- c("extinction_rate.", "speciation_rate.")
  idxs <- list(
    c(1,2),
    c(1,3),
    c(1,4),
    c(2,3),
    c(2,4),
    c(3,4)
  )
  prettylabel <- c("speciation_rate." = "Speciation rate",
                   "extinction_rate." = "Extinction rate")
  res <- list(); j = 1
  for (idx in idxs){
    idx1 <- idx[1]
    idx2 <- idx[2]
    
    ks_rates <- list()
    for (rate_name in rate_names){
      ks_episodes <- list()
      
      run1 <- dfs[[idx1]] %>%
        dplyr::select(starts_with(rate_name))
      run2 <- dfs[[idx2]] %>%
        dplyr::select(starts_with(rate_name))
      
      for (i in seq_along(colnames(run1))){
        cname <- colnames(run1)[i]
        ks_res <- ks.test(run1[[cname]], run2[[cname]])
        ks_df <- tibble("D" = ks_res$statistic,
                        "p" = ks_res$p.value,
                        "comparison" = paste0("Run ", idx1, " vs Run ", idx2),
                        "episode" = i,
                        "rate_name" = prettylabel[rate_name])
        ks_episodes[[i]] <- ks_df
      }
      ks_rates[[rate_name]] <- bind_rows(ks_episodes)
    }
    res[[j]] <- bind_rows(ks_rates)
    
    j <- j  + 1 
  }
  ks <- bind_rows(res)
  return(ks)
}

make_ksplot <- function(ks){
  ksplot <- ks %>% ggplot(aes(x = episode, y = D, color = rate_name)) +
    geom_line() +
    theme_classic() +
    ylim(c(0.0, 0.09)) +
    scale_x_reverse() +
    geom_hline(aes(yintercept = 0.09), linetype = 2, color = "red") +
    scale_color_manual(values = c("orange", "black")) +
    ylab("Kolmogorov-Smirnov test (D statistic)") +
    theme(legend.position = c(0.50, 0.9),
          legend.title = element_blank(),
          legend.background = element_blank()) +
    facet_wrap(vars(comparison))
  return(ksplot)
}

###########
treepaths <- c("~/projects/diversify/trees/condamine_etal_2019_Picidae.tre",
               "~/projects/diversify/trees/condamine_etal_2019_Tyrannidae.tre",
               "~/projects/diversify/trees/condamine_etal_2019_Psittacidae.tre")

trees <- lapply(treepaths, read.tree)

logpaths <- c("~/projects/diversify/output/HSMRF_condamine_etal_2019_Picidae2_run_*.log",
              "~/projects/diversify/output/HSMRF_condamine_etal_2019_Tyrannidae_run_*.log",
              "~/projects/diversify/output/HSMRF_condamine_etal_2019_Psittacidae_run_*.log")
lfpaths <- lapply(logpaths, Sys.glob)

outnames <- c("Picidae", "Tyrannidae", "Psittacidae")

## Picidae
ps <- list()
for (i in seq_along(outnames)){
  print(".")
  tree <- trees[[i]]
  fpaths <- lfpaths[[i]]
  ksdata <- make_ksdata(fpaths)
  ksplot <- make_ksplot(ksdata)
  dfs <- lapply(fpaths, function(path) tibble(read.table(path, header = TRUE, sep = "\t")))
  df_quantiles <- make_quantiles(dfs, tree)
  p_rates <- make_facetplot(df_quantiles, limit = c(0, 0.5))
  
  p <- ksplot | p_rates
  ps[[outnames[i]]] <- p
  
  fname <- paste0("figures/suppmat/convergence/", outnames[i], "_ksplot.pdf")
  ggsave(fname, p, height = 140, width = 250, units = "mm")
}

###################
