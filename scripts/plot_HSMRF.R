library(ape)
library(dplyr)
library(ggthemes)

#setwd("~/projects/diversify")

my_quantiles <- function(x){
  res <- quantile(x, probs = c(0.025, 0.5, 0.975))
  names(res) <- c("lower", "median", "upper")
  return(res)
}

extract_rate <- function(item, df, times){
  df2 <- df %>% select(starts_with(paste0(item, ".")))

  df3 <- as.data.frame(t(apply(df2, 2, my_quantiles)))
  df3$time <- times
  rownames(df3) <- NULL
  df3$item <- item
  
  return(df3)
}

foo <- function(DS){
  # Get file paths for logs
  paths <- Sys.glob(paste0("output/HSMRF_", DS, "*.log"))
  
  if(length(paths) > 0){
    # Read the log files
    dfs <- lapply(paths, function(e) read.table(e, header = TRUE, sep = "\t"))

    
    # Concatenate the runs
    df <- do.call(rbind, dfs)
    
    # Read the tree and compute tree height, interval times
    phy <- read.tree(paste0("trees/", DS, ".tre"))
    th <- max(node.depth.edgelength(phy))
    n_intervals <- 100
    times <- seq(0, th, length.out = n_intervals)
    
    ## add netdiv
    sp <- df %>% select(starts_with(paste0("speciation_rate", ".")))
    ex <- df %>% select(starts_with(paste0("extinction_rate", ".")))
    netdiv <- sp - ex; colnames(netdiv) <- paste0("netdiv", ".", 1:100, ".")
    df <- bind_cols(df, netdiv)
    
    # Compute median and quantiles for the posterior
    items <- c("extinction_rate", "speciation_rate", "netdiv")
    
    df2 <- do.call(rbind, lapply(items, function(item) extract_rate(item, df, times)))
    
    return(df2)
  }else{
    return(NULL)
  }
}

DATASETS <- sapply(Sys.glob("trees/*.tre"), function(e) gsub(".tre", "", basename(e)), USE.NAMES = FALSE)
#DATASETS <- c(DATASETS, "condamine_etal_2019_Picidae2")

## Borrow plot code from RG [Cant install here since it requires R >= 4.0.0]
source("scripts/geomstepribbon.R")

# Loop over dataset names
for (DS in DATASETS){
  print(DS) 
  
  plotdata <- foo(DS)
  
  if(!is.null(plotdata)){
    p <- plotdata %>%
      ggplot2::ggplot(ggplot2::aes(time, median, color = item))  +
      ggplot2::theme_bw() +
      ggplot2::geom_step(ggplot2::aes(time, median),
                         direction = "vh") +
      geom_stepribbon(ggplot2::aes(x = time,
                                   ymin = lower,
                                   ymax = upper,
                                   fill = item),
                      direction = "vh",
                      alpha = 0.4,
                      color = NA) +
      ggplot2::scale_x_reverse() +
      ggplot2::xlab("Age (Ma)") +
      ggplot2::ylab("Rate") +
      ggplot2::theme(legend.title = ggplot2::element_blank(),
                     legend.position = "none",
                     panel.grid.major = ggplot2::element_blank(),
                     panel.grid.minor = ggplot2::element_blank(),
                     strip.background = ggplot2::element_blank()) +
      facet_wrap(dplyr::vars(item))
    
    ggsave(paste0("figures/rate_plots/HSMRF_", DS, "_rates.pdf"))
    write.table(plotdata, paste0("figures/plotdata/HSMRF_", DS, "_rates.csv"), sep = ";", row.names = FALSE)
  }
}
