models2df <- function(models, ...){
  l <- list()
  
  for (i in seq_along(models)){
    df <- model2df(models[[i]], ...)
    df$name <- names(models)[[i]]
    l[[i]] <- df
  }
  res <- bind_rows(l)
  
  return(res)
}

spaghetti_plots <- function(cg, columns = c("Speciation", "Extinction")){
  cg_df <- models2df(cg) %>%
    filter(rate %in% columns) %>%
    mutate(rate = factor(rate, levels = columns))
  
  df4 <- cg_df %>% 
    filter(startsWith(name, "reference"))
  
  df5 <- cg_df %>% 
    filter(!startsWith(name, "reference"))
  
  df5$category <- sapply(strsplit(df5$name, "\\."), function(e) e[[1]])
  df5$category <- factor(df5$category, levels = names(ref_foos))
  df5$index <- stringr::str_match_all(df5$name, "[0-9]+") %>% 
    unlist() %>% 
    as.numeric() %>%
    `+`(ifelse(df5$rate == "Speciation", 5, 0)) %>%
    as.factor()
  
  
  cols <- c(colorspace::sequential_hcl(8, "orange")[1:5],
            colorspace::sequential_hcl(8, "blues")[1:5])
  p <- df5 %>%
    ggplot(aes(x = Time, y = value, color = index)) +
    geom_line() +
    geom_line(data=df4, colour="black") +
    facet_grid(category ~ rate, labeller = labeller(category = lf), scales = "free") + 
    theme_classic() +
    scale_color_manual(values = cols) +
    scale_x_reverse() +
    labs(y = "Rate (events per lineage per time)", x = "Time (Ma)") +
    theme(legend.position = "none")
  return(p)  
}


f <- function(m, name){
  times <- m$times
  df <- tibble("times" = times,
               "Speciation rate" = m$lambda(times),
               "Extinction rate" = m$mu(times),
               "Model" = name)
  return(df)
}


read_dataset <- function(dataset_name){
  df <- tibble(read.csv(paste0("figures/plotdata/HSMRF_", dataset_name, "_rates.csv"), sep = ";"))
  
  phy <- read.tree(paste0("trees/", dataset_name, ".tre"))
  height <- max(node.depth.edgelength(phy))
  
  speciation <- dplyr::filter(df, item == "speciation_rate")
  extinction <- dplyr::filter(df, item == "extinction_rate")
  netdiv <- dplyr::filter(df, item == "netdiv")
  
  lambda <- approxfun(speciation$time, speciation$median)
  mu <- approxfun(extinction$time, extinction$median)
  times <- seq(from = min(speciation$time), to = max(speciation$time), length.out = 500)
  
  reference <- create.model(lambda, mu, times = times)
  
  res <- list(reference = reference, 
              speciation = speciation,
              extinction = extinction,
              netdiv = netdiv,
              height = height, 
              phy = phy)
  return(res)
}

logistic_shift <- function(min, max, steepness, x0){
  res <- function(x) {
    ((max-min) / (1 + 1.0*steepness^(-1*(x - x0)))) + min
  }
  return(res)
}

alternative_lambdas <- function(dataset){
  reference <- dataset$reference
  speciation <- dataset$speciation
  extinction <- dataset$extinction
  height <- dataset$height
  
  
  models <- list()
  
  ## SET UP MODELS
  
  lambda0 <- reference$lambda(0.0)
  
  models[["constant"]] <- congruent.models(reference, lambdas = function(t) lambda0)
  
  ## LINEARLY INCREASING
  fold_changes <- c(1.4, 1.8, 2.2, 2.6)
  bs <- sapply(fold_changes, function(fold_change) lambda0 *(1 - fold_change)/(fold_change*height))
  linear_lambda_inc <- sapply(bs, function(b) function(t) lambda0 + b * t)
  
  models[["linear"]] <- congruent.models(reference, lambdas = linear_lambda_inc)
  
  ## LINEAR DECREASING
  linear_lambda_dec <- sapply(bs, function(b) function(t) lambda0 - b * t)
  models[["linear_dec"]] <- congruent.models(reference, lambdas = linear_lambda_dec)
  
  ## EXPONENTIAL
  bs <- sapply(fold_changes, function(fold_change) -log(1/fold_change)/height)
  exp_lambdas <- sapply(bs, function(b) function(t) lambda0 * exp(-b*t))
  models[["exp1"]] <- congruent.models(reference, lambdas = exp_lambdas)
  
  ## two-epoch
  twoepoch_lambdas_inc <- sapply(fold_changes, function(fold_change) function(t) ifelse(t <= 0.2*height, lambda0, lambda0/fold_change))
  models[["two_inc"]] <- congruent.models(reference, lambdas = twoepoch_lambdas_inc)
  
  ##two-epoch decrease
  twoepoch_lambdas_dec <- sapply(fold_changes, function(fold_change) function(t) ifelse(t <= 0.2*height, lambda0, lambda0*fold_change))
  models[["two_dec"]] <- congruent.models(reference, lambdas = twoepoch_lambdas_dec)
  
  ##three-epoch
  threeepoch_lambdas <- sapply(fold_changes, function(fold_change) function(t) ifelse(t <= 0.45*height, lambda0, ifelse(t <= 0.55*height, lambda0*fold_change, lambda0)))
  models[["three"]] <- congruent.models(reference, lambdas = threeepoch_lambdas)
  
  class(models) <- c("list", "ACDCsets")
  
  return(models)
}


alternative_mus <- function(dataset, addconstant){
  reference <- dataset$reference
  speciation <- dataset$speciation
  extinction <- dataset$extinction
  height <- dataset$height
  
  models = list()
  class(models) <- c("list", "ACDCsets")
  
  ## CONSTANT MODELS
  mus = c(0.1, 0.2, 0.3, 0.4)
  fmus <- sapply(mus, function(mu) function(t) mu + addconstant)
  models[["constants"]] <- congruent.models(reference, mus = fmus)
  
  ## EXPONENTIAL
  bs <- c(0.05, 0.07, 0.1)
  fmus <- sapply(bs, function(b) function(t) 0.3*exp(-b*t) + addconstant)
  models[["exp1"]] <- congruent.models(reference, mus = fmus)
  
  ## Exponential decreasing
  fmus <- sapply(bs, function(b) function(t) 0.3*exp(-b*(height-t)) + addconstant)
  models[["exp2"]] <- congruent.models(reference, mus = fmus)
  
  ## bell curve middle
  bs <- c(0.1, 0.2, 0.3, 0.4)
  gaus <- sapply(bs, function(b) function(t) b * exp((-(t - height*0.5)^2) / 6) + 0.1 + addconstant)
  models[["modal1"]] <- congruent.models(reference, mus = gaus)
  
  ## shifted modal
  gaus <- sapply(bs, function(b) function(t) b * exp((-(t - height*0.2)^2) / 6) + 0.1 + addconstant)
  models[["modal2"]] <- congruent.models(reference, mus = gaus)
  
  ## shifted modal at lower mu0
  gaus <- sapply(bs, function(b) function(t) b * exp((-(t - height*0.2)^2) / 6) + 0.005 + addconstant)
  models[["modal3"]] <- congruent.models(reference, mus = gaus)
  
  
  
  class(models) <- c("list", "ACDCsets")
  
  return(models)
}

rateplots <- function(m2, 
                      ylabs = c("Increasing", "Decreasing", "Modal"),
                      rate_names = c("Extinction rate", "Speciation rate"),
                      set_names = c("exp1", "exp2", "modal3"),
                      ylim = c(0, 0.6)
){
  
  res <- list()
  
  
  ## Iter over dec/inc/modal
  for (i in seq_along(set_names)){
    ps <- list()
    plotdata <- do.call(bind_rows, mapply(f, m2[[set_names[i]]], names(m2[[set_names[i]]]), SIMPLIFY = FALSE)) %>%
      pivot_longer(`Speciation rate`:`Extinction rate`, names_to = "rate", values_to = "value")
    
    cols <- list("Speciation rate" = c(head(colorspace::sequential_hcl(palette = "Blues", n = length(unique(plotdata$Model))+2), n = -3), "black"),
                 "Extinction rate" = c(head(colorspace::sequential_hcl(palette = "Oranges", n = length(unique(plotdata$Model))+2), n = -3), "black"))
    
    ## Iter over lambda, mu
    for (j in 1:2){
      rate_name <- rate_names[j]
      df <- plotdata %>% 
        dplyr::filter(rate == rate_name)
      
      p <- df %>% 
        ggplot(aes(x = times, y = value, color = Model)) +
        geom_line() +
        scale_x_reverse() +
        theme_classic() +
        ylab(ylabs[i]) +
        scale_color_manual(values = cols[[rate_name]]) +
        theme(
          strip.background = element_blank(),
          strip.text.x = element_blank(),
          legend.position =  "none",
          axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          axis.line.y = element_blank(),
          axis.ticks.y = element_blank()
        )
      
      if(!is.na(head(ylim, n = 1))){
        p <- p + coord_cartesian(ylim = ylim)
      }
      
      
      if (j == 1){
        p <- p + 
          theme(
            axis.title.y = element_text(angle = 90),
            axis.text.y = element_text(),
            axis.line.y = element_line(colour = "black"),
            axis.ticks.y = element_line(colour = "black")
          )
      }
      
      if (i == length(set_names)){ ## check if last row
        p <- p + 
          theme(
            axis.title.x = element_text(),
            axis.text.x = element_text()
          ) +
          xlab("Time (Ma)")
      }
      
      ps[[rate_name]] <- p
    }
    res[[set_names[i]]] <- ps
  }
  return(res)
}