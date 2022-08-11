## Suppose we set up a reference model, with
## 
## \lambda(t) = sigmoidally increasing
## \mu(t) = 0.28

source("scripts/hypothetical_models_frommu.R")
plot(cgs_frommu$up$reference)

## And conversely,
## \lambda(t) = 0.28
## \mu(t) = sigmoidally increasing

source("scripts/hypothetical_models_fromlambda.R")
plot(cgs_fromlambda$up$reference)

######

up_models_frommu <- models2df(cgs_frommu$up)
up_models_fromlambda <- models2df(cgs_fromlambda$up)

var_sp <- up_models_frommu %>%
  group_by(Time) %>%
  summarize("v" = var(value))
var_sp$rate <- "Speciation"

var_ex <- up_models_fromlambda %>%
  group_by(Time) %>%
  summarize("v" = var(value))
var_ex$rate <- "Extinction"
var_df <- bind_rows(var_sp, var_ex)

p1 <- ggplot(var_df, aes(x = Time, y = v, color = rate)) +
  geom_line() +
  theme_classic() +
  scale_x_reverse() +
  scale_y_log10() +
  ylab("Var") +
  theme(legend.position = c(0.2, 0.6),
        legend.background = element_blank(),
        legend.title = element_blank()) +
  scale_color_manual(values = c("orange", "black"))

p1

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



