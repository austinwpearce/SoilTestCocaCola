library(tidyverse)
library(nlraa)
library(minpack.lm)

corr <- readxl::read_xlsx("datasets.xlsx") |> 
    #filter(dataset == "Olsen P") |> 
    mutate(sig_response = as_factor(sig_response_),
           responsive = if_else(sig_response == "Responsive", 1, 0))


theme_set(theme_bw())

corr |> 
    ggplot(aes(x = stv, y = ry_fitmax, color = sig_response)) +
    geom_point(size = 2) +
    facet_wrap(vars(dataset), scales = "free_x")

# Typical correlation
corr |> 
    ggplot(aes(x = stv, y = ry_max)) +
    geom_point(size = 2, aes(color = sig_response)) +
    geom_line(stat = "smooth",
              method = "nls",
              formula = y ~ SSquadp3xs(x, a, b, j),
              se = FALSE) +
    facet_wrap(vars(dataset), scales = "free_x")

# Alternative correlation
corr |> 
    ggplot(aes(x = stv, y = responsive)) + 
    geom_point(size = 2, aes(color = sig_response)) + 
    stat_smooth(method = "glm",
                method.args = list(family = "binomial"),
                se=FALSE) +
    # geom_vline(
    #     #xintercept = c(38, 47, 56)
    #     xintercept = c(83, 107, 132)
    #     #xintercept = c(3, 4.8, 12.7)
    #     ) +
    facet_wrap(vars(dataset), scales = "free_x")

# model?
model <- glm(sig_response ~ stv, data = corr, family = "binomial")

model
summary(model)
confint(model)

exp(coef(model))
