# This testing file includes a soil test correlation dataset from `agridat`
# for testing the lin_plateau, quad_plateau, mitscherlich, and alcc functions.

library(tidyverse)
# these packages are for completing code examples after the alcc stuff
library(agridat) # for obtaining a testing dataset
library(nlraa) # for self-starting functions and predicted intervals
library(minpack.lm) # for nlsLM, a robust backup to nls
library(nlstools) # for residuals plots
library(modelr) # for the r-squared and rmse
library(devtools)

# Load correlation functions
base_url <- "https://raw.githubusercontent.com/austinwpearce/SoilTestCocaCola/main/"
# linear plateau function
source_url(str_c(base_url, "lin_plateau.R"))

# quadratic plateau function
source_url(str_c(base_url, "quad_plateau.R"))

# mitscherlich function
source_url(str_c(base_url, "mitscherlich.R"))

# ALCC function
source_url(str_c(base_url, "alcc.R"))
source_url(str_c(base_url, "alcc_plot.R"))



# =============================================================================
# small n
crop <- tibble(x = seq(1, 4, 0.5),
               y = c(20, 70, 90, 95, 98, 92, 100))

lin_plateau(crop, plot = TRUE)
quad_plateau(crop, plot = TRUE)
mitscherlich(crop, plot = TRUE)
mitscherlich(crop, plot = TRUE, band = TRUE)
alcc_plot(crop, x, y, sufficiency = 3.5)

nls(y ~ mb(x, asym, b, c), data = crop, start = c(asym = 3.5, b = 7, c = -0.9))
nls(y ~ SSasymp(x, a, b, c), data = crop)

# relative cotton yield vs soil test potassium
head(agridat::cate.potassium)
crop <- agridat::cate.potassium %>% 
    rename(x = potassium,
           y = yield)

lin_plateau(crop, plot = TRUE)
quad_plateau(crop, plot = TRUE)
#AgroReg::quadratic.plateau(trat = crop$x, resp = crop$y)

# While the plotting function could possibly be combined with the previous
# function, keeping them separate is simpler

# =============================================================================
cotton <- tibble(stk = agridat::cate.potassium$potassium,
                 ry = agridat::cate.potassium$yield, 
                 dataset = "cotton")

cotton %>% 
    mutate(x = stk, y = ry) %>% 
    lin_plateau(plot = TRUE)

cotton %>% 
    mutate(x = stk, y = ry) %>% 
    quad_plateau(plot = TRUE)

cotton %>% 
    mutate(x = stk, y = ry) %>% 
    mitscherlich(plot = TRUE)

# ALCC

# Relative yield can be a ratio or percentage
# The ALCC method requires RY ratio values between 0-1,
# and percentage values between 0-100%

plot(cotton$ry ~ cotton$stk) |> abline(h = 100)
count(cotton, stk > 100) # 3 site-years exceeded 100

# Round down to 100 or let function do it
# cotton <- cotton %>% 
#     # cap RY at 100
#     mutate(ry = if_else(ry > 100, 100, ry))

# Create new dataset from correlation data
alcc(cotton, stk, ry)

alcc(cotton, stk, ry, sma = FALSE) %>% 
    ggplot(aes(xt_centered, yt)) +
    geom_smooth(method = "lm") +
    geom_point(size = 2, alpha = 0.5) +
    geom_vline(xintercept = 0) +
    geom_hline(aes(yintercept = intercept))

# Alternatively, if you have many groups/datasets in one table
multiple <- cotton %>% 
    group_by(dataset) %>%
    group_modify(~ alcc(data = .x, stk, ry)) %>% 
    ungroup()

multiple

##### PLOT #####
# for a single dataset

alcc_plot(cotton, stk, ry, sufficiency = 90)

# alternatively, continue using group_by + group_map framework for analyzing multiple datasets seamlessly
multiple %>% 
    group_by(dataset) %>%
    group_map(~ alcc_plot(data = .x, stk, ry))


##### CREATE simplified results table for export #####