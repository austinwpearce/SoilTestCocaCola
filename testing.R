# This testing file includes a soil test correlation dataset from `agridat`
# for testing the lin_plateau, quad_plateau, mitscherlich, and alcc functions.

library(ggplot2)
library(devtools)
library(stringr)
library(dplyr)
# these packages are for completing code examples after the alcc stuff
library(agridat) # for obtaining a testing dataset
library(nlraa) # for self-starting functions and predicted intervals
library(minpack.lm) # for nlsLM, a robust backup to nls
library(nlstools) # for residuals plots
library(modelr) # for the r-squared and rmse
library(soiltestcorr)

# Load correlation functions
base_url <- "https://raw.githubusercontent.com/austinwpearce/SoilTestCocaCola/main/"
# linear plateau function
source_url(str_c(base_url, "lin_plateau.R"))

# quadratic plateau function
source_url(str_c(base_url, "quad_plateau.R"))

# mitscherlich function
source_url(str_c(base_url, "mitscherlich.R"))

# mitscherlich_1000 (forces asymptote to 100 and intercept to 0)
source_url(str_c(base_url, "mitscherlich_1000.R"))

# ALCC function
source_url(str_c(base_url, "alcc.R"))
source_url(str_c(base_url, "alcc_plot.R"))

# =============================================================================
corr_data <- tibble(
    stk = agridat::cate.potassium$potassium,
    ry = agridat::cate.potassium$yield, 
    dataset = "cotton")

corr_data <- tibble(
    x = agridat::cate.potassium$potassium,
    y = agridat::cate.potassium$yield, 
    dataset = "cotton")


cate_nelson(corr_data, stk, ry)


# Note that the X variable is an integer! This is handled within the function

lin_plateau(corr_data)
# must specify the ST and RY columns
lin_plateau(corr_data, stk, ry)
lin_plateau(stv =  corr_data$stk, ry = corr_data$ry)
lin_plateau(corr_data, stk, ry, plot = TRUE)
lin_plateau(corr_data, stk, ry, plot = TRUE, extrapolate = TRUE)

quad_plateau(corr_data)
# must specify the ST and RY columns
quad_plateau(corr_data, stk, ry)
quad_plateau(stv =  corr_data$stk, ry = corr_data$ry)
quad_plateau(corr_data, stk, ry, plot = TRUE)
quad_plateau(corr_data, stk, ry, plot = TRUE, band = TRUE)

mitscherlich(corr_data)
# must specify the ST and RY columns
mitscherlich(corr_data, stk, ry)
mitscherlich(stv = corr_data$stk, ry = corr_data$ry)
mitscherlich(corr_data, x, y, plot = TRUE)

# ALCC

# Relative yield can be a ratio or percentage
# The ALCC method requires RY ratio values between 0-1,
# and percentage values between 0-100%

plot(corr_data$ry ~ corr_data$stk) |> abline(h = 100)
count(corr_data, stk > 100) # 3 site-years exceeded 100

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

#### All in one ####
all_in_one(cotton, stk, ry)


# =============================================================================
# small n
crop <- tibble(x = seq(0, 50, 5),
               y = quadp3xs(x, 0, 10, 10) + rnorm(11, 0, 4))

crop <- tibble(x = seq(5, 50, 10),
               y = c(20, 70, 90, 95, 98, 92, 100))

lin_plateau(crop, plot = TRUE)
quad_plateau(crop, plot = TRUE)
mitscherlich(crop, x, y, plot = TRUE)
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