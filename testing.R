# This testing file includes a soil test correlation dataset from the agridat
# package as an example for testing the lin_plateau, quad_plateau, mitscherlich,
# and ALCC functions.

library(tidyverse)
library(devtools)
theme_set(theme_classic())

# Load correlation functions
# linear plateau function
source_url("https://raw.githubusercontent.com/austinwpearce/SoilTestCocaCola/main/lin_plateau.R")

# quadratic plateau function
source_url("https://raw.githubusercontent.com/austinwpearce/SoilTestCocaCola/main/quad_plateau.R")

# ALCC function
source_url("https://raw.githubusercontent.com/austinwpearce/SoilTestCocaCola/main/alcc.R")


# =============================================================================
# small n
crop <- tibble(x = c(1, 2, 3, 4, 5),
               y = c(1, 2, 4, 4, 3))

lin_plateau(crop, plot = TRUE)
quad_plateau(crop, plot = TRUE)
mitscherlich(crop, plot = TRUE)

nls(y ~ mb(x, asym, b, c), data = crop, start = c(asym = 3.5, b = 7, c = -0.9))
nls(y ~ SSasymp(x, a, b, c), data = crop)

# relative cotton yield vs soil test potassium
head(agridat::cate.potassium)
crop <- agridat::cate.potassium %>% 
    rename(x = potassium,
           y = yield)

lin_plateau(crop, plot = TRUE)
quad_plateau(crop, plot = TRUE)
AgroReg::linear.plateau(trat = crop$x, resp = crop$y)

# While the plotting function could possibly be combined with the previous
# function, keeping them separate is simpler

# ALCC ========================================================================

##### DATA IMPORT #####
cotton <- tibble(stk = agridat::cate.potassium$potassium,
                 ry = agridat::cate.potassium$yield, 
                 dataset = "cotton")

# Relative yield can be a ratio or percentage
# The ALCC method requires RY ratio values between 0-1,
# and percentage values between 0-100%

plot(cotton$ry ~ cotton$stk) |> abline(h = 100)
count(cotton, stk > 100) # 3 site-years exceeded 100

cotton <- cotton %>% 
    # cap RY at 100
    mutate(ry = if_else(ry > 100, 100, ry))

plot(cotton$ry ~ cotton$stk) |> abline(h = 100)

# Create new dataset from correlation data
alcc_results <- alcc(cotton, stk, ry)

alcc_results

alcc_results %>% 
    ggplot(aes(xt_centered, yt)) +
    geom_point(size = 2, alpha = 0.5) +
    geom_smooth(method = "lm") +
    geom_vline(xintercept = 0) +
    geom_hline(yintercept = alcc_results$intercept)

# Alternatively, if you have many groups/datasets in one table
alcc_results <- cotton %>% 
    group_by(dataset) %>%
    group_modify(~ alcc(data = .x, stk, ry))

##### PLOT #####
# for a single dataset
alcc_plot(cotton, stk, ry)

# alternatively, continue using group_by + group_map framework for analyzing multiple datasets seamlessly
alcc_results %>% 
    group_by(dataset) %>%
    group_map(~ alcc_plot(data = .x, sufficiency = 95))


##### CREATE simplified results table for export #####