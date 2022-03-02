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

