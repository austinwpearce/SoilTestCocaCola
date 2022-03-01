# This testing file includes a soil test correlation dataset from the agridat
# package as an example for testing the lin_plateau, quad_plateau, mitscherlich,
# and ALCC functions.

library(tidyverse)
theme_set(theme_classic())

# relative cotton yield vs soil test potassium
agridat::cate.potassium
plot(agridat::cate.potassium$yield ~ agridat::cate.potassium$potassium)
#
agridat::engelstad.nitro # QP
plot(agridat::engelstad.nitro$yield ~ agridat::engelstad.nitro$nitro)
#
agridat::gartner.corn

agridat::gartner.corn %>% 
    ggplot(aes(long, lat, color = mass)) +
    geom_point(size = 7, alpha = 0.3) +
    scale_color_viridis_c(option = "mako")
#
agridat::lasrosas.corn

agridat::lasrosas.corn %>%
    ggplot(aes(long, lat, color = yield)) +
    geom_point(size = 3, alpha = 0.5) +
    facet_wrap(vars(year)) +
    scale_color_viridis_c(option = "mako")
#
agridat::sinclair.clover
agridat::sinclair.clover %>% 
    ggplot(aes(P, yield, color = yield)) +
    geom_point(size = 2, alpha = 1) +
    scale_color_viridis_c(option = "mako")


cotton <- tibble(stk = agridat::cate.potassium$potassium,
                 ry = agridat::cate.potassium$yield)