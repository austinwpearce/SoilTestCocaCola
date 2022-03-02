#
head(agridat::engelstad.nitro) # QP
crop <- agridat::engelstad.nitro %>%
    rename(x = nitro, y = yield) %>% 
    group_by(x) %>% 
    summarise(y = mean(y)) %>% 
    ungroup()

crop %>% 
    ggplot(aes(x, y, color = y)) +
    geom_point(size = 2, alpha = 1) +
    scale_color_viridis_c(option = "mako")

lin_plateau(crop, plot = TRUE)

#
head(agridat::sinclair.clover)
crop <- tibble(x = agridat::sinclair.clover$P,
               y = agridat::sinclair.clover$yield) %>% 
    group_by(x) %>% 
    summarise(y = mean(y))

crop %>% 
    ggplot(aes(x, y, color = y)) +
    geom_point(size = 2, alpha = 1) +
    scale_color_viridis_c(option = "mako")

lin_plateau(crop, plot = TRUE)

#
head(agridat::lasrosas.corn)
agridat::lasrosas.corn %>%
    ggplot(aes(long, lat, color = yield)) +
    geom_point(size = 3, alpha = 0.5) +
    facet_wrap(vars(year)) +
    scale_color_viridis_c(option = "mako")
crop <- tibble(agridat::lasrosas.corn) %>% 
    filter(year == 1999) %>% 
    group_by(year, rep, nitro) %>% 
    summarise(y = mean(yield)) %>% 
    rename(x = nitro) %>% 
    ungroup()
crop %>%
    ggplot(aes(x, y, color = y)) +
    geom_point(size = 3, alpha = 0.5) +
    scale_color_viridis_c(option = "mako")

lin_plateau(crop, plot = TRUE)

#
head(agridat::gartner.corn)

agridat::gartner.corn %>% 
    ggplot(aes(x = long, y = lat, color = mass)) +
    geom_point(size = 7, alpha = 0.3) +
    scale_color_viridis_c(option = "mako")
