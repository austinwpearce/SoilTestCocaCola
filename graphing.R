# My own ggplot2 theme
sanserif <- c("Roboto Condensed")

theme_gc <- 
    theme_minimal(#base_family = sanserif,
        base_size = 14) +
    ggplot2::theme(
        plot.background = NULL,
        plot.margin = margin(t = 2, r = 2, b = 2, l = 2, unit = "pt"),
        panel.grid = element_line(color = "#F4F4F4"),
        panel.spacing = unit(2, "lines"),
        panel.border = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.title.y = element_text(
            hjust = 1,
            margin = margin(t = 0, r = 10, b = 0, l = 0, unit = "pt")),
        axis.title.x = element_text(
            hjust = 0,
            margin = margin(t = 10, r = 0, b = 0, l = 0, unit = "pt")),
        #axis.text = element_text(family = sanserif),
        legend.title.align = 0,
        legend.key.height = unit(x = 5, units = "mm"),
        legend.justification = c(1, 1))

ggplot2::theme_set(theme_gc)
