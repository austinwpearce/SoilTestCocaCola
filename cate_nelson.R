## Cate-Nelson ===========================================
library(rcompanion)

red <- "#CE1141"
gold <- "#EAAA00"
blue <- "#13274F"
black <- "#000000"

cate_nelson <-
    function(data,
             trend = "positive",
             verbose = FALSE,
             plot = FALSE) {
        par(mfrow = c(2, 3))
        cn <- cateNelson(
            x          = data$x,
            y          = data$y,
            plotit     = FALSE,
            hollow     = TRUE,
            verbose    = verbose,
            progress = FALSE,
            xlab       = "stv",
            ylab       = "ry",
            trend      = trend,
            clx        = 1,
            cly        = 1,
            xthreshold = 0.10,
            ythreshold = 0.15
        )
        minx <- min(data$x)
        maxx <- max(data$x)
        rangex <- maxx - minx
        maxy <- max(data$y)
        cstv <- cn$CLx
        cry <- cn$CLy
        
        if (plot == FALSE) {
            # from Steve Culman
            dat <- data.frame(x = data$x,
                              y = data$y)
            
            dat <- dat[order(dat$x), ] # Sort the data by x
            x <- dat$x
            y <- dat$y
            
            # Create a data.frame to store the results
            
            out <- data.frame(
                x = NA,
                mean1 = NA,
                css1 = NA,
                mean2 = NA,
                css2 = NA,
                r2 = NA
            )
            
            css <- function(x) {
                var(x) * (length(x) - 1)
            }
            tcss <- css(y) # Total corrected sum of squares
            
            for (i in 2:(length(y) - 2)) {
                y1 <- y[1:i]
                y2 <- y[-(1:i)]
                out[i, 'x'] <- x[i]
                out[i, 'mean1'] <- mean(y1)
                out[i, 'mean2'] <- mean(y2)
                out[i, 'css1'] <- css1 <- css(y1)
                out[i, 'css2'] <- css2 <- css(y2)
                out[i, 'r2'] <- (tcss - (css1 + css2)) / tcss
                
            }
            
            cn_rsq <- as_tibble(out) %>%
                arrange(desc(r2)) %>%
                slice_max(r2)
            
            return(tibble(cstv = round(cstv, 1),
                          ry = round(cry, 1)) %>%
                bind_cols(cn_rsq))
        } else {
            
            data %>%
                ggplot(aes(x, y)) +
                geom_vline(xintercept = cstv,
                           alpha = 1,
                           color = blue) +
                geom_hline(yintercept = cry,
                           alpha = 1,
                           color = blue) +
                geom_point(size = 3, alpha = 0.5) +
                geom_rug(alpha = 0.2, length = unit(2, "pt")) +
                annotate(
                    "text",
                    label = paste("CSTV =", round(cstv, 1), "ppm"),
                    x = cstv,
                    y = 0,
                    angle = 90,
                    hjust = 0,
                    vjust = 1.5,
                    alpha = 0.5
                ) +
                annotate(
                    "text",
                    label = paste0("RY = ", round(cry, 1), "%"),
                    x = max(data$x),
                    y = cry,
                    angle = 0,
                    hjust = 1,
                    vjust = 1.5,
                    alpha = 0.5
                ) +
                scale_x_continuous(
                    breaks = seq(0, maxx * 2, by = if_else(
                        condition = rangex >= 300,
                        true = 30,
                        false = if_else(
                            condition = rangex >= 100,
                            true = 20,
                            false = if_else(
                                condition = rangex >= 50,
                                true = 5,
                                false = 2))))) +
                scale_y_continuous(limits = c(0, maxy),
                                   breaks = seq(0, maxy * 2, 10)) +
                labs(
                    x = "Soil test value (mg/kg)",
                    y = "Relative yield (%)",
                    caption = paste("Each point is a site. n =", nrow(data))
                ) +
                theme(legend.position =  "none")
        }
    }
