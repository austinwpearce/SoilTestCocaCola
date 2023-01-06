## Cate-Nelson ===========================================
library(rcompanion)

red <- "#CE1141"
gold <- "#EAAA00"
blue <- "#13274F"
black <- "#000000"



cate_nelson <-
    function(data,
             x,
             y,
             trend = "positive",
             verbose = FALSE,
             plot = FALSE,
             sufficiency = NULL) {
        
        if (missing(x)) {
            stop("Please specify the explanatory variable name (e.g. soil test concentration) using the `x` argument")
        }
        
        if (missing(y)) {
            stop("Please specify the response variable name (e.g. relative yield) using the `y` argmuent")
        }
        
        # Re-define x and y from STV and RY (tip to AC)
        x <- rlang::eval_tidy(data = data, rlang::quo({{x}}) )
        
        y <- rlang::eval_tidy(data = data, rlang::quo({{y}}) )
        
        if (max(y) < 2) {
            stop("The reponse variable appears to not be on a percentage scale.
             If so, please multiply it by 100.")
        }
        
        corr_data <- dplyr::tibble(x = as.numeric(x), 
                                   y = as.numeric(y))
        
        if (nrow(corr_data) < 4) {
            stop("Too few distinct input values to fit LP. Try at least 4.")
        }
        
        par(mfrow = c(2, 3))
        
        minx <- min(corr_data$x)
        maxx <- max(corr_data$x)
        rangex <- maxx - minx
        maxy <- max(corr_data$y)

        if (is.null(sufficiency)){
        cn_free <- rcompanion::cateNelson(
            x          = corr_data$x,
            y          = corr_data$y,
            plotit     = FALSE,
            verbose    = verbose,
            progress = FALSE,
            xlab       = "STV",
            ylab       = "RY",
            trend      = trend,
            clx        = 1,
            cly        = 1,
            xthreshold = 0.10,
            ythreshold = 0.15
        )
        
        cstv <- cn_free$CLx
        cry <- cn_free$CLy
        }
        
        if(!is.null(sufficiency)) {
            cn_fixed <- rcompanion::cateNelsonFixedY(
                x = corr_data$x,
                y = corr_data$y,
                cly = sufficiency,
                trend = trend,
                plotit = FALSE,
                xlab = "STV",
                ylab = "RY")
        
            cstv <- cn_fixed$Critx[[1]]
            cry <- cn_fixed$Crity[[1]]
            
        }
        
            # from Steve Culman
            dat <- data.frame(x = corr_data$x,
                              y = corr_data$y)
            
            dat <- corr_data
            
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
            
            cn_table <- tibble(cstv = round(cstv, 1),
                               cry = round(cry, 1)) %>%
                bind_cols(cn_rsq)
            
            if (plot == FALSE){
            
            return(cn_table)
                
            } else {
            
            cn_plot <- corr_data %>%
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
                    x = maxx,
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
                    caption = paste("Each point is a site. n =", nrow(corr_data))
                ) +
                theme(legend.position =  "none")
            
            return(cn_plot)
        }
    }


## Cate-Nelson

Follows the function from the `rcompanion` package. The corresponding function in the `soiltestcorr` package is `cate_nelson_1965()` but this function makes an error and needs fixing.

Warning: DO NOT use the highest R^2^ value obtained from the intermediate Cate-Nelson table, but rather the R^2^ value associated with the CSTV at 95% RY. The absolute highest R^2^ value may be associated instead with another RY level not 95%. I just caught that error in my previous work.

```{r}
get_cstv_cn <- function(data, target){
    cstv <- cateNelsonFixedY(
        x = data$x, y = data$y, 
        cly = target,
        plotit = FALSE,
        trend = "positive")[1, 1]
    return(cstv)
}

get_r2_cn <- function(data, cstv){
    x <- data$x
    y <- data$y
    
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
    
    r2 <- filter(out, x == cstv) %>% 
        slice_max(r2) %>% 
        pull(r2) %>% 
        round(2)
    
    return(r2)
    
}
```

```{r}
cn_models <- corr_nested %>% 
    mutate(method = "CN",
           CSTV = map2(data, 95, get_cstv_cn) %>% as.numeric(),
           AIC = NA,
           r2 = map2(data, CSTV, get_r2_cn) %>% as.numeric(),
           suff = 95)

cn_models
```
