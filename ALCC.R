# Soil-test correlation using Arcsine-Log Calibration Curve (ALCC)
# Author: Austin Pearce
# ================================================================

##### All packages needed #####
library(tidyverse) # For data import, wrangling, and plotting

##### Visualization settings #####
red <- "#CE1141"
blue <- "#13274F"

# optional plot theme setting (for me mostly)
theme_set(
    theme_minimal(base_family = "serif", base_size = 14) +
        theme(
            plot.background = NULL,
            plot.margin = margin(
                t = 2,
                r = 2,
                b = 2,
                l = 2,
                unit = "pt"
            ),
            panel.grid = element_line(color = "#F4F4F4"),
            panel.spacing = unit(2, "lines"),
            panel.border = element_blank(),
            axis.line = element_blank(),
            axis.ticks = element_blank(),
            axis.title.y = element_text(
                hjust = 1,
                margin = margin(
                    t = 0,
                    r = 10,
                    b = 0,
                    l = 0,
                    unit = "pt"
                )
            ),
            axis.title.x = element_text(
                hjust = 0,
                margin = margin(
                    t = 10,
                    r = 0,
                    b = 0,
                    l = 0,
                    unit = "pt"
                )
            ),
            axis.text = element_text(family = "serif"),
            legend.title.align = 0,
            legend.key.height = unit(x = 5, units = "mm"),
            legend.justification = c(1, 1)
            #legend.position = c(1, 1)
        )
)


##### Two fancy functions #####
# Function alcc_sma() creates new variables on existing dataset
alcc_sma <- function(data, sufficiency = 95){
    steps_1_4 <- data %>% 
        mutate(
            # get sample size
            n = n(),
            # Step 1 Transform (remember, soil test goes to Y and RY goes to X!)
            y = log(stv),
            x = asin(sqrt(ry / 100)),
            # Step 2 Center
            sufficiency = sufficiency,
            x_adjust = asin(sqrt(sufficiency / 100)),
            x_center = x - x_adjust,
            # Step 3 Correlation
            pearson = cor(x_center, y, method = "pearson"),
            t_stat  = (pearson * sqrt(n - 2)) / sqrt(1 - pearson ^ 2),
            pvalue  = pt(t_stat, df = n - 1, lower.tail = FALSE),
            # Step 4 Means
            mean_x  = mean(x_center),
            mean_y  = mean(y)
        )
    
    # Step 5
    ols_center <- lm(y ~ x_center, data = steps_1_4)
    
    steps_5_9 <- steps_1_4 %>% 
        mutate(intercept = coef(ols_center)[[1]],
               slope = coef(ols_center)[[2]],
               # Step 6 Rotate the regression (SMA)
               # slope must come first
               slope_sma = slope / pearson,
               intercept_sma = mean_y - slope_sma * mean_x,
               # Step 7 Estimate Critical Soil Test Concentration
               cstv = exp(intercept_sma),
               # Step 8 Estimate the confidence interval
               pred_y = intercept_sma + slope_sma * x_center,
               mse    = sum((y - pred_y) ^ 2) / (n - 2),
               ssx    = var(x_center) * (n - 1),
               se     = sqrt(mse * ((1 / n) + ((
                   mean_x ^ 2
               ) / ssx))),
               lower_cl = exp(intercept_sma - se * qt(1 - 0.05 / 2, df = n - 2)),
               upper_cl = exp(intercept_sma + se * qt(1 - 0.05 / 2, df = n - 2)),
               # Step 9 Back-transform
               stv_sma = exp(pred_y),
               ry_sma = 100 * (sin(x_adjust + ((pred_y - intercept_sma) / slope_sma
               ))) ^ 2)
    
    return(steps_5_9)
}

# Function alcc_plot() creates a scatter plot of soil test correlation data
alcc_plot <- function(data, sufficiency = 95) {
    output <- alcc_sma(data, sufficiency = sufficiency)
    # for plot annotations
    cstv <- (unique(output$cstv))
    cstv <- if_else(cstv < 10, round(cstv, 1), round(cstv, 0))
    
    lower_cl <- (unique(output$lower_cl))
    lower_cl <- if_else(lower_cl < 10, round(lower_cl, 1), round(lower_cl, 0))
    
    upper_cl <- (unique(output$upper_cl))
    upper_cl <- if_else(upper_cl < 10, round(upper_cl, 1), round(upper_cl, 0))
    
    sufficiency <- unique(output$sufficiency)
    
    pearson <- round(unique(output$pearson), 2)
    
    # ggplot style
    output %>%
        ggplot(aes(stv, ry)) +
        geom_point(size = 2, alpha = 0.2) +
        geom_vline(xintercept = lower_cl,
                   alpha = 0.5,
                   linetype = 3) +
        geom_vline(xintercept = upper_cl,
                   alpha = 0.5,
                   linetype = 3) +
        # twice the CSTV
        geom_vline(xintercept = cstv * 2,
                   alpha = 0.5,
                   linetype = 2) +
        # the predicted back-transformed regression line
        geom_line(aes(x = stv_sma,
                      y = ry_sma),
                  color = red) +
        # the cstv at X% sufficiency
        geom_vline(xintercept = cstv, color = blue, alpha = 1) +
        geom_hline(yintercept = sufficiency, alpha = 0.2) +
        annotate(
            "text",
            label = paste0("CSTV = ", cstv, " ppm",
                          "\n95% CI [", lower_cl, " - ", upper_cl, "]",
                          "\nPearson correlation = ", pearson),
            x = max(output$stv),
            y = min(output$ry),
            family = "serif",
            hjust = 1,
            vjust = 0,
            alpha = 1
        ) +
        scale_x_continuous(breaks =
                               seq(0, 1000, if_else(
                                   max(output$stv) >= 200, 40,
                                   if_else(max(output$stv) >= 100, 20,
                                           if_else(max(output$stv) >= 50, 10, 5))
                               ))) +
        scale_y_continuous(breaks = seq(0, 105, 10)) +
        labs(
            x = expression(Soil ~ Test ~ Value ~ (mg ~ kg ^ -1)),
            y = "Relative yield (%)",
            title = "Critical soil test value with back-transformed correlation curve from Arcsine-Log Calibration Curve (ALCC)",
            caption = paste0(
                "Results determined for ",
                sufficiency,
                "% sufficiency.\nLower and upper 95% confidence limits shown by vertical dotted lines"
            )
        )
}

# While the plotting function could possibly be combined with the previous
# function, keeping them separate is simpler

##### DATA IMPORT #####
corr_data <-
    readxl::read_xlsx(path = "data_alcc_2022_01_26.xlsx") %>%
    rename(stv = soil_test, ry = relative_yield) %>%
    select(dataset, stv, ry, definition, everything())

corr_data %>% glimpse()

# Relative yield can be a ratio or percentage
# The ALCC method requires RY ratio values between 0-1,
# and percentage values between 0-100%

count(corr_data, ry > 100)

# Create new dataset from correlation data
# group_by and group_modify work together to analyze multiple datasets
alcc_data <- corr_data %>% 
    filter(dataset %in% c("correndo", "dyson"), trial_quality == "A") %>% 
    group_by(dataset, level, definition, nutrient, crop, method) %>%
    group_modify(~ alcc_sma(data = .x))

alcc_data

# alternatively, analyze a single dataset
single <- corr_data %>% 
    filter(dataset %in% c("VA-soybean-K"),
           trial_quality == "A",
           definition == "MAX"
           #definition == "FITMAX-CAP",
           #stv < 120
           )

single <- corr_data %>% 
    filter(dataset %in% c("dyson"),
           #trial_quality == "A"
           #stv < 120
    )

alcc_results <- alcc_sma(single, sufficiency = 95)

alcc_results

alcc_results %>% 
    ggplot(aes(x_center, y)) +
    geom_point(size = 2, alpha = 0.5) +
    geom_smooth(method = "lm") +
    geom_vline(xintercept = 0) +
    geom_hline(yintercept = alcc_results$intercept)

##### PLOT #####
# for a single dataset

alcc_plot(single, sufficiency = 90)

# alternatively, continue using group_by + group_map framework for analyzing multiple datasets seamlessly
single %>% 
    group_by(dataset, level, definition, nutrient, crop, method) %>%
    group_map(~ alcc_plot(data = .x, sufficiency = 90))


##### CREATE simplified results table for export #####
alcc_results %>%
    select(dataset,
           level,
           crop,
           nutrient,
           method,
           definition,
           sufficiency,
           cstv_adj = cstv,
           lower_cl,
           upper_cl,
           pvalue,
           pearson) %>%
    distinct(across(everything()))
    mutate(model = "ALCC")

##### Compare to Linear-plateau and Mitscherlich-Bray with bootstrap CI #####
library(nlraa)
library(minpack.lm) # in case nls doesn't converge

single %>% 
    ggplot(aes(stv, ry)) +
    geom_point(size = 2, alpha = 0.5) +
    # geom_line(stat = "smooth",
    #           method = "nls",
    #           formula = y ~ SSasymp(x, a, b, c),
    #           se = FALSE,
    #           color = red) +
    geom_line(stat = "smooth",
              method = "nls",
              formula = y ~ SSlinp(x, a, b, c),
              se = FALSE,
              color = blue) +
    scale_x_continuous(breaks = seq(0, 1000, 10))

### Linear-plateau
nls_lp <- nls(ry ~ SSlinp(stv, a, b, jp),
              data = single)

find_jp <- function(model) coef(model)[[3]]

find_jp(nls_lp)

# just going with simpler resampling of residuals though not ideal
boot_lp <- boot_nls(object = nls_lp, f = find_jp, data = single)

boot::boot.ci(boot_lp, type = "perc")
# Join point at 20.5 with 95% CI from 


### Mitscherlich-Bray
mb <- function(x, asym, b, c) {
    asym - b * exp(c * x)
}

nls_mitsch <- nls(ry ~ mb(stv, asym, b, c), 
                  data = single,
                  start = c(asym = 95, b = 50, c = -0.1))

find_cstv <- function(model){
    asym <- coef(model)[[1]]
    b <- coef(model)[[2]]
    c <- coef(model)[[3]]
    cstv <- log((asym - (asym * 95 / 100)) / b) / c
    return(cstv)
}

find_cstv(nls_mitsch)

boot_mitsch <- boot_nls(object = nls_mitsch, f = find_cstv, R = 1000, data = single)

boot::boot.ci(boot_mitsch, type = "perc")
