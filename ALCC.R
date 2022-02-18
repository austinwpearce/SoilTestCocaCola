#' Soil-test correlation using Arcsine-Log Calibration Curve (ALCC)
#' Author: Austin Pearce
#' Last update: 2022-02-18
# ================================================================

# package libraries needed (won't just work in base R)
library(agridat) # for obtaining a testing dataset
library(tidyverse) # a suite of packages for wrangling and plotting
library(nlraa) # for self-starting functions and predicted intervals
library(minpack.lm) # for nlsLM, a robust backup to nls
library(nlstools) # for residuals plots
library(modelr) # for the r-squared and rmse

# For now, the dataset needs `stv` (Soil Test Value) and `ry` (Relative Yield)
# as the starting names for the X and Y variables
# So prior to using this function, the X and Y columns that represent stv and ry
# will need to be renamed to stv and ry
# for example: %>% rename/mutate(stv = x, ry = y) %>%

# =============================================================================
theme_set(
    theme_minimal(base_size = 14) +
        theme(
            plot.background = NULL,
            plot.margin = margin(t = 2, r = 10, b = 2, l = 2, unit = "pt"),
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
            axis.text = element_text(),
            legend.title.align = 0,
            legend.key.height = unit(x = 5, units = "mm"),
            legend.justification = c(1, 1)
            #legend.position = c(1, 1)
        ))

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
                   alpha = 0.2,
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
            hjust = 1,
            vjust = 0,
            alpha = 0.8
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
cotton <- tibble(x = agridat::cate.potassium$potassium,
                 y = agridat::cate.potassium$yield, 
                 dataset = "cotton")

# Relative yield can be a ratio or percentage
# The ALCC method requires RY ratio values between 0-1,
# and percentage values between 0-100%

plot(cotton$y ~ cotton$x) |> abline(h = 100)
count(cotton, y > 100) # 3 site-years exceeded 100

cotton <- cotton %>% 
    # cap RY at 100
    mutate(stv = x, ry = if_else(y > 100, 100, y))

plot(cotton$ry ~ cotton$stv) |> abline(h = 100)

# Create new dataset from correlation data
alcc_results <- alcc_sma(cotton)

alcc_results

alcc_results %>% 
    ggplot(aes(x_center, y)) +
    geom_point(size = 2, alpha = 0.5) +
    geom_smooth(method = "lm") +
    geom_vline(xintercept = 0) +
    geom_hline(yintercept = alcc_results$intercept)

# Alternatively, if you have many groups/datasets in one table
alcc_results <- cotton %>% 
    group_by(dataset) %>%
    group_modify(~ alcc_sma(data = .x))

##### PLOT #####
# for a single dataset
alcc_plot(alcc_results, sufficiency = 95)

# alternatively, continue using group_by + group_map framework for analyzing multiple datasets seamlessly
alcc_results %>% 
    group_by(dataset) %>%
    group_map(~ alcc_plot(data = .x, sufficiency = 95))


##### CREATE simplified results table for export #####
alcc_results %>%
    select(dataset,
           sufficiency,
           cstv_adj = cstv,
           lower_cl,
           upper_cl,
           pvalue,
           pearson) %>%
    distinct(across(everything()))
    mutate(model = "ALCC")