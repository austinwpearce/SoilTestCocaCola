#' Soil-test correlation using Arcsine-Log Calibration Curve (ALCC)
#' Author: Austin Pearce
#' Last update: 2022-02-18
# ================================================================
# 
# package libraries needed (won't just work in base R)
require(tidyverse) # a suite of packages for wrangling and plotting
# these packages are for completing code examples after the alcc stuff
require(agridat) # for obtaining a testing dataset
library(nlraa) # for self-starting functions and predicted intervals
library(minpack.lm) # for nlsLM, a robust backup to nls
library(nlstools) # for residuals plots
library(modelr) # for the r-squared and rmse

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

#' Doesn't work without installing other packages just work in base R
#'
#' @name alcc
#' @author Austin Pearce
#' @author Dyson and Conyers 2013
#' @author Correndo et al. 2017
#' @param data a data frame with XY data
#' @param x
#' @param y
#' @param sma choose if Standardized Major Axis is used
#' @param sufficiency choose at which RY value to get CSTV and estimate CI
#' @param summary choose if full table or just summary should be returned
#' Could potentially add an argument that checks if data is percentage or ratio
##### Two fancy functions #####
# Function alcc_sma() creates new variables on existing dataset
alcc <- function(data,
                 soil_test,
                 ry,
                 sma = TRUE,
                 sufficiency = 90,
                 summary = FALSE) {
    if (nrow(data) < 8) {
        stop("Too few distinct input values to perform ALCC-SMA.
             Try at least 8.")
    }
    
    if (missing(soil_test)) {
        stop("Specify the soil test variable (e.g., soil_test = STK)")
    }
    
    if (missing(ry)) {
        stop("Enter name of relative yield (%) variable (e.g., ry = RY)")
    }
    
    # enquo let's the user enter their own column names for x and y variables
    x <- enquo(soil_test)
    y <- enquo(ry)
    
    input <- data %>% 
        select(!!x, !!y)
    
    if (max(input[,2]) > 100) {
        warning("One or more original RY values exceeded 100%.\nAll RY values greater than 100% have been capped to 100% for arcsine transformation.", call. = FALSE)
    }
    
    steps_1_4 <- data %>%
        mutate(
            # RY values greater than 100 are capped at 100
            ry_cap = if_else(!!y > 100, 100, !!y),
            model = case_when(sma == TRUE ~ "MALCC",
                              sma == FALSE ~ "ALCC"),
            # get sample size
            n = n(),
            # Step 1 Transform (t for "transformed" added to x and y)
            # for ALCC, soil test goes to Y-axis and RY goes to X-axis
            xt = asin(sqrt(ry_cap / 100)),
            yt = log(!!x),
            # Step 2 Center
            sufficiency = sufficiency,
            adjust_by = asin(sqrt(sufficiency / 100)),
            xt_centered = xt - adjust_by,
            # Step 3 Correlation
            pearson = cor(xt_centered, yt, method = "pearson"),
            t_stat  = (pearson * sqrt(n - 2)) / sqrt(1 - pearson ^ 2),
            pvalue  = pt(t_stat, df = n - 1, lower.tail = FALSE),
            # Step 4 Means
            mean_xt  = mean(xt_centered),
            mean_yt  = mean(yt)
        )

    # Step 5 Fit linear model to transformed, centered data
    ols_center <- lm(yt ~ xt_centered, data = steps_1_4)

    steps_5_9 <- steps_1_4 %>%
        mutate(
            intercept = coef(ols_center)[[1]],
            slope = coef(ols_center)[[2]],
            # Step 6 Rotate the regression (SMA)
            # slope must come first
            slope = case_when(sma == TRUE ~ slope / pearson,
                              sma == FALSE ~ slope),
            intercept = case_when(sma == TRUE ~ mean_yt - slope * mean_xt,
                                  sma == FALSE ~ intercept),
            # Step 7 Estimate Critical Soil Test Concentration
            cstv = exp(intercept),
            # Step 8 Estimate the confidence interval
            pred_yt = intercept + slope * xt_centered,
            mse    = sum((yt - pred_yt) ^ 2) / (n - 2),
            ssx    = var(xt_centered) * (n - 1),
            se     = sqrt(mse * ((1 / n) + ((mean_xt ^ 2) / ssx))),
            lower_cl = exp(intercept - se * qt(1 - 0.05 / 2, df = n - 2)),
            upper_cl = exp(intercept + se * qt(1 - 0.05 / 2, df = n - 2)),
            # Step 9 Back-transform
            fitted_stv = exp(pred_yt),
            fitted_ry = 100 * (sin(
                adjust_by + ((pred_yt - intercept) / slope))) ^ 2
        ) %>%
        # 'dataset' might be problematic, not defined in scope
        select(model, dataset, sufficiency, cstv, lower_cl, upper_cl, 
               fitted_stv, fitted_ry, pvalue, pearson, everything())
    
    alcc_summary <- steps_5_9 %>%
        select(
            model,
            dataset,
            sufficiency,
            cstv,
            lower_cl,
            upper_cl,
            pvalue,
            pearson) %>% 
    distinct(across(everything()))

    if(summary == TRUE) {
        return(alcc_summary)
    } else {
        return(steps_5_9)
    }
}

# Function alcc_plot() creates a scatter plot of soil test correlation data
alcc_plot <- function(data, soil_test, ry, sma = TRUE, sufficiency = 90) {
    
    if (missing(soil_test)) {
        stop("Specify the soil test variable (e.g., soil_test = STK)")
    }
    
    if (missing(ry)) {
        stop("Enter name of relative yield (%) variable (e.g., ry = RY)")
    }
    
    # enquo let's the user enter their own column names for x and y variables
    x <- enquo(soil_test)
    y <- enquo(ry)
    # pass to alcc function
    output <- alcc(data, soil_test = !!x, ry = !!y,
                   sma = sma,
                   sufficiency = sufficiency)
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
        ggplot(aes(!!x, !!y)) +
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
        # fitted values from back-transformed regression line from alcc()
        geom_line(aes(x = fitted_stv,
                      y = fitted_ry),
                  color = "#CC0000") +
        # the cstv at X% sufficiency
        geom_vline(xintercept = cstv, color = "#4156A1", alpha = 1) +
        geom_hline(yintercept = sufficiency, alpha = 0.2) +
        annotate(
            "text",
            label = paste0("CSTV = ", cstv, " ppm",
                          "\n95% CI [", lower_cl, " - ", upper_cl, "]",
                          "\nPearson correlation = ", pearson),
            x = cstv * 2,
            y = 0,
            hjust = 0,
            vjust = 0,
            alpha = 0.8
        ) +
        # scale_x_continuous(breaks =
        #                        seq(0, 1000, if_else(
        #                            max(output$stv) >= 200, 40,
        #                            if_else(max(output$stv) >= 100, 20,
        #                                    if_else(max(output$stv) >= 50, 10, 5))
        #                        ))) +
        scale_y_continuous(breaks = seq(0, 105, 10)) +
        labs(
            x = expression(Soil ~ Test ~ Value ~ (mg ~ kg ^ -1)),
            y = "Relative yield (%)",
            caption = paste0(
                "Results determined for ",
                sufficiency,
                "% sufficiency.\nLower and upper 95% confidence limits shown by vertical dotted lines"
            )
        )
}
