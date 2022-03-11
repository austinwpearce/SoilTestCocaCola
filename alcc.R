#' @title modified Arcsine-Log Calibration Curve, or ALCC
#' @author Austin Pearce
#' Last update: 2022-03-10
#' @references Correndo et al. 2017
#' @references Dyson and Conyers 2013
#' @name alcc
#' @description perform ALCC/MALCC for soil test correlation
#' @description creates new variables on existing dataset
#' @param data a data frame with XY data
#' @param soil_test column for soil test values
#' @param ry column for relative yield values 0-100%
#' @param sma choose if Standardized Major Axis (SMA) is used for Modified ALCC (MALCC)
#' @param sufficiency choose at which RY value to get CSTV
#' @param confidence choose at which confidence level to estimate CI of CSTV
#' @param summary choose if full table or just summary should be returned
# =============================================================================
# Could potentially add an argument that checks if data is percentage or ratio


# package libraries needed (won't just work in base R)
library(tidyverse) # a suite of packages for wrangling and plotting

# =============================================================================
# ALCC
# =============================================================================
alcc_core <- function(data,
                      soil_test,
                      ry,
                      sma,
                      sufficiency,
                      confidence,
                      summary) {
    
    if (nrow(data) < 8) {
        stop("Too few distinct input values. Try at least 8.")
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
            stv = !!x,
            # RY values greater than 100 are capped at 100
            ry_cap = if_else(!!y > 100, 100, !!y),
            model = case_when(sma == TRUE ~ "MALCC",
                              sma == FALSE ~ "ALCC"),
            # get sample size
            n = n(),
            # Step 1 Transform (t for "transformed" added to x and y)
            # for ALCC, soil test goes to Y-axis and RY goes to X-axis
            xt = asin(sqrt(ry_cap / 100)),
            yt = log(stv),
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
            pred_yt  = intercept + slope * xt_centered,
            mse      = sum((yt - pred_yt) ^ 2) / (n - 2),
            ssx      = var(xt_centered) * (n - 1),
            se       = sqrt(mse * ((1 / n) + ((mean_xt ^ 2) / ssx))),
            lower_cl = exp(intercept - se * qt(1 - (1 - confidence / 100) / 2,
                                               df = n - 2)),
            upper_cl = exp(intercept + se * qt(1 - (1 - confidence / 100) / 2,
                                               df = n - 2)),
            # Step 9 Back-transform
            fitted_stv = exp(pred_yt),
            fitted_ry = 100 * (sin(
                adjust_by + ((pred_yt - intercept) / slope))) ^ 2
        ) %>%
        # 'dataset' might be problematic, not defined in scope
        select(model, 
               #dataset,
               sufficiency, cstv, lower_cl, upper_cl, 
               fitted_stv, fitted_ry, pvalue, pearson, everything())
}

alcc <- function(data,
                 soil_test,
                 ry,
                 sma = TRUE,
                 sufficiency = 90,
                 confidence = 95,
                 remove2x = FALSE,
                 summary = FALSE){
    
    if (missing(soil_test)) {
        stop("Specify the soil test variable (e.g., soil_test = STK)")
    }
    
    if (missing(ry)) {
        stop("Enter name of relative yield (%) variable (e.g., ry = RY)")
    }
    
    # enquo let's the user enter their own column names for x and y variables
    x <- enquo(soil_test)
    y <- enquo(ry)
    
    alcc_table <- alcc_core(data,
                            soil_test = !!x,
                            ry = !!y,
                            sma = sma,
                            sufficiency = sufficiency,
                            confidence = confidence)
    
    if (remove2x == TRUE) {
        cstv_2x <- (unique(alcc_table$cstv)) * 2
        # redo alcc() with data greater than twice the CSTV removed
        alcc_table <- alcc_table %>% 
            filter(stv <= cstv_2x) %>% 
            alcc(soil_test = stv,
                 ry = ry_cap,
                 sma = sma,
                 sufficiency = sufficiency,
                 confidence = confidence) %>% 
            mutate(remove2x = "TRUE")
    } else {
        alcc_table <- alcc_table %>% 
            mutate(remove2x = "FALSE")
    }
    
    if(summary == TRUE) {
        return(
            alcc_table %>%
                select(
                    model,
                    sufficiency,
                    cstv,
                    lower_cl,
                    upper_cl,
                    pvalue,
                    pearson,
                    remove2x) %>% 
                distinct(across(everything()))
        )
    } else {
        return(alcc_table)
    }
    
}
