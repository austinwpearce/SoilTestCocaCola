#' The following function fits a Mitscherlich-type exponential model
#' May also be known as exponential, asymptotic, exponential rise to the max
#' It is designed for soil test correlation data 
#' This function can provide results in a table format or as a plot
#' Author: Austin Pearce
#' Last updated: 2022-05-31
#'
#' @name exponential
#' @param data a data frame with XY data
#' @param stv column for soil test values
#' @param ry column for relative yield
#' @param percent_of_max if wanting to find the X value for a point along the
#' quadratic portion at certain Y value
#' @param resid choose whether to create residuals plots
#' @param plot choose whether to create correlation plot rather than table
#' @param extrapolate choose whether the fitted line extends to X = 0
#' no effect if plot = FALSE
#' @export

# packages/dependencies needed
library(dplyr) # a suite of packages for wrangling and plotting
library(rlang) # evaluate column names for STV and RY (tip to AC)
library(nlraa) # for self-starting functions and predicted intervals
library(minpack.lm) # for nlsLM, a robust backup to nls
library(nlstools) # for residuals plots
library(modelr) # for the r-squared and rmse
library(ggplot2) # plots

# Colors for plot later on

red <- "#CE1141"
gold <- "#EAAA00"
blue <- "#13274F"
black <- "#000000"

# ========================================================

# "Mitscherlich" type model with three parameters
# CSTV is evaluated at Y some % of asymptote
# y = a + (b - a) * e^(-e^(c) * x)
# following the form used in SSasymp()
# alternatively we could use the mit() function
# a = horizontal asymptote (maximum yield potential)
# b = intercept when soil test is 0
# intercept is theoretical as soil never quite reaches 0
# c = curvature, natural log of the rate constant (ought to be negative in nls)

fun_exp <- function(x, a, b, c) {
    a + (b - a) * exp(-exp(c) * x)
}

get_exp_cstv <- function(model, pct_of_max = 95, target = NULL){
    a <- coef(model)[[1]]
    b <- coef(model)[[2]]
    c <- coef(model)[[3]]
    if(is.null(target)){
        cstv <- log(((a * pct_of_max / 100) - a)/(b - a))/(-exp(c))
    } else {
        cstv <- log((target - a)/(b - a))/(-exp(c))
    }
    return(cstv)
}

# SSasymp uses this formula, so for those interested in not providing starting values, that could be an option

exponential <- function(data = NULL,
                        stv,
                        ry,
                        percent_of_max = 95,
                        resid = FALSE,
                        plot = FALSE,
                        extrapolate = FALSE) {
    
    if (missing(stv)) {
        stop("Please specify the variable name for soil test concentrations using the `stv` argument")
    }
    
    if (missing(ry)) {
        stop("Please specify the variable name for relative yields using the `ry` argmuent")
    }
    
    if (percent_of_max < 1){
        warning("Please express percent on scale 0 - 100 rather than 0 - 1 (e.g. 95 = 95%)")
    }
    
    # Re-define x and y from STV and RY (tip to AC)
    x <- rlang::eval_tidy(data = data, rlang::quo({{stv}}) )
    
    y <- rlang::eval_tidy(data = data, rlang::quo({{ry}}) )
    
    if (max(y) < 2) {
        stop("The reponse variable appears to not be on a percentage scale.
             If so, please multiply it by 100.")
    }
    
    corr_data <- dplyr::tibble(x = as.numeric(x), 
                               y = as.numeric(y))
    
    if (nrow(corr_data) < 4) {
        stop("Too few distinct input values to fit LP. Try at least 4.")
    }
    
    minx <- min(corr_data$x) 
    maxx <- max(corr_data$x)
    rangex <- maxx - minx
    miny <- min(corr_data$y)
    maxy <- max(corr_data$y)
    start_c <- -(rangex) / (maxy - miny) / 2
    
    # build the model/fit ==================================================
    # even though the functions are selfStarting, providing starting values
    # increases the chance the SS functions converge on something reasonable
    # starting values shown are based on whether asymptote or origin are forced
    
    corr_model <- try(
        nlsLM(
        formula = y ~ fun_exp(x, a, b, c),
        data = corr_data,
        start = list(a = maxy, b = miny, c = start_c),
        upper = c(a = Inf, b = maxy, c = -1e-7), # force c to be negative is theoretical
        lower = c(a = miny, b = -Inf, c = -100) # at c = -100 the line would almost be flat
        )
        )
    
    if (inherits(corr_model, "try-error")) {
        stop("Mitscherlich-type exponential model could not be fit. Consider another model.")
    }
    
    # How did the model do overall?
    AIC      <- nlraa::IC_tab(corr_model)[3] %>% round()
    AICc     <- nlraa::IC_tab(corr_model)[4] %>% round()
    rmse     <- round(modelr::rmse(corr_model, corr_data), 2)
    rsquared <- round(modelr::rsquare(corr_model, corr_data), 2)
    
    # get model coefficients
    a <- coef(corr_model)[[1]]
    b <- coef(corr_model)[[2]]
    c <- exp(coef(corr_model)[[3]]) # equation based on natural log
    # derived values
    ry_pom <- a * percent_of_max / 100
    cx <- get_exp_cstv(corr_model, pct_of_max = percent_of_max)
    cstv <- round(cx, 0)
    
    # this equation is modified from original notation for c 
    equation <- paste0(round(a, 1), " + (", round(b,1), " - ", round(a, 1),
                       ") * e^(-", round(c, 3), "x)")
    
    # Table output =================================================
    if (plot == FALSE) {
        {
            if (resid == TRUE)
                plot(nlstools::nlsResiduals(corr_model), which = 0)
        }
        tibble(
            asymptote = round(a, 1),
            intercept = round(b, 2),
            rate_constant = round(c, 2),
            equation,
            cstv,
            ry_pom = round(ry_pom, 1),
            percent_of_max,
            AIC,
            AICc,
            rmse,
            rsquared
        )
    } else {
        # Residual plots and normality
        {
            if (resid == TRUE)
                plot(nlstools::nlsResiduals(corr_model), which = 0)
        }
        
        # To get fitted line from corr_model
        pred_y <- dplyr::tibble(x = seq(
            from = if_else(extrapolate == TRUE, 0, minx),
            to = maxx, by = 0.1)) %>%
            modelr::gather_predictions(corr_model)
        
        # ggplot of correlation
        exp_plot <- corr_data %>% 
            ggplot(aes(x, y)) +
            {
                if (extrapolate == TRUE)
                    geom_vline(xintercept = 0, alpha = 0.2)
            } +
            geom_vline(xintercept = cx,
                       alpha = 1,
                       color = blue) +
            geom_hline(yintercept = ry_pom,
                       alpha = 0.2) +
            # fitted line
            geom_line(data = pred_y,
                      aes(x, pred),
                      color = red) +
            geom_point(size = 2, alpha = 0.5) +
            geom_rug(alpha = 0.2, length = unit(2, "pt")) +
            scale_y_continuous(limits = c(0, maxy),
                               breaks = seq(0, maxy * 2, 10)) +
            scale_x_continuous(
                breaks = seq(0, maxx * 2, by = if_else(
                    condition = rangex >= 200,
                    true = 20,
                    false = if_else(
                        condition = rangex >= 100,
                        true = 10,
                        false = if_else(
                            condition = rangex >= 50,
                            true = 5,
                            false = 2))))) +
            annotate(
                "text",
                label = paste("CSTV =", cstv, "ppm"),
                x = cx,
                y = 0,
                angle = 90,
                hjust = 0,
                vjust = 1.5,
                alpha = 0.5
            ) +
            annotate(
                "text",
                label = paste0(percent_of_max, "% of asymptote = ",
                               round(ry_pom, 1), "% RY"),
                x = maxx,
                y = ry_pom,
                alpha = 0.5,
                vjust = 1.5,
                hjust = 1
            )+
            annotate(
                "text",
                alpha = 0.5,
                label = paste0("y = ", equation,
                               "\nAIC, AICc = ", AIC, ", ",AICc,
                               "\nRMSE = ", rmse,
                               "\nR-squared = ", rsquared),
                x = maxx,
                y = 0,
                vjust = 0,
                hjust = 1
            ) +
            labs(x = "Soil test value (mg/kg)",
                 y = "Relative yield (%)",
                 caption = paste("Each point is a site. n =", nrow(corr_data)))
        
        return(exp_plot)
    }
    
}
