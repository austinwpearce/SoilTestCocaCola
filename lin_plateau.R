#' The following function fits a linear plateau model
#' It was designed with soil test correlation data in mind
#' This function can provide results in a table format or as a plot
#' Consider using the soiltestcorr package which is better maintained and less experimental
#' Author: Austin Pearce
#' Last updated: 2022-06-01
#'
#' @name lin_plateau
#' @param data a data frame with XY data
#' @param x column for soil test values
#' @param y column for response (e.g. relative yield)
#' @param force100 force model to plateau at 100% RY
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

# =============================================================================
# supporting functions
# Linear plateau model
# y = if{x <= cx, a + b1x; a + b * cx}
# a = intercept
# b = slope
# cx = critical X value = join point = Critical Soil Test Value (CSTV)

lp <- function(x, a, b, cx){
    if_else(x < cx,
            a + b * x,
            a + b * cx)
}

# Linear-plateau (plateau = 100)
lp100 <- function(x, a, b){
    cx <- (100 - a) / b
    if_else(x < cx, 
            a + b * x,
            a + b * cx)
}

# or use nlraa::SSlinp for self-starting function

get_lp_plateau <-function(lp_model){
    p <- coef(lp_model)[[1]] + coef(lp_model)[[2]] * coef(lp_model)[[3]]
    return(round(p, 1))
}

get_lp_cstv <- function(model, pct_of_max = 100, target = NULL){
    a <- coef(model)[[1]]
    b <- coef(model)[[2]]
    cx <- coef(model)[[3]]
    plateau <- a + b * cx
    if(is.null(target)){
        cstv <- round(((plateau * pct_of_max / 100) - a) / b, 1)
    } else {
        cstv <- if_else(target > plateau,
                        Inf,
                        round((target - a) / b, 1))
    }
    
    return(cstv)
}

# =============================================================================

lin_plateau <- function(data = NULL,
                        x,
                        y,
                        force100 = FALSE,
                        resid = FALSE,
                        plot = FALSE,
                        extrapolate = FALSE) {
    
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
    
    minx <- min(corr_data$x)
    meanx <- mean(corr_data$x)
    maxx <- max(corr_data$x)
    rangex <- maxx - minx
    miny <- min(corr_data$y)
    maxy <- max(corr_data$y)
    
    # build the model/fit =====
    # starting values (sv)
    # even though the functions are selfStarting, providing starting values
    # increases the chance the SS functions converge on something reasonable
    sv <- list(a = miny, b = 1, cx = meanx)
    
    # even though there is a risk that nlsLM results in a false convergence, this risk is likely low
    if (force100 == FALSE) {
        corr_model <- try(minpack.lm::nlsLM(y ~ SSlinp(x, a, b, cx),
                                        data = corr_data,
                                        start = sv,
                                        upper = c(a = maxy, b = Inf, cx = Inf),
                                        lower = c(a = -Inf, b = 0, cx = minx)),
                         silent = TRUE)
    } else {
        
        corr_model <- try(minpack.lm::nlsLM(y ~ lp100(x, a, b),
                                        data = corr_data,
                                        start = c(a = miny, b = 1),
                                        upper = c(a = 100, b = Inf),
                                        lower = c(a = -Inf, b = 0)),
                         silent = TRUE)
    }

    if (inherits(corr_model, "try-error")) {
        stop("LP model could not be fit with nlsLM.
             Consider another model.")
    } else {
        corr_model <- corr_model
    }
    
    # How did the model do overall?
    AICc     <- nlraa::IC_tab(corr_model)[4] %>% round(1)
    rmse     <- round(modelr::rmse(corr_model, corr_data), 2)
    rsquared <- round(modelr::rsquare(corr_model, corr_data), 2)
    
    #booted <- nlraa::boot_nls(corr_model, data = corr_data)
    # confint(booted, type = "perc", level = 0.95)
    
    # get model coefficients
    a <- coef(corr_model)[[1]]
    b <- coef(corr_model)[[2]]
    
    if (force100 == FALSE){
        cx <- coef(corr_model)[[3]]
    } else {
        cx <- (100 - a) / b
    }
    
    if (force100 == FALSE){
        plateau <- a + b * cx
    } else {
        plateau <- 100
    }
    
    cstv <- round(cx, 1)
    
    # makes an exact line with clean bend rather than relying on predictions
    lp_line <- dplyr::tibble(
        x = c(if_else(extrapolate == TRUE, 0, minx), cx, maxx),
        y = c(if_else(extrapolate == TRUE, a, a + b * minx), plateau, plateau))
    
    equation <- paste0(round(a, 1), " + ",
                       round(b, 2), "x")
    
    # Table output =================================================
    if (plot == FALSE) {
        {
            if (resid == TRUE)
                plot(nlstools::nlsResiduals(corr_model), which = 0)
        }
        dplyr::tibble(
            intercept = round(a, 2),
            slope = round(b, 2),
            equation,
            cstv,
            plateau = round(plateau, 1),
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
        
        ## ggplot of correlation
        lp_plot <- corr_data %>%
            ggplot(aes(x, y)) +
            {
                if (extrapolate == TRUE)
                    geom_vline(xintercept = 0, alpha = 0.2)
            } +
            geom_vline(xintercept = cx,
                       alpha = 1,
                       color = blue) +
            geom_hline(yintercept = plateau,
                       alpha = 0.2) +
            # fitted line
            geom_line(data = lp_line,
                      aes(x = x, y = y),
                      color = red) +
            geom_point(size = 2, alpha = 0.5) +
            geom_rug(alpha = 0.2, length = unit(2, "pt")) +
            scale_y_continuous(
                # start from 0 helps show the overall response
                limits = c(0, maxy),
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
            annotate("text",
                     label = paste("CSTV =", cstv, "ppm"),
                     x = cstv,
                     y = 0,
                     angle = 90,
                     hjust = 0,
                     vjust = 1.5,
                     alpha = 0.5) +
            annotate("text",
                     alpha = 0.5,
                     label = paste0("Plateau = ", round(plateau, 1), "%"),
                     x = maxx,
                     y = plateau,
                     hjust = 1,
                     vjust = 1.5) +
            annotate("text",
                     alpha = 0.5,
                     label = paste0("y = ", equation,
                                    "\nAICc = ",AICc,
                                    "\nRMSE = ", rmse,
                                    "\nR-squared = ", rsquared),
                     x = maxx,
                     y = 0,
                     vjust = 0,
                     hjust = 1) +
            labs(
                x = "Soil test value (mg/kg)",
                y = "Relative yield (%)",
                caption = paste("Each point is a site. n =", nrow(corr_data))
            )
        
        return(lp_plot)
    }
    
}

lin_plateau(corr_data, x, y, plot = TRUE)
lin_plateau(corr_data, x, y, plot = TRUE, force100 = TRUE)
