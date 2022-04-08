#' The following function fits a linear plateau model
#' It is designed for soil test correlation data 
#' This function can provide results in a table format or as a plot
#' Author: Austin Pearce
#' Last updated: 2022-04-07
#'
#' @name lin_plateau
#' @param data a data frame with XY data
#' @param resid choose whether to create residuals plots
#' @param plot choose whether to create correlation plot rather than table
#' @param band choose whether the correlation plot displays confidence band
#' no effect if plot = FALSE
#' @export

# packages/dependencies needed
library(dplyr) # a suite of packages for wrangling and plotting
library(rlang) # evaluate column names for STV and RY (tip to AC)
library(nlraa) # for self-starting functions and predicted intervals
library(minpack.lm) # for nlsLM, a robust backup to nls
library(nlstools) # for residuals plots
library(modelr) # for the r-squared and rmse
library(ggplot2)

# Colors for plot later on
red <- "#CE1141"
gold <- "#EAAA00"
blue <- "#13274F"
black <- "#000000"

# =============================================================================
# Linear plateau model
# y = if{x < cx, b0 + b1*x; b0 + b1 * cx}
# b0 = intercept
# b1 = slope
# cx = critical X value = join point = Critical Soil Test Value (CSTV)

# =============================================================================
# lin_plateau function

lin_plateau <- function(data = NULL,
                        stv, # column for soil test values
                        ry,  # column relative yield
                        resid = FALSE,
                        plot = FALSE,
                        band = FALSE) {
    
    if (missing(stv)) {
        stop("Please specify the variable name for soil test concentrations using the `stv` argument")
    }
    
    if (missing(ry)) {
        stop("Please specify the variable name for relative yields using the `ry` argmuent")
    }
    
    # Re-define x and y from STV and RY (tip to AC)
    x <- rlang::eval_tidy(data = data, rlang::quo({{stv}}) )
    
    y <- rlang::eval_tidy(data = data, rlang::quo({{ry}}) )
    
    corr_data <- dplyr::tibble(x = x, y = y)
    
    if (nrow(corr_data) < 4) {
        stop("Too few distinct input values to fit LP. Try at least 4.")
    }
    
    minx <- min(corr_data$x)
    meanx <- mean(corr_data$x)
    maxx <- max(corr_data$x)
    miny <- min(corr_data$y)
    maxy <- max(corr_data$y)
    
    # build the model/fit =====
    # starting values
    sv <- list(b0 = miny, b1 = 1, cx = meanx)
    
    nls_model <- try(
        nls(y ~ SSlinp(x, b0, b1, cx),
            data = corr_data, start = sv)
        )
    
    if (inherits(nls_model, "try-error")) {
        corr_model <- try(
            minpack.lm::nlsLM(y ~ SSlinp(x, b0, b1, cx),
                              data = corr_data,
                              start = sv)
            )
    } else {
         corr_model <- nls_model
    }

    if (inherits(corr_model, "try-error")) {
        stop("LP model could not be fit with nls or nlsLM.
             Consider another model.")
    } else {
        corr_model <- corr_model
    }
    
    # Find p-value and pseudo R-squared
    AIC      <- nlraa::IC_tab(corr_model)[3] %>% round()
    AICc     <- nlraa::IC_tab(corr_model)[4] %>% round()
    rsquared <- round(modelr::rsquare(corr_model, corr_data), 2)
    rmse     <- round(modelr::rmse(corr_model, corr_data), 2)
    
    #booted <- nlraa::boot_nls(corr_model, data = corr_data)
    # confint(booted, type = "perc", level = 0.95)
    
    # get model coefficients
    b0 <- coef(corr_model)[[1]]
    b1 <- coef(corr_model)[[2]]
    cx <- coef(corr_model)[[3]]
    plateau <- b0 + b1 * cx
    cstv <- round(cx, 0)
    
    # CSTV at defined % of max/plateau (use percent_of_max argument)
    ##### currently ignoring this #####
    # 
    # cstv_adj <-
    #     round(((plateau * percent_of_max / 100) - b0) / b1, 0)
    # 
    # cstv_90ry <- if_else(
    #     condition = plateau >= 90,
    #     true = round((90 - b0) / b1, 0),
    #     false = NULL
    # )
    
    # makes an exact line rather than relying on predictions
    lp_line <- dplyr::tibble(x = c(minx, cx, maxx),
                             y = c(b0 + b1 * minx, plateau, plateau))
    
    equation <- paste0(round(b0, 1), " + ",
                       round(b1, 2), "x")
    
    # Table output =================================================
    if (plot == FALSE) {
        {
            if (resid == TRUE)
                plot(nlstools::nlsResiduals(corr_model), which = 0)
        }
        dplyr::tibble(
            intercept = round(b0, 2),
            slope = round(b1, 2),
            equation,
            cstv,
            plateau = round(plateau, 1),
            AIC,
            AICc,
            rmse,
            rsquared
            # cstv_adj,
            # percent_of_max,
            # cstv_90ry = if_else(cstv_90ry > 0,
            #                     cstv_90ry,
            #                     NULL),
            # at 90% RY
        )
    } else {
        # Residual plots and normality
        {
            if (resid == TRUE)
                plot(nlstools::nlsResiduals(corr_model), which = 0)
        }
        
        {
            if (band == TRUE)
                corr_pred <- dplyr::tibble(x = seq(minx, maxx, 10))
                pred_band <- nlraa::predict2_nls(
                    object = corr_model,
                    newdata = corr_pred,
                    interval = "confidence",
                    level = 0.95) %>%
                    as_tibble() %>% 
                    bind_cols(corr_pred)
        }
        
        # Alternative way to get model predictions for the fitted line
        # but the lp_line is cleaner
        # pred_y <- tibble(x = seq(minx, maxx, 0.1)) %>%
        #     modelr::gather_predictions(corr_model)
        
        lp_plot <- corr_data %>%
            ggplot(aes(x, y)) +
            {
                if (band == TRUE)
                    geom_ribbon(data = pred_band,
                                aes(x = x,
                                    y = Estimate,
                                    ymin = Q2.5,
                                    ymax = Q97.5),
                                #inherit.aes = FALSE,
                                alpha = 0.1)
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
                    condition = maxx >= 300,
                    true = 30,
                    false = if_else(
                        condition = maxx >= 100,
                        true = 20,
                        false = if_else(
                            condition = maxx >= 50,
                            true = 10,
                            false = 5))))) +
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
                                    "\nAIC, AICc = ", AIC, ", ",AICc,
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


fit_SSlinp <- function(data) {
    # nlraa model
    fit <- nlsLM(
        formula = y ~ SSlinp(x, b0, b1, cx),
        data = data,
        control = nls.lm.control(maxiter = 500),
        upper = c(
            b0 = max(data$y),
            b1 = 10000,
            cx = max(data$x)
        ),
        lower = c(
            b0 = -10000,
            b1 = 0,
            cx = min(data$x)
        )
    )
    
    return(fit)
}
