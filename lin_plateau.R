# package libraries needed
library(tidyverse) # a suite of packages for wrangling and plotting
library(nlraa) # for self-starting functions and predicted intervals
library(minpack.lm) # for nlsLM, a robust backup to nls
library(nlstools) # for residuals plots
library(modelr) # for the r-squared and rmse

# For now, everything is written in terms of x and y
# So prior to using this function, the columns for soil test value and RY, for
# example, will need to be renamed to x and y
# for example: %>% rename/mutate(x = stv, y = ry) %>%

# =============================================================================

#' The following function fits a linear plateau model to soil test
#' correlation XY data and provides results in a table format or as a plot
#'
#' Last updated: 2022-02-18
#'
#' This function is essentially a wrapper that uses other packages' functions
#' Won't just work in base R
#'
#' @name lin_plateau
#' @param data a data frame with XY data
#' @param resid choose whether to create residuals plots
#' @param plot choose whether to create correlation plot rather than table
#' @param band choose whether the correlation plot displays confidence band
#' no effect if plot = FALSE
#' @export


# Linear plateau model
# y = b0 + b1x
# b0 = intercept
# b1 = slope
# jp = join point = critical concentration

# =============================================================================
# lin_plateau function

lin_plateau <- function(data,
                        resid = FALSE,
                        plot = FALSE,
                        band = FALSE) {
    minx <- min(data$x)
    maxx <- max(data$x)
    miny <- min(data$y)
    maxy <- max(data$y)
    
    # build the model/fit =====
    corr_model <- try(nls(y ~ SSlinp(x, b0, b1, jp), data = data))
    
    if (inherits(corr_model, "try-error")) {
        corr_model <-
            try(minpack.lm::nlsLM(y ~ SSlinp(x, b0, b1, jp), data = data))
    } else {
        corr_model <- corr_model
    }
    
    # Find p-value and pseudo R-squared
    AIC      <- round(AIC(corr_model), 0) 
    rsquared <- round(modelr::rsquare(corr_model, data), 2)
    rmse     <- round(modelr::rmse(corr_model, data), 2)
    
    # booted <- nlraa::boot_nls(corr_model, data = data)
    # confint(booted, type = "perc", level = 0.95)
    
    # get model coefficients
    b0 <- coef(corr_model)[[1]]
    b1 <- coef(corr_model)[[2]]
    jp <- coef(corr_model)[[3]]
    join_point <- round(jp, 0)
    plateau <- b0 + b1 * jp
    
    # CSTV at defined % of max/plateau (use percent_of_max argument)
    ##### currently ignoring this #####
    # 
    # cstv_adj <-
    #     round(((plateau * percent_of_max / 100) - b0) / b1, 0)
    # 
    # cstv_95ry <- if_else(
    #     condition = plateau >= 95,
    #     true = round((95 - b0) / b1, 0),
    #     false = NULL
    # )
    # 
    # cstv_90ry <- if_else(
    #     condition = plateau >= 90,
    #     true = round((90 - b0) / b1, 0),
    #     false = NULL
    # )
    
    # have to make a line because the SSlinp doesn't plot right
    lp_line <- tibble(x = c(minx, jp, maxx),
                      y = c(b0 + b1 * minx, plateau, plateau))
    
    equation <- paste0(round(b0, 1), " + ",
                       round(b1, 2), "x")
    
    # Table output =====
    if (plot == FALSE) {
        {
            if (resid == TRUE)
                plot(nlstools::nlsResiduals(corr_model), which = 0)
        }
        tibble(
            intercept = b0,
            slope = b1,
            equation,
            join_point,
            plateau = round(plateau, 1),
            AIC,
            rmse,
            rsquared
            # cstv_adj,
            # percent_of_max,
            # cstv_95ry = if_else(cstv_95ry > 0,
            #                     cstv_95ry,
            #                     NULL),
            # # at 95% RY
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
        predicted <- 
            nlraa::predict2_nls(corr_model,
                                newdata = data,
                                interval = "confidence") %>%
            as_tibble() %>%
            bind_cols(data)
        
        lp_plot <- predicted %>%
            ggplot(aes(x, y)) +
            {
                if (band == TRUE)
                    geom_ribbon(aes(
                        y = Estimate,
                        ymin = Q2.5,
                        ymax = Q97.5
                    ),
                    alpha = 0.05)
            } +
            geom_vline(xintercept = jp,
                       alpha = 0.5,
                       linetype = 3) +
            geom_hline(yintercept = plateau,
                       alpha = 0.5,
                       linetype = 3) +
            geom_line(data = lp_line,
                      aes(x = x, y = y),
                      color = "#CC0000") +
            geom_point(size = 3, alpha = 0.5) +
            geom_rug(alpha = 0.2, length = unit(2, "pt")) +
            scale_y_continuous(limits = c(0, maxy),
                               breaks = seq(0, maxy * 2, 10)) +
            annotate(
                "text",
                label = paste("CSTV =", join_point, "ppm"),
                x = join_point,
                y = 0,
                angle = 90,
                hjust = 0,
                vjust = 1.5,
                alpha = 0.5
            ) +
            annotate(
                "text",
                label = paste0("Plateau = ", round(plateau, 1), "%"),
                x = maxx,
                y = plateau,
                alpha = 0.5,
                hjust = 1,
                vjust = 1.5
            ) +
            annotate(
                "text",
                alpha = 0.5,
                label = paste0(
                    "y = ",
                    equation,
                    "\nAIC = ",
                    AIC,
                    "\nRMSE = ",
                    rmse,
                    "\nR-squared = ",
                    rsquared
                ),
                x = maxx,
                y = 0,
                vjust = 0,
                hjust = 1
            ) +
            labs(
                x = "Soil test value (mg/kg)",
                y = "Relative yield (%)",
                caption = paste("Each point is a site. n =", nrow(data))
            )
        
        return(lp_plot)
    }
    
}


fit_SSlinp <- function(data) {
    # nlraa model
    fit <- nlsLM(
        formula = y ~ SSlinp(x, b0, b1, jp),
        data = data,
        control = nls.lm.control(maxiter = 500),
        upper = c(
            b0 = max(data$y),
            b1 = 10000,
            jp = max(data$x)
        ),
        lower = c(
            b0 = -10000,
            b1 = 0,
            jp = min(data$x)
        )
    )
    
    return(fit)
}
