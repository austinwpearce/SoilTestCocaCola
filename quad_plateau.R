#' The following function fits a quadratic plateau model
#' It is designed for soil test correlation data 
#' This function can provide results in a table format or as a plot
#' Author: Austin Pearce
#' Last updated: 2022-04-08
#'
#' @name quad_plateau
#' @param data a data frame with XY data
#' @param stv column for soil test values
#' @param ry column for relative yield
#' @param resid choose whether to create residuals plots
#' @param plot choose whether to create correlation plot rather than table
#' @param band choose whether the correlation plot displays confidence band
#' no effect if plot = FALSE
#' @param percent_of_max if wanting to find the X value for a point along the
#' quadratic portion at certain Y value
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
# The QP model and parameters
# y = if{x <= cx, b0 + b1x + b2x^2; b0 + b1*cx + b2*cx^2}
# b0 = intercept
# b1 = slope
# b2 = quadratic term = -0.5 * b1 / cx
# cx = critical X value = join point = Critical Soil Test Value (CSTV) 
# cx = -0.5 * b1 / b2

quad_plateau <- function(data = NULL,
                         stv,
                         ry,
                         resid = FALSE,
                         plot = FALSE,
                         band = FALSE,
                         percent_of_max = 95) {
    
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
    # even though the functions are selfStarting, providing starting values
    # increases the chance the SS functions converge on something reasonable
    # starting values (sv)
    sv <- list(b0 = miny, b1 = 1, cx = meanx)
    
    nls_model <- try(
        nls(y ~ SSquadp3xs(x, b0, b1, cx),
            data = corr_data,
            start = sv))
    # SSquadp3() also an option, especially if setting bounds on b2 param
    
    if (inherits(nls_model, "try-error")) {
        corr_model <- try(
            minpack.lm::nlsLM(y ~ SSquadp3xs(x, b0, b1, cx),
                              data = corr_data,
                              start = sv))
    } else {
        corr_model <- nls_model
    }
    
    if (inherits(corr_model, "try-error")) {
        stop("QP model could not be fit with nls or nlsLM.
             Consider another model.")
    } else {
        corr_model <- corr_model
    }
    
    # How did the model do overall?
    AIC      <- nlraa::IC_tab(corr_model)[3] %>% round()
    AICc     <- nlraa::IC_tab(corr_model)[4] %>% round()
    rmse     <- round(modelr::rmse(corr_model, corr_data), 2)
    rsquared <- round(modelr::rsquare(corr_model, corr_data), 2)
    
    # get model coefficients
    b0 <- coef(corr_model)[[1]]
    b1 <- coef(corr_model)[[2]]
    cx <- coef(corr_model)[[3]]
    # derived values
    b2 <- -0.5 * b1 / cx
    plateau <- b0 + (b1 * cx) + (b2 * cx * cx)
    cstv <- round(cx, 0)
    
    equation <- paste0(round(b0, 1), " + ",
                       round(b1, 2), "x + ",
                       round(b2, 3), "x^2")
    
    # CSTV at defined % of max/plateau
    # To find an X value at a given Y less than predicted plateau
    newplateau <- plateau * percent_of_max / 100
    discriminant <- (b1 ^ 2) - (4 * (b0 - newplateau) * b2)
    cstv_pom <- (-b1 + sqrt(discriminant)) / (2 * b2)
    
    # Printouts
    if (plot == FALSE) {
        {
            if (resid == TRUE)
                plot(nlsResiduals(corr_model), which = 0)
        }
        tibble(
            intercept = round(b0, 2),
            slope = round(b1, 2),
            curve = round(b2, 4),
            equation,
            cstv,
            plateau = round(plateau, 1),
            AIC,
            AICc,
            rmse,
            rsquared,
            cstv_pom = round(cstv_pom, 0),
            percent_of_max
        )
    } else {
        # Residual plots and normality
        {
            if (resid == TRUE)
                plot(nlsResiduals(corr_model), which = 0)
        }
        
        {
            if (band == TRUE)
                conf_band <- nlraa::predict2_nls(
                    object = corr_model,
                    newdata = corr_data,
                    interval = "confidence",
                    level = 0.95) %>%
                dplyr::as_tibble() %>%
                dplyr::bind_cols(corr_data)
        }
        
        # To get fitted line from corr_model
        pred_y <- dplyr::tibble(x = seq(minx, maxx, 0.1)) %>%
            modelr::gather_predictions(corr_model)
        
        # ggplot of correlation
        qp_plot <- corr_data %>%
            ggplot(aes(x, y)) +
            {
                if (band == TRUE)
                    geom_ribbon(data = conf_band,
                                aes(x = x,
                                    y = Estimate,
                                    ymin = Q2.5,
                                    ymax = Q97.5),
                                alpha = 0.1)
            } +
            geom_vline(xintercept = cx,
                       alpha = 1,
                       color = blue) +
            geom_hline(yintercept = plateau,
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
                    condition = maxx >= 300,
                    true = 30,
                    false = if_else(
                        condition = maxx >= 100,
                        true = 20,
                        false = if_else(
                            condition = maxx >= 50,
                            true = 10,
                            false = 5))))) +
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
                label = paste0("y = ", equation,
                               "\nAIC, AICc = ", AIC, ", ",AICc,
                               "\nRMSE = ", rmse,
                               "\nR-squared = ", rsquared),
                x = maxx,
                y = 0,
                vjust = 0,
                hjust = 1
            ) +
            labs(
                x = "Soil test value (mg/kg)",
                y = "Relative yield (%)",
                caption = paste("Each point is a site. n =", nrow(corr_data))
            )
        
        return(qp_plot)
    }
    
}

# =============================================================================
# preferred theme for ggplot

theme_set(
    theme_minimal(base_size = 14) +
        theme(
            plot.background = NULL,
            plot.margin = margin(
                t = 2,
                r = 10,
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
            axis.text = element_text(),
            legend.title.align = 0,
            legend.key.height = unit(x = 5, units = "mm"),
            legend.justification = c(1, 1)
            #legend.position = c(1, 1)
        )
)

# =============================================================================
# other functions for fitting nls model only
# 
# qp <- function(x, b0, b1, cx) {
#     b2 <- -0.5 * b1 / cx
#     if_else(
#         condition = x < cx,
#         true  = b0 + (b1 * x) + (b2 * x * x),
#         false = b0 + (b1 * cx) + (b2 * cx * cx)
#     )
# }
# 
# fit_qp <- function(data) {
#     start <- lm(y ~ poly(x, 2, raw = TRUE), data = data)
#     # nls model
#     fit <- nlsLM(
#         formula = y ~ qp(x, b0, b1, cx),
#         data = data,
#         start = list(
#             b0 = start$coef[[1]],
#             b1 = start$coef[[2]],
#             cx = mean(data$x)
#         ),
#         control = nls.lm.control(maxiter = 500)
#     )
#     return(fit)
# }
# 
# fit_SSquadp3xs <- function(data) {
#     # nlraa model
#     fit <- nlsLM(
#         formula = y ~ SSquadp3xs(x, b0, b1, cx),
#         data = data,
#         control = nls.lm.control(maxiter = 500)
#         # upper = c(b0 = max(data$y),
#         #           b1 = 1000,
#         #           b2 = 0),
#         # lower = c(b0 = -1000,
#         #           b1 = 0,
#         #           b2 = -Inf)
#     )
#     return(fit)
# }