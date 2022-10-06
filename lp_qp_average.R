#' This experimental function fits two models: a linear- and quadratic-plateau
#' and then simply averages them, without bootstrapping
#' It was designed with soil test correlation data in mind
#' This function can provide results in a table format or as a plot
#' Consider using the soiltestcorr package
#' Author: Austin Pearce
#' Last updated: 2022-10-06
#'
#' @name lp_qp_average
#' @param data a data frame with XY data
#' @param x column for soil test values
#' @param y column for response (e.g. relative yield)
#' @param confint estimate a 95% confidence interval by bootstrap
#' @param boot_R number of bootstrap replicates
#' @param resid choose whether to create residuals plots
#' @param plot choose whether to create correlation plot rather than table
#' @param extrapolate choose whether the fitted line extends to X = 0
#' no effect if plot = FALSE
#' @export

# packages/dependencies needed
library(dplyr) # a suite of packages for wrangling and plotting
library(tidyr) # tidying functions
library(purrr) # map functions
library(rlang) # evaluate column names for STV and RY (tip to AC)
library(nlraa) # for self-starting functions and predicted intervals
library(minpack.lm) # for nlsLM, a robust backup to nls
library(nlstools) # for residuals plots
library(rsample)
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

get_plateau_lp <-function(lp_model){
    p <- coef(lp_model)[[1]] + coef(lp_model)[[2]] * coef(lp_model)[[3]]
    return(round(p, 1))
}

# =============================================================================

lp_qp_average <- function(data = NULL,
                          x,
                          y,
                          weighted = TRUE,
                          increment = 1,
                          # confint = FALSE,
                          # boot_R = 500,
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
    lp_model <<- try(minpack.lm::nlsLM(y ~ SSlinp(x, a, b, cx),
                                      data = corr_data,
                                      start = sv,
                                      upper = c(a = maxy, b = Inf, cx = 10 * maxx),
                                      lower = c(a = -Inf, b = 0, cx = minx)),
                    silent = TRUE)
    
    qp_model <<- try(minpack.lm::nlsLM(y ~ SSquadp3xs(x, a, b, cx),
                                      data = corr_data,
                                      start = sv,
                                      upper = c(a = maxy, b = Inf, cx = 10 * maxx),
                                      lower = c(a = -Inf, b = 0, cx = minx)),
                    silent = TRUE)
    
    if (inherits(lp_model, "try-error") | inherits(qp_model, "try-error")) {
        stop("One of the two model could not be fit with nlsLM.
             Consider another model.")
    }
    
    lp_cstv <- coef(lp_model)[[3]]
    qp_cstv <- coef(qp_model)[[3]]
    
    ictab <<- IC_tab(lp_model, qp_model) %>% as_tibble()
    
    if (weighted == TRUE) {
        lp_weight <- ictab %>% filter(model == "lp_model") %>% pull(weight)
        qp_weight <- ictab %>% filter(model == "qp_model") %>% pull(weight)
    } else {
        lp_weight <- 0.5
        qp_weight <- 0.5
    }
    
    cstv <- round((lp_cstv * lp_weight) + (qp_cstv * qp_weight))
        
    # predicted values
    pred_dat <- tibble(x = seq(0, maxx, increment))
    
    if (weighted == TRUE) {
        pred_curve <- predict_nls(lp_model,
                                  qp_model,
                                  interval = "conf",
                                  newdata = ndat,
                                  criteria = "AIC") %>%
            as_tibble() %>%
            bind_cols(ndat)
    } else {
        pred_curve <- predict_nls(lp_model,
                             qp_model,
                             interval = "conf",
                             newdata = ndat,
                             weights = c(0.5, 0.5)) %>%
        as_tibble() %>% 
        bind_cols(ndat)
    }
    
    # # 95% Bootstrap confidence intervals
    # if (confint == TRUE & force100 == FALSE) {
    #     fit_LP <- function(split) {
    #         fit <- nlsLM(formula = y ~ lp(x, a, b, cx),
    #                      data = analysis(split),
    #                      start = as.list(coef(corr_model)))
    #         
    #         return(fit)
    #     }
    #     set.seed(911)
    #     
    #     boot_ci <- corr_data %>%
    #         bootstraps(times = boot_R) %>% 
    #         mutate(models = map(splits, possibly(fit_LP, otherwise = NULL)),
    #                coefs = map(models, tidy)) %>% 
    #         int_pctl(coefs)
    #     
    #     lcl <- boot_ci$.lower[3]
    #     ucl <- boot_ci$.upper[3]
    # }
    
    # Table output =================================================
    ## ggplot of correlation
        avg_plot <- corr_data %>%
            ggplot(aes(x, y)) +
            {
                if (extrapolate == TRUE)
                    geom_vline(xintercept = 0, alpha = 0.2)
            } +
            geom_vline(xintercept = cstv,
                       alpha = 1,
                       color = blue) +
            # {
            #     if(confint == TRUE)
            # geom_vline(xintercept = c(lcl, ucl),
            #            alpha = 0.8,
            #            color = blue,
            #            linetype = 3)
            # } +
            # fitted line
            geom_line(data = pred_curve,
                      aes(x = x, y = Estimate),
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
            # {
            #     if(confint == TRUE)
            # annotate("text",
            #          label = paste("LCL =", round(lcl,1), "ppm"),
            #          x = lcl,
            #          y = 0,
            #          angle = 90,
            #          hjust = 0,
            #          vjust = -0.5,
            #          alpha = 0.5)
            # } + {
            #     if(confint == TRUE)
            # annotate("text",
            #          label = paste("UCL =", round(ucl,1), "ppm"),
            #          x = ucl,
            #          y = 0,
            #          angle = 90,
            #          hjust = 0,
            #          vjust = 1.5,
            #          alpha = 0.5)
            # } +
            # annotate("text",
            #          alpha = 0.5,
            #          label = paste0("Plateau = ", round(plateau, 1), "%"),
            #          x = maxx,
            #          y = plateau,
            #          hjust = 1,
            #          vjust = 1.5) +
            # annotate("text",
            #          alpha = 0.5,
            #          label = paste0("y = ", equation,
            #                         "\nAICc = ",AICc,
            #                         "\nRMSE = ", rmse,
            #                         "\nR-squared = ", rsquared),
            #          x = maxx,
            #          y = 0,
            #          vjust = 0,
            #          hjust = 1) +
            labs(
                x = "Soil test value (mg/kg)",
                y = "Relative yield (%)",
                caption = paste("Each point is a site. n =", nrow(corr_data))
            )
        
        return(avg_plot)
    
}

