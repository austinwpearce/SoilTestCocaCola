library(tidyverse) # yes this loads lots of packages but these should be installed anyways
library(broom)
library(modelr)
library(minpack.lm) # for nlsLM which is more robust
library(nlraa) # for self-starting functions


# Update has everything in terms of x and y
# So prior to using this function, the columns STV and RY, for example, will
# need to be renamed to x and y
# for example: %>% rename(x = stv, y = ry) %>% 

# The MB NLS Model with three parameters (asymptote, "intercept", curvatuve)
# must be evaluated at Y less than asymptote for cstv
# ========================================================
# y = a - b + exp(-c*x)
# a = asymptote
# b = difference to intercept, where a - b = intercept I believe
# c = curvature (ought to be negative)

mb <- function(x, asym, b, c) {
    asym - b * exp(c * x)
}


fit_mb <- function(data, ...){
    start <- lm(y ~ poly(x, 2, raw = TRUE), data = data)
    fit <- nlsLM(formula = y ~ mb(x, asym, b, c),
                 data = data,
                 start = list(asym = mean(data$y),
                              b = start$coef[[2]],
                              c = start$coef[[3]]),
                 control = nls.lm.control(maxiter = 500),
                 #upper = c(asym = Inf, b = Inf, c = Inf),
                 #lower = c(asym = min(data$y), b = -Inf, c = -Inf),
                 ...)
    return(fit)
}

## Quadratic-plateau Bundled Function with Plotting
# ========================================================
# only works if x = "stv" and y = "ry"
# returns either table OR plot
mitscherlich <- function(data,
                   #confidence = 95,
                   alpha = 0.05,
                   percent_of_max = 95,
                   plot = FALSE) {
    
    minx <- min(data$x)
    maxx <- max(data$x)
    miny <- min(data$y)
    maxy <- max(data$y)
    
    # now get the main model with b2
    corr_model <- fit_mb(data)
    
    # Find p-value and pseudo R-squared
    aic <- round(AIC(corr_model), 0)
    rsquared <- round(modelr::rsquare(corr_model, data), 2)
    rmse <- round(modelr::rmse(corr_model, data), 2)
    
    # get model coefficients
    asym <- tidy(corr_model)$estimate[1]
    b <- tidy(corr_model)$estimate[2]
    c <- tidy(corr_model)$estimate[3]
    new_asym <- asym * percent_of_max / 100
    cstv_adj <- log((asym - new_asym) / b) / c
    
    # cstv_95ry <- if_else(condition = asym >= 95, 
    #                      true = round(log((asym - 95) / b) / c, 0),
    #                      false = NULL
    # )
    # 
    # cstv_90ry <- if_else(condition = asym >= 90, 
    #                      true = round(log((asym - 90) / b) / c, 0),
    #                      false = NULL
    # )
    
    equation <- paste0(round(asym, 1), " - ",
                       round(b, 2), "e^",
                       round(c, 3), "x")
    
    # get error for each parameter
    se_a <- round(tidy(corr_model)$std.error[1], 2)
    se_b <- round(tidy(corr_model)$std.error[2], 2)
    se_c <- round(tidy(corr_model)$std.error[3], 3)
    
    # Printouts
    if (plot == FALSE) {
        tibble(
            equation,
            asymptote = round(asym, 0),
            percent_of_max,
            #new_ry = round(new_asym, 0),
            cstv_adj = round(cstv_adj, 0),
            # cstv_95ry = if_else(cstv_95ry > 0,
            #                     cstv_95ry,
            #                     NULL), # at 95% RY
            # cstv_90ry = if_else(cstv_90ry > 0,
            #                     cstv_90ry,
            #                     NULL), # at 90% RY
            AIC = aic,
            rsquared,
            rmse,
            model_error
        )
    } else {
        # Residual plots and normality
        #nls_resid <- nlstools::nlsResiduals(corr_model)
        #plot(nls_resid, which = 0)
        predicted <- nlraa::predict_nls(corr_model,
                                 newdata = data,
                                 interval = "confidence") %>%
            as_tibble() %>%
            bind_cols(data)
        
        # ggplot of data
        predicted %>%
            plottr_d(x, y) +
            # if you want to show the confidence bands
            # geom_ribbon(aes(y = Estimate,
            #                 ymin = Q2.5,
            #                 ymax = Q97.5),
            #             alpha = 0.2) +
            ggpp::geom_vhlines(xintercept = cstv_adj,
                               yintercept = new_asym,
                               linetype = 3,
                               alpha = 0.5) +
            geom_line(
                stat="smooth",
                method = "nlsLM",
                formula = y ~ mb(x, asym, b, c),
                method.args = list(start = as.list(coef(corr_model))),
                se = FALSE,
                color = "#CC0000"
            ) +
            annotate(
                "text",
                label = paste("CSTV =", round(cstv_adj,0), "ppm"),
                x = cstv_adj,
                y = miny,
                angle = 90,
                family = rc,
                hjust = 0,
                vjust = 1.5,
                alpha = 0.5
            ) +
            annotate(
                "text",
                label = paste0(percent_of_max, "% of asymptote = ", round(new_asym,1), "% RY"),
                x = maxx,
                y = new_asym,
                alpha = 0.5,
                family = rc,
                hjust = 1,
                vjust = 1.5
            )+
            annotate(
                "text",
                label = paste0("y = ", equation,
                               "\nAIC = ", aic, "  RMSE = ", rmse,
                              "\nR-squared = ", rsquared
                              #"\nCSTV at 95 and 90% RY = ", cstv_95ry," & ", cstv_90ry, " ppm"
                              ),
                x = max(data$x) * 0.7,
                y = min(data$y) + 5,
                hjust = 0,
                family = rc
            ) +
            labs(x = label_stv,
                 y = label_ry0,
                 caption = paste(caption_site, nrow(data)))
    }
    
}
