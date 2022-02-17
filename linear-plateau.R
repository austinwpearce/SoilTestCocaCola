# Last updated 2021-11-03

library(tidyverse)
library(broom)
library(modelr)
library(minpack.lm)
library(nlstools)
library(nlraa)
# Update has everything in terms of x and y
# So prior to using this function, the columns STV and RY, for example, will
# need to be renamed to x and y
# for example: %>% rename(x = stv, y = ry) %>% 


# Linear-plus-plateau model
# y = b0 + b1x
# intercept (b0)
# slope (b1)
# jp = join point = critical concentration

# fitting function that finds starting values based on linear regression
# this results in fast results, much faster than nls_multstart or SSlinp
# but the trade-off is that on bootstraps the results will be all over
# the place so the fit_lp/lin_plateau functions work fine if you just need a
# quick model.
# nlraa::SSlinp is actually pretty good, using a fancy algorithm and all

fit_SSlinp <- function(data) {
    # nlraa model
    fit <- nlsLM(
        formula = y ~ SSlinp(x, b0, b1, jp),
        data = data,
        control = nls.lm.control(maxiter = 500),
        upper = c(b0 = max(data$y), b1 = 10000, jp = max(data$x)),
        lower = c(b0 = -10000, b1 = 0, jp = min(data$x)))
    
    return(fit)
}

# Linear-plateau ==========================================================
lin_plateau <- function(data,
                        confidence = 95,
                        alpha = 0.05,
                        percent_of_max = 95,
                        plot = FALSE,
                        resid = FALSE,
                        band = FALSE) {
    minx <- min(data$x)
    maxx <- max(data$x)
    miny <- min(data$y)
    maxy <- max(data$y)
    
    # build the model/fit =====
    corr_model <<- fit_SSlinp(data)
    
    # Find p-value and pseudo R-squared
    rsquared <- round(modelr::rsquare(corr_model, data), 2)
    rmse <- round(modelr::rmse(corr_model, data), 2)
    aic <- round(AIC(corr_model), 0)
    
    # WALD IS PROBLEMATIC ACCORDING TO FERNANDEZ AND MULLA AND MIGUEZ
    # # get wald confidence interval of join point at 100% 
    # ci <- nlstools::confint2(corr_model, level = confidence / 100)
    # lower_cl <- round(ci[3, 1], 0)
    # upper_cl <- round(ci[3, 2], 0)
    
    # booted <- nlraa::boot_nls(corr_model, data = data)
    # confint(booted, type = "perc", level = 0.95)
    
    # get model coefficients
    b0 <- tidy(corr_model)$estimate[1]
    b1 <- tidy(corr_model)$estimate[2]
    jp <- tidy(corr_model)$estimate[3]
    
    # get error for each paramater
    se_b0 <- round(tidy(corr_model)$std.error[1], 2)
    se_b1 <- round(tidy(corr_model)$std.error[2], 2)
    
    # X% Sufficiency ====
    plateau <- b0 + b1 * jp
    
    cstv_adj <- round(((plateau * percent_of_max / 100) - b0) / b1, 0)
    
    cstv_95ry <- if_else(condition = plateau >= 95, 
                    true = round((95 - b0) / b1, 0),
                    false = NULL
                    )
    
    cstv_90ry <- if_else(condition = plateau >= 90, 
                       true = round((90 - b0) / b1, 0),
                       false = NULL
    )
    
    # have to make a line because the SSlinp doesn't plot right
    lp_line <- tibble(x = c(minx, jp, maxx),
                      y = c(b0 + b1 * minx, plateau, plateau))
    
    equation <- paste0(round(b0, 1), " + ",
                       round(b1, 2), "x")
    
    join_point <- round(jp, 0)
    
    # Table output =====
    if (plot == FALSE) {
        tibble(
            intercept = b0,
            slope = b1,
            equation,
            join_point,
            plateau = round(plateau, 1),
            cstv_adj,
            percent_of_max,
            cstv_95ry = if_else(cstv_95ry > 0,
                                cstv_95ry,
                                NULL), # at 95% RY
            cstv_90ry = if_else(cstv_90ry > 0,
                                cstv_90ry,
                                NULL), # at 90% RY
            AIC = aic,
            rsquared,
            rmse,
            se_b0,
            se_b1
        )
    } else {
        # Residual plots and normality
        {if (resid == TRUE) 
            plot(nlstools::nlsResiduals(corr_model), which = 0)}
        predicted <- nlraa::predict2_nls(corr_model,
                                 newdata = data,
                                 interval = "confidence") %>%
            as_tibble() %>% 
            bind_cols(data)
        
        predicted %>%
            plottr_d(x, y, size = 3) +
            {if (band == TRUE) 
                geom_ribbon(aes(y = Estimate,
                                ymin = Q2.5,
                                ymax = Q97.5),
                            alpha = 0.05)} +
            ggpp::geom_vhlines(xintercept = jp,
                               yintercept = plateau,
                               alpha = 0.5,
                               linetype = 3) +
            annotate(
                "text",
                label = paste("CSTV =", join_point, "ppm"),
                x = join_point,
                y = miny,
                angle = 90,
                family = sanserif,
                hjust = 0,
                vjust = 1.5,
                alpha = 0.5
                ) +
            annotate(
                "text",
                label = paste0("Plateau = ", round(plateau,1), "%"),
                x = maxx,
                y = plateau,
                alpha = 0.5,
                family = sanserif,
                hjust = 1,
                vjust = 1.5
            ) +
            annotate(
                "text",
                alpha = 0.5,
                label = paste0("y = ", equation,
                    "\nAIC = ", aic, "  RMSE = ", rmse,
                    "\nR-squared = ", rsquared,
                    "\nCSTV at 95 and 90% RY = ", cstv_95ry," & ", cstv_90ry, " ppm"),
                x = maxx * 0.7,
                y = quantile(data$y, 0.05),
                hjust = 0,
                family = sanserif) +
            # annotate("rect",
            #          xmin = lower_cl,
            #          xmax = upper_cl,
            #          ymin = -Inf,
            #          ymax = Inf,
            #          alpha = 0.05) +
            geom_line(data = lp_line,
                      aes(x = x, y = y),
                      color = "#CC0000") +
            labs(x = label_stv,
                 y = label_ry0,
                 caption = paste(caption_site, nrow(data)))
    }
    
}

