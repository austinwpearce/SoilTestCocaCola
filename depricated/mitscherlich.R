library(tidyverse) # yes this loads lots of packages but these should be installed anyways
library(modelr)
library(minpack.lm) # for nlsLM which is more robust
library(nlraa) # for self-starting functions


# Update has everything in terms of x and y
# So prior to using this function, the columns STV and RY, for example, will
# need to be renamed to x and y
# for example: %>% rename(x = stv, y = ry) %>% 

red <- "#CE1141"
gold <- "#EAAA00"
blue <- "#13274F"
black <- "#000000"

# The MB NLS Model with three parameters (asymptote, "intercept", curvatuve)
# must be evaluated at Y less than asymptote for cstv
# ========================================================
# y = a - b + exp(-c*x)
# a = asymptote
# b = difference to intercept, where a - b = intercept I believe
# c = curvature (ought to be negative)

mb <- function(x, a, b, c) {
    a - b * exp(c * x)
}


fit_mb <- function(data){
    start <- lm(y ~ poly(x, 2, raw = TRUE), data = data)
    fit <- nls(formula = y ~ mb(x, a, b, c),
               data = data,
               start = list(a = max(data$y),
                            b = mean(data$y),
                            c = start$coef[[3]]),
                 control = nls.lm.control(maxiter = 500))
    return(fit)
}

fit_mbLM <- function(data){
    start <- lm(y ~ poly(x, 2, raw = TRUE), data = data)
    fit <- nlsLM(formula = y ~ mb(x, a, b, c),
                 data = data,
                 start = list(a = max(data$y),
                              b = mean(data$y),
                              c = start$coef[[3]]),
                 control = nls.lm.control(maxiter = 500)
                 #upper = c(asym = Inf, b = Inf, c = Inf),
                 #lower = c(asym = min(data$y), b = -Inf, c = -Inf),
                 )
    return(fit)
}

## Bundled Function with Plotting
# ========================================================
# only works if x = "stv" and y = "ry"
# returns either table OR plot
mitscherlich <- function(data,
                         percent_of_max = 95,
                         resid = FALSE,
                         plot = FALSE,
                         band = FALSE) {
    
    if (nrow(data) < 4) {
        stop("Too few distinct input values to fit QP. Try at least 4.")
    }
    
    if ("x" %in% colnames(data) == FALSE) {
        stop("Rename the explanatory variable 'x'")
    }
    
    if ("y" %in% colnames(data) == FALSE) {
        stop("Rename the response variable 'y'")
    }
    
    minx <- min(data$x)
    maxx <- max(data$x)
    miny <- min(data$y)
    maxy <- max(data$y)
    
    # build the model/fit =====
    nls_model <- try(fit_mb(data))
    
    if (inherits(nls_model, "try-error")) {
        corr_model <-
            try(fit_mbLM(data))
    } else {
        corr_model <- nls_model
    }
    
    if (inherits(corr_model, "try-error")) {
        stop("Mitscherlich model could not be fit with nls or nlsLM.
             Try something else.")
    } else {
        corr_model <- corr_model
    }
    
    # Find p-value and pseudo R-squared
    AIC <- round(AIC(corr_model), 0)
    rsquared <- round(modelr::rsquare(corr_model, data), 2)
    rmse <- round(modelr::rmse(corr_model, data), 2)
    
    # booted <- nlraa::boot_nls(corr_model, data = data)
    # confint(booted, type = "perc", level = 0.95)
    
    # get model coefficients
    a <- coef(corr_model)[[1]]
    b <- coef(corr_model)[[2]]
    c <- coef(corr_model)[[3]]
    new_asym <- a * percent_of_max / 100
    cstv <- log((a - new_asym) / b) / c
    
    # cstv_90ry <- if_else(condition = asym >= 90, 
    #                      true = round(log((asym - 90) / b) / c, 0),
    #                      false = NULL
    # )
    
    equation <- paste0(round(a, 1), " - ",
                       round(b, 2), "e^",
                       round(c, 3), "x")
    
    # # get error for each parameter
    # se_a <- round(tidy(corr_model)$std.error[1], 2)
    # se_b <- round(tidy(corr_model)$std.error[2], 2)
    # se_c <- round(tidy(corr_model)$std.error[3], 3)
    
    # Table output =================================================
    if (plot == FALSE) {
        {
            if (resid == TRUE)
                plot(nlstools::nlsResiduals(corr_model), which = 0)
        }
        tibble(
            asymptote = round(a, 1),
            b = round(b, 2),
            c = round(c, 3),
            equation,
            percent_of_max,
            #new_ry = round(new_asym, 0),
            cstv = round(cstv, 0),
            # cstv_90ry = if_else(cstv_90ry > 0,
            #                     cstv_90ry,
            #                     NULL), # at 90% RY
            AIC,
            rmse,
            rsquared
        )
    } else {
        # Residual plots and normality
        {
            if (resid == TRUE)
                plot(nlstools::nlsResiduals(corr_model), which = 0)
        }
        
        {
            if (band == TRUE)
                predicted <- 
                # changed to predict_nls instead of predict2_nls
                nlraa::predict_nls(corr_model,
                                   newdata = data,
                                   interval = "confidence") %>%
                as_tibble() %>%
                bind_cols(data)
        }
        
        {
            if (band == FALSE)
                predicted <- data
        }
        
        mit_plot <- predicted %>%
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
            geom_vline(xintercept = cstv,
                       alpha = 1,
                       color = blue) +
            geom_hline(yintercept = new_asym,
                       alpha = 0.2) +
            geom_line(
                stat="smooth",
                method = "nls",
                formula = y ~ mb(x, a, b, c),
                method.args = list(start = as.list(coef(corr_model))),
                se = FALSE,
                color = "#CC0000"
            ) +
            geom_point(size = 3, alpha = 0.5) +
            geom_rug(alpha = 0.2, length = unit(2, "pt")) +
            scale_y_continuous(limits = c(0, maxy),
                               breaks = seq(0, maxy * 2, 10)) +
            annotate(
                "text",
                label = paste("CSTV =", round(cstv,0), "ppm"),
                x = cstv,
                y = 0,
                angle = 90,
                hjust = 0,
                vjust = 1.5,
                alpha = 0.5
            ) +
            annotate(
                "text",
                alpha = 0.5,
                label = paste0(percent_of_max, "% of asymptote = ",
                               round(new_asym,1), "% RY"),
                x = maxx,
                y = new_asym,
                vjust = 1.5,
                hjust = 1
            )+
            annotate(
                "text",
                alpha = 0.5,
                label = paste0(
                    "y = ", equation,
                    "\nAIC = ", AIC,
                    "\nRMSE = ", rmse,
                    "\nR-squared = ", rsquared
                    #"\nCSTV at 95 and 90% RY = ", cstv_95ry," & ", cstv_90ry, " ppm"
                    ),
                x = maxx,
                y = 0,
                vjust = 0,
                hjust = 1
            ) +
            labs(x = "Soil test value (mg/kg)",
                 y = "Relative yield (%)",
                 caption = paste("Each point is a site. n =", nrow(data)))
        
        return(mit_plot)
    }
    
}

