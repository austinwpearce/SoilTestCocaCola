# This function uses SSasym, a self-starting function.
# These work better than guessing at starting values

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

# The "Mitscherlich" nls model with three parameters
# must be evaluated at Y less than asymptote for cstv
# ========================================================
# y = a + (a - b) * e^(-e^(c) * x)
# a = horizontal asymptote
# b = intercept when X is 0
# c = curvature, natural log of the rate constant (ought to be negative)

mit <- function(x, a, b, c) {
    a + (b - a) * exp(-exp(c) * x)
}

# Alternative is to use SSasym
fit_mit <- function(data){
    fit <- nlsLM(formula = y ~ mit(x, a, b, c),
                 data = data,
                 start = list(a = max(data$y),
                              b = min(data$y),
                              c = -1),
                 control = nls.lm.control(maxiter = 500),
                 upper = c(a = Inf, b = max(data$y), c = 0),
                 lower = c(a = min(data$y), b = -Inf, c = -Inf)
                 )
    return(fit)
}

## "Mitscherlich" Bundled Function with Plotting
# ========================================================
# only works if x = "stv" and y = "ry"
# returns either table OR plot
mitscherlich <- function(data,
                         force_origin = FALSE,
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
    corr_model <- try(nlsLM(
        formula = y ~ SSasymp(x, a, b, c),
        data = data,
        start = list(
            a = maxy,
            b = if_else(force_origin == TRUE, 0, miny),
            c = -0.5
        ),
        upper = c(
            a = Inf,
            b = if_else(force_origin == TRUE, 0, maxy),
            c = 0),
        lower = c(
            a = miny,
            b = if_else(force_origin == TRUE, 0, -Inf),
            c = -Inf
        )
    )
    )
    
    if (inherits(corr_model, "try-error")) {
        stop("Mitscherlich model could not be fit. Try something else.")
    } # else {
    #     corr_model <- corr_model
    # }
    
    # Find p-value and pseudo R-squared
    AIC <- round(AIC(corr_model), 0)
    rsquared <- round(modelr::rsquare(corr_model, data), 2)
    rmse <- round(modelr::rmse(corr_model, data), 2)
    
    # booted <- nlraa::boot_nls(corr_model, data = data)
    # confint(booted, type = "perc", level = 0.95)
    
    # get model coefficients
    a <- coef(corr_model)[[1]]
    b <- coef(corr_model)[[2]]
    c <- exp(coef(corr_model)[[3]])
    # new values
    ry_cstv <- a * percent_of_max / 100
    cstv <- log((ry_cstv - a) / (b - a)) / -c
    
    equation <- paste0(round(a, 1), " + (", round(b,1), " - ", round(a, 1),
                       ") * e^(-", round(c, 3), "x)")
    
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
            intercept = round(b, 2),
            rate_constant = round(c, 2),
            equation,
            percent_of_max,
            ry_cstv = round(ry_cstv, 0),
            cstv = round(cstv, 0),
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
                pred_band <- 
                # changed to predict_nls instead of predict2_nls
                nlraa::predict_nls(corr_model,
                                   newdata = data,
                                   interval = "confidence") %>%
                as_tibble()
        }
        
        pred_y <- tibble(x = seq(minx, maxx, 0.1)) %>%
            gather_predictions(corr_model)
        
        mit_plot <- data %>% 
            ggplot(aes(x, y)) +
            {
                if (band == TRUE)
                    geom_ribbon(data = pred_band,
                                aes(y = Estimate,
                                    ymin = Q2.5,
                                    ymax = Q97.5),
                                alpha = 0.05)
            } +
            # fitted line
            geom_line(data = pred_y, aes(x, pred),
                      color = red) +
            geom_vline(xintercept = cstv,
                       alpha = 1,
                       color = blue) +
            geom_hline(yintercept = ry_cstv,
                       alpha = 0.2) +
            geom_point(size = 3, alpha = 0.5) +
            geom_rug(alpha = 0.2, length = unit(2, "pt")) +
            scale_y_continuous(limits = c(0, maxy),
                               breaks = seq(0, maxy * 2, 10)) +
            annotate(
                "text",
                label = paste("CSTV =", round(cstv,1), "ppm"),
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
                               round(ry_cstv,1), "% RY"),
                x = maxx,
                y = ry_cstv,
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

