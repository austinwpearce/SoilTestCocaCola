library(tidyverse) # yes this loads lots of packages but these should be loaded anyways
library(broom)  # for tidying models into dataframes
library(modelr) # for easy rsquare calculation
library(minpack.lm) # for nlsLM which is more robust and fails less than 
library(nlstools)
library(nlraa) # for self-starting functions


# Update has everything in terms of x and y
# So prior to using this function, the columns STV and RY, for example, will
# need to be renamed to x and y
# for example: %>% rename(x = stv, y = ry) %>% 

# I'm not sure how to calculate the significance of the B2 term without including it in the model. But I would also like to see the confidence interval for the join point at 100%. So I'll use the nlraa SSquadp3 for a B2 model, but then I'll run either my fit_qp bootstrap resamples for confidence intervals. For LP models I can get them from confint2() but not with QP models.

# Also you cannot make a 4 parameter version without the fit looking like a quad + plateau where the plateau is not joined to the quad max.

# The QP NLS Model with three parameters (intercept, slope, join point)
# ========================================================
# y = b0 + b1x + b2x^2
# b0 = intercept
# b1 = slope
# b2 = quadratic term
# jp = join point = critical concentration = -0.5 * b1 / b2

qp <- function(x, b0, b1, jp) {
    b2 <- -0.5 * b1 / jp
    if_else(condition = x < jp,
            true  = b0 + (b1 * x) + (b2 * x * x),
            false = b0 + (b1 * jp) + (b2 * jp * jp))
}

fit_qp <- function(data) {
    start <- lm(y ~ poly(x, 2, raw = TRUE), data = data)
    # nls model
    fit <- nlsLM(formula = y ~ qp(x, b0, b1, jp),
                 data = data,
                 start = list(b0 = start$coef[[1]],
                              b1 = start$coef[[2]],
                              jp = mean(data$x)),
                 control = nls.lm.control(maxiter = 500)
                 )
    return(fit)
}

fit_SSquadp3xs <- function(data) {
    # nlraa model
    fit <- nlsLM(formula = y ~ SSquadp3xs(x, b0, b1, jp),
                 data = data,
                 control = nls.lm.control(maxiter = 500)
                 # upper = c(b0 = max(data$y),
                 #           b1 = 1000,
                 #           b2 = 0),
                 # lower = c(b0 = -1000,
                 #           b1 = 0,
                 #           b2 = -Inf)
                 )
    return(fit)
}

## Quadratic-plateau Bundled Function with Plotting
# ========================================================
# only works if x = "stv" and y = "ry"
# returns either table OR plot
quad_plateau <- function(data,
                         alpha = 0.05,
                         confidence = 95,
                         percent_of_max = 95,
                         plot = FALSE,
                         band = FALSE) {
    minx <- min(data$x)
    maxx <- max(data$x)
    miny <- min(data$y)
    maxy <- max(data$y)
    
    # now get the main model with b2
    corr_model <- fit_SSquadp3xs(data)
    
    # Find p-value and pseudo R-squared
    aic <- round(AIC(corr_model),0)
    rmse <- round(modelr::rmse(corr_model, data), 2)
    rsquared <- round(modelr::rsquare(corr_model, data), 2)
    
    # get model coefficients
    b0 <- tidy(corr_model)$estimate[1]
    b1 <- tidy(corr_model)$estimate[2]
    jp <- tidy(corr_model)$estimate[3]
    b2 <- -0.5 * b1 / jp
    join_point <- round(jp, 0)
    
    # get error for each parameter
    se_b0 <- round(tidy(corr_model)$std.error[1], 2)
    se_b1 <- round(tidy(corr_model)$std.error[2], 2)
    se_jp <- round(tidy(corr_model)$std.error[3], 3)
    
    # Plateau and equation
    plateau <- b0 + (b1 * jp) + (b2 * jp * jp)
    
    equation <- paste0(round(b0, 1), " + ",
                       round(b1, 2), "x + ",
                       round(b2, 3), "x^2")
    
    # CSTV at defined % of max/plateau
    adjust_cstv <- function(b0, b1, b2) {
        jp <- -0.5 * b1 / b2
        newplateau <- (b0 + (b1 * jp) + (b2 * jp * jp)) * percent_of_max / 100
        newb0 <- b0 - newplateau
        discriminant <- (b1 ^ 2) - (4 * newb0 * b2)
        
        if (discriminant < 0) {
            return(NA)
        }
        else if (discriminant > 0) {
            cstv_adj <- (-b1 + sqrt(discriminant)) / (2 * b2)
            
            return(cstv_adj)
        }
    }
    
    cstv_adj <- round(adjust_cstv(b0, b1, b2), 0)
    
    # Printouts
    if (plot == FALSE) {
        tibble(
            intercept = b0,
            slope = b1,
            curve = b2,
            equation,
            join_point,
            plateau = round(plateau, 1),
            cstv_adj,
            percent_of_max,
            AIC = aic,
            rsquared,
            rmse,
            se_b0,
            se_b1,
            se_jp
        )
    } else {
        # Residual plots and normality
        # nls_resid <- nlstools::nlsResiduals(corr_model)
        # plot(nls_resid, which = 0)
        predicted <- nlraa::predict_nls(corr_model,
                                 newdata = data,
                                 interval = "confidence") %>%
            as_tibble() %>%
            bind_cols(data)

        # ggplot of data
        predicted %>%
            plottr_d(x, y, size = 3) +
            {if (band == TRUE) 
                geom_ribbon(aes(y = Estimate,
                                ymin = Q2.5,
                                ymax = Q97.5),
                            alpha = 0.05)}
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
            geom_line(
                stat="smooth",
                method = "nlsLM",
                formula = y ~ SSquadp3xs(x, b0, b1, jp),
                method.args = list(start = as.list(coef(corr_model))),
                se = FALSE,
                color = "#CC0000"
            ) +
            annotate("text",
                     alpha = 0.5,
                label = paste0("y = ", equation,
                               "\nAIC = ", aic, "  RMSE = ", rmse,
                               "\nR-squared = ", rsquared),
                x = maxx * 0.7,
                y = quantile(data$y, 0.05),
                hjust = 0,
                family = sanserif
            ) +
            labs(x = label_stv,
                 y = label_ry0,
                 caption = paste(caption_site, nrow(data)))
    }
    
}

