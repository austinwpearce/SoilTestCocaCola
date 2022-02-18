#' This function fits a quadratic plateau model to soil test correlation XY data
#' and either provides the results in a table format or as a plot
#'
#' Last updated: 2022-02-17
#'
#' This function is essentially a wrapper that uses other packages' functions
#' Won't just work in base R
#'
#' @name quad_plateau
#' @param data a data frame with XY data
#' @param resid choose whether to create residuals plots
#' @param plot choose whether to create correlation plot rather than table
#' @param band choose whether the correlation plot displays confidence band
#' no effect if plot = FALSE
#' @param percent_of_max if wanting to find the X value for a point along the
#' quadratic portion at certain Y value
#' @export

#' The QP model and parameters
#' y = b0 + b1x + b2x^2
#' b0 = intercept
#' b1 = slope
#' b2 = quadratic term
#' jp = join point = critical concentration = -0.5 * b1 / b2

# =============================================================================
# package libraries needed
library(tidyverse) # a suite of packages for wrangling and plotting
library(nlraa) # for self-starting functions and predicted intervals
library(minpack.lm) # for nlsLM, a robust backup to nls
library(nlstools) # for residuals plots
library(modelr) # for the r-squared and rmse

# For now, everything is written in terms of x and y
# So prior to using this function, the columns for soil test value and RY, for # example, will need to be renamed to x and y
# for example: %>% rename(x = stv, y = ry) %>%

# =============================================================================
# quad_plateau function

quad_plateau <- function(data,
                         resid = FALSE,
                         plot = FALSE,
                         band = FALSE,
                         percent_of_max = 95) {
    minx <- min(data$x)
    maxx <- max(data$x)
    miny <- min(data$y)
    maxy <- max(data$y)
    
    # uses the fit_qp function defined earlier
    corr_model <-
        try(nls(y ~ SSquadp3xs(x, b0, b1, jp), data = data))
    
    if (inherits(corr_model, "try-error")) {
        corr_model <-
            try(minpack.lm::nlsLM(y ~ SSquadp3xs(x, b0, b1, jp), data = data))
    } else {
        corr_model <- corr_model
    }
    
    # How did the model do overall?
    AIC      <- round(AIC(corr_model), 0)
    rsquared <- round(modelr::rsquare(corr_model, data), 2)
    rmse     <- round(modelr::rmse(corr_model, data), 2)
    
    # get model coefficients
    b0 <- coef(corr_model)[[1]]
    b1 <- coef(corr_model)[[2]]
    # join point of segmented regression = critical soil test value
    jp <- coef(corr_model)[[3]]
    join_point <- round(jp, 0)
    b2 <- -0.5 * b1 / jp
    
    plateau <- round(b0 + (b1 * jp) + (b2 * jp * jp), 1)
    
    equation <- paste0(round(b0, 1), " + ",
                       round(b1, 2), "x + ",
                       round(b2, 3), "x^2")
    
    # CSTV at defined % of max/plateau
    adjust_cstv <- function(b0, b1, b2) {
        jp <- -0.5 * b1 / b2
        newplateau <-
            (b0 + (b1 * jp) + (b2 * jp * jp)) * percent_of_max / 100
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
    
    cstv_adj <- adjust_cstv(b0, b1, b2)
    
    # Printouts
    if (plot == FALSE) {
        {
            if (resid == TRUE)
                plot(nlsResiduals(corr_model), which = 0)
        }
        tibble(
            intercept = b0,
            slope = b1,
            curve = b2,
            equation,
            join_point,
            plateau = round(plateau, 1),
            AIC,
            rmse,
            rsquared,
            adjusted_cstv = round(cstv_adj, 0),
            percent_of_max
        )
    } else {
        # Residual plots and normality
        {
            if (resid == TRUE)
                plot(nlsResiduals(corr_model), which = 0)
        }
        predicted <- nlraa::predict2_nls(corr_model,
                                         newdata = data,
                                         interval = "confidence") %>%
            as_tibble() %>%
            bind_cols(data)
        
        # ggplot of data
        qp_plot <- predicted %>%
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
            geom_line(
                stat = "smooth",
                method = "nls",
                formula = y ~ SSquadp3xs(x, b0, b1, jp),
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

qp <- function(x, b0, b1, jp) {
    b2 <- -0.5 * b1 / jp
    if_else(
        condition = x < jp,
        true  = b0 + (b1 * x) + (b2 * x * x),
        false = b0 + (b1 * jp) + (b2 * jp * jp)
    )
}

fit_qp <- function(data) {
    start <- lm(y ~ poly(x, 2, raw = TRUE), data = data)
    # nls model
    fit <- nlsLM(
        formula = y ~ qp(x, b0, b1, jp),
        data = data,
        start = list(
            b0 = start$coef[[1]],
            b1 = start$coef[[2]],
            jp = mean(data$x)
        ),
        control = nls.lm.control(maxiter = 500)
    )
    return(fit)
}

fit_SSquadp3xs <- function(data) {
    # nlraa model
    fit <- nlsLM(
        formula = y ~ SSquadp3xs(x, b0, b1, jp),
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