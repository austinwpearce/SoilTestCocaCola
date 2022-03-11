#' @title modified Arcsine-Log Calibration Curve, or ALCC
#' @author Austin Pearce
#' Last update: 2022-03-10
#' @references Correndo et al. 2017
#' @references Dyson and Conyers 2013
#' @description creates new variables on existing dataset
#' @name alcc_plot
#' @description perform ALCC and make plot without data table output
#' @param data a data frame with XY data
#' @param soil_test column for soil test values
#' @param ry column for relative yield values 0-100%
#' @param sma choose if Standardized Major Axis (SMA) is used for Modified ALCC (MALCC)
#' @param sufficiency choose at which RY value to get CSTV
#' @param confidence choose at which confidence level to estimate CI of CSTV
#' @param remove2x if TRUE, redo alcc() with data greater than twice the CSTV removed
# =============================================================================
# Could potentially add an argument that checks if data is percentage or ratio


# package libraries needed (won't just work in base R)
library(tidyverse) # a suite of packages for wrangling and plotting

# import alcc() function
source_url("https://raw.githubusercontent.com/austinwpearce/SoilTestCocaCola/main/alcc.R")


# =============================================================================
red <- "#CE1141"
blue <- "#13274F"

theme_set(
    theme_minimal(base_size = 14) +
        theme(
            plot.background = NULL,
            plot.margin = margin(t = 2, r = 10, b = 2, l = 2, unit = "pt"),
            panel.grid = element_line(color = "#F9F9F9"),
            panel.spacing = unit(2, "lines"),
            panel.border = element_blank(),
            #panel.grid.minor.y = element_blank(),
            axis.line = element_blank(),
            axis.ticks = element_blank(),
            axis.title.y = element_text(
                hjust = 1,
                margin = margin(t = 0, r = 10, b = 0, l = 0, unit = "pt")),
            axis.title.x = element_text(
                hjust = 0,
                margin = margin(t = 10, r = 0, b = 0, l = 0, unit = "pt")),
            axis.text = element_text(),
            legend.title.align = 0,
            legend.key.height = unit(x = 5, units = "mm"),
            legend.justification = c(1, 1)
            #legend.position = c(1, 1)
        ))


# =============================================================================
# ALCC
# =============================================================================

# Function alcc_plot() creates a scatter plot of soil test correlation data
# arguments are passed to alcc()
alcc_plot <- function(data,
                      soil_test,
                      ry,
                      sma = TRUE,
                      sufficiency = 90,
                      confidence = 95,
                      remove2x = FALSE) {
    
    if (missing(soil_test)) {
        stop("Specify the soil test variable (e.g., soil_test = STK)")
    }
    
    if (missing(ry)) {
        stop("Enter name of relative yield (%) variable (e.g., ry = RY)")
    }
    
    # enquo let's the user enter their own column names for x and y variables
    x <- enquo(soil_test)
    y <- enquo(ry)
    # pass to alcc function
    alcc_table <- alcc(data,
                       soil_test = !!x,
                       ry = !!y,
                       sma = sma,
                       sufficiency = sufficiency,
                       confidence = confidence)
    
    if (remove2x == TRUE) {
        cstv_2x <- (unique(alcc_table$cstv)) * 2
        # redo alcc() with data greater than twice the CSTV removed
        output <- alcc_table %>% 
            filter(stv <= cstv_2x) %>% 
            alcc(soil_test = stv,
                 ry = ry_cap,
                 sma = sma,
                 sufficiency = sufficiency,
                 confidence = confidence) %>% 
            mutate(remove2x = "TRUE")
    } else {
        output <- alcc_table %>% 
            mutate(remove2x = "FALSE")
    }
    
    # for plot annotations
    
    maxx <- max(alcc_table$stv)
    
    cstv <- (unique(output$cstv))
    cstv <- if_else(cstv < 10, round(cstv, 1), round(cstv, 0))
    
    lower_cl <- (unique(output$lower_cl))
    lower_cl <- if_else(lower_cl < 10, round(lower_cl, 1), round(lower_cl, 0))
    
    upper_cl <- (unique(output$upper_cl))
    upper_cl <- if_else(upper_cl < 10, round(upper_cl, 1), round(upper_cl, 0))
    
    #sufficiency <- unique(output$sufficiency)
    
    pearson <- round(unique(output$pearson), 2)
    
    # plot includes all data, even if remove2x == TRUE, but curve will follow remove2x
    alcc_table %>%
        ggplot(aes(stv, ry_cap)) +
        geom_point(size = 2, alpha = 0.5) +
        geom_vline(xintercept = lower_cl,
                   alpha = 0.5,
                   linetype = 3) +
        geom_vline(xintercept = upper_cl,
                   alpha = 0.5,
                   linetype = 3) +
        # twice the CSTV
        # geom_vline(xintercept = cstv * 2,
        #            alpha = 0.2,
        #            linetype = 2) +
        # fitted values from back-transformed regression line from alcc()
        geom_line(aes(x = fitted_stv,
                      y = fitted_ry),
                  color = red) +
        # the cstv at X% sufficiency
        geom_vline(xintercept = cstv, color = blue, alpha = 1) +
        geom_hline(yintercept = sufficiency, alpha = 0.2) +
        annotate(
            "text",
            label = paste0(
                "CSTV = ", cstv, " ppm at ", sufficiency, "% RY",
                "\n", confidence,"% CI [", lower_cl, " - ", upper_cl, "]",
                "\nPearson correlation = ", pearson),
            x = maxx,
            y = 0,
            hjust = 1,
            vjust = 0,
            alpha = 0.8
        ) +
        # scale_x_continuous(breaks =
        #                        seq(0, 1000, if_else(
        #                            max(output$stv) >= 200, 40,
        #                            if_else(max(output$stv) >= 100, 20,
        #                                    if_else(max(output$stv) >= 50, 10, 5))
        #                        ))) +
        scale_y_continuous(breaks = seq(0, 105, 10)) +
        labs(
            subtitle = if_else(sma == TRUE,
                            "Modified Arcsine-Log Calibration Curve with SMA",
                            "Arcsine-Log Calibration Curve"),
            x = expression(Soil ~ Test ~ Value ~ (mg ~ kg ^ -1)),
            y = "Relative yield (%)",
            caption = paste0(
                if_else(remove2x == TRUE,
                        "Data with soil test values more than twice the initial CSTV excluded from model.",
                        "All data included in model."),
                "\nVertical dotted lines show lower and upper confidence limits"
            )
        )
}
