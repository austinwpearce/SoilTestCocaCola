#' @title modified Arcsine-Log Calibration Curve, or ALCC
#' @author Austin Pearce
#' Last update: 2022-03-10
#' @references Correndo et al. 2017
#' @references Dyson and Conyers 2013
#' @description creates new variables on existing dataset
#' @name alcc_plot
#' @description perform ALCC-SMA and make plot without data table output
#' @param data a data frame with XY data
#' @param soil_test column for soil test values
#' @param ry column for relative yield values 0-100%
#' @param sufficiency choose at which RY value to get CSTV
#' @param confidence choose at which confidence level to estimate CI of CSTV
#' @param remove2x if TRUE, redo alcc_sma() with data greater than twice the CSTV removed
# =============================================================================
# Could potentially add an argument that checks if data is percentage or ratio


# package libraries needed (won't just work in base R)
library(tidyverse) # a suite of packages for wrangling and plotting

# import alcc() function
source_url("https://raw.githubusercontent.com/austinwpearce/SoilTestCocaCola/main/alcc_sma.R")


# =============================================================================
red <- "#CE1141"
gold <- "#EAAA00"
blue <- "#13274F"
black <- "#000000"

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
            axis.text = element_text()
            #legend.title.align = 0,
            #legend.key.height = unit(x = 5, units = "mm"),
            #legend.justification = c(1, 1),
            #legend.title = element_blank(),
            #legend.position = c(1, 0.5)
        ))


# =============================================================================
# ALCC
# =============================================================================

# Function alcc_plot() creates a scatter plot of soil test correlation data
# arguments are passed to alcc_sma()
alcc_plot <- function(data,
                      soil_test,
                      ry,
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
    
    input <- data %>% 
        transmute(stv = !!x,
                  # RY values greater than 100 are capped at 100
                  ry_cap = if_else(!!y > 100, 100, !!y))
    
    # pass to alcc_sma function
    output <- alcc_sma(data,
                       soil_test = !!x,
                       ry = !!y,
                       sufficiency = sufficiency,
                       confidence = confidence,
                       remove2x = remove2x)
    
    # for plot annotations
    
    maxx <- max(input$stv)
    
    cstv <- unique(output$cstv)
    cstv <- if_else(cstv < 10, round(cstv, 1), round(cstv, 0))
    cstv_100 <- unique(output$cstv_100)
    cstv90_2x <- unique(output$cstv90_2x)
    
    lower_cl <- (unique(output$lower_cl))
    lower_cl <- if_else(lower_cl < 10, round(lower_cl, 1), round(lower_cl, 0))
    
    upper_cl <- (unique(output$upper_cl))
    upper_cl <- if_else(upper_cl < 10, round(upper_cl, 1), round(upper_cl, 0))
    
    #sufficiency <- unique(output$sufficiency)
    
    pearson <- round(unique(output$pearson), 2)

    # plot includes all data, even if remove2x == TRUE, but curve will follow remove2x
    output %>%
        ggplot() +
        geom_point(data = input,
                   aes(stv, ry_cap),
                   size = 2, alpha = 0.7,
                   color = case_when(
                               remove2x == FALSE & input$stv > cstv_100 ~ red,
                               remove2x == FALSE & input$stv > cstv90_2x ~ gold,
                               remove2x == TRUE & input$stv > cstv90_2x ~ red,
                               TRUE ~ black)) +
        geom_vline(xintercept = lower_cl,
                   alpha = 0.5,
                   linetype = 3) +
        geom_vline(xintercept = upper_cl,
                   alpha = 0.5,
                   linetype = 3) +
        geom_point(aes(x = fitted_stv,
                      y = fitted_ry),
                  color = red) +
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
        scale_x_continuous(breaks =
                               seq(0, 10000, by = if_else(
                                   maxx >= 200, 40,
                                   if_else(maxx >= 100, 20,
                                           if_else(maxx >= 50, 10, 5))
                               ))) +
        scale_y_continuous(breaks = seq(0, 105, 10)) +
        labs(
            subtitle = "Modified Arcsine-Log Calibration Curve with SMA",
            x = expression(Soil ~ Test ~ Value ~ (mg ~ kg ^ -1)),
            y = "Relative yield (%)",
            caption = paste0(
                if_else(remove2x == TRUE,
                        "Red points > CSTV90 * 2 excluded from model",
                        "Red points > CSTV100 may be outliers, but were included. Yellow points > CSTV90 * 2 not excluded."),
                "\nVertical dotted lines show lower and upper confidence limits"
            )
        )
}

# Attempt to color code the points in the model
# input <- input %>% 
#     mutate(removed = case_when(
#         remove2x == FALSE & input$stv > cstv_100 ~ "Excluded > CSTV100",
#         remove2x == FALSE & input$stv > cstv90_2x ~ "Included > CSTV90*2",
#         remove2x == TRUE & input$stv > cstv90_2x ~ "Excluded > CSTV90*2",
#         TRUE ~ "Included < CSTV90*2"))
# { 
#     if(remove2x == TRUE)
#         scale_color_manual(values = c(red, black)) 
#     else scale_color_manual(values = c(red, black, gold))
# } +
