library(agridat)
library(rcompanion)
library(dplyr)

df <- data.frame(stv = agridat::cate.potassium$potassium,
                 ry = agridat::cate.potassium$yield)

df$ry <- ifelse(df$ry > 100, 100, df$ry)

df <- filter(df, stv < 140)

sufficiency <- 95

cn_fixed <- rcompanion::cateNelsonFixedY(
    x = df$stv,
    y = df$ry,
    cly = sufficiency,
    trend = "positive",
    plotit = TRUE,
    xlab = "STV",
    ylab = "RY")

cstv <- cn_fixed$Critx[[1]]

cstv

df <- df[order(df$stv), ] # Sort the data by x
x <- df$stv
y <- df$ry

# Create a data.frame to store the results

out <- data.frame(
    x = NA,
    mean1 = NA,
    css1 = NA,
    mean2 = NA,
    css2 = NA,
    r2 = NA
)

css <- function(x) {
    var(x) * (length(x) - 1)
}
tcss <- css(y) # Total corrected sum of squares

for (i in 2:(length(y) - 2)) {
    y1 <- y[1:i]
    y2 <- y[-(1:i)]
    out[i, 'x'] <- x[i]
    out[i, 'mean1'] <- mean(y1)
    out[i, 'mean2'] <- mean(y2)
    out[i, 'css1'] <- css1 <- css(y1)
    out[i, 'css2'] <- css2 <- css(y2)
    out[i, 'r2'] <- (tcss - (css1 + css2)) / tcss
    
}

cn_rsquared <- as_tibble(out) %>%
    arrange(desc(r2)) %>% 
    filter(x == cstv)

cn_rsquared$r2

# alternatively, though the CSTV value is slightly different
soiltestcorr::cate_nelson_1965(df, stv, ry, 95)$R2
