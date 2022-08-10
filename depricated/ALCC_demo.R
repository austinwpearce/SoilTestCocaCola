library(ggplot2)
library(agridat)
# ALCC demo

df <- data.frame(stv = agridat::cate.potassium$potassium,
                 ry = agridat::cate.potassium$yield)

qplot(data = df, stv, ry, geom = "point")

# limit to 100 for ALCC
df$ry <- ifelse(df$ry > 100, 100, df$ry)

qplot(data = df, stv, ry, geom = "point")

# get sample size
n <- nrow(df)
# Step 1 Transform (t for "transformed" added to x and y)
# for ALCC, soil test goes to Y-axis and RY goes to X-axis
xt <- asin(sqrt(df$ry / 100))
yt <- log(df$stv)
# Step 2 Center
sufficiency <- 95
adjust_by <- asin(sqrt(sufficiency / 100))
xt_centered <- xt - adjust_by
# Step 3 Correlation
pearson <- cor(xt_centered, yt, method = "pearson")
t_stat  <- (pearson * sqrt(n - 2)) / sqrt(1 - pearson ^ 2)
pvalue  <- pt(t_stat, df = n - 1, lower.tail = FALSE)
# Step 4 Means
mean_xt  <- mean(xt_centered)
mean_yt  <- mean(yt)

# Step 5 Fit linear model to transformed, centered data
fit_lm <- lm(yt ~ xt_centered)

intercept <- coef(fit_lm)[[1]]
slope     <- coef(fit_lm)[[2]]
# Step 6 Rotate the regression (SMA)
# slope must come first
slope     <- slope / pearson
intercept <- mean_yt - slope * mean_xt
# Step 7 Estimate Critical Soil Test Concentration
cstv <- exp(intercept)
cstv
# Step 8 Estimate the confidence interval
pred_yt  <- intercept + slope * xt_centered
mse      <- sum((yt - pred_yt) ^ 2) / (n - 2)
ssx      <- var(xt_centered) * (n - 1)
confidence <- 95
se       <- sqrt(mse * ((1 / n) + ((mean_xt ^ 2) / ssx)))
lower_cl <-
    exp(intercept - se * qt(1 - (1 - confidence / 100) / 2,
                            df = n - 2))
upper_cl <- exp(intercept + se * qt(1 - (1 - confidence / 100) / 2,
                                  df = n - 2))
# Step 9 Back-transform
# New RY values to create smoother curve to 0
tmp <- seq(0, 100, by = 0.2)
tmp2 <- asin(sqrt(tmp/100)) - asin(sqrt(sufficiency/100))
new_pred_y <- intercept + slope * tmp2
fitted_stv <- exp(new_pred_y)
fitted_ry <- 100 * (sin(adjust_by + ((new_pred_y - intercept) / slope))) ^ 2

ggplot() +
    geom_point(aes(df$stv, df$ry)) +
    geom_line(aes(fitted_stv, fitted_ry)) +
    geom_hline(yintercept = sufficiency, lty = 2) +
    geom_vline(xintercept = c(lower_cl, upper_cl), lty = 3) +
    geom_vline(xintercept = cstv, color = "#CC0000") +
    geom_point(aes(cstv, sufficiency), color = "#CC0000", size = 3)
