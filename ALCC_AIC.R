library(tidyverse)
library(AICcmodavg)
# ALCC demo

df <- data.frame(
    stv = c(179, 52, 43, 17, 16, 22, 83, 38, 28, 28, 50, 94, 54, 196, 74, 83, 43,
          40, 90, 98, 76, 65, 56, 76, 108, 27, 26, 60, 174, 177, 92, 146, 71, 81,
          80, 66, 35, 30, 51, 50, 51), 
    ry = c(92.39102, 97.34460, 98.83400, 21.51055, 31.69724, 42.37560, 96.12094,
          69.81919, 58.44371, 66.13142, 92.36598, 85.84545, 97.95544, 97.67861,
          98.66029, 95.53441, 61.64348, 82.68268, 92.35338, 95.93326, 96.47150,
          85.64553, 100.00000, 87.57085, 100.00000, 56.42384, 60.79521, 97.39357,
          100.00000, 98.98271, 84.51251, 96.53179, 99.39140, 99.67509, 97.22519,
          88.28658, 80.46175, 74.29626, 91.46719,  96.81733,  99.84994))

df_2 <- data.frame(stv = c(1, 1, 11, 11, 21, 21, 31, 31),
                 ry = c(95, 100, 95, 100, 85, 90, 85, 90))

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

AICc(fit_lm) # can we use AIC here?

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

# 8.1 Create residuals
## residuals based on transformed soil test values
alcc_res <- yt - pred_yt
## or do we use residuals based on back-transformed soil test values?
alcc_res2 <- exp(yt) - exp(pred_yt) 

# 8.2 Estimate CI
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

# logLik

res <- alcc_res ## Replace this line with extraction of residuals from ALCC

N <- length(res)

w <- rep_len(1, N)

zw <- w == 0

N <- sum(!zw)

val <- -N * (log(2 * pi) + 1 - log(N) - sum(log(w + zw))/N + log(sum(res^2)))/2

val # value is negative

# AIC