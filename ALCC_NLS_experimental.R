library(tidyverse)
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

df$ry <- ifelse(df$ry > 100, 100, df$ry)

ry_prop <- df$ry/100

ry_transformed <- asin(sqrt(ry_prop))
stv_transformed <- log(df$stv)

var(ry_prop)
var(ry_transformed)

plot(ry_transformed ~ stv_transformed)

df <- df %>% filter(stv < 300)

minx <- min(df$stv)

fit_alcc <- function(x, a, b){
    100 * sin(a*log(b * x))^2
}

fit <- nls(ry ~ fit_alcc(stv, a, b),
           data = df,
           start = c(a = 0.3, b = 0.3),
           # lower = c(a = 0.1, b = 0.1),
           # upper = c(a = 2, b = 30)
)

fit_alcc <- function(x, a){
    100 * sin(a*log(x))^2
}


fit <- nls(ry ~ fit_alcc(stv, a),
           data = df,
           start = c(a = 0.3),
           # lower = c(a = 0.1),
           # upper = c(a = 2)
           )

fit

fit_alcc <- function(x, a){
    100 * sin((a*x))^2
}


fit <- nls(ry ~ fit_alcc(stv, a),
           data = df,
           start = c(a = 0.3),
           # lower = c(a = 0.1),
           # upper = c(a = 2)
)

fit

tmp_df <- tibble(stv = seq(minx, 200, 0.1))

pred <- modelr::add_predictions(tmp_df, fit)

maxx <- filter(pred, pred > 99.9 & stv > minx) %>% pull(stv) %>% mean()

cstv <- filter(pred, pred > 94.9 & pred < 95.1 & stv < maxx) %>% pull(stv) %>% mean()

cstv

ggplot(data = df) +
    geom_point(aes(stv, ry)) +
    geom_line(data = pred, aes(stv, pred)) +
    geom_hline(yintercept = 95) +
    geom_vline(xintercept = cstv, color = "red")

soiltestcorr::mod_alcc(df, ry, stv, target = 95, plot = TRUE)

AIC(fit)
modelr::rsquare(fit, df)

# compare to alcc_curve

alcc_curve <-
    # same as soiltestcorr::mod_alcc but with finer resolution
    soiltestcorr::mod_alcc(df, ry, stv,
           target = 95, tidy = TRUE) %>% 
    unnest(Curve) %>% 
    select(ends_with("fitted")) %>% 
    rename(x = STV.fitted,
           ALCC = RY.fitted) %>% 
    arrange(x)

ggplot(data = df) +
    geom_point(aes(stv, ry)) +
    geom_line(data = pred, aes(stv, pred)) +
    geom_line(data = alcc_curve, aes(x, ALCC), color = "red3") +
    geom_hline(yintercept = 95) +
    geom_vline(xintercept = cstv, color = "blue3")
