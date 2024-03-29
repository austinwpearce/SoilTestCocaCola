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

# df <- df %>% filter(stv < 300)

minx <- min(df$stv)

fit_alcc <- function(x, a, b){
    100 * sin((-a / b) + (log(x) / b))^2
}

fit <- nls(ry ~ fit_alcc(stv, a, b),
           data = df,
           algorithm = "port",
           start = c(a = 0.5, b = 2.5),
           lower = c(a = -10, b = 1.5),
           upper = c(a = 10, b = 10)
)

fit_alcc <- function(y, a, b){
    exp(a + b * asin(sqrt(y / 100)))
}

fit <- nls(stv ~ fit_alcc(ry, a, b),
           data = df,
           start = c(a = 0.5, b = 2.5))

fit

tmp_df <- tibble(stv = seq(minx, 200, 0.1))

tmp_df <- tibble(ry = seq(0, 100, 1))

pred <- modelr::add_predictions(tmp_df, fit)

maxx <- filter(pred, pred > 99.9 & stv > minx) %>% pull(stv) %>% mean()

cstv <- filter(pred, pred > 94.9 & pred < 95.1 & stv < maxx) %>% pull(stv) %>% mean()

cstv

tmp_df <- subset(df, stv < maxx)

pred <- predict(fit, newdata = tmp_df)

resid_ry <- pred %>%
    filter(stv < maxx & stv %in% df$stv) %>% 
    left_join(df) %>% 
    mutate(resid = ry - pred)

ggplot(data = df) +
    geom_point(aes(stv, ry)) +
    geom_line(data = pred, aes(stv, pred)) +
    geom_point(data = resid_ry,
               aes(stv, pred), color = "red") +
    geom_hline(yintercept = 95) +
    geom_vline(xintercept = cstv, color = "red")


soiltestcorr::mod_alcc(df %>% filter(stv < 140), ry, stv, target = 95, plot = TRUE)

AIC(fit)
modelr::rsquare(fit, df)

library(nlraa)
fit2 <- nls(ry ~ SSlinp(stv, a, b, j), data = df)
AIC(fit2)

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

# AIC

# logLik

res <- resid_ry$resid ## Replace this line with extraction of residuals from ALCC

N <- length(res)

w <- rep_len(1, N)

zw <- w == 0

N <- sum(!zw)

val <- -N * (log(2 * pi) + 1 - log(N) - sum(log(w + zw))/N + log(sum(res^2)))/2

val # value is negative

# AIC


alccfun <- function(y, a, b){
    c <- asin(sqrt(0.95))
    exp(b * asin(sqrt(y/100)) + (a - b * c))
}

tmp <- m1k_data %>% filter(x < 153)

mod <- nls(x ~ alccfun(y, a, b), data = tmp,
           start = c(a = 4, b = 2))


ggplot()+
    geom_point(aes(m1k_data$x, y = m1k_data$y, color = m1k_data$sig_response)) +
    geom_line(aes(x = predict(mod), y = tmp$y)) +
    geom_vline(xintercept = 67.7) + geom_hline(yintercept = 95)


alccfun <- function(x, a, b){
    c <- asin(sqrt(0.95))
    100 * (sin((log(x)/b) - ((a - b * c) / b)))^2
}

mod <- nls(y ~ alccfun(x, a, b), data = tmp,
           start = c(a = 4, b = 2))
# not quite, but close

ggplot()+
    geom_point(aes(m1k_data$x, y = m1k_data$y, color = m1k_data$sig_response)) +
    geom_line(aes(x = predict(mod), y = tmp$y)) +
    geom_vline(xintercept = 67.7) + geom_hline(yintercept = 95)
