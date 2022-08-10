library(ggplot2)
library(agridat)
# ALCC demo

df <- data.frame(stv = agridat::cate.potassium$potassium,
                 ry = agridat::cate.potassium$yield)
df$ry <- ifelse(df$ry > 100, 100, df$ry)

df <- df %>% filter(stv < 300)

minx <- min(df$stv)

fit_alcc <- function(x, a, b){
    100 * sin(a*log(b * x))^2
}

fit <- nls(ry ~ fit_alcc(stv, a, b), data = df, start = c(a = 1, b = 0.1))

fit

tmp_df <- tibble(stv = seq(1, 200, 0.1))

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
