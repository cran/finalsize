## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  message = FALSE,
  warning = FALSE,
  fig.height = 4,
  fig.width = 5,
  dpi = 150
)

## ----setup, message=FALSE, warning=FALSE, class.source = 'fold-hide'----------
# load finalsize
library(finalsize)
library(ggplot2)

## -----------------------------------------------------------------------------
# future time period considered
years_range <- seq(0, 20)

# proportion losing immunity per year
annual_wane <- 0.05

# proportion newly susceptible per year
birth_rate <- 14 / 1000

# define r0 values
r0_mean <- 2.0
r0_samples <- rnorm(n = 100, mean = r0_mean, sd = 0.3) # for first outbreak
r0_samples2 <- rnorm(n = 100, mean = r0_mean, sd = 0.3) # for second outbreak

## -----------------------------------------------------------------------------
# run final size for each value and collate estimates
# use `Map()` to iterate over values and an index to add to the final data
final_size_est <- Map(
  function(r0, i) {
    fs <- final_size(r0)$p_infected # run final size model, get estimate only
    data.frame(fs = fs, i = i, r0 = r0) # data.frame for each value considered
  },
  r0 = r0_samples,
  i = seq_along(r0_samples)
)

# use Reduce to combine the list of data.frames into a single data.frame
final_size_df <- Reduce(rbind, final_size_est)

## -----------------------------------------------------------------------------
# define changes in susceptibility over time from waning
# and population turnover (assuming) rectangular age distribution
relative_immunity <- (1 - annual_wane)^(years_range) *
  pmax(1 - birth_rate * years_range, 0) # use pmax in case of longer year range

## -----------------------------------------------------------------------------
# calculate effective R for each value and collate estimates
r_eff_est <- Map(
  function(r0, fs) {
    r_eff <- (1 - fs * relative_immunity) * r0

    # data.frame for each value considered
    data.frame(r_eff = r_eff, yr = years_range)
  },
  r0 = r0_samples2,
  fs = final_size_df$fs
)

# use Reduce to combine the list of data.frames into a single data.frame
r_eff_df <- Reduce(rbind, r_eff_est)

## ----class.source = 'fold-hide', class.source = 'fold-hide', fig.cap="Predicted value of the effective reproduction number $R$ following an initial epidemic of an emerging infection in year 0. We assume 5% of the immune population becomes susceptible again each year following the epidemic, and there is an influx of new susceptibles via a birth rate of 14/1000 per year. We also assume $R_0$ has a mean of 2, with a standard deviation around this estimate of 0.3."----
# plot results
ggplot(r_eff_df) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "gray") +
  geom_vline(xintercept = 0, colour = "red", linetype = "dashed") +
  annotate(
    geom = "text",
    x = 0, y = 2, angle = 90, vjust = "outward",
    colour = "red",
    label = "Initial epidemic at year 0"
  ) +
  stat_summary(
    aes(
      yr, r_eff
    ),
    fun = mean,
    fun.min = function(x) {
      quantile(x, 0.025)
    },
    fun.max = function(x) {
      quantile(x, 0.975)
    }
  ) +
  scale_y_continuous(
    limits = c(0, 3)
  ) +
  theme_classic() +
  coord_cartesian(
    expand = TRUE
  ) +
  labs(
    x = "Years after initial epidemic",
    y = "Effective R (mean and 95% CI)"
  )

