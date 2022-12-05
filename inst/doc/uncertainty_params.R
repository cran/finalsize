## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  message = FALSE,
  warning = FALSE,
  dpi = 300
)
options(rmarkdown.html_vignette.check_title = FALSE)

## ----setup, message=FALSE, warning=FALSE, class.source = 'fold-hide'----------
library(finalsize)

# load necessary packages
library(data.table)
library(ggplot2)

## -----------------------------------------------------------------------------
# get UK polymod data
polymod <- socialmixr::polymod
contact_data <- socialmixr::contact_matrix(
  polymod,
  countries = "United Kingdom",
  age.limits = c(0, 5, 18, 40, 65),
  symmetric = TRUE
)

# get the contact matrix and demography data
contact_matrix <- t(contact_data$matrix)
demography_vector <- contact_data$demography$population

# scale the contact matrix so the largest eigenvalue is 1.0
contact_matrix <- contact_matrix / max(eigen(contact_matrix)$values)

# divide each row of the contact matrix by the corresponding demography
contact_matrix <- contact_matrix / demography_vector

n_demo_grps <- length(demography_vector)

## -----------------------------------------------------------------------------
# mean R0 is 2.0
r0_mean <- 2.0

## -----------------------------------------------------------------------------
# susceptibility is uniform
susc_uniform <- matrix(
  data = 1,
  nrow = n_demo_grps,
  ncol = 1L
)

# p_susceptibility is uniform
p_susc_uniform <- susc_uniform

## -----------------------------------------------------------------------------
# create an R0 samples vector
r0_samples <- rnorm(n = 1000, mean = r0_mean, sd = 0.1)

## -----------------------------------------------------------------------------
# run final size on each sample with the same data
final_size_data <- Map(
  r0_samples, seq_along(r0_samples),
  f = function(r0, i) {
    # the i-th final size estimate
    fs <- final_size(
      r0 = r0,
      contact_matrix = contact_matrix,
      demography_vector = demography_vector,
      susceptibility = susc_uniform,
      p_susceptibility = p_susc_uniform
    )

    fs$replicate <- i
    fs$r0_estimate <- r0
    fs
  }
)

# combine data
final_size_data <- rbindlist(final_size_data)

# order age groups
final_size_data$demo_grp <- factor(
  final_size_data$demo_grp,
  levels = contact_data$demography$age.group
)

## ----class.source = 'fold-hide', class.source = 'fold-hide', fig.cap="Estimated ranges of the final size of a hypothetical SIR epidemic in age groups of the U.K. population, when the $R_0$ is estimated to be 2.0, with a standard deviation around this estimate of 0.1. In this example, relatively low uncertainty in $R_0$ estimates can also lead to uncertainty in the estimated final size of the epidemic. Points represent means, while ranges extend between the 5th and 95th percentiles.", fig.width=5, fig.height=4----
ggplot(final_size_data) +
  stat_summary(
    aes(
      demo_grp, p_infected
    ),
    fun = mean,
    fun.min = function(x) {
      quantile(x, 0.05)
    },
    fun.max = function(x) {
      quantile(x, 0.95)
    }
  ) +
  scale_y_continuous(
    labels = scales::percent,
    limits = c(0.25, 1)
  ) +
  theme_classic() +
  theme(
    legend.position = "top",
    legend.key.height = unit(2, "mm"),
    legend.title = ggtext::element_markdown(
      vjust = 1
    )
  ) +
  coord_cartesian(
    expand = TRUE
  ) +
  labs(
    x = "Age group",
    y = "% Infected"
  )

