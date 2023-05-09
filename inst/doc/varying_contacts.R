## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  message = FALSE,
  warning = FALSE,
  dpi = 300
)

## -----------------------------------------------------------------------------
# load finalsize
library(finalsize)

# load necessary packages
if (!require("socialmixr")) install.packages("socialmixr")
if (!require("ggplot2")) install.packages("ggplot2")

library(ggplot2)

## -----------------------------------------------------------------------------
# define r0 as 1.5
r0 <- 1.5

## -----------------------------------------------------------------------------
# get UK polymod data
polymod <- socialmixr::polymod
contact_data <- socialmixr::contact_matrix(
  polymod,
  countries = "United Kingdom",
  age.limits = c(0, 5, 18, 40, 65),
  symmetric = TRUE
)

# view the elements of the contact data list
# the contact matrix
contact_data$matrix

# the demography data
contact_data$demography

# get the contact matrix and demography data
contact_matrix <- t(contact_data$matrix)
demography_vector <- contact_data$demography$population
demography_data <- contact_data$demography

# scale the contact matrix so the largest eigenvalue is 1.0
# this is to ensure that the overall epidemic dynamics correctly reflect
# the assumed value of R0
contact_matrix <- contact_matrix / max(Re(eigen(contact_matrix)$values))

# divide each row of the contact matrix by the corresponding demography
# this reflects the assumption that each individual in group {j} make contacts
# at random with individuals in group {i}
contact_matrix <- contact_matrix / demography_vector

n_demo_grps <- length(demography_vector)

## -----------------------------------------------------------------------------
# all individuals are equally and highly susceptible
n_susc_groups <- 1L
susc_guess <- 1.0

## -----------------------------------------------------------------------------
susc_uniform <- matrix(
  data = susc_guess,
  nrow = n_demo_grps,
  ncol = n_susc_groups
)

## -----------------------------------------------------------------------------
p_susc_uniform <- matrix(
  data = 1.0,
  nrow = n_demo_grps,
  ncol = n_susc_groups
)

## -----------------------------------------------------------------------------
# calculate final size
final_size_data <- final_size(
  r0 = r0,
  contact_matrix = contact_matrix,
  demography_vector = demography_vector,
  susceptibility = susc_uniform,
  p_susceptibility = p_susc_uniform
)

# view the output data frame
final_size_data

## ----class.source = 'fold.hide'-----------------------------------------------
# order demographic groups as factors
final_size_data$demo_grp <- factor(
  final_size_data$demo_grp,
  levels = demography_data$age.group
)

## ----fig.cap ="Final size of an SIR epidemic in each age group. The final size is the cumulative number of infections in each age group over the course of the epidemic, expressed as a proportion of the respective age group.", fig.width=5, fig.height=4, class.source = 'fold-hide'----
# plot data
ggplot(final_size_data) +
  geom_col(
    aes(
      demo_grp, p_infected
    ),
    colour = "black", fill = "grey"
  ) +
  scale_y_continuous(
    labels = scales::percent,
    limits = c(0, 1)
  ) +
  expand_limits(
    x = c(0.5, nrow(final_size_data) + 0.5)
  ) +
  theme_classic() +
  coord_cartesian(
    expand = FALSE
  ) +
  labs(
    x = "Age group",
    y = "% Infected"
  )

## -----------------------------------------------------------------------------
# prepare demography data
demography_data <- contact_data$demography

# merge final size counts with demography vector
final_size_data <- merge(
  final_size_data,
  demography_data,
  by.x = "demo_grp",
  by.y = "age.group"
)

# reset age group order
final_size_data$demo_grp <- factor(
  final_size_data$demo_grp,
  levels = contact_data$demography$age.group
)

# multiply counts with proportion infected
final_size_data$n_infected <- final_size_data$p_infected *
  final_size_data$population

## ----fig.cap="Final size of an epidemic outbreak in a population, for different values of infection $R_0$. Converting the final size proportions in each age group to counts shows that individuals aged 18 -- 64 make up the bulk of cases in this scenario. This may be attributed to this being both the largest age range in the analysis (more years in this range than any other), and because more people fall into this wide range than others. Contrast this figure with the one above, in which similar _proportions_ of each age group are infected.", fig.width=5, fig.height=4, class.source = 'fold-hide'----
ggplot(final_size_data) +
  geom_col(
    aes(
      x = demo_grp, y = n_infected
    ),
    fill = "grey", col = "black"
  ) +
  expand_limits(
    x = c(0.5, nrow(final_size_data) + 0.5)
  ) +
  scale_y_continuous(
    labels = scales::comma_format(
      scale = 1e-6, suffix = "M"
    ),
    limits = c(0, 15e6)
  ) +
  theme_classic() +
  coord_cartesian(
    expand = FALSE
  ) +
  labs(
    x = "Age group",
    y = "Number infected (millions)"
  )

