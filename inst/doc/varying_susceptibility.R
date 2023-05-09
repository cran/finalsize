## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  message = FALSE,
  warning = FALSE,
  dpi = 300
)

## ----setup, message=FALSE, warning=FALSE, class.source = 'fold-hide'----------
# load finalsize
library(finalsize)

# load necessary packages
if (!require("socialmixr")) install.packages("socialmixr")
if (!require("ggplot2")) install.packages("ggplot2")
if (!require("colorspace")) install.packages("colorspace")

library(ggplot2)
library(colorspace)

## ----class.source = 'fold-hide'-----------------------------------------------
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
contact_matrix <- contact_matrix / max(Re(eigen(contact_matrix)$values))

# divide each row of the contact matrix by the corresponding demography
contact_matrix <- contact_matrix / demography_vector

n_demo_grps <- length(demography_vector)

## ----class.source = 'fold-hide'-----------------------------------------------
r0 <- 1.5

## -----------------------------------------------------------------------------
# susceptibility is higher for the old
susc_variable <- matrix(
  data = c(0.75, 0.8, 0.85, 0.9, 1.0)
)
n_susc_groups <- 1L

## -----------------------------------------------------------------------------
p_susc_uniform <- matrix(
  data = 1.0,
  nrow = n_demo_grps,
  ncol = n_susc_groups
)

## -----------------------------------------------------------------------------
# calculate the effective R0 using `r_eff()`
r_eff(
  r0 = r0,
  contact_matrix = contact_matrix,
  demography_vector = demography_vector,
  susceptibility = susc_variable,
  p_susceptibility = p_susc_uniform
)

## -----------------------------------------------------------------------------
# run final_size with default solvers and control options
# final size with heterogeneous susceptibility
final_size_heterog <- final_size(
  r0 = r0,
  contact_matrix = contact_matrix,
  demography_vector = demography_vector,
  susceptibility = susc_variable,
  p_susceptibility = p_susc_uniform
)

## -----------------------------------------------------------------------------
# prepare uniform susceptibility matrix
susc_uniform <- matrix(1.0, nrow = n_demo_grps, ncol = n_susc_groups)

# run final size with uniform susceptibility
final_size_uniform <- final_size(
  r0 = r0,
  contact_matrix = contact_matrix,
  demography_vector = demography_vector,
  susceptibility = susc_uniform,
  p_susceptibility = p_susc_uniform
)

## -----------------------------------------------------------------------------
# assign scenario name and join data
final_size_heterog$scenario <- "heterogeneous"
final_size_uniform$scenario <- "uniform"

# join dataframes
final_size_data <- rbind(
  final_size_heterog,
  final_size_uniform
)

# prepare age group order
final_size_data$demo_grp <- factor(
  final_size_data$demo_grp,
  levels = contact_data$demography$age.group
)

# examine the combined data
final_size_data

## ----class.source = 'fold-hide', fig.cap="Final sizes of epidemics in populations wherein susceptibility to the infection is either uniform (green), or heterogeneous (purple), with older individuals more susceptible to the infection.", fig.width=5, fig.height=4----
ggplot(final_size_data) +
  geom_col(
    aes(
      x = demo_grp, y = p_infected,
      fill = scenario
    ),
    col = "black",
    position = position_dodge(
      width = 0.75
    )
  ) +
  expand_limits(
    x = c(0.5, length(unique(final_size_data$demo_grp)) + 0.5)
  ) +
  scale_fill_discrete_qualitative(
    palette = "Cold",
    name = "Population susceptibility",
    labels = c("Heterogeneous", "Uniform")
  ) +
  scale_y_continuous(
    labels = scales::percent,
    limits = c(0, 1)
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
    expand = FALSE
  ) +
  labs(
    x = "Age group",
    y = "% Infected"
  )

## -----------------------------------------------------------------------------
# immunisation effect
immunisation_effect <- 0.25

## -----------------------------------------------------------------------------
# model an immunised group with a 25% lower susceptibility
susc_immunised <- cbind(
  susc_variable,
  susc_variable * (1 - immunisation_effect)
)

# assign names to groups
colnames(susc_immunised) <- c("Un-immunised", "Immunised")
n_risk_groups <- ncol(susc_immunised)

## -----------------------------------------------------------------------------
# immunisation rate is uniform across age groups
immunisation_rate <- rep(0.5, n_demo_grps)

# add a second column to p_susceptibility
p_susc_immunised <- cbind(
  susceptible = p_susc_uniform - immunisation_rate,
  immunised = immunisation_rate
)

## -----------------------------------------------------------------------------
# we run final size over all r0 values
final_size_immunised <- final_size(
  r0 = r0,
  contact_matrix = contact_matrix,
  demography_vector = demography_vector,
  susceptibility = susc_immunised,
  p_susceptibility = p_susc_immunised
)

## -----------------------------------------------------------------------------
# add scenario identifier
final_size_immunised$scenario <- "immunisation"

# prepare age group order
final_size_heterog$demo_grp <- factor(
  final_size_heterog$demo_grp,
  levels = contact_data$demography$age.group
)

final_size_immunised$demo_grp <- factor(
  final_size_immunised$demo_grp,
  levels = contact_data$demography$age.group
)

# examine the data
final_size_immunised

## ----fig.cap="Final size of an SIR epidemic with $R_0$ = 1.5, in a population wherein 50% of each age group is immunised against the infection. The immunisation is assumed to reduce the initial susceptibility of each age group by 25%. This leads to both within- and between-group heterogeneity in susceptibility. Vaccinating even 50% of each age group can substantially reduce the epidemic final size in comparison with a scenario in which there is no immunisation (grey). Note that the final sizes in this figure are all below 50%.", fig.width=5, fig.height=4, class.source = 'fold-hide'----
ggplot(final_size_immunised) +
  geom_col(
    data = final_size_heterog,
    aes(
      x = demo_grp, y = p_infected,
      fill = "baseline",
      colour = "baseline"
    ),
    width = 0.75,
    show.legend = TRUE
  ) +
  geom_col(
    aes(
      x = demo_grp, y = p_infected,
      fill = susc_grp
    ),
    col = "black",
    position = position_dodge()
  ) +
  facet_grid(
    cols = vars(scenario),
    labeller = labeller(
      scenario = c(
        heterogeneous = "Between groups only",
        immunisation = "Within & between groups"
      )
    )
  ) +
  expand_limits(
    x = c(0.5, length(unique(final_size_immunised$demo_grp)) + 0.5)
  ) +
  scale_fill_discrete_qualitative(
    palette = "Dynamic",
    rev = TRUE,
    limits = c("Immunised", "Un-immunised"),
    name = "Immunisation scenario",
    na.value = "lightgrey"
  ) +
  scale_colour_manual(
    values = "black",
    name = "No immunisation",
    labels = "Susceptibility homogeneous\nwithin groups"
  ) +
  scale_y_continuous(
    labels = scales::percent,
    limits = c(0, 0.5)
  ) +
  theme_classic() +
  theme(
    legend.position = "bottom",
    legend.key.height = unit(2, "mm"),
    legend.title = ggtext::element_markdown(
      vjust = 1
    ),
    strip.background = element_blank(),
    strip.text = element_text(
      face = "bold",
      size = 11
    )
  ) +
  guides(
    colour = guide_legend(
      override.aes = list(fill = "lightgrey"),
      title.position = "top",
      order = 1
    ),
    fill = guide_legend(
      nrow = 2,
      title.position = "top",
      order = 2
    )
  ) +
  coord_cartesian(
    expand = FALSE
  ) +
  labs(
    x = "Age group",
    y = "% Infected",
    title = "Heterogeneous susceptibility",
    fill = "Immunisation\nscenario"
  )

## -----------------------------------------------------------------------------
# define r0
r0 <- 1.5

# define UK population size and prepare contact matrix
uk_pop <- 67 * 1e6
contact_matrix <- matrix(1.0) / uk_pop

# define susceptibility matrix
susceptibility <- matrix(c(1.0, 0.7), nrow = 1, ncol = 2)

# define p_susceptibility
p_susceptibility <- matrix(c(0.7, 0.3), nrow = 1, ncol = 2)

# running final_size()
final_size(
  r0 = r0,
  demography_vector = uk_pop,
  contact_matrix = contact_matrix,
  susceptibility = susceptibility,
  p_susceptibility = p_susceptibility
)

