## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  message = FALSE,
  warning = FALSE,
  dpi = 300
)

## ----setup, message=FALSE, warning=FALSE, class.source = 'fold-hide'----------
library(finalsize)

# load necessary packages
library(data.table)
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
contact_matrix <- contact_matrix / max(eigen(contact_matrix)$values)

# divide each row of the contact matrix by the corresponding demography
contact_matrix <- contact_matrix / demography_vector

n_demo_grps <- length(demography_vector)

## ----class.source = 'fold-hide'-----------------------------------------------
r0 <- 2.0

## -----------------------------------------------------------------------------
# susceptibility is higher for the old
susc_variable <- matrix(
  data = c(0.2, 0.5, 0.6, 0.9, 1.0)
)
n_susc_groups <- 1L

## -----------------------------------------------------------------------------
p_susc_uniform <- matrix(
  data = 1.0,
  nrow = n_demo_grps,
  ncol = n_susc_groups
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
# immunisation increases with age between 20% (infants) and 90% (65+)
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

## ----fig.cap="Final size of an SIR epidemic with $R_0$ = 2.0, in a population wherein 50% of each age group is immunised against the infection. The immunisation is assumed to reduce the initial susceptibility of each age group by 25%. This leads to both within- and between-group heterogeneity in susceptibility. Vaccinating even 50% of each age group can substantially reduce the epidemic final size in comparison with a scenario in which there is no immunisation (grey). Note that the final sizes in this figure are all below 50%.", fig.width=5, fig.height=4, class.source = 'fold-hide'----
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
  scale_fill_discrete_qualitative(
    palette = "Dynamic",
    rev = TRUE,
    limits = c("Immunised", "Un-immunised"),
    name = "Immunisation scenario",
    na.value = "lightgrey"
  ) +
  scale_colour_manual(
    values = "black",
    name = NULL,
    labels = "Susceptibility heterogeneous\nbetween groups,\nno immunisation"
  ) +
  scale_y_continuous(
    labels = scales::percent,
    limits = c(0, 0.5)
  ) +
  theme_classic() +
  theme(
    legend.position = "top",
    legend.key.height = unit(2, "mm"),
    legend.title = ggtext::element_markdown(
      vjust = 1
    )
  ) +
  guides(
    colour = guide_legend(
      override.aes = list(fill = "lightgrey")
    ),
    fill = guide_legend(
      nrow = 2,
      title.position = "top"
    )
  ) +
  coord_cartesian(
    expand = FALSE
  ) +
  labs(
    x = "Age group",
    y = "% Infected",
    fill = "Immunisation\nscenario"
  )

