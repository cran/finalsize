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
# load necessary packages
if (!require("dplyr")) install.packages("dplyr")
if (!require("tibble")) install.packages("tibble")

library(dplyr)
library(tibble)

## -----------------------------------------------------------------------------
# susceptibility matrix
susceptibility <- tibble(
  age_group = c("[0,5)", "[5,18)", "[18,40)", "[40,65)", "65+"),
  susceptible = c(0.8, 0.8, 0.8, 0.8, 0.8)
) %>%
  column_to_rownames(var = "age_group") %>%
  as.matrix()

susceptibility

## -----------------------------------------------------------------------------
# demography-in-susceptibility matrix
p_susceptibility <- tibble(
  age_group = c("[0,5)", "[5,18)", "[18,40)", "[40,65)", "65+"),
  susceptible = c(1.0, 1.0, 1.0, 1.0, 1.0)
) %>%
  column_to_rownames(var = "age_group") %>%
  as.matrix()

p_susceptibility

## -----------------------------------------------------------------------------
# susceptibility matrix
susceptibility <- tibble(
  age_group = c("[0,5)", "[5,18)", "[18,40)", "[40,65)", "65+"),
  susceptible = c(0.2, 0.5, 0.6, 0.9, 1.0)
) %>%
  column_to_rownames(var = "age_group") %>%
  as.matrix()

susceptibility

# demography-in-susceptibility matrix
p_susceptibility <- tibble(
  age_group = c("[0,5)", "[5,18)", "[18,40)", "[40,65)", "65+"),
  susceptible = c(1.0, 1.0, 1.0, 1.0, 1.0)
) %>%
  column_to_rownames(var = "age_group") %>%
  as.matrix()

p_susceptibility

## -----------------------------------------------------------------------------
immunization_effect <- 0.25

# susceptibility matrix
susceptibility <- tibble(
  age_group = c("[0,5)", "[5,18)", "[18,40)", "[40,65)", "65+"),
  unimmunised = c(0.2, 0.5, 0.6, 0.9, 1.0)
) %>%
  mutate(immunised = unimmunised * (1 - immunization_effect)) %>%
  column_to_rownames(var = "age_group") %>%
  as.matrix()

susceptibility

# demography-in-susceptibility matrix
p_susceptibility <- tibble(
  age_group = c("[0,5)", "[5,18)", "[18,40)", "[40,65)", "65+"),
  immunised = c(0.4, 0.4, 0.4, 0.4, 0.4)
) %>%
  mutate(unimmunised = 1 - immunised) %>%
  column_to_rownames(var = "age_group") %>%
  as.matrix()

p_susceptibility

## -----------------------------------------------------------------------------
immunization_effect <- 0.25

# susceptibility matrix
susceptibility <- tibble(
  age_group = c("[0,5)", "[5,18)", "[18,40)", "[40,65)", "65+"),
  unimmunised = c(0.2, 0.5, 0.6, 0.9, 1.0)
) %>%
  mutate(immunised = unimmunised * (1 - immunization_effect)) %>%
  column_to_rownames(var = "age_group") %>%
  as.matrix()

susceptibility

# demography-in-susceptibility matrix
p_susceptibility <- tibble(
  age_group = c("[0,5)", "[5,18)", "[18,40)", "[40,65)", "65+"),
  immunised = c(0.2, 0.4, 0.6, 0.7, 0.9)
) %>%
  mutate(unimmunised = 1 - immunised) %>%
  column_to_rownames(var = "age_group") %>%
  as.matrix()

p_susceptibility

## -----------------------------------------------------------------------------
immunization_effect <- 0.25

# susceptibility matrix
susceptibility <- tibble(
  age_group = c("[0,5)", "[5,18)", "[18,40)", "[40,65)", "65+"),
  susceptible = c(1.0, 1.0, 1.0, 1.0, 1.0),
  unimmunised = c(0.2, 0.5, 0.6, 0.9, 1.0)
) %>%
  mutate(immunised = unimmunised * (1 - immunization_effect)) %>%
  column_to_rownames(var = "age_group") %>%
  as.matrix()

susceptibility

# demography-in-susceptibility matrix
p_susceptibility <- tibble(
  age_group = c("[0,5)", "[5,18)", "[18,40)", "[40,65)", "65+"),
  susceptible = c(0.1, 0.1, 0.1, 0.1, 0.1),
  immunised = c(0.4, 0.4, 0.4, 0.4, 0.4)
) %>%
  mutate(unimmunised = 1 - immunised - susceptible) %>%
  column_to_rownames(var = "age_group") %>%
  as.matrix()

p_susceptibility

