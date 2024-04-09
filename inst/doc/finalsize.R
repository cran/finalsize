## ----include = FALSE----------------------------------------------------------
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

## -----------------------------------------------------------------------------
# define r0 as 1.5
r0 <- 1.5

## -----------------------------------------------------------------------------
# get UK population size
uk_pop <- 67 * 1e6

## -----------------------------------------------------------------------------
# prepare contact matrix
contact_matrix <- matrix(1.0) / uk_pop

## -----------------------------------------------------------------------------
# all individuals are fully susceptible
susceptibility <- matrix(1.0)

## -----------------------------------------------------------------------------
# all individuals are in the single, high-susceptibility group
p_susceptibility <- matrix(1.0)

## -----------------------------------------------------------------------------
# calculate final size
final_size_data <- final_size(
  r0 = r0,
  contact_matrix = contact_matrix,
  demography_vector = uk_pop,
  susceptibility = susceptibility,
  p_susceptibility = p_susceptibility
)

# view the output data frame
final_size_data

## -----------------------------------------------------------------------------
final_size(r0)

