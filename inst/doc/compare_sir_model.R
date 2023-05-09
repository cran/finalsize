## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  message = FALSE,
  warning = FALSE,
  dpi = 300
)

## ----setup, class.source = 'fold-hide'----------------------------------------
library(finalsize)

## ----class.source = 'fold-hide'-----------------------------------------------
# Replicate the contact matrix
repmat <- function(X, m, n) {
  mx <- dim(X)[1]
  nx <- dim(X)[2]
  matrix(t(matrix(X, mx, nx * n)), mx * m, nx * n, byrow = TRUE)
}

# Calculate the proportion of individuals infected
final_size_r <- function(r0, contact_matrix,
                         demography_vector, susceptibility, p_susceptibility) {
  # check p_susceptibility
  stopifnot(all(rowSums(p_susceptibility) == 1))

  p_susc <- as.vector(p_susceptibility)
  susc <- as.vector(susceptibility)
  demo_v <- rep(demography_vector, ncol(p_susceptibility)) * p_susc
  nsteps <- 10000

  # epidemiological parameters
  h <- 0.1
  beta <- r0
  gamma <- 1

  # initial conditions
  I <- 1e-8 * demo_v
  S <- demo_v - I
  R <- 0

  # prepare contact matrix
  m <- repmat(contact_matrix, ncol(p_susceptibility), ncol(p_susceptibility))

  # multiply contact matrix by susceptibility
  m <- m * susc

  # iterate over compartmental transitions
  for (i in seq_len(nsteps)) {
    force <- (m %*% I) * S * h * beta
    S <- S - force
    g <- I * gamma * h
    I <- I + force - g
    R <- R + g
  }

  # return proportion of individuals recovered
  as.vector(R / demo_v)
}

## -----------------------------------------------------------------------------
r0 <- 2.0

## -----------------------------------------------------------------------------
# estimated population size for the UK is 67 million
population_size <- 67e6

# estimate using finalsize::final_size
final_size_finalsize <- final_size(
  r0 = r0,
  contact_matrix = matrix(1) / population_size,
  demography_vector = population_size,
  susceptibility = matrix(1),
  p_susceptibility = matrix(1),
  solver = "newton"
)

# estimate from SIR model
final_size_sir <- final_size_r(
  r0 = r0,
  contact_matrix = 0.99 * matrix(1) / population_size,
  demography_vector = population_size,
  susceptibility = matrix(1),
  p_susceptibility = matrix(1)
)

# View the estimates
final_size_finalsize

final_size_sir

## -----------------------------------------------------------------------------
# function to check error in final size estimates
final_size_error <- function(x, r0, tolerance = 1e-6) {
  error <- abs(x - (1.0 - exp(-r0 * x)))
  if (any(error > tolerance)) {
    print("Final size estimate error is greater than tolerance!")
  }
  # return error
  error
}

## -----------------------------------------------------------------------------
# error for estimate using finalsize::final_size()
final_size_error(final_size_finalsize$p_infected, r0)

# error for estimate from SIR model
final_size_error(final_size_sir, r0)

## ----class.source = 'fold-hide'-----------------------------------------------
# Load example POLYMOD data included with the package
data(polymod_uk)

# Define contact matrix (entry {ij} is contacts in group i
# reported by group j)
contact_matrix <- polymod_uk$contact_matrix

# Define population in each age group
demography_vector <- polymod_uk$demography_vector

## -----------------------------------------------------------------------------
contact_matrix

demography_vector

## -----------------------------------------------------------------------------
# Define susceptibility of each group
susceptibility <- matrix(
  data = c(1.0, 0.5, 0.5),
  nrow = length(demography_vector),
  ncol = 1
)

# Assume uniform susceptibility within age groups
p_susceptibility <- matrix(
  data = 1.0,
  nrow = length(demography_vector),
  ncol = 1
)

## -----------------------------------------------------------------------------
# from {finalsize} using finalsize::final_size
final_size(
  r0 = r0,
  contact_matrix = contact_matrix,
  demography_vector = demography_vector,
  susceptibility = susceptibility,
  p_susceptibility = p_susceptibility,
  solver = "newton"
)

# using SIR model
final_size_r(
  r0 = r0,
  contact_matrix = contact_matrix * 0.99,
  demography_vector = demography_vector,
  susceptibility = susceptibility,
  p_susceptibility = p_susceptibility
)

## -----------------------------------------------------------------------------
# Define susceptibility of each group
# Here
susceptibility <- matrix(
  c(1.0, 0.1),
  nrow = length(demography_vector), ncol = 2,
  byrow = TRUE
)

# view the susceptibility matrix
susceptibility

# Assume variation in susceptibility within age groups
# A higher proportion of 0 -- 20s are in the lower susceptibility group
p_susceptibility <- matrix(1, nrow = length(demography_vector), ncol = 2)
p_susceptibility[, 2] <- c(0.6, 0.5, 0.3)
p_susceptibility[, 1] <- 1 - p_susceptibility[, 2]

# view p_susceptibility
p_susceptibility

## -----------------------------------------------------------------------------
# estimate group-specific final sizes using finalsize::final_size
final_size(
  r0 = r0,
  contact_matrix = contact_matrix,
  demography_vector = demography_vector,
  susceptibility = susceptibility,
  p_susceptibility = p_susceptibility
)

# estimate group-specific final sizes from the SIR model
final_size_r(
  r0 = r0,
  contact_matrix = contact_matrix,
  demography_vector = demography_vector,
  susceptibility = susceptibility,
  p_susceptibility = p_susceptibility
)

## -----------------------------------------------------------------------------
# Define susceptibility of each group
susceptibility <- matrix(1, nrow = length(demography_vector), ncol = 2)
susceptibility[, 2] <- 0.0

# view susceptibility
susceptibility

# Assume that some proportion are completely immune
# complete immunity is more common among younger individuals
p_susceptibility <- matrix(1, nrow = length(demography_vector), ncol = 2)
p_susceptibility[, 2] <- c(0.6, 0.5, 0.3)
p_susceptibility[, 1] <- 1 - p_susceptibility[, 2]

# view p_susceptibility
p_susceptibility

## -----------------------------------------------------------------------------
final_size(
  r0 = r0,
  contact_matrix = contact_matrix,
  demography_vector = demography_vector,
  susceptibility = susceptibility,
  p_susceptibility = p_susceptibility
)

final_size_r(
  r0 = r0,
  contact_matrix = contact_matrix,
  demography_vector = demography_vector,
  susceptibility = susceptibility,
  p_susceptibility = p_susceptibility
)

