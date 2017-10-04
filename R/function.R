#' @title ode_functions
#' @description The set of ODE functions describing the minimal network for the central nervous system development.
#' @param x_old A named vector containing gene expression for 13 genes in the minimal network.
#' @param dt The time step used in the simulation
#' @param mx the number of lags to include
#' @param n the number of lags to include
#' @param N the number of lags to include
#' @param a the number of lags to include
#' @param a_s the number of lags to include
#' @param a_e the number of lags to include
#' @param k the number of lags to include
#' @param eta the number of lags to include
#' @param eta_m the number of lags to include
#' @param mature_mu the number of lags to include
#' @param eta_b the number of lags to include
#' @return the predicted value of y
#' @keywords predict
#' @seealso \code{\link{simulate_cns}}
#' @examples
#' x_old <- c(2, rep(0, 12))
#' names(x_old) <- c("Pax6","Mash1","Brn2","Zic1","Tuj1","Hes5","Scl","Olig2","Stat3","Myt1L","Aldh1L","Sox8","Mature")
# ode_functions(x_old, dt = 0.01)
#'
ode_functions <- function(x_old, dt, mx = 10,
                          N_end = 40, n = 4,
                          N = 13, a = 4, a_s = 2.2, a_e = 6,
                          k = 1, eta = 0.25, eta_m = 0.125,
                          mature_mu = 0, eta_b = 0.1) {
  dx <- c() # vector for the data

  dx["Pax6"] <- a_s * 1 / (1 + eta^n*(x_old["Tuj1"]+x_old["Aldh1L"]+x_old["Olig2"])^n * x_old["Mature"]^n ) - k * x_old["Pax6"]
  dx["Mash1"] <- a * ( x_old["Pax6"]^n ) / (1 + x_old["Pax6"]^n + x_old["Hes5"]^n) - k * x_old["Mash1"]
  dx["Brn2"] <- a * ( x_old["Mash1"]^n ) / (1 + x_old["Mash1"]^n) - k * x_old["Brn2"]
  dx["Zic1"] <- a * ( x_old["Mash1"]^n ) / (1 + x_old["Mash1"]^n) - k * x_old["Zic1"]
  dx["Tuj1"] <- a_e * ( x_old["Brn2"]^n + x_old["Zic1"]^n + x_old["Myt1L"]^n ) / (1 + x_old["Brn2"]^n + x_old["Zic1"]^n + x_old["Myt1L"]^n) - k*x_old["Tuj1"]
  dx["Hes5"] <- a * ( x_old["Pax6"]^n ) / (1 +  x_old["Pax6"]^n + x_old["Mash1"]^n)  - k * x_old["Hes5"]
  dx["Scl"] <- a_e * ( eta^n * x_old["Hes5"]^n ) / (1 + eta ^ n * x_old["Hes5"]^n + x_old["Olig2"]^n) - k * x_old["Scl"]
  dx["Olig2"] <- a_e * ( (eta * x_old["Hes5"])^n ) / (1 + (x_old["Scl"])^n + (eta * x_old["Hes5"])^n  ) - k * x_old["Olig2"]
  dx["Stat3"] <- a * ( eta^n * x_old["Hes5"]^n * x_old["Scl"]^n ) / (1 + eta^n * x_old["Hes5"]^n * x_old["Scl"]^n ) - k * x_old["Stat3"]
  dx["Myt1L"] <- a * ( x_old["Olig2"]^n ) / (1 + x_old["Olig2"]^n) - k * x_old["Myt1L"]
  dx["Aldh1L"] <- a_e * ( x_old["Stat3"]^n ) / (1 + x_old["Stat3"]^n) - k * x_old["Aldh1L"]
  dx["Sox8"] <- a * ( eta_m^n * x_old["Olig2"]^n ) / (1 + eta_m^n * x_old["Olig2"]^n) - k * x_old["Sox8"]
  dx["Mature"] <- mature_mu * (1 - x_old["Mature"] / mx)

  dx <- dx * dt + rnorm(N, mean = 0, sd = .01)

  dx["Mature"] <- 0

  return(x_old + dx)
}

#' @title simulate_cns
#' @description This function simulates the differentiation process of central nervous system and generating a time series
#' where different cells differentiate into different cell lineages (neuron, astrocyte or oligodendrocytes).
#' @param steps Number of total steps for the simulation.
#' @param cells Numer of cells you want to simulate.
#' @param N_end The end time point for your simulation.
#' @param mx parameter mx in the ODE model
#' @param n parameter mx in the ODE model
#' @param N parameter mx in the ODE model
#' @param a parameter mx in the ODE model
#' @param a_s parameter mx in the ODE model
#' @param a_e parameter mx in the ODE model
#' @param k parameter mx in the ODE model
#' @param eta parameter mx in the ODE model
#' @param eta_m parameter mx in the ODE model
#' @param mature_mu parameter mx in the ODE model
#' @param eta_b parameter mx in the ODE model
#' @return the predicted value of y
#' @keywords predict
#' @seealso \code{\link{ode_functions}}
#' @examples
#' simulate_cns()
#'
simulate_cns <- function(steps = 1000,
                         cells = 25,
                         N_end = 40,
                         mx = 10, n = 4, N = 13, a = 4, a_s = 2.2,
                         a_e = 6, k = 1, eta = 0.25, eta_m = 0.125,
                         mature_mu = 0, eta_b = 0.1) {

  mu <- matrix(0, nrow = steps, ncol = N - 1)
  dt <- N_end / steps

  cell_simulate <- array(0, dim = c(N, steps + 1, cells))
  f_em <- matrix(0, nrow = N, ncol = steps + 1)

  f_init <- rep(0, N)
  names(f_init) <- c("Pax6","Mash1","Brn2","Zic1","Tuj1","Hes5","Scl","Olig2","Stat3","Myt1L","Aldh1L","Sox8","Mature")
  f_init[1] <- 2

  D <- array(0, dim = c(N - 1, N - 1, steps))
  f_em[, 1] <- f_init

  for(cell in 1:cells) {
    f_temp <- f_init

    for(step in 1:steps) {
      f_temp <- ode_functions(f_temp, dt,
                              mx = mx, n = n, N = N, a = a, a_s = a_s,
                              a_e = a_e, k = k, eta = eta, eta_m = eta_m,
                              mature_mu = mature_mu, eta_b = eta_b)
      f_temp[f_temp < 0] <- 0
      f_em[, step + 1] <- f_temp
    }
    cell_simulate[, , cell] <- f_em
  }

  return(cell_simulate)
}

