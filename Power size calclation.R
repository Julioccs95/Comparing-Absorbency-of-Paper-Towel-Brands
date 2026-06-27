rcbd_power_delta <- function(blocks, #b
                             treatments, #v
                             reps, #S in the book treatments
                             Delta,
                             sigma2,
                             alpha = 0.05) {
  
  v <- treatments
  b <- blocks
  r <- reps 
  
  #v1 in the book numerator degrees of freedom
  df1 <- 1
  
  #n=bvs
  N <- v * b * r
  df_error <- N-v
  
  # Noncentrality parameter for pairwise Tukey comparison
  lambda <- (b * r * Delta^2) / (2 * sigma2)
  
  # Tukey critical value, with the number of means being treatment, and df_error
  qcrit <- qtukey(1 - alpha, nmeans = v, df = df_error)
  
  # Equivalent F critical value for Tukey
  fcrit <- (qcrit^2) / 2
  
  # Power using noncentral F
  power <- 1 - pf(fcrit, df1 = df1, df2 = df_error, ncp = lambda)
  
  return(power)
}











blocks <- 2
treatments <- 3
target_power <- 0.80

Delta <- 2
sigma2 <- 2.20

reps <- 2:12
reps
results <- data.frame(
  reps_per_block_treatment = reps,
  samples_per_block = treatments * reps,
  total_sample_size = blocks * treatments * reps,
  #Take each value of reps and store results in power
  power = sapply(reps, function(r) {
    rcbd_power_delta(
      blocks = blocks,
      treatments = treatments,
      reps = r,
      Delta = Delta,
      sigma2 = sigma2,
      alpha = 0.05
    )
  })
)

results




sample(18)











#---------------------------------------------------------------------------------------
#Use Sheffe to take into account all uncertain other possibilities
#Delta divided t-1
#
Tukey_sample_size <- data.frame()

for (r in 2:50) {
  
  alpha <- 0.05
  t <- 3
  b <- 2
  variance <- 2.20
  
  # Pairwise contrast coefficients, e.g. treatment 1 - treatment 2
  sum_c_square <- (1)^2 + (-1)^2
  
  
  # Total sample size
  n <- t * b * r
  
  # RCBD error degrees of freedom
  df <- n - b - t + 1
  
  # Tukey studentized range critical value
  q <- qtukey(1 - alpha, nmeans = t, df = df)
  
  # Tukey multiplier
  w <- q / sqrt(2)
  
  # Margin of error
  ME <- w * sqrt(variance * sum_c_square / (b * r))
  
  # Full confidence interval width
  L <- 2 * ME
  
  Tukey_sample_size <- rbind(
    Tukey_sample_size,
    data.frame(
      r = r,
      total_n = n,
      df = df,
      q = q,
      w = w,
      ME = ME,
      L = L
    )
  )
}

Tukey_sample_size