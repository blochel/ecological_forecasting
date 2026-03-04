# Load required package
library(deSolve)

# Seasonal water depth function (in cm)
water_depth <- function(t) {
  # Sinusoidal cycle: mean 30cm, amplitude 25cm, period 365 days
  30 + 25 * sin(2 * pi * t / 365)
}

# Bird abundance function - now depends on water depth AND small fish
bird_abundance <- function(t, S, bird_response, half_sat) {
  W <- water_depth(t)
  
  # Birds arrive when water is shallow
  if (W > 40) {
    B_water <- 10
  } else {
    B_water <- 10 + 90 * (40 - W) / 40
  }
  
  # Additional birds attracted by small fish density
  B_fish <- bird_response * S / (half_sat + S)
  
  return(B_water + B_fish)
}

# ODE model
fish_model <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    
    # Current water depth
    W <- water_depth(t)
    
    # Current bird abundance (depends on water and small fish)
    B <- bird_abundance(t, S, bird_response, half_sat)
    
    # Large fish removal/emigration
    if (W < 10) {
      L <- 0  # Complete removal below threshold
    } else {
      # Emigration rate increases as depth approaches 10cm
      emigration_rate <- e_max * exp(-e_shape * (W - 10))
    }
    
    # Bird predation increases with both bird abundance and concentration
    concentration <- ifelse(W > 10, 50 / W, 5)
    bird_predation <- b * concentration * B
    
    # Small fish dynamics
    dS <- r_s * S * (1 - S/K_s) -           # Logistic growth
      bird_predation * S -               # Bird predation
      a * S * L                          # Large fish predation
    
    # Large fish dynamics
    if (W < 10) {
      dL <- -L  # Rapid removal
    } else {
      dL <- c * a * S * L -                  # Growth from eating small fish
        d * L -                          # Natural mortality
        emigration_rate * L              # Emigration
    }
    
    return(list(c(dS, dL)))
  })
}

# Parameters
params <- c(
  r_s = 0.5,         # Small fish growth rate
  K_s = 1000,        # Small fish carrying capacity
  b = 0.001,         # Bird predation base rate
  a = 0.005,        # Large fish predation rate # a = 0.0005 -> stable fish
  c = 0.3,           # Conversion efficiency
  d = 0.1,           # Large fish mortality
  e_max = 0.5,       # Maximum emigration rate
  e_shape = 0.3,     # Emigration response shape
  bird_response = 50, # Maximum additional birds from fish availability
  half_sat = 200     # Half-saturation for bird numerical response
)

# Initial conditions
initial_state <- c(S = 500, L = 50)

# Time vector (2 years)
times <- seq(0, 730, by = 1)

# Solve ODE
output <- ode(y = initial_state, times = times, func = fish_model, parms = params)

# Calculate bird abundance for all times (now depends on small fish too)
birds <- sapply(1:length(times), function(i) {
  bird_abundance(times[i], output[i, 2], params["bird_response"], params["half_sat"])
})

# Calculate water for all times
water <- sapply(times, water_depth)

# Create single plot with all variables
par(mar = c(5, 5, 4, 5))

# Plot small fish
plot(output[, 1], output[, 2], type = "l", col = "blue", lwd = 2,
     xlab = "Time (days)", ylab = "Small fish / Large fish / Bird abundance",
     ylim = c(0, max(output[, 2], output[, 3], birds) * 1.1))

# Add large fish
lines(output[, 1], output[, 3], col = "red", lwd = 2)

# Add birds
lines(times, birds, col = "darkgreen", lwd = 2)

# Add water depth on secondary axis
par(new = TRUE)
plot(times, water, type = "l", col = "cyan", lwd = 2,
     xaxt = "n", yaxt = "n", xlab = "", ylab = "")
axis(4)
mtext("Water depth (cm)", side = 4, line = 3)
abline(h = 10, lty = 2, col = "gray")

# Add legend
legend("topleft", legend = c("Small fish", "Large fish", "Birds", "Water depth"),
       col = c("blue", "red", "darkgreen", "cyan"), lwd = 2, bty = "n")













