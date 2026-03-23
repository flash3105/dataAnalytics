# Compute theoretical ACF for AR(2): phi1 = 1.5, phi2 = -0.75
acf_values <- ARMAacf(ar = c(1.5, -0.75), lag.max = 24)

# Plot
plot(
  x    = 0:24,
  y    = acf_values,
  type = "h",                          # vertical lines (spike plot)
  lwd  = 2,
  col  = ifelse(acf_values >= 0, "steelblue", "green"),
  xlab = "Lag (k)",
  ylab = expression(rho(k)),
  main = expression("Theoretical ACF of AR(2): " ~
                      X[t] == 1.5*X[t-1] - 0.75*X[t-2] + epsilon[t]),
  ylim = c(-1, 1),
  xaxt = "n"
)

# Axes and reference line
axis(1, at = 0:24)
abline(h = 0, lty = 2, col = "gray40")

# Mark lag 0
points(0, 1, pch = 19, col = "black")


#  Question 1(b) - Population Spectrum of AR(2) Process
#  phi1 = 1.5, phi2 = -0.75, sigma^2 = 1

# --- Parameters ---
phi1 <- 1.5
phi2 <- -0.75
sig2 <- 1

# --- Grid: omega in [0, pi] in steps of pi/24 ---
k     <- 0:24
omega <- k * pi / 24

# --- Population spectrum formula ---
denom   <- 1 + phi1^2 + phi2^2 -
  2 * phi1 * (1 - phi2) * cos(omega) -
  2 * phi2  * cos(2 * omega)
f_omega <- sig2 / (2 * pi * denom)

# --- Print table of values ---
cat("=== Spectral Values ===\n")
cat(sprintf("%-5s %-12s %-12s\n", "k", "omega", "f(omega)"))
cat(rep("-", 32), "\n", sep = "")
for (i in seq_along(k)) {
  cat(sprintf("%-5d %-12.4f %-12.5f\n", k[i], omega[i], f_omega[i]))
}

# --- Cyclical frequency ---
discriminant <- phi1^2 + 4 * phi2
cat("\n=== Cyclical Frequency ===\n")
cat("Discriminant (phi1^2 + 4*phi2):", discriminant, "\n")

if (discriminant < 0) {
  cat("=> Complex roots: damped sinusoidal ACF and spectral peak\n\n")
  R     <- sqrt(-phi2)
  Theta <- acos(phi1 / (2 * R))
  T_cyc <- 2 * pi / Theta
  cat("R  (damping factor) :", round(R,     4), "\n")
  cat("omega* (radians)    :", round(Theta, 4), "= pi /", round(pi / Theta, 2), "\n")
  cat("Period T (lags)     :", round(T_cyc, 4), "\n")
}

# ============================================================
#  Plot
# ============================================================
png("plot1b.png", width = 800, height = 500, res = 120)

plot(
  x    = omega,
  y    = f_omega,
  type = "l",
  lwd  = 2.5,
  col  = "steelblue",
  xlab = expression(omega ~ "(radians)"),
  ylab = expression(f(omega)),
  main = expression("Population Spectrum  —  AR(2):  " ~
                      phi[1] == 1.5 ~ "," ~ phi[2] == -0.75 ~ "," ~ sigma^2 == 1),
  xaxt = "n",
  ylim = c(0, max(f_omega) * 1.12)
)

# x-axis with pi fractions
axis(1,
     at     = c(0, pi/6, pi/4, pi/3, pi/2, 2*pi/3, 3*pi/4, 5*pi/6, pi),
     labels = expression(0, pi/6, pi/4, pi/3, pi/2,
                         2*pi/3, 3*pi/4, 5*pi/6, pi))

# Grid and zero line
abline(h   = 0,     lty = 2, col = "gray70")
abline(v   = Theta, lty = 2, lwd = 2, col = "tomato")

# Dot at spectral peak
f_star <- sig2 / (2 * pi *
                    (1 + phi1^2 + phi2^2
                     - 2 * phi1 * (1 - phi2) * cos(Theta)
                     - 2 * phi2  * cos(2 * Theta)))
points(Theta, f_star, pch = 19, cex = 1.6, col = "tomato")

# Annotation
text(Theta + 0.08, f_star,
     labels = expression(omega^"*" == pi/6),
     col = "tomato", adj = 0, cex = 0.9)
text(Theta + 0.08, f_star * 0.88,
     labels = paste0("f(pi/6) = ", round(f_star, 3)),
     col = "tomato", adj = 0, cex = 0.8)
text(Theta + 0.08, f_star * 0.76,
     labels = "T = 12 lags",
     col = "tomato", adj = 0, cex = 0.8)

# Legend
legend("topright",
       legend = c(expression(f(omega)),
                  expression(omega^"*" == pi/6 ~ " (cyclical frequency)")),
       col    = c("steelblue", "tomato"),
       lty    = c(1, 2),
       lwd    = c(2.5, 2),
       bty    = "n",
       cex    = 0.85)

dev.off()
cat("\nPlot saved to plot1b.png\n")
