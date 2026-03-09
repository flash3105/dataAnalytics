# ============================================================
# TRA880 Assignment 1 - Full R Code
# 2 Broke Girls (2BG) Cupcake Business Time Series
# ============================================================

# ---- READ DATA ----
data <- read.table("2BG.txt", na.strings = ".")
colnames(data) <- c("Season", "Episode", "Amount")
head(data)
nrow(data)

# ---- PART (a): MISSING VALUE IMPUTATION ----
which(is.na(data$Amount))

# Linear interpolation for missing values
data$Amount[117] <- (data$Amount[116] + data$Amount[118]) / 2
data$Amount[124] <- (data$Amount[123] + data$Amount[125]) / 2

# Verify no more NAs
which(is.na(data$Amount))

# ---- EXPLORATORY PLOTS ----
# Time plot
plot(data$Episode, data$Amount,
     type = "o",
     pch = 16,
     col = "blue",
     lwd = 2,
     xlab = "Episode",
     ylab = "Final Tally ($)",
     main = "Final Tally Over Episodes (2 Broke Girls)",
     cex.main = 1.2)
grid()

# Boxplot
boxplot(data$Amount)

# Outlier detection via IQR
Q1      <- quantile(data$Amount, 0.25, na.rm = TRUE)
Q3      <- quantile(data$Amount, 0.75, na.rm = TRUE)
IQR_val <- IQR(data$Amount, na.rm = TRUE)
upper   <- Q3 + 1.5 * IQR_val
lower   <- Q1 - 1.5 * IQR_val
data[data$Amount > upper | data$Amount < lower, ]

# Time plot with outliers highlighted
plot(data$Episode, data$Amount,
     type = "l",
     lwd = 2,
     col = "navy",
     xlab = "Episode",
     ylab = "Final Tally ($)",
     main = "Final Tally Over Episodes")
points(data$Episode[data$Amount > upper],
       data$Amount[data$Amount > upper],
       col = "red",
       pch = 19)
grid()

# ---- PART (b): BOX-COX TRANSFORMATION ----
library(forecast)
library(e1071)

# Convert to time series object
w <- ts(data$Amount, start = 1, frequency = 1)

# Check minimum value
min(w)
sum(w <= 0)

# Shift to make all values strictly positive
shift   <- abs(min(w)) + 1     # shift = 15
w_shifted <- w + shift

# Estimate lambda
lambda <- BoxCox.lambda(w_shifted)
lambda  # report this value

# Apply Box-Cox transformation
w_transformed <- BoxCox(w_shifted, lambda)

# Compare skewness before and after
skewness(w)              # original
skewness(w_transformed)  # transformed

# Visual comparison
par(mfrow = c(1, 2))
hist(w, main = "Original", col = "lightblue", xlab = "Amount")
hist(w_transformed, main = "Box-Cox Transformed",
     col = "lightgreen", xlab = "Transformed Amount")

# ---- TIME PLOT OF TRANSFORMED SERIES ----
par(mfrow = c(1, 1))
plot(w_transformed,
     type = "o",
     pch = 16,
     col = "darkgreen",
     lwd = 2,
     xlab = "Episode",
     ylab = "Transformed Amount",
     main = "Box-Cox Transformed Final Tally (2 Broke Girls)",
     cex.main = 1.2)
abline(h = mean(w_transformed), col = "red", lty = 2, lwd = 1.5)
grid()

# ---- PART (c): ARMA(p,q) IDENTIFICATION ----
# ACF and PACF plots
par(mfrow = c(1, 2))
acf(w_transformed,
    main = "ACF - Transformed Series",
    col = "darkgreen",
    lwd = 2)
pacf(w_transformed,
     main = "PACF - Transformed Series",
     col = "darkblue",
     lwd = 2)

# Automated AIC model selection
p_max <- 3
q_max <- 3

aic_matrix <- matrix(NA, nrow = p_max + 1, ncol = q_max + 1,
                     dimnames = list(paste0("p=", 0:p_max),
                                     paste0("q=", 0:q_max)))

for (p in 0:p_max) {
  for (q in 0:q_max) {
    tryCatch({
      fit <- arima(w_transformed, order = c(p, 0, q), method = "ML")
      aic_matrix[p+1, q+1] <- AIC(fit)
    }, error = function(e) NA)
  }
}

# Display AIC table
print(round(aic_matrix, 2))

# Identify best model
best <- which(aic_matrix == min(aic_matrix, na.rm = TRUE), arr.ind = TRUE)
cat("Best model: ARMA(", best[1]-1, ",", best[2]-1, ") with AIC =",
    round(min(aic_matrix, na.rm = TRUE), 2))

# ---- PART (d): FIT MA(2) MODEL ----
ma2_fit <- arima(w_transformed, order = c(0, 0, 2), method = "ML")
ma2_fit  # view coefficients and standard errors

# Residual diagnostics
par(mfrow = c(2, 2))

# 1. Residuals over time
plot(residuals(ma2_fit),
     type = "o",
     col = "darkblue",
     main = "Residuals over Time",
     ylab = "Residuals")
abline(h = 0, col = "red", lty = 2)

# 2. ACF of residuals
acf(residuals(ma2_fit), main = "ACF of Residuals")

# 3. Histogram of residuals
hist(residuals(ma2_fit),
     col = "lightblue",
     main = "Histogram of Residuals",
     xlab = "Residuals",
     prob = TRUE)
curve(dnorm(x,
            mean = mean(residuals(ma2_fit)),
            sd   = sd(residuals(ma2_fit))),
      add = TRUE, col = "red", lwd = 2)

# 4. QQ plot
qqnorm(residuals(ma2_fit))
qqline(residuals(ma2_fit), col = "red", lwd = 2)

# Formal tests
Box.test(residuals(ma2_fit), lag = 10, type = "Ljung-Box")
shapiro.test(residuals(ma2_fit))

