#################################
# Load required libraries
# - ggplot2: visualization
# - gamlss: fitting flexible regression models
# - gamlss.dist: access to distribution families (e.g., Gamma)
#################################

library(ggplot2)
library(gamlss)
library(gamlss.dist)

#################################
# Simulate data
# N    = sample size
# age  = uniform random ages between 5 and 15
# ability = latent ability term (normally distributed)
# mu(age, ability) = expected RT; decreases with age/ability (faster with development/skill)
# rt ~ Gamma(mean = mu, CV = sigma_true), generated via shape/scale
#################################

set.seed(1)
N <- 1500

age <- round(runif(N, 5, 15), 1)
tau <- 1
ability <- rnorm(N, 0, tau)

# Mean RT: strictly positive; lower for older/higher-ability
mu <- exp(2.2 - 0.4 * log(age - 4) - 0.5 * ability)  # seconds (arbitrary scale)

sigma_true <- 0.30                     # coefficient of variation (CV) for Gamma
shape <- 1 / sigma_true^2
scale <- mu * sigma_true^2             # ensures mean = shape*scale = mu

rt <- rgamma(N, shape = shape, scale = scale)

dx <- data.frame(
  id = 1:N,
  age = age,
  ability = ability,
  rt = rt
)

# Visualize raw data with scatterplot and smooth trend
ggplot(dx, aes(x = age, y = rt)) +
  geom_smooth() +
  geom_point(alpha = .2) + 
  theme(text=element_text(size=18)) +
  scale_x_continuous(breaks=1:100)

#################################
# Fit a Gamma regression model
# - Response: positive continuous RT
# - Predictor: smooth monotonic function of age (pbm ensures decreasing trend)
# - Family: GA (Gamma), appropriate for skewed positive outcomes
#################################

fit <- gamlss(
  rt ~ pbm(age, mono = "down"),
  family = GA,
  data   = dx
)
summary(fit)

plot(fit)
#################################
# Example of deriving age-specific norms
# - Predict model parameters for a given age (e.g., 8 years)
# - Use qGA to obtain the RT corresponding to a percentile (e.g., 15th)
#   Note: For RT, lower values typically indicate better/faster performance.
#################################

pA <- predictAll(fit, newdata = data.frame(age = 5:30))
plot(5:30,pA$mu)
qGA(p = seq(0.05,0.95,0.05), mu = pA$mu, sigma = pA$sigma)
pGA(q = 8, mu = pA$mu, sigma = pA$sigma)

#################################



