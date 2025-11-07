# predglmm

Marginal predictions (expected values on the response scale) for generalized 
linear mixed models.

Note: This package has undergone only minimal testing and comes with no warranty!

## Overview

When fitting GLMMs, predictions that account for random effects uncertainty are 
often needed. This package computes marginal expectations by integrating over 
the random effects distribution, using either analytical expressions (when available) 
or Gaussian-Hermite quadrature.

## Installation

```r
# Install from GitHub
# devtools::install_github("koekvall/predglmm")
```

## Usage

```r
library(lme4)
library(glmmTMB)
library(predglmm)

# Poisson GLMM (closed form expressions for predictions)
fit_tmb <- glmmTMB(count ~ mined + (1|site),
  family=poisson, data=Salamanders)
fit_lme4 <- glmer(count ~ mined + (1|site),
  family=poisson, data=Salamanders)
  
# Get marginal predictions on the response scale
pred_tmb <- pred_glmmtmb(fit_tmb)
pred_lme4 <- pred_glmer(fit_lme4)

# Compare predictions
pred_tmb[1:10]
pred_lme4[1:10]

# Predictions from usual prediction functions which do not give marginal means:
predict(fit_tmb, type = "response")[1:10] # REs evaluated at predicted values
predict(fit_tmb, type = "response", re.form = NA)[1:10] # REs "set to zero"


# Logistic GLMM example (using quadrature)
fit_tmb <- glmmTMB(I(count > 2) ~ mined + (1|site),
  family=binomial, data=Salamanders)

fit_lme4 <- glmer(I(count > 2) ~ mined + (1|site),
  family=binomial, data=Salamanders)
  
# Get marginal predictions on the response scale
pred_tmb <- pred_glmmtmb(fit_tmb)
pred_lme4 <- pred_glmer(fit_lme4)

# Compare predictions
pred_tmb[1:10]
pred_lme4[1:10]

# Predictions from usual prediction functions which do not give marginal means:
predict(fit_tmb, type = "response")[1:10] # REs evaluated at predicted values
predict(fit_tmb, type = "response", re.form = NA)[1:10] # REs "set to zero"
```

Works with models from:
- **lme4**: `pred_glmer()` for `glmer()` and `lmer()` fits
- **glmmTMB**: `pred_glmmtmb()` for `glmmTMB()` fits

## Supported Link Functions

**Analytical** (fast): `identity`, `log`, `sqrt`  
**Numerical** (via quadrature): any custom inverse link function

## How It Works

For a GLMM with linear predictor η = Xβ and random effects b ~ N(0, Ψ):

- **Linear models**: Returns E[Y] = Xβ (no integration needed)
- **Nonlinear models**: Computes E[g⁻¹(η + Zb)] by integrating over b

For supported links, uses analytical formulas (e.g., E[exp(η + Zb)] = exp(η + σ²/2) for log link). Otherwise, uses Gaussian-Hermite quadrature.

## License

MIT
