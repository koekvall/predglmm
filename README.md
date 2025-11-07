# predglmm

Marginal predictions (expected values on the response scale) for generalized 
linear mixed models.

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
library(predglmm)

# Fit a Poisson GLMM
fit <- glmer(y ~ x + (1|group), data = mydata, family = poisson())

# Get marginal predictions on the response scale
predictions <- pred_glmer(fit)
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
