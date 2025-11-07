# predglmm Copilot Instructions

## Project Overview
This is an R package for computing marginal predictions (expected values) from generalized linear mixed models (GLMMs). The package integrates over random effects using Gaussian quadrature to calculate predictions on the natural scale.

## Architecture

### Three-Layer Prediction System
1. **`pred_base()`**: Core computation engine that handles the mathematical integration
   - Uses Gaussian-Hermite quadrature via `mvQuad` package for numerical integration
   - Supports custom inverse link functions OR named links ("identity", "log", "sqrt")
   - Returns marginal expectations by integrating over random effects
   
2. **`pred_glmer()`**: Adapter for `lme4::glmer` fitted models
   - Extracts fixed effects (`X`, `beta`) and random effects structure (`Z`, `Lambda`) using `lme4::getME()`
   - Computes conditional variance `sigma` from random effects: `sigma = sqrt(rowSums((Z %*% Lambda)^2))`
   - Delegates to `pred_base()` for integration

3. **`pred_glmmtmb()`**: Adapter for `glmmTMB::glmmTMB` fitted models
   - Reconstructs random effects using `lme4::mkReTrms()` with `reorder.terms = FALSE` for glmmTMB compatibility
   - Manually fills sparse precision matrix `Lambdat` from variance components
   - Uses `Matrix::forceSymmetric()` to ensure symmetry after filling upper triangle

## Key Dependencies & Integration Points

- **mvQuad**: Gaussian quadrature nodes/weights (`createNIGrid`, `getNodes`, `getWeights`)
- **lme4**: Model extraction (`getME`, `mkReTrms`, `findbars`) used by BOTH adapters
- **glmmTMB**: Variance extraction (`VarCorr`) and model structure
- **Matrix**: Sparse matrix operations for random effects covariance

## Critical Conventions

### Function Parameters
- **Inconsistency Alert**: `pred_base()` roxygen docs don't match implementation:
  - Docs say `(X, Beta, sigma, link, num_nodes)` 
  - Code implements `(eta, sigma, link_name, num_nodes, inv_link)`
  - This needs reconciliation when editing

### Random Effects Variance Calculation
- Standard deviation of random effects (`sigma`) is calculated as row sums of squared random effects design:
  ```r
  H <- Z %*% Lambda  # Random effects contribution per observation
  sigma <- sqrt(rowSums(H^2))  # Standard deviation for each observation
  ```

### Link Function Handling
- For analytical solutions (log, sqrt, identity): use `link_name` parameter
- For custom/complex links: pass function via `inv_link` and use quadrature
- Warning system alerts if both or neither are provided

## Development Workflow

### Documentation with roxygen2
- Use `roxygen2::roxygenize()` or RStudio Build → Document to regenerate `.Rd` files
- Only `pred_base()` is currently exported (`@export` tag)
- NAMESPACE is auto-generated - don't edit manually

### Package Building
```r
# In R console
devtools::document()      # Update documentation
devtools::load_all()      # Test changes without install
devtools::check()         # R CMD check
devtools::install()       # Install locally
```

### Testing Approach
- No formal test suite currently exists
- When adding tests, verify predictions match for:
  - Identity link (should equal linear predictor)
  - Log link with analytical formula vs quadrature
  - Models with/without random effects (sigma = 0 edge case)

## Common Patterns

### Adding New Link Functions
To support a new link with analytical marginal expectation:
1. Add branch in `pred_base()` after line 47
2. Compute marginal expectation analytically if possible (e.g., `exp(eta + 0.5*sigma^2)` for log)
3. Fall back to quadrature via `inv_link` for complex cases

### Supporting New GLMM Packages
Follow the adapter pattern from `pred_glmer()` and `pred_glmmtmb()`:
1. Extract linear predictor `eta = X %*% beta`
2. Compute per-observation random effects SD `sigma`
3. Call `pred_base(eta, sigma, link_name = ...)`

## Known Issues
- Documentation mismatch between roxygen and implementation for `pred_base()`
- Only `pred_base()` is exported; adapter functions are package-internal
- No dependency declarations in DESCRIPTION file (needs Imports: lme4, glmmTMB, mvQuad, Matrix)
