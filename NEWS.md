# predglmm 0.2.0

* New `boot_pred()`: parametric bootstrap confidence intervals for marginal
  predictions and for functions of them, for glmmTMB and lme4 fits. The
  statistic receives the data the model was fit to with a column `pred` of
  marginal predictions added. Refits on the boundary (singular fits) are
  kept and counted, not dropped.
* `pred_nlme()` is removed. For linear `nlme::lme` fits,
  `predict(fit, level = 0)` gives the marginal prediction.
* Offsets are included in the linear predictors of `pred_glmer()` and
  `pred_glmmtmb()` (conditional and zero-inflation components).
* `pred_glmmtmb()` supports rank-deficient designs.
* `pred_base()` errors on an unsupported `link_name` without `inv_link`
  instead of warning and returning `NA`.

# predglmm 0.1.0

* First version: `pred_base()`, `pred_glmer()`, `pred_glmmtmb()`,
  `pred_nlme()`.
