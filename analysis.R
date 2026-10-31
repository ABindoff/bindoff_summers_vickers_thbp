# Cognitive trajectories in the Tasmanian Healthy Brain Project: main analysis
#
# Fits the four cognitive-trajectory models reported in the manuscript and prints the
# coefficients and evidence ratios reported there, beside the published values.
#
# DATA. If data/episodic.rds, working.rds, executive.rds and language.rds are present
# (the study data, supplied under agreement), they are used unchanged and the fits are
# identical to the manuscript's. Otherwise the synthetic data in data/synthetic/ are used:
# the analysis is the same but every number will differ from the manuscript.
#
# Run from this folder: source("analysis.R"). Each model is 4 chains x 10,000 iterations.
#
# The fitted models are left in `fits` (fits$episodic, fits$working, fits$executive,
# fits$language) for inspection, e.g.
#   model_results(fits$episodic)                     # full report with convergence checks
#   plot(fits$episodic)                              # trace plots
#   pp_check(fits$episodic)                          # posterior predictive check
#   posterior::summarise_draws(fits$episodic$draws)  # R-hat, ESS for every parameter

library(smoothbp)

# smoothbp 0.2.8 fixed an error in the random-effects sampler that affects every
# random-intercept fit made with earlier 0.2.x releases.
stopifnot("smoothbp >= 0.2.8 is required" = packageVersion("smoothbp") >= "0.2.8")
options(smoothbp.cores = parallel::detectCores())

# ---------------------------------------------------------------------------
# Data: one data frame per domain, one row per participant-wave.
#   value                 composite score (z; higher is better in every domain)
#   years_since_baseline  elapsed time since the baseline assessment, whole years
#   *_s                   z-scored copies of continuous baseline covariates
# ---------------------------------------------------------------------------
DOMAINS <- c(episodic  = "Episodic memory",    working  = "Working memory",
             executive = "Executive function", language = "Language processing")

real_files <- file.path("data", paste0(names(DOMAINS), ".rds"))

if (all(file.exists(real_files))) {
  data_source <- "study data"
  domains <- lapply(real_files, readRDS)
} else {
  data_source <- "SYNTHETIC data: results will differ from the manuscript"
  domains <- lapply(file.path("data/synthetic", c("epis.csv", "wm.csv", "ef.csv", "lang.csv")),
                    function(f) transform(read.csv(f), IDCode = factor(IDCode)))
}
names(domains) <- names(DOMAINS)
message("Running on the ", data_source)

# ---------------------------------------------------------------------------
# Model (see Methods). Smooth piecewise regression with a transition point omega:
#   mu = b0 + b1 * d + delta1 * d * logistic(rho * d),   d = years_since_baseline - omega
# b1 is the pre-transition slope and b1 + delta1 the post-transition slope (z per year).
# Group enters b0, b1 and delta1; baseline covariates enter b0; random intercept per
# participant; rho fixed at 3. Identical specification in all four domains.
# ---------------------------------------------------------------------------
THBP_PRIORS <- smoothbp_priors(
  b0      = prior_normal(0, 2),
  b1      = prior_normal(0, 2),
  deltas  = prior_normal(0, 1),
  omega   = list(list("(Intercept)" = prior_normal(3, 2, lb = 1, ub = 6))),
  rho     = prior_normal(2, 2, lb = 1),   # not used: rho is fixed below
  sigma   = prior_invgamma(shape = 2, scale = 1),
  sigma_u = prior_invgamma(shape = 2, scale = 1)
)

fit_domain <- function(dat, seed = 42L) {
  smoothbp(
    formula = value ~ years_since_baseline,
    b0      = ~ 1 + Group + education_total_s + Age_t1_s + Sex +
                hads_anx_baseline + hads_dep_baseline +
                wtar_fsiq_baseline_s + drs_tot_raw_baseline_s +
                (1 | IDCode),
    b1      = ~ 1 + Group,
    deltas  = list(~ 1 + Group),
    omega   = list(~ 1),
    rho     = list(fixed(3)),
    data    = dat,
    priors  = THBP_PRIORS,
    chains  = 4L, iter = 10000L, warmup = 5000L,
    target_accept = 0.99, seed = seed
  )
}

fits <- lapply(domains, fit_domain)

# ---------------------------------------------------------------------------
# Results. Each quantity is a posterior mean [95% credible interval] of one parameter or
# a sum of parameters; slopes and contrasts are z per year. Group is coded Control/Experimental ("Education" in the manuscript).
# For the group contrasts the hypothesis is Education > Control: P(H) is the proportion of
# draws above zero and ER = P(H) / (1 - P(H)). ERs above 1000 rest on fewer than 20 of the
# 20,000 draws and are shown as "> 1000".
# ---------------------------------------------------------------------------
QUANTITIES <- list(
  list(label = "beta1, control (pre-transition slope)",       vars = "b1_(Intercept)",                                     digits = 2, er = FALSE),
  list(label = "beta1 + delta1, control (post-transition)",  vars = c("b1_(Intercept)", "delta1_(Intercept)"),             digits = 2, er = FALSE),
  list(label = "delta1 (change in slope)",                    vars = "delta1_(Intercept)",                                 digits = 2, er = FALSE),
  list(label = "omega (transition point, years)",             vars = "omega1_(Intercept)",                                 digits = 2, er = FALSE),
  list(label = "Group contrast: beta1",                       vars = "b1_GroupExperimental",                               digits = 3, er = TRUE),
  list(label = "Group contrast: beta1 + delta1",              vars = c("b1_GroupExperimental", "delta1_GroupExperimental"), digits = 3, er = TRUE)
)

# Published values (manuscript Results and Table 2; Supplementary Note 4 for the
# working-memory and executive-function control slopes), in the order of QUANTITIES.
MANUSCRIPT <- list(
  episodic  = c("0.10 [0.05, 0.15]", "-0.02 [-0.06, 0.01]", "-0.12 [-0.19, -0.05]", "2.85 [2.52, 3.18]",
                "0.069 [0.018, 0.122]; ER 228.9 [1.00]", "0.027 [-0.013, 0.067]; ER 10.2 [0.91]"),
  working   = c("0.03 [-0.02, 0.09]", "not reported", "-0.06 [-0.13, 0.02]", "2.49 [1.31, 3.59]",
                "0.015 [-0.043, 0.082]; ER 2.3 [0.69]", "-0.005 [-0.040, 0.030]; ER 0.6 [0.38]"),
  executive = c("0.08 [-0.01, 0.24]", "not reported", "-0.09 [-0.26, 0.04]", "2.78 [1.07, 5.65]",
                "0.005 [-0.131, 0.085]; ER 1.8 [0.64]", "0.010 [-0.085, 0.064]; ER 2.3 [0.70]"),
  language  = c("0.06 [-0.01, 0.18]", "-0.08 [-0.19, -0.04]", "-0.14 [-0.27, -0.04]", "2.23 [1.06, 5.17]",
                "0.027 [-0.065, 0.127]; ER 2.9 [0.74]", "0.071 [0.027, 0.181]; ER > 1000 [1.00]")
)

draws_of <- function(m, vars)
  rowSums(as.data.frame(posterior::as_draws_df(m$draws))[, vars, drop = FALSE])

summarise_quantity <- function(m, q) {
  x   <- draws_of(m, q$vars)
  out <- sprintf("%.*f [%.*f, %.*f]", q$digits, mean(x), q$digits, quantile(x, .025),
                 q$digits, quantile(x, .975))
  if (q$er) {
    p  <- mean(x > 0)
    er <- p / (1 - p)
    out <- sprintf("%s; ER %s [%.2f]", out, if (er > 1000) "> 1000" else sprintf("%.1f", er), p)
  }
  out
}

results <- lapply(names(DOMAINS), function(d)
  data.frame(Quantity   = vapply(QUANTITIES, `[[`, "", "label"),
             Manuscript = MANUSCRIPT[[d]],
             This_run   = vapply(QUANTITIES, summarise_quantity, "", m = fits[[d]])))
names(results) <- names(DOMAINS)

cat("\nData:", data_source, "\n")
cat("Participants:", length(unique(unlist(lapply(domains, function(d) as.character(d$IDCode))))),
    "| observations:", paste(sprintf("%s %d", DOMAINS, sapply(domains, nrow)), collapse = ", "), "\n\n")
op <- options(width = max(getOption("width"), 130))
for (d in names(DOMAINS)) {
  cat("
", DOMAINS[[d]], "
", sep = "")
  print(results[[d]], row.names = FALSE, right = FALSE)
}
options(op)
