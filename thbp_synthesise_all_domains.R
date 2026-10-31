# Joint synthetic release of the THBP cognitive domains (language, episodic memory,
# working memory, executive function).
#
# One synthesis, one release: baseline characteristics are drawn once and all four
# domains are synthesised conditional on them, so the released files cannot be
# averaged or linked against each other to recover real records.
#
# Input: four long data frames, one row per person per assessment, sharing IDCode,
# years_since_baseline and the baseline covariates, each with its domain score in `value`.
# Output: data/synthetic/<domain>.csv, long format, ready for the analysis code.
#
# The real data never leave this script's session. Do not commit anything except the
# synthetic CSVs and this script. `models = FALSE` keeps the fitted trees (which contain
# real records) out of the synds object.

library(synthpop)
library(dplyr)
library(tidyr)

# ---------------------------------------------------------------------------
# 1. Inputs: edit the object names and covariate names to match yours
# ---------------------------------------------------------------------------
domains <- list(lang = language, epis = episodic, wm = working, ef = executive)

baseline_vars <- c("Group", "Sex", "Age_t1", "education_total", "wtar_fsiq_baseline",
                   "drs_tot_raw_baseline", "hads_anx_baseline", "hads_dep_baseline")
keys <- c("Group", "Sex", "Age_t1", "education_total")   # what an outsider could know

# Column that identifies an assessment occasion uniquely within a person and is the SAME
# occasion across domains (e.g. a wave number or assessment date). years_since_baseline
# cannot do this: it is in whole years, so two sessions in the same year tie.
VISIT_KEY <- "session"

SEED     <- 2026L
MINBUCKET <- 10L
out_dir  <- "data/synthetic"

# ---------------------------------------------------------------------------
# 2. Visits: the union of assessment times across domains, per person
#    (a person may lack a score in one domain at a visit; that stays NA)
# ---------------------------------------------------------------------------
long_all <- bind_rows(lapply(names(domains), function(d)
  domains[[d]] |> ungroup() |>
    transmute(IDCode, occ = .data[[VISIT_KEY]], years_since_baseline, value) |>
    mutate(domain = d)))
stopifnot(!anyDuplicated(long_all[c("IDCode", "domain", "occ")]))   # one score per domain per occasion

# One row per occasion. The same occasion must carry the same years value in every domain.
occ_years <- long_all |> distinct(IDCode, occ, years_since_baseline)
stopifnot(!anyDuplicated(occ_years[c("IDCode", "occ")]))

visits <- occ_years |>
  arrange(IDCode, occ) |>
  group_by(IDCode) |>
  mutate(visit = row_number(),
         gap   = years_since_baseline - lag(years_since_baseline)) |>   # can be 0 for same-year sessions
  ungroup()
K <- max(visits$visit)

# Baseline covariates: one row per person, taken from whichever domain table has them
base <- bind_rows(lapply(domains, function(x) x |> ungroup() |> select(IDCode, all_of(baseline_vars)))) |>
  distinct(IDCode, .keep_all = TRUE)
stopifnot(!anyDuplicated(base$IDCode))

# ---------------------------------------------------------------------------
# 3. Wide file: baseline, then per visit: gap_k, lang_k, epis_k, wm_k, ef_k
# ---------------------------------------------------------------------------
scores_wide <- long_all |>
  left_join(visits |> select(IDCode, occ, visit), by = c("IDCode", "occ")) |>
  select(IDCode, visit, domain, value) |>
  pivot_wider(names_from = c(domain, visit), values_from = value, names_sep = "_")
gaps_wide <- visits |>
  filter(visit > 1) |>
  select(IDCode, visit, gap) |>
  pivot_wider(names_from = visit, values_from = gap, names_prefix = "gap_")

W_id <- base |> left_join(scores_wide, by = "IDCode") |> left_join(gaps_wide, by = "IDCode")

seq_vars <- c(baseline_vars,
              paste0(names(domains), "_1"),
              unlist(lapply(2:K, function(k) c(paste0("gap_", k), paste0(names(domains), "_", k)))))
for (v in setdiff(seq_vars, names(W_id))) W_id[[v]] <- NA_real_   # a domain absent at some visit index
W <- as.data.frame(W_id[, seq_vars])
cat("Wide file:", nrow(W), "people x", ncol(W), "columns\n")

# ---------------------------------------------------------------------------
# 4. Rules: dropout is monotone, and there are no scores at a visit that didn't happen
# ---------------------------------------------------------------------------
rules <- rvalues <- list()
for (k in 2:K) {
  g <- paste0("gap_", k)
  if (k > 2) { rules[[g]] <- paste0("is.na(gap_", k - 1, ")"); rvalues[[g]] <- NA }
  for (d in names(domains)) {
    v <- paste0(d, "_", k)
    rules[[v]] <- paste0("is.na(", g, ")"); rvalues[[v]] <- NA
  }
}

# ---------------------------------------------------------------------------
# 5. Predictor matrix: each domain score is predicted from baseline, the visit gaps,
#    and that domain's own earlier scores only. This keeps the models from learning
#    individual patterns across ~60 columns. Gaps may use everything before them.
# ---------------------------------------------------------------------------
P <- matrix(0L, ncol(W), ncol(W), dimnames = list(seq_vars, seq_vars))
for (i in seq_along(seq_vars)) {
  v <- seq_vars[i]
  earlier <- seq_vars[seq_len(i - 1)]
  if (i == 1) next
  if (v %in% baseline_vars || grepl("^gap_", v)) {
    P[v, earlier] <- 1L
  } else {
    d <- sub("_\\d+$", "", v)
    own <- earlier[grepl(paste0("^", d, "_\\d+$"), earlier)]
    P[v, intersect(earlier, c(baseline_vars, grep("^gap_", earlier, value = TRUE), own))] <- 1L
  }
}

# ---------------------------------------------------------------------------
# 6. Synthesise once
# ---------------------------------------------------------------------------
sds <- syn(W,
           method           = c("sample", rep("cart", ncol(W) - 1)),
           visit.sequence   = seq_along(W),
           predictor.matrix = P,
           rules = rules, rvalues = rvalues,
           smoothing        = "spline",
           cart.minbucket   = MINBUCKET,
           models           = FALSE,
           seed             = SEED)

# ---------------------------------------------------------------------------
# 7. Utility and disclosure, on the joint file
# ---------------------------------------------------------------------------
print(compare(sds, W, vars = c(baseline_vars, paste0(names(domains), "_1"))))
print(utility.gen(sds, W))

# Exact keys (no banding of age or education); scores banded into 5 groups as targets
md <- multi.disclosure(sds, W, keys = keys, targets = paste0(names(domains), "_1"),
                       ngroups_targets = 5)
print(md)

sds_safe <- sdc(sds, W, keys = keys, rm.replicated.uniques = TRUE)
# Re-check the file that will actually be released
md_safe <- multi.disclosure(sds_safe, W, keys = keys, targets = paste0(names(domains), "_1"),
                            ngroups_targets = 5)
print(md_safe)

# ---------------------------------------------------------------------------
# 8. Back to one long file per domain, with derived variables rebuilt
# ---------------------------------------------------------------------------
to_long <- function(S, d) {
  S <- S |> mutate(IDCode = row_number())                       # new IDs, unrelated to real ones
  gaps <- S |> select(IDCode, starts_with("gap_")) |>
    pivot_longer(-IDCode, names_to = "visit", names_prefix = "gap_", values_to = "gap") |>
    mutate(visit = as.integer(visit)) |>
    bind_rows(tibble(IDCode = S$IDCode, visit = 1L, gap = 0)) |>
    filter(!is.na(gap)) |>
    arrange(IDCode, visit) |>
    group_by(IDCode) |> mutate(years_since_baseline = cumsum(gap)) |> ungroup()
  sc <- S |> select(IDCode, all_of(baseline_vars), matches(paste0("^", d, "_\\d+$"))) |>
    pivot_longer(matches(paste0("^", d, "_\\d+$")), names_to = "visit",
                 names_prefix = paste0(d, "_"), values_to = "value") |>
    mutate(visit = as.integer(visit)) |>
    filter(!is.na(value))
  sc |> inner_join(gaps |> select(IDCode, visit, years_since_baseline), by = c("IDCode", "visit")) |>
    arrange(IDCode, visit) |>                                   # visit order, not years (years can tie)
    group_by(IDCode) |>
    mutate(session = row_number()) |>                           # completed assessments on THIS domain
    ungroup() |>
    mutate(Age_v = Age_t1 + years_since_baseline,
           log_session = log(session),
           Age_v_s  = (Age_v  - mean(Age_v))  / sd(Age_v),
           Age_t1_s = (Age_t1 - mean(Age_t1)) / sd(Age_t1),
           education_total_s      = (education_total      - mean(education_total))      / sd(education_total),
           wtar_fsiq_baseline_s   = (wtar_fsiq_baseline   - mean(wtar_fsiq_baseline))   / sd(wtar_fsiq_baseline),
           drs_tot_raw_baseline_s = (drs_tot_raw_baseline - mean(drs_tot_raw_baseline)) / sd(drs_tot_raw_baseline)) |>
    select(-visit)
}

S <- sds_safe$syn
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
for (d in names(domains)) {
  syn_d <- to_long(S, d)
  write.csv(syn_d, file.path(out_dir, paste0(d, ".csv")), row.names = FALSE)
  cat(sprintf("%-4s: %d synthetic people, %d rows\n", d, length(unique(syn_d$IDCode)), nrow(syn_d)))
}

# Record what was done, for the README
writeLines(c(
  sprintf("synthpop %s, seed %d, method cart (first variable sampled), smoothing spline, cart.minbucket %d",
          as.character(packageVersion("synthpop")), SEED, MINBUCKET),
  sprintf("Joint synthesis of %s from one wide file (%d columns); replicated uniques on keys %s removed.",
          paste(names(domains), collapse = ", "), ncol(W), paste(keys, collapse = ", ")),
  "These data are synthetic and will not reproduce the published estimates."),
  file.path(out_dir, "SYNTHESIS_NOTES.txt"))
