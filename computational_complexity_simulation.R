###############################################################################
# Computational complexity simulation for TARP on simulated Arabidopsis-like
# biallelic genotype data.
#
# Reference implementation: uploaded arabidopsis_tarp_full_study.Rmd
#   screen_set()      : adaptive delta and Bernoulli RIS probabilities
#   tarp_omega()      : three-point sparse projection distribution
#   nig_mspe()        : NIG posterior-mean prediction on compressed features
#   agg_tarp_pred()   : B random projections; m ~ U[2 log(p), 3 n_train / 4]
#
# This is a standalone simulation. It does NOT load the real Arabidopsis files,
# alter the main R Markdown, or replace TARP with a different algorithm.
#
# Run from ~/arabidopsis_project (or set PROJECT_DIR below), then source this
# file in RStudio. QUICK_RUN is TRUE by default for the first end-to-end check.
###############################################################################

## ----------------------------- USER SETTINGS -------------------------------
QUICK_RUN <- FALSE

# Full-study design. QUICK_RUN overrides these with a small test design.
FULL_N <- c(200L, 400L)
FULL_P <- c(2000L, 10000L, 50000L)
FULL_RHO <- c(0, 0.30, 0.60, 0.85)   # latent within-block Gaussian correlation
FULL_LAYOUTS <- c("distributed", "clustered")
FULL_REPS <- 3L

# TARP settings copied from the R Markdown implementation.
TARP_B <- 100L
TARP_PSI <- 0.25
BLOCK_SIZE <- 50L
N_CAUSAL <- 12L

# Comparator settings. GBLUP is intentionally fixed-lambda, not tuned GBLUP.
SIS_D_MAX <- 2000L
PCA_K <- 100L
GBLUP_LAMBDA_FIXED <- 1.0

# Project root. The recommended location is ~/arabidopsis_project.
PROJECT_DIR <- path.expand("~/arabidopsis_project")
if (!dir.exists(PROJECT_DIR)) PROJECT_DIR <- getwd()

## -------------------------- CONFIGURE THE RUN -------------------------------
if (!dir.exists(PROJECT_DIR)) {
  stop("Project directory does not exist: ", PROJECT_DIR,
       "\nSet PROJECT_DIR near the top of this script.")
}
setwd(PROJECT_DIR)

if (!requireNamespace("Matrix", quietly = TRUE)) {
  stop("Package 'Matrix' is required. Install it with install.packages('Matrix').")
}
if (requireNamespace("RhpcBLASctl", quietly = TRUE)) {
  try(RhpcBLASctl::blas_set_num_threads(1L), silent = TRUE)
}

if (isTRUE(QUICK_RUN)) {
  N_GRID <- 100L
  P_GRID <- 600L
  RHO_GRID <- c(0, 0.60)
  LAYOUT_GRID <- c("distributed", "clustered")
  N_REPS <- 1L
  B_USE <- 3L
  K_PCA_USE <- 20L
  N_CAUSAL_USE <- 8L
  SIS_D_MAX_USE <- 40L
} else {
  N_GRID <- FULL_N
  P_GRID <- FULL_P
  RHO_GRID <- FULL_RHO
  LAYOUT_GRID <- FULL_LAYOUTS
  N_REPS <- FULL_REPS
  B_USE <- TARP_B
  K_PCA_USE <- PCA_K
  N_CAUSAL_USE <- N_CAUSAL
  SIS_D_MAX_USE <- SIS_D_MAX
}

if (any(N_GRID < 40L) || any(P_GRID < 100L)) {
  stop("Use n >= 40 and p >= 100 so the training split and TARP dimension range are valid.")
}
if (any(RHO_GRID < 0 | RHO_GRID >= 1)) {
  stop("Each latent block-correlation rho must be in [0, 1).")
}
if (TARP_PSI <= 0 || TARP_PSI >= 0.5) {
  stop("TARP_PSI must be strictly between 0 and 0.5.")
}

OUT_DIR <- file.path(PROJECT_DIR, "computational_study")
RESULT_DIR <- file.path(OUT_DIR, "results")
CHECKPOINT_DIR <- file.path(OUT_DIR, "checkpoints")
LOG_DIR <- file.path(OUT_DIR, "logs")
TEMP_DIR <- file.path(OUT_DIR, "temporary")
for (d in c(OUT_DIR, RESULT_DIR, CHECKPOINT_DIR, LOG_DIR, TEMP_DIR)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

RUN_LOG <- file.path(LOG_DIR, "run.log")
FAILURE_LOG <- file.path(LOG_DIR, "failures.csv")
MANIFEST_FILE <- file.path(RESULT_DIR, "case_manifest.csv")
RAW_FILE <- file.path(RESULT_DIR, "raw_results.csv")
SUMMARY_FILE <- file.path(RESULT_DIR, "case_summary.csv")
META_FILE <- file.path(RESULT_DIR, "benchmark_metadata.txt")
RDS_FILE <- file.path(RESULT_DIR, "benchmark_results.rds")

## ----------------------------- UTILITIES ------------------------------------
now_string <- function() format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")

log_message <- function(...) {
  msg <- paste0("[", now_string(), "] ", paste0(..., collapse = ""))
  cat(msg, "\n")
  cat(msg, "\n", file = RUN_LOG, append = TRUE)
  invisible(msg)
}

atomic_save_rds <- function(object, path) {
  tmp <- tempfile(pattern = "checkpoint_", tmpdir = TEMP_DIR, fileext = ".rds")
  on.exit(if (file.exists(tmp)) unlink(tmp), add = TRUE)
  saveRDS(object, tmp)
  if (file.exists(path)) unlink(path)
  if (!file.rename(tmp, path)) {
    stop("Could not move temporary RDS file into place: ", path)
  }
  invisible(path)
}

atomic_write_csv <- function(x, path) {
  tmp <- tempfile(pattern = "table_", tmpdir = TEMP_DIR, fileext = ".csv")
  on.exit(if (file.exists(tmp)) unlink(tmp), add = TRUE)
  utils::write.csv(x, tmp, row.names = FALSE, na = "")
  if (file.exists(path)) unlink(path)
  if (!file.rename(tmp, path)) {
    stop("Could not move temporary CSV file into place: ", path)
  }
  invisible(path)
}

append_failure <- function(case_id, message) {
  row <- data.frame(timestamp = now_string(), case_id = case_id,
                    error = gsub("[\r\n]+", " ", message),
                    stringsAsFactors = FALSE)
  utils::write.table(row, FAILURE_LOG, sep = ",", row.names = FALSE,
                     col.names = !file.exists(FAILURE_LOG), append = file.exists(FAILURE_LOG),
                     qmethod = "double")
}

## ------------------------ SIMULATED GENOTYPES --------------------------------
# Each SNP is coded 0/2, matching the homozygous biallelic coding described in
# the Arabidopsis R Markdown. Within a block, latent Gaussian variables share a
# common factor. rho is the latent-normal correlation parameter; the realized
# Pearson correlation between binary genotype columns is generally different.
simulate_genotypes <- function(n, p, rho, block_size = BLOCK_SIZE) {
  maf <- stats::runif(p, min = 0.05, max = 0.50)
  X <- matrix(0, nrow = n, ncol = p)
  blocks <- split(seq_len(p), ceiling(seq_len(p) / block_size))

  for (idx in blocks) {
    nb <- length(idx)
    shared <- stats::rnorm(n)
    eps <- matrix(stats::rnorm(n * nb), nrow = n, ncol = nb)
    latent <- if (rho == 0) eps else {
      sqrt(rho) * shared + sqrt(1 - rho) * eps
    }
    thresholds <- stats::qnorm(1 - maf[idx])
    X[, idx] <- 2 * (sweep(latent, 2L, thresholds, FUN = ">"))
  }
  storage.mode(X) <- "double"
  X
}

choose_causal_snps <- function(p, n_causal, layout, block_size = BLOCK_SIZE) {
  blocks <- split(seq_len(p), ceiling(seq_len(p) / block_size))
  eligible <- which(lengths(blocks) >= n_causal)
  if (layout == "clustered" && length(eligible) > 0L) {
    b <- eligible[sample.int(length(eligible), 1L)]
    b_members <- blocks[[b]]
    return(b_members[sample.int(length(b_members), n_causal, replace = FALSE)])
  }

  if (layout == "distributed") {
    eligible <- which(lengths(blocks) >= 1L)
    use_n <- min(n_causal, length(eligible))
    chosen_blocks <- eligible[sample.int(length(eligible), use_n, replace = FALSE)]
    return(vapply(chosen_blocks, function(b) {
      members <- blocks[[b]]
      members[sample.int(length(members), 1L)]
    }, integer(1)))
  }

  # Fallback when p is smaller than one complete block of causal SNPs.
  sample.int(p, n_causal, replace = FALSE)
}

empirical_block_correlation <- function(X, block_size = BLOCK_SIZE,
                                        max_blocks = 20L, markers_per_block = 8L) {
  blocks <- split(seq_len(ncol(X)), ceiling(seq_len(ncol(X)) / block_size))
  eligible <- which(lengths(blocks) >= 2L)
  if (length(eligible) == 0L) return(NA_real_)
  n_use <- min(max_blocks, length(eligible))
  chosen <- eligible[sample.int(length(eligible), n_use, replace = FALSE)]
  vals <- numeric(0)
  for (b in chosen) {
    members <- blocks[[b]]
    k <- min(markers_per_block, length(members))
    selected <- members[sample.int(length(members), k, replace = FALSE)]
    cc <- suppressWarnings(stats::cor(X[, selected, drop = FALSE]))
    vv <- cc[upper.tri(cc)]
    vals <- c(vals, vv[is.finite(vv)])
  }
  if (length(vals) == 0L) NA_real_ else mean(vals)
}

simulate_response <- function(X, causal_snps, noise_sd = 1) {
  s <- length(causal_snps)
  signs <- sample(c(-1, 1), s, replace = TRUE)
  beta <- signs / sqrt(s)
  Xc <- X[, causal_snps, drop = FALSE]
  cm <- colMeans(Xc)
  cs <- apply(Xc, 2L, stats::sd)
  cs[!is.finite(cs) | cs == 0] <- 1
  Xcz <- sweep(sweep(Xc, 2L, cm, "-"), 2L, cs, "/")
  signal <- as.numeric(Xcz %*% beta)
  signal_sd <- stats::sd(signal)
  if (!is.finite(signal_sd) || signal_sd <= 0) {
    stop("The simulated causal signal has zero variance.")
  }
  signal <- signal / signal_sd
  y <- signal + stats::rnorm(nrow(X), sd = noise_sd)
  list(y = as.numeric(y), signal = signal, beta = beta)
}

## ------------------------ EXACT TARP BUILDING BLOCKS -------------------------
# Equivalent to the R Markdown screen_set() correlation calculation.
# The t-statistic-to-r transformation in screen_set() is algebraically the
# ordinary Pearson correlation for a one-predictor regression with an intercept.
marginal_correlations <- function(X_train, y_train) {
  y_centered <- y_train - mean(y_train)
  n_train <- nrow(X_train)
  x_sum <- colSums(X_train)
  x_ss <- colSums(X_train * X_train) - (x_sum * x_sum) / n_train
  xy <- as.numeric(crossprod(X_train, y_centered)) -
    x_sum * sum(y_centered) / n_train
  denom <- sqrt(pmax(x_ss, 0) * sum(y_centered^2))
  r <- xy / denom
  r[!is.finite(r)] <- 0
  pmax(-1, pmin(1, r))
}

ris_probabilities <- function(r, p_total, n_train) {
  delta <- max(0, (1 + log(p_total / n_train)) / 2)
  score <- abs(r)^delta
  max_score <- max(score, na.rm = TRUE)
  if (!is.finite(max_score) || max_score <= 0) {
    q <- rep(0, length(r))
    q[which.max(abs(r))] <- 1
  } else {
    q <- score / max_score
    q[!is.finite(q)] <- 0
  }
  list(q = q, delta = delta)
}

# Sparse-storage implementation of the exact three-point distribution in
# tarp_omega(): +1/sqrt(2*psi) with probability psi, the negative value with
# probability psi, and zero with probability 1-2*psi. The number of nonzeros
# is Binomial(p_gamma*m, 2*psi), with uniformly sampled positions and signs.
tarp_omega_sparse <- function(p_gamma, m, psi = TARP_PSI) {
  total <- as.double(p_gamma) * as.double(m)
  if (p_gamma < 1L || m < 1L || total > .Machine$integer.max) {
    stop("Invalid or too-large TARP projection dimensions.")
  }
  nnz <- stats::rbinom(1L, size = as.integer(total), prob = 2 * psi)
  if (nnz == 0L) {
    return(Matrix::sparseMatrix(i = integer(0), j = integer(0), x = numeric(0),
                                dims = c(p_gamma, m)))
  }
  positions <- sample.int(as.integer(total), size = nnz, replace = FALSE)
  i <- ((positions - 1L) %% p_gamma) + 1L
  j <- ((positions - 1L) %/% p_gamma) + 1L
  val <- 1 / sqrt(2 * psi)
  signs <- ifelse(stats::runif(nnz) < 0.5, -val, val)
  Matrix::sparseMatrix(i = i, j = j, x = signs, dims = c(p_gamma, m))
}

# Exact NIG posterior-mean calculation from nig_mspe(), with predictions returned.
nig_predict <- function(Z_tr, Z_te, y_centered, y_mean) {
  m <- ncol(Z_tr)
  mu <- solve(diag(m) + crossprod(Z_tr), crossprod(Z_tr, y_centered))
  as.numeric(Z_te %*% mu) + y_mean
}

# Equivalent to agg_tarp_pred(), except Omega is stored as a sparse Matrix.
# Every draw uses m ~ Uniform[2 log(p_total), 3*n_train/4] and predictions are
# averaged over B independent projection matrices.
agg_tarp_predict <- function(X_tr, X_te, y_centered, y_mean, p_total,
                             B = 100L, psi = TARP_PSI) {
  p_gamma <- ncol(X_tr)
  n_train <- nrow(X_tr)
  m_range <- c(2 * log(p_total), 3 * n_train / 4)
  if (m_range[2] < m_range[1]) {
    stop("Invalid projection-dimension interval: 3*n_train/4 < 2*log(p_total).")
  }
  pred_sum <- numeric(nrow(X_te))
  used_dims <- integer(B)
  for (b in seq_len(B)) {
    m_b <- max(1L, as.integer(round(stats::runif(1L, m_range[1], m_range[2]))))
    used_dims[b] <- m_b
    Om <- tarp_omega_sparse(p_gamma, m_b, psi)
    Z_tr <- as.matrix(X_tr %*% Om)
    Z_te <- as.matrix(X_te %*% Om)
    pred_sum <- pred_sum + nig_predict(Z_tr, Z_te, y_centered, y_mean)
  }
  list(pred = pred_sum / B, projection_dims = used_dims)
}

## -------------------------- MODEL PREPARATION -------------------------------
standardize_train_test <- function(X_train, X_test) {
  cm <- colMeans(X_train)
  xc <- sweep(X_train, 2L, cm, "-")
  cs <- sqrt(colSums(xc * xc) / max(1, nrow(X_train) - 1L))
  cs[!is.finite(cs) | cs == 0] <- 1
  list(train = sweep(xc, 2L, cs, "/"),
       test = sweep(sweep(X_test, 2L, cm, "-"), 2L, cs, "/"))
}

nig_mspe_from_X <- function(X_tr, X_te, y_train, y_test) {
  y_mean <- mean(y_train)
  pred <- nig_predict(X_tr, X_te, y_train - y_mean, y_mean)
  mean((y_test - pred)^2)
}

## ---------------------------- COMPARATORS -----------------------------------
# SIS: deterministic top-d absolute marginal correlations, followed by the same
# closed-form NIG posterior-mean regression used by the R Markdown TARP fit.
run_sis_nig <- function(X_train_raw, X_test_raw, y_train, y_test, r,
                        d_max = SIS_D_MAX_USE) {
  n_train <- nrow(X_train_raw)
  d_use <- min(as.integer(d_max), max(10L, floor(3 * n_train / 4)), ncol(X_train_raw))
  t_rank <- system.time({
    ord <- order(abs(r), decreasing = TRUE)
    selected <- ord[seq_len(d_use)]
  })[["elapsed"]]
  t_std <- system.time({
    xs <- standardize_train_test(X_train_raw[, selected, drop = FALSE],
                                 X_test_raw[, selected, drop = FALSE])
  })[["elapsed"]]
  t_fit <- system.time({
    mspe <- nig_mspe_from_X(xs$train, xs$test, y_train, y_test)
  })[["elapsed"]]
  list(mspe = mspe, elapsed = t_rank + t_std + t_fit,
       time_rank = t_rank, time_standardize = t_std, time_fit = t_fit,
       n_selected = d_use)
}

# Exact PCA using the n_train x n_train Gram matrix, as in gram_pca_basis().
run_exact_gram_pca <- function(X_train, X_test, y_train, y_test, k_max = PCA_K) {
  y_mean <- mean(y_train)
  y_centered <- y_train - y_mean
  t_eigen <- system.time({
    gram <- tcrossprod(X_train)
    eg <- eigen(gram, symmetric = TRUE)
    vals <- pmax(eg$values, 0)
    positive <- which(vals > max(vals, 1) * 1e-10)
    k <- min(as.integer(k_max), nrow(X_train) - 1L, length(positive))
    if (k < 1L) stop("Exact Gram PCA found no positive components.")
    U <- eg$vectors[, seq_len(k), drop = FALSE]
    lam <- vals[seq_len(k)]
    V <- sweep(crossprod(X_train, U), 2L, sqrt(lam), "/")
    Z_tr <- sweep(U, 2L, sqrt(lam), "*")
    Z_te <- X_test %*% V
  })[["elapsed"]]
  t_fit <- system.time({
    pred <- nig_predict(Z_tr, Z_te, y_centered, y_mean)
    mspe <- mean((y_test - pred)^2)
  })[["elapsed"]]
  list(mspe = mspe, elapsed = t_eigen + t_fit, time_pca = t_eigen,
       time_fit = t_fit, k = k)
}

# Fixed-lambda GBLUP computational reference. The kernel is K = XX'/p after
# training-standardization; lambda is deliberately not cross-validated/tuned.
run_gblup_fixed <- function(X_train, X_test, y_train, y_test,
                            lambda = GBLUP_LAMBDA_FIXED) {
  y_mean <- mean(y_train)
  y_centered <- y_train - y_mean
  p <- ncol(X_train)
  t_fit <- system.time({
    K_tr <- tcrossprod(X_train) / p
    K_te <- (X_test %*% t(X_train)) / p
    alpha <- solve(K_tr + diag(lambda, nrow(K_tr)), y_centered)
    pred <- as.numeric(K_te %*% alpha) + y_mean
    mspe <- mean((y_test - pred)^2)
  })[["elapsed"]]
  list(mspe = mspe, elapsed = t_fit, time_fit = t_fit, lambda = lambda)
}

## -------------------------- ONE SIMULATION CASE -----------------------------
run_one_case <- function(n, p, rho, layout, replicate, case_index) {
  data_seed <- as.integer(100000L + case_index * 1000L + replicate)
  set.seed(data_seed)
  X <- simulate_genotypes(n, p, rho, BLOCK_SIZE)
  empirical_r <- empirical_block_correlation(X, BLOCK_SIZE)
  causal <- choose_causal_snps(p, min(N_CAUSAL_USE, p), layout, BLOCK_SIZE)
  sim_y <- simulate_response(X, causal, noise_sd = 1)
  y <- sim_y$y

  set.seed(data_seed + 101L)
  n_train <- round(0.8 * n)
  train_rows <- sample.int(n, size = n_train, replace = FALSE)
  test_rows <- setdiff(seq_len(n), train_rows)
  X_train_raw <- X[train_rows, , drop = FALSE]
  X_test_raw <- X[test_rows, , drop = FALSE]
  y_train <- y[train_rows]
  y_test <- y[test_rows]
  rm(X, sim_y)
  invisible(gc(FALSE))
  y_mean <- mean(y_train)
  baseline_mspe <- mean((y_test - y_mean)^2)

  # The marginal scores are computed once and re-used by RIS and SIS.
  t_score <- system.time({
    r <- marginal_correlations(X_train_raw, y_train)
  })[["elapsed"]]

  # TARP/RIS: preserve the exact adaptive delta and q_j from screen_set().
  set.seed(data_seed + 202L)
  t_ris <- system.time({
    ris <- ris_probabilities(r, p_total = p, n_train = n_train)
    included <- which(stats::rbinom(p, size = 1L, prob = ris$q) == 1L)
    if (length(included) == 0L) included <- which.max(abs(r))
  })[["elapsed"]]
  if (length(included) < 1L) stop("RIS returned no predictors after fallback.")

  t_std_tarp <- system.time({
    X_tarp <- standardize_train_test(X_train_raw[, included, drop = FALSE],
                                     X_test_raw[, included, drop = FALSE])
  })[["elapsed"]]
  set.seed(data_seed + 303L)
  t_agg <- system.time({
    tarp <- agg_tarp_predict(X_tarp$train, X_tarp$test,
                             y_train - y_mean, y_mean, p_total = p,
                             B = B_USE, psi = TARP_PSI)
  })[["elapsed"]]
  tarp_mspe <- mean((y_test - tarp$pred)^2)
  rm(X_tarp)

  # SIS comparator; uses the same training-only marginal scores.
  sis <- run_sis_nig(X_train_raw, X_test_raw, y_train, y_test, r,
                     d_max = SIS_D_MAX_USE)

  # Full-p standardization is measured as part of both full-design comparators.
  t_std_full <- system.time({
    X_full <- standardize_train_test(X_train_raw, X_test_raw)
  })[["elapsed"]]
  rm(X_train_raw, X_test_raw, r)
  invisible(gc(FALSE))
  pca <- run_exact_gram_pca(X_full$train, X_full$test, y_train, y_test,
                            k_max = K_PCA_USE)
  gb <- run_gblup_fixed(X_full$train, X_full$test, y_train, y_test,
                        lambda = GBLUP_LAMBDA_FIXED)

  scenario <- if (rho == 0) "independent" else paste0("block_rho_", sprintf("%.2f", rho))
  common <- data.frame(
    scenario = scenario, rho_latent = rho, signal_layout = layout,
    n = n, p = p, replicate = replicate, data_seed = data_seed,
    n_train = n_train, n_test = length(test_rows),
    n_causal = length(causal), baseline_mspe = baseline_mspe,
    empirical_within_block_r = empirical_r, stringsAsFactors = FALSE
  )

  mkrow <- function(method, mspe, elapsed, time_score = NA_real_, time_select = NA_real_,
                    time_standardize = NA_real_, time_fit = NA_real_,
                    n_selected = NA_integer_, n_projection = NA_integer_,
                    mean_projection_dim = NA_real_, delta = NA_real_,
                    lambda_fixed = NA_real_, k_pca = NA_integer_) {
    cbind(common, data.frame(
      method = method, mspe = as.numeric(mspe), elapsed_seconds = as.numeric(elapsed),
      time_score_seconds = as.numeric(time_score), time_select_seconds = as.numeric(time_select),
      time_standardize_seconds = as.numeric(time_standardize), time_fit_seconds = as.numeric(time_fit),
      n_selected = as.integer(n_selected), n_projections = as.integer(n_projection),
      mean_projection_dim = as.numeric(mean_projection_dim), delta = as.numeric(delta),
      lambda_fixed = as.numeric(lambda_fixed), k_pca = as.integer(k_pca),
      stringsAsFactors = FALSE
    ))
  }

  rows <- rbind(
    mkrow("TARP_RIS_RP_aggregated", tarp_mspe,
          t_score + t_ris + t_std_tarp + t_agg,
          time_score = t_score, time_select = t_ris,
          time_standardize = t_std_tarp, time_fit = t_agg,
          n_selected = length(included), n_projection = B_USE,
          mean_projection_dim = mean(tarp$projection_dims), delta = ris$delta),
    mkrow("SIS_top_d_plus_NIG", sis$mspe,
          t_score + sis$elapsed, time_score = t_score,
          time_select = sis$time_rank,
          time_standardize = sis$time_standardize, time_fit = sis$time_fit,
          n_selected = sis$n_selected),
    mkrow("Exact_Gram_PCA_plus_NIG", pca$mspe,
          t_std_full + pca$elapsed, time_standardize = t_std_full,
          time_fit = pca$time_pca + pca$time_fit, k_pca = pca$k),
    mkrow("GBLUP_fixed_lambda", gb$mspe,
          t_std_full + gb$elapsed, time_standardize = t_std_full,
          time_fit = gb$time_fit, lambda_fixed = gb$lambda),
    mkrow("Mean_only_baseline", baseline_mspe, 0)
  )

  list(rows = rows, diagnostics = list(
    causal_snps = causal, n_ris_selected = length(included),
    ris_expected_selected = sum(ris$q), delta = ris$delta,
    projection_dims = tarp$projection_dims,
    projection_nonzero_fraction_expected = 2 * TARP_PSI,
    block_size = BLOCK_SIZE, rho_latent = rho
  ))
}

## ------------------------- DESIGN AND CHECKPOINTS ----------------------------
design <- expand.grid(
  n = N_GRID, p = P_GRID, rho = RHO_GRID,
  signal_layout = LAYOUT_GRID, replicate = seq_len(N_REPS),
  KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE
)
design <- design[order(design$rho, design$signal_layout, design$n,
                        design$p, design$replicate), , drop = FALSE]
rownames(design) <- NULL

design$scenario <- ifelse(design$rho == 0, "independent",
                           paste0("block_rho_", sprintf("%.2f", design$rho)))
design$case_id <- sprintf("%s_%s_n%05d_p%07d_rep%03d",
                          design$scenario, design$signal_layout,
                          design$n, design$p, design$replicate)
design$checkpoint <- file.path(CHECKPOINT_DIR, paste0(design$case_id, ".rds"))
design$status <- "pending"
design$started_at <- NA_character_
design$finished_at <- NA_character_
design$error <- NA_character_

manifest <- design

save_manifest <- function() atomic_write_csv(manifest, MANIFEST_FILE)

# Reconcile the manifest with already-completed checkpoints before starting.
for (i in seq_len(nrow(manifest))) {
  cp <- manifest$checkpoint[i]
  if (file.exists(cp)) {
    ok <- tryCatch({
      z <- readRDS(cp)
      is.list(z) && is.data.frame(z$rows) && nrow(z$rows) > 0L
    }, error = function(e) FALSE)
    if (isTRUE(ok)) {
      manifest$status[i] <- "completed"
      manifest$finished_at[i] <- now_string()
    } else {
      manifest$status[i] <- "pending"
      log_message("Invalid checkpoint will be recomputed: ", manifest$case_id[i])
    }
  }
}
save_manifest()

## ----------------------------- RUN THE STUDY --------------------------------
log_message("Starting computational study. QUICK_RUN=", QUICK_RUN,
            "; cases=", nrow(manifest), "; B=", B_USE,
            "; project=", PROJECT_DIR)
log_message("TARP settings: adaptive delta=max(0,(1+log(p/n_train))/2), psi=",
            TARP_PSI, "; m in [2 log(p), 3*n_train/4].")

for (i in seq_len(nrow(manifest))) {
  if (manifest$status[i] == "completed") next

  case_id <- manifest$case_id[i]
  manifest$status[i] <- "running"
  manifest$started_at[i] <- now_string()
  manifest$error[i] <- NA_character_
  save_manifest()

  log_message("Case ", i, "/", nrow(manifest), " starting: ", case_id)
  case_start <- Sys.time()

  result <- tryCatch({
    ans <- run_one_case(
      n = as.integer(manifest$n[i]), p = as.integer(manifest$p[i]),
      rho = as.numeric(manifest$rho[i]), layout = manifest$signal_layout[i],
      replicate = as.integer(manifest$replicate[i]), case_index = i
    )
    ans$rows$case_id <- case_id
    checkpoint_object <- list(
      case_id = case_id, completed_at = now_string(), rows = ans$rows,
      diagnostics = ans$diagnostics
    )
    atomic_save_rds(checkpoint_object, manifest$checkpoint[i])
    list(ok = TRUE, error = NA_character_)
  }, error = function(e) {
    list(ok = FALSE, error = conditionMessage(e))
  })

  elapsed_case <- as.numeric(difftime(Sys.time(), case_start, units = "secs"))
  if (isTRUE(result$ok)) {
    manifest$status[i] <- "completed"
    manifest$finished_at[i] <- now_string()
    manifest$error[i] <- NA_character_
    log_message("Case completed: ", case_id, " (", round(elapsed_case, 2), " sec)")
  } else {
    manifest$status[i] <- "failed"
    manifest$finished_at[i] <- now_string()
    manifest$error[i] <- result$error
    append_failure(case_id, result$error)
    log_message("CASE FAILED: ", case_id, " | ", result$error)
  }
  save_manifest()
}

## -------------------------- REBUILD PERSISTENT OUTPUTS -----------------------
completed_rows <- list()
completed_diagnostics <- list()
for (i in seq_len(nrow(manifest))) {
  cp <- manifest$checkpoint[i]
  if (!file.exists(cp)) next
  z <- tryCatch(readRDS(cp), error = function(e) NULL)
  if (is.null(z) || !is.data.frame(z$rows)) next
  completed_rows[[length(completed_rows) + 1L]] <- z$rows
  completed_diagnostics[[length(completed_diagnostics) + 1L]] <- z$diagnostics
}

if (length(completed_rows) > 0L) {
  raw_results <- do.call(rbind, completed_rows)
  rownames(raw_results) <- NULL
  group_cols <- c("scenario", "rho_latent", "signal_layout", "n", "p", "method")
  group_key <- do.call(paste, c(raw_results[group_cols], sep = "\r"))
  groups <- split(raw_results, group_key)
  case_summary <- do.call(rbind, lapply(groups, function(z) {
    data.frame(
      scenario = z$scenario[1], rho_latent = z$rho_latent[1],
      signal_layout = z$signal_layout[1], n = z$n[1], p = z$p[1],
      method = z$method[1], n_replicates = nrow(z),
      mean_MSPE = mean(z$mspe, na.rm = TRUE),
      sd_MSPE = if (nrow(z) > 1L) stats::sd(z$mspe, na.rm = TRUE) else NA_real_,
      se_MSPE = if (nrow(z) > 1L) stats::sd(z$mspe, na.rm = TRUE) / sqrt(nrow(z)) else NA_real_,
      mean_elapsed_seconds = mean(z$elapsed_seconds, na.rm = TRUE),
      median_elapsed_seconds = stats::median(z$elapsed_seconds, na.rm = TRUE),
      mean_n_selected = if (all(is.na(z$n_selected))) NA_real_ else mean(z$n_selected, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  }))
  rownames(case_summary) <- NULL
} else {
  raw_results <- data.frame()
  case_summary <- data.frame()
}

atomic_write_csv(raw_results, RAW_FILE)
atomic_write_csv(case_summary, SUMMARY_FILE)
atomic_write_csv(manifest, MANIFEST_FILE)

config <- list(
  quick_run = QUICK_RUN, n_grid = N_GRID, p_grid = P_GRID,
  rho_grid = RHO_GRID, signal_layouts = LAYOUT_GRID, replicates = N_REPS,
  B = B_USE, psi = TARP_PSI, block_size = BLOCK_SIZE,
  n_causal = N_CAUSAL_USE, sis_d_max = SIS_D_MAX_USE,
  pca_k = K_PCA_USE, gblup_fixed_lambda = GBLUP_LAMBDA_FIXED,
  delta_rule = "max(0, (1 + log(p / n_train)) / 2)",
  projection_dimension = "round(Uniform(2*log(p), 3*n_train/4))",
  genotype_coding = "biallelic homozygous calls coded 0/2",
  dependence = "Gaussian one-factor blocks; rho is latent-normal correlation",
  tarp_model = "RIS + three-point random projection + NIG posterior mean"
)

atomic_save_rds(list(
  raw_results = raw_results, case_summary = case_summary,
  case_manifest = manifest, config = config,
  diagnostics = completed_diagnostics
), RDS_FILE)

metadata_lines <- c(
  "Arabidopsis TARP computational complexity simulation",
  paste("Completed at:", now_string()),
  paste("Project directory:", PROJECT_DIR),
  paste("R version:", R.version.string),
  paste("Platform:", R.version$platform),
  paste("QUICK_RUN:", QUICK_RUN),
  paste("Cases in design:", nrow(manifest)),
  paste("Completed cases:", sum(manifest$status == "completed")),
  paste("Failed cases:", sum(manifest$status == "failed")),
  paste("Pending cases:", sum(manifest$status == "pending")),
  paste("n grid:", paste(N_GRID, collapse = ", ")),
  paste("p grid:", paste(P_GRID, collapse = ", ")),
  paste("latent rho grid:", paste(RHO_GRID, collapse = ", ")),
  paste("signal layouts:", paste(LAYOUT_GRID, collapse = ", ")),
  paste("replicates per design point:", N_REPS),
  paste("TARP aggregation B:", B_USE),
  paste("TARP psi:", TARP_PSI),
  paste("RIS delta:", config$delta_rule),
  paste("TARP projection dimension:", config$projection_dimension),
  "NIG posterior mean: solve(I + crossprod(Z_tr), crossprod(Z_tr, y_centered)).",
  "GBLUP comparator: fixed lambda only; not nested-tuned or optimised.",
  "SIS comparator: top-d absolute training marginal correlations + the same NIG posterior-mean fit.",
  "Exact PCA comparator: eigen-decomposition of the n_train x n_train Gram matrix + NIG fit.",
  "Projection matrix is stored as a sparse Matrix. At psi=0.25, its expected nonzero fraction is 0.50; sparse storage is exact but is not guaranteed to use less memory than dense storage.",
  "rho is the latent Gaussian within-block correlation parameter; observed binary-genotype Pearson correlations are attenuated and depend on allele frequencies.",
  "All screening, standardisation and model fitting use training rows only."
)
writeLines(metadata_lines, META_FILE)

# Remove only temporary staging files. Persistent checkpoints, logs and results
# are deliberately left untouched so that a later source() resumes safely.
tmp_files <- list.files(TEMP_DIR, full.names = TRUE, all.files = TRUE, no.. = TRUE)
if (length(tmp_files)) unlink(tmp_files, recursive = TRUE, force = TRUE)

log_message("Study run finished. Completed=", sum(manifest$status == "completed"),
            "; failed=", sum(manifest$status == "failed"),
            "; pending=", sum(manifest$status == "pending"))
log_message("Results: ", RESULT_DIR)

if (any(manifest$status != "completed")) {
  warning("Some cases are not complete. See case_manifest.csv and logs/failures.csv; source the script again to resume.")
} else {
  message("All simulation cases completed successfully.")
}

###############################################################################
# END
###############################################################################
