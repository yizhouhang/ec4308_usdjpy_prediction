library(HDeconometrics)
library(sandwich)
library(hdm)
library(tidyverse)

setwd("/Users/yizhouhang/Documents/Y4S1/EC4308/ec4308 project/src")
df_1m  <- read_csv("../data/lag1m.csv",  show_col_types = FALSE)
df_3m  <- read_csv("../data/lag3m.csv",  show_col_types = FALSE)
df_6m  <- read_csv("../data/lag6m.csv",  show_col_types = FALSE)
df_12m <- read_csv("../data/lag12m.csv", show_col_types = FALSE)

split_date <- as.Date("2012-01-01")

rmse <- function(pred, truth) {
  pred  <- as.numeric(pred); truth <- as.numeric(truth)
  stopifnot(length(pred) == length(truth))
  sqrt(mean((truth - pred)^2))
}
r2 <- function(pred, y) {
  pred <- as.numeric(pred); y <- as.numeric(y)
  stopifnot(length(pred) == length(y))
  1 - sum((y - pred)^2) / sum((y - mean(y))^2)
}

# OLS summary
ols_quick_summary <- function(fit, newX_te, y_tr_used, y_te) {
  s       <- summary(fit)
  pred_tr <- as.numeric(fitted(fit))                 # exactly rows lm() used
  pred_te <- as.numeric(predict(fit, newdata = newX_te))
  
  cat("\n=== OLS (compact) ===\n")
  cat(sprintf("R^2 (train) = %.3f   Adj R^2 = %.3f   Residual SE = %.5f\n",
              s$r.squared, s$adj.r.squared, s$sigma))
  cat(sprintf("F-stat = %.3f (p = %.4g)\n",
              s$fstatistic[1],
              pf(s$fstatistic[1], s$fstatistic[2], s$fstatistic[3], lower.tail = FALSE)))
  cat(sprintf("AIC = %.1f   BIC = %.1f   LogLik = %.1f\n",
              AIC(fit), BIC(fit), as.numeric(logLik(fit))))
  cat(sprintf("Train: RMSE = %.5f   R^2 = %.3f\n", rmse(pred_tr, y_tr_used), r2(pred_tr, y_tr_used)))
  cat(sprintf("Test : RMSE = %.5f   R^2 = %.3f\n",  rmse(pred_te, y_te),    r2(pred_te, y_te)))
}

prep_xy <- function(df, date_col = "date", y_col = "y", split_date) {
  
  # assume df is already sorted and has no NA
  
  feat_cols <- setdiff(names(df), c(date_col, y_col))
  
  tr_idx <- df[[date_col]] < split_date
  te_idx <- !tr_idx
  
  y_tr <- df[[y_col]][tr_idx]
  y_te <- df[[y_col]][te_idx]
  
  X_tr <- as.matrix(df[tr_idx, feat_cols, drop = FALSE])
  X_te <- as.matrix(df[te_idx, feat_cols, drop = FALSE])
  
  # zero-variance guard based on train only
  sd_tr <- apply(X_tr, 2, sd)
  keep  <- which(!is.na(sd_tr) & sd_tr > 0)
  X_tr  <- X_tr[, keep, drop = FALSE]
  X_te  <- X_te[, keep, drop = FALSE]
  
  # scale using train stats only (leak-free)
  mu <- colMeans(X_tr)
  sd <- apply(X_tr, 2, sd)
  scale_with <- function(M, center, scale) sweep(sweep(M, 2, center, "-"), 2, scale, "/")
  X_tr_s <- scale_with(X_tr, mu, sd)
  X_te_s <- scale_with(X_te, mu, sd)
  
  # safe column names
  safe_names <- make.names(colnames(X_tr_s), unique = TRUE)
  colnames(X_tr_s) <- safe_names
  colnames(X_te_s) <- safe_names
  
  list(
    X_tr_s = X_tr_s, X_te_s = X_te_s,
    y_tr = y_tr, y_te = y_te,
    df_tr = as.data.frame(X_tr_s, check.names = FALSE),
    df_te = as.data.frame(X_te_s, check.names = FALSE)
  )
}


# =========================================
# One-horizon pipeline (quiet, informative)
# =========================================
run_pipeline <- function(df, horizon_label) {
  cat("\n==============================\n")
  cat("Horizon:", horizon_label, "\n")
  
  # map target to y and drop original target column
  stopifnot("USDJPY_logreturn" %in% names(df))
  df <- df %>%
    mutate(y = USDJPY_logreturn) %>%
    select(date, y, everything(), -USDJPY_logreturn)
  
  D <- prep_xy(df, date_col = "date", y_col = "y", split_date = split_date)
  X_tr_s <- D$X_tr_s; X_te_s <- D$X_te_s
  y_tr   <- D$y_tr;   y_te   <- D$y_te
  df_tr  <- D$df_tr;  df_te  <- D$df_te
  
  # ---------------- ROLLING *WITHIN TRAIN* ----------------
  window_length <- 120   # e.g. 10 years of monthly data
  if (length(y_tr) > window_length) {
    roll_preds <- rep(NA, length(y_tr) - window_length)
    roll_nonzero_count <- rep(NA, length(y_tr) - window_length)
    
    for (i in seq(window_length, length(y_tr)-1)) {
      idx <- (i - window_length + 1):i
      X_roll <- X_tr_s[idx, , drop = FALSE]
      y_roll <- y_tr[idx]
      
      fit_roll <- ic.glmnet(X_roll, y_roll, family = "gaussian", alpha = 1, standardize = FALSE)
      
      b <- as.numeric(coef(fit_roll))
      beta0 <- b[1]
      beta  <- b[-1]
      
      # one-step-ahead prediction
      roll_preds[i - window_length + 1] <- sum(X_tr_s[i+1, ] * beta) + beta0
      
      # count # nonzero variables selected
      roll_nonzero_count[i - window_length + 1] <- sum(beta != 0)
    }
    
    # compute rolling train-OOS RMSE
    roll_y_true <- y_tr[(window_length+1):length(y_tr)]
    rolling_rmse <- rmse(roll_preds, roll_y_true)
    
    cat(sprintf("\nRolling LASSO (train-only) — Window=%d months\n", window_length))
    cat(sprintf("Rolling RMSE (within train OOS): %.5f\n", rolling_rmse))
    cat(sprintf("Avg # variables selected: %.2f\n", mean(roll_nonzero_count)))
    
  } else {
    cat("\nRolling skipped (train too short for chosen window).\n")
  }
  
  
  # baseline: 0 return (RW for returns)
  pred_rw_te <- rep(0, length(y_te))
  cat(sprintf("RW baseline — Test RMSE: %.5f  R^2: %.3f\n",
              rmse(pred_rw_te, y_te), r2(pred_rw_te, y_te)))
  
  # ------------------- OLS -------------------
  if (ncol(X_tr_s) < length(y_tr)) {
    fit_ols <- lm(y_tr ~ ., data = cbind.data.frame(y_tr = y_tr, df_tr), model = TRUE)
    ols_quick_summary(
      fit_ols,
      newX_te   = df_te,
      y_tr_used = model.response(model.frame(fit_ols)),
      y_te      = y_te
    )
  } else {
    cat("OLS skipped (p >= n on train).\n")
  }
  
  # ----- LASSO (IC) -----
  fit_lasso <- ic.glmnet(X_tr_s, y_tr, family = "gaussian", alpha = 1, standardize = FALSE)
  
  # grab coefficients as a numeric vector [beta0, beta_1..p]
  b_las  <- as.numeric(coef(fit_lasso))
  beta0  <- b_las[1]
  beta   <- b_las[-1]
  
  # manual prediction to guarantee matching lengths
  pred_te_lasso <- as.numeric(drop(X_te_s %*% beta + beta0))
  
  sel_idx   <- which(beta != 0)
  sel_names <- if (length(sel_idx)) colnames(X_tr_s)[sel_idx] else character(0)
  
  cat(sprintf("\nLASSO (IC) — Test RMSE: %.5f  R^2: %.3f  (selected = %d)\n",
              rmse(pred_te_lasso, y_te), r2(pred_te_lasso, y_te), length(sel_idx)))
  if (length(sel_names)) {
    cat("Selected variable(s):", paste(sel_names, collapse = ", "), "\n")
  } else {
    cat("Selected variable(s): <none>\n")
  }
  
  # ----- Elastic Net (alpha = 0.5) -----
  fit_en <- ic.glmnet(X_tr_s, y_tr, family = "gaussian", alpha = 0.5, standardize = FALSE)
  
  b_en   <- as.numeric(coef(fit_en))
  beta0e <- b_en[1]
  betae  <- b_en[-1]
  
  pred_te_en <- as.numeric(drop(X_te_s %*% betae + beta0e))
  en_sel_idx <- which(betae != 0)
  
  cat(sprintf("Elastic Net (alpha=0.5) — Test RMSE: %.5f  R^2: %.3f  (selected = %d)\n",
              rmse(pred_te_en, y_te), r2(pred_te_en, y_te), length(en_sel_idx)))
  
  
  # ------------------- Post-LASSO -------------------
  if (length(sel_idx) > 0 && length(sel_idx) < length(y_tr)) {
    fit_post <- lm(y_tr ~ ., data = cbind.data.frame(y_tr = y_tr, df_tr[, sel_idx, drop = FALSE]))
    cat("\nPost-LASSO (OLS on selected vars)\n")
    ols_quick_summary(
      fit_post,
      newX_te   = df_te[, sel_idx, drop = FALSE],
      y_tr_used = model.response(model.frame(fit_post)),
      y_te      = y_te
    )
  } else {
    cat("\nPost-LASSO skipped (no vars or too many vars selected).\n")
  }
  
  invisible(NULL)
}

# ===========================
# Run all horizons
# ===========================
# df_1m$USDJPY_logreturn etc. must already be the h-ahead return aligned at time t.
run_pipeline(df_1m,  "1-month ahead")
run_pipeline(df_3m,  "3-month ahead")
run_pipeline(df_6m,  "6-month ahead")
run_pipeline(df_12m, "12-month ahead")



# ---- Rolling LASSO sparsity diagnostic (train-only) ----
# For each rolling window (length = window), fit LASSO (IC) on y ~ X and
# record how many coefficients are non-zero. Uses only dates < split_date.
roll_lasso_counts <- function(df, split_date, 
                              target_col = "USDJPY_logreturn",
                              date_col   = "date",
                              window     = 120,
                              alpha      = 1) {
  stopifnot(all(c(target_col, date_col) %in% names(df)))
  df[[date_col]] <- as.Date(df[[date_col]])
  
  # Train-only slice (no leakage to test)
  df_tr <- df[df[[date_col]] < split_date, , drop = FALSE]
  
  # Map target -> y and keep features (everything else except date & y)
  df_tr <- df_tr |>
    dplyr::mutate(y = .data[[target_col]]) |>
    dplyr::select(dplyr::all_of(c(date_col, "y")), dplyr::everything(), -dplyr::all_of(target_col))
  
  feat_cols <- setdiff(names(df_tr), c(date_col, "y"))
  # If you *know* there are no NAs you can skip this; keeping it is safer:
  df_tr <- tidyr::drop_na(df_tr, dplyr::all_of(c("y", feat_cols)))
  
  # Matrices
  dates <- df_tr[[date_col]]
  y     <- df_tr$y
  X     <- as.matrix(df_tr[, feat_cols, drop = FALSE])
  
  n <- nrow(X)
  if (n <= window + 1) stop("Not enough train rows for the chosen window.")
  
  out <- vector("list", length = n - window)  # each element -> one forecast origin
  k_names <- colnames(X)
  
  for (t in (window + 1):n) {
    # Window: (t-window) ... (t-1)  -> forecast origin at dates[t]
    idx    <- (t - window):(t - 1)
    X_win  <- X[idx, , drop = FALSE]
    y_win  <- y[idx]
    
    # Scale-within-window (no leakage)
    sd_win <- apply(X_win, 2, sd)
    keep   <- which(!is.na(sd_win) & sd_win > 0)
    if (length(keep) == 0) {
      out[[t - window]] <- data.frame(date = dates[t], n_selected = 0L)
      next
    }
    X_win  <- X_win[, keep, drop = FALSE]
    mu     <- colMeans(X_win)
    sdv    <- apply(X_win, 2, sd)
    scale_with <- function(M, center, scale) sweep(sweep(M, 2, center, "-"), 2, scale, "/")
    X_win_s <- scale_with(X_win, mu, sdv)
    
    # Fit LASSO by IC (HDeconometrics::ic.glmnet)
    fit <- HDeconometrics::ic.glmnet(X_win_s, y_win, family = "gaussian",
                                     alpha = alpha, standardize = FALSE)
    b   <- as.numeric(coef(fit))  # [beta0, beta_1..p]
    nz  <- sum(b[-1] != 0)
    
    out[[t - window]] <- data.frame(date = dates[t], n_selected = as.integer(nz))
  }
  
  dplyr::bind_rows(out)
}

# ---- Plot helper ----
plot_roll_sparsity <- function(cnt_df, title = "Rolling LASSO sparsity (train only)") {
  stopifnot(all(c("date", "n_selected") %in% names(cnt_df)))
  ggplot2::ggplot(cnt_df, ggplot2::aes(x = date, y = n_selected)) +
    ggplot2::geom_line() +
    ggplot2::geom_smooth(method = "loess", span = 0.2, se = FALSE) +
    ggplot2::labs(x = "Forecast origin (within train)", y = "# non-zero coefficients", title = title) +
    ggplot2::theme_minimal()
}

# ===== Run it for each horizon (example: 1m) =====
# Assumes you already have df_1m, df_3m, df_6m, df_12m and split_date in memory.
cnt_1m <- roll_lasso_counts(df_1m, split_date, window = 120, alpha = 1)
print(plot_roll_sparsity(cnt_1m, "Rolling LASSO sparsity — 1m horizon, window=120"))

# You can do the same for other horizons:
# cnt_3m <- roll_lasso_counts(df_3m, split_date, 120, 1); print(plot_roll_sparsity(cnt_3m, "3m horizon"))
# cnt_6m <- roll_lasso_counts(df_6m, split_date, 120, 1); print(plot_roll_sparsity(cnt_6m, "6m horizon"))
# cnt_12m <- roll_lasso_counts(df_12m, split_date, 120, 1); print(plot_roll_sparsity(cnt_12m, "12m horizon"))


