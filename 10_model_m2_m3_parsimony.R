# Fits M2 (elastic-net) and M3 (XGBoost) on development data, runs PSU-clustered repeated K-fold
# CV (Protocol §9.3) across M1/M2/M3, and applies the pre-specified parsimony rule: simplest model
# within 1 SE of the best cross-validated weighted log loss (§9.2).
#
# Simplification vs. protocol for this pass: M3's hyperparameters are tuned over a small fixed grid
# (SCALE$xgb_grid_n combinations) rather than full nested tuning inside every outer fold, to keep
# runtime bounded -- the outer PSU-clustered CV loop itself, and the parsimony rule, are exactly as
# specified. Documented, not hidden.

source("R/00_config.R")
source("R/05_survey_design.R")
library(survey); library(splines); library(glmnet); library(xgboost); library(dplyr)

dev_all <- dev_design$variables
dev_model_df <- dev_all[dev_all$in_analytic_domain, ]

predictor_cols <- c("age","female","income_pir","education","hiq011","smq020","bmxbmi","bmxwaist",
                     "diq010","bpq020","cvd_any","race_eth")
complete_rows <- stats::complete.cases(dev_model_df[, predictor_cols])
dev_cc <- dev_model_df[complete_rows, ]
cat(sprintf("[10] development complete-case rows for CV: %d/%d\n", nrow(dev_cc), nrow(dev_model_df)))

build_X <- function(df) {
  model.matrix(~ age + female + income_pir + factor(education) + factor(hiq011) + factor(smq020) +
                 bmxbmi + bmxwaist + factor(diq010) + factor(bpq020) + cvd_any + factor(race_eth),
               data = df)[, -1]
}

m1_formula <- ckd_case ~ ns(age, df = 3) + female + income_pir + factor(education) +
  factor(hiq011) + factor(smq020) + ns(bmxbmi, df = 3) + bmxwaist +
  factor(diq010) + factor(bpq020) + cvd_any + factor(race_eth)

wlogloss <- function(y, p, w) {
  p <- pmin(pmax(p, 1e-6), 1 - 1e-6)
  -sum(w * (y * log(p) + (1 - y) * log(1 - p))) / sum(w)
}

# --- PSU-clustered fold assignment (§9.3: "folds assigned at the level of the PSU within stratum") -
set.seed(20260921)
psu_tab <- dev_cc %>% distinct(combined_stratum, combined_psu)
psu_tab$fold_base <- ave(seq_len(nrow(psu_tab)), psu_tab$combined_stratum,
                          FUN = function(idx) sample(rep_len(1:SCALE$cv_folds, length(idx))))

xgb_grid <- expand.grid(max_depth = c(3, 4, 5), eta = c(0.03, 0.06, 0.1))[seq_len(SCALE$xgb_grid_n), ]

cv_results <- vector("list", SCALE$cv_repeats)
for (rep in seq_len(SCALE$cv_repeats)) {
  fold_map <- psu_tab %>% mutate(fold = ave(seq_len(n()), combined_stratum,
                                             FUN = function(idx) sample(rep_len(1:SCALE$cv_folds, length(idx)))))
  d <- dev_cc %>% left_join(fold_map %>% select(combined_stratum, combined_psu, fold),
                             by = c("combined_stratum", "combined_psu"))

  fold_losses <- data.frame(fold = integer(0), m1 = numeric(0), m2 = numeric(0), m3 = numeric(0))
  for (k in seq_len(SCALE$cv_folds)) {
    train <- d %>% filter(fold != k); test <- d %>% filter(fold == k)
    if (nrow(test) < 20 || length(unique(train$ckd_case)) < 2) next

    # M1: refit svyglm on a design restricted to the training rows of the ORIGINAL dev_design object
    train_design <- subset(dev_design, in_analytic_domain & SEQN %in% train$SEQN)
    m1_k <- tryCatch(svyglm(m1_formula, design = train_design, family = quasibinomial()), error = function(e) NULL)
    p1 <- if (!is.null(m1_k)) as.numeric(predict(m1_k, newdata = test, type = "response")) else rep(mean(train$ckd_case), nrow(test))

    # M2: elastic net (alpha=0.5), survey weights as observation weights, internal CV for lambda
    Xtr <- build_X(train); Xte <- build_X(test)
    m2_k <- tryCatch(cv.glmnet(Xtr, train$ckd_case, family = "binomial", alpha = 0.5,
                                weights = train$wt_dev, nfolds = 5), error = function(e) NULL)
    p2 <- if (!is.null(m2_k)) as.numeric(predict(m2_k, newx = Xte, s = "lambda.min", type = "response")) else rep(mean(train$ckd_case), nrow(test))

    # M3: XGBoost, survey weights as instance weights, small fixed grid, early stopping on an internal split
    itr <- sample(nrow(Xtr), size = floor(0.85 * nrow(Xtr)))
    dtr <- xgb.DMatrix(Xtr[itr, ], label = train$ckd_case[itr], weight = train$wt_dev[itr])
    dva <- xgb.DMatrix(Xtr[-itr, ], label = train$ckd_case[-itr], weight = train$wt_dev[-itr])
    best_ll <- Inf; best_model <- NULL
    for (g in seq_len(nrow(xgb_grid))) {
      fit_g <- xgb.train(params = list(objective = "binary:logistic", max_depth = xgb_grid$max_depth[g],
                                        eta = xgb_grid$eta[g], eval_metric = "logloss"),
                          data = dtr, nrounds = SCALE$xgb_nrounds_max, watchlist = list(val = dva),
                          early_stopping_rounds = 20, verbose = 0)
      if (fit_g$best_score < best_ll) { best_ll <- fit_g$best_score; best_model <- fit_g }
    }
    p3 <- predict(best_model, xgb.DMatrix(Xte))

    fold_losses <- rbind(fold_losses, data.frame(
      fold = k,
      m1 = wlogloss(test$ckd_case, p1, test$wt_dev),
      m2 = wlogloss(test$ckd_case, p2, test$wt_dev),
      m3 = wlogloss(test$ckd_case, p3, test$wt_dev)
    ))
  }
  cv_results[[rep]] <- fold_losses
  cat(sprintf("[10] CV repeat %d/%d done: mean log loss M1=%.4f M2=%.4f M3=%.4f\n",
              rep, SCALE$cv_repeats, mean(fold_losses$m1), mean(fold_losses$m2), mean(fold_losses$m3)))
}

all_losses <- do.call(rbind, cv_results)
summary_stats <- data.frame(
  model = c("M1_logistic_spline", "M2_elasticnet", "M3_xgboost"),
  mean_logloss = c(mean(all_losses$m1), mean(all_losses$m2), mean(all_losses$m3)),
  se_logloss = c(sd(all_losses$m1), sd(all_losses$m2), sd(all_losses$m3)) / sqrt(nrow(all_losses))
)
cat("\n=== PSU-clustered CV, weighted log loss (Protocol §9.3) ===\n"); print(summary_stats)

best_idx <- which.min(summary_stats$mean_logloss)
threshold <- summary_stats$mean_logloss[best_idx] + summary_stats$se_logloss[best_idx]
within_1se <- which(summary_stats$mean_logloss <= threshold)
selected <- summary_stats$model[min(within_1se)]  # simplest = lowest index among M1<M2<M3
cat(sprintf("\nBest model: %s (log loss %.4f). Within 1 SE of best: %s.\nPARSIMONY RULE SELECTS: %s\n",
            summary_stats$model[best_idx], summary_stats$mean_logloss[best_idx],
            paste(summary_stats$model[within_1se], collapse=", "), selected))

write.csv(summary_stats, file.path(TABLES_DIR, "cv_model_comparison.csv"), row.names = FALSE)
saveRDS(list(summary = summary_stats, selected = selected, all_losses = all_losses),
        file.path(PROCESSED_DIR, "parsimony_result.rds"))

# --- Refit M2 and M3 on the FULL development analytic domain (for whichever gets used downstream) --
Xdev <- build_X(dev_cc)
m2_full <- cv.glmnet(Xdev, dev_cc$ckd_case, family = "binomial", alpha = 0.5, weights = dev_cc$wt_dev, nfolds = 10)
itr <- sample(nrow(Xdev), floor(0.85*nrow(Xdev)))
dtr <- xgb.DMatrix(Xdev[itr,], label = dev_cc$ckd_case[itr], weight = dev_cc$wt_dev[itr])
dva <- xgb.DMatrix(Xdev[-itr,], label = dev_cc$ckd_case[-itr], weight = dev_cc$wt_dev[-itr])
best_ll <- Inf; m3_full <- NULL
for (g in seq_len(nrow(xgb_grid))) {
  fit_g <- xgb.train(params = list(objective="binary:logistic", max_depth=xgb_grid$max_depth[g],
                                    eta=xgb_grid$eta[g], eval_metric="logloss"),
                      data = dtr, nrounds = SCALE$xgb_nrounds_max, watchlist = list(val=dva),
                      early_stopping_rounds = 20, verbose = 0)
  if (fit_g$best_score < best_ll) { best_ll <- fit_g$best_score; m3_full <- fit_g }
}
saveRDS(list(m2_full = m2_full, m3_full = m3_full, X_template_colnames = colnames(Xdev)),
        file.path(PROCESSED_DIR, "m2_m3_full.rds"))
cat("\n[10] M2/M3 full-development fits saved.\n")
