#!/usr/bin/env Rscript
# =============================================================================
# 00_setup_survival.R — CSF解析用 定数・共変量・設定
# Causal Survival Forest + RMST 解析
# =============================================================================

suppressPackageStartupMessages({
  library(grf)
  library(readr)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(survival)
  library(here)
})

# パス設定
# code_release/ をルートとする自己完結構成（公開版）。
# here::here() は code_release/ 直下（.Rprofile / renv.lock のある階層）を指す。
# 原データ・生成物は本リポジトリに含まれない（README の "Data availability" 参照）。
PROJECT_DIR  <- here::here()
SURV_DIR     <- PROJECT_DIR
SURV_CSV     <- file.path(SURV_DIR, "outputs", "joint0405_survival.csv")
OUT_DIR      <- file.path(SURV_DIR, "outputs")
TABLE_DIR    <- file.path(OUT_DIR, "tables")
FIGURE_DIR   <- file.path(OUT_DIR, "figures")
LOG_DIR      <- file.path(OUT_DIR, "logs")

dir.create(TABLE_DIR,  showWarnings = FALSE, recursive = TRUE)
dir.create(FIGURE_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(LOG_DIR,    showWarnings = FALSE, recursive = TRUE)

# =============================================================================
# 共変量定義
# =============================================================================

COVARIATES_INTEGRATED <- c(
  "TRIAL_ID", "AGE", "bmi", "MENOPAUSE_AGE",
  "TUG_0", "ONE_LEG_STANDING_0",
  "TRACP_5B_0", "EGFR_0",
  "LDL_0", "HDL_0", "TG_0", "TOTAL_CHOLESTEROL_0",
  "PREVALENT_VF_COUNT", "PREVALENT_VF_MAX_GR", "PREVALENT_PROXIMAL_FEMOR_FX",
  "COMORBIDITY_DM", "COMORBIDITY_HL", "COMORBIDITY_HT",
  "PRIOR_BP_TREATMENT"
)

# 骨折歴除外（感度分析用）
COVARIATES_INTEGRATED_NOHX <- setdiff(
  COVARIATES_INTEGRATED,
  c("PREVALENT_VF_COUNT", "PREVALENT_VF_MAX_GR")
)

# =============================================================================
# CSF ハイパーパラメータ
# =============================================================================
HORIZON       <- 2.0        # 24ヶ月 = 2.0年
CSF_NUM_TREES <- 2000      # grf default
CSF_MIN_NODE  <- 5
CSF_SEED      <- 42
CSF_TARGET    <- "RMST"     # "RMST" or "survival.probability"

# GATE 分位数
GATE_PROBS <- c(0, 0.25, 0.5, 0.75, 1.0)

# 方向性:
#   CSF/RMST では benefit = 無骨折期間が長い = 正の CATE
#   ※ binary CF とは逆（binary では benefit = 骨折抑制 = 負のCATE）
CSF_BENEFIT_DIR <- "pos_cate"  # AUTOC ranking_score = +tau_hat

# =============================================================================
# アウトカム定義
# =============================================================================
SURV_OUTCOMES <- list(
  list(code = "E1", time_col = "T_VF",          event_col = "D_VF",
       label = "Vertebral Fracture (RMST diff)",
       priority = 1),
  list(code = "E2", time_col = "NVF_PERSON_YEAR", event_col = "NVF_BINARY",
       label = "Non-Vertebral Fracture (RMST diff)",
       priority = 2),
  list(code = "E3", time_col = "T_E3",          event_col = "D_E3",
       label = "Composite Fracture VF+NVF (RMST diff)",
       priority = 3),
  list(code = "E4", time_col = "T_E4",          event_col = "D_E4",
       label = "Hip/Femoral Fracture (exploratory)",
       priority = 4)
)

# =============================================================================
# 出力パス
# =============================================================================
table_path_surv  <- function(code, suffix)
  file.path(TABLE_DIR,  paste0("CSF_table_", code, "_", suffix, ".csv"))
figure_path_surv <- function(code, suffix)
  file.path(FIGURE_DIR, paste0("CSF_fig_",   code, "_", suffix, ".png"))

cat(sprintf(
  "CSF Setup loaded. Covariates: %d (NOHX: %d), horizon=%.1fyr, target=%s\n",
  length(COVARIATES_INTEGRATED), length(COVARIATES_INTEGRATED_NOHX),
  HORIZON, CSF_TARGET
))
