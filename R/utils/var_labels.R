# =============================================================================
# var_labels.R — 共変量の「生カラム名 → 読者向け表示ラベル」辞書
# -----------------------------------------------------------------------------
# 図表に生のデータ項目名（COMORBIDITY_HT 等）が出ないよう、論文掲載用の
# 自然な英語ラベルに置換するための一元管理ファイル。
# manuscript.md / dist/manuscript.docx の本文・表ラベルと表記を揃えること。
#
# 使い方:
#   source("R/utils/var_labels.R")
#   relabel_vars(c("COMORBIDITY_HT", "TUG_0"))  # -> c("Hypertension", "TUG test")
# =============================================================================

VAR_LABELS <- c(
  # 試験・人口統計
  "TRIAL_ID"                    = "Trial (JOINT-04/05)",
  "AGE"                         = "Age",
  "bmi"                         = "BMI",
  "BMI"                         = "BMI",
  "MENOPAUSE_AGE"               = "Age at menopause",
  # 身体機能
  "TUG_0"                       = "TUG test",
  "ONE_LEG_STANDING_0"          = "One-leg standing time",
  # 骨代謝・腎機能
  "TRACP_5B_0"                  = "TRACP-5b",
  "EGFR_0"                      = "eGFR",
  # 脂質
  "LDL_0"                       = "LDL cholesterol",
  "HDL_0"                       = "HDL cholesterol",
  "TG_0"                        = "Triglycerides",
  "TOTAL_CHOLESTEROL_0"         = "Total cholesterol",
  # 骨折歴
  "PREVALENT_VF_COUNT"          = "Prevalent vertebral fractures (n)",
  "PREVALENT_VF_MAX_GR"         = "Prevalent vertebral fracture grade",
  "PREVALENT_PROXIMAL_FEMOR_FX" = "Prior proximal-femur fracture",
  # 併存症
  "COMORBIDITY_DM"              = "Diabetes mellitus",
  "COMORBIDITY_HL"             = "Dyslipidemia",
  "COMORBIDITY_HT"              = "Hypertension",
  # 前治療
  "PRIOR_BP_TREATMENT"          = "Prior bisphosphonate treatment",
  # 切片
  "(Intercept)"                 = "(Intercept)"
)

# ベクトルを表示ラベルへ変換（辞書に無い名前はそのまま返す）
relabel_vars <- function(x) {
  out <- VAR_LABELS[as.character(x)]
  out[is.na(out)] <- as.character(x)[is.na(out)]
  unname(out)
}
