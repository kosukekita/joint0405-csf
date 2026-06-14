#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
JOINT-04/05 統合データ前処理 — Causal Survival Forest 用
==========================================================
入力 : data/joint0405_raw.xlsx   （★ 原データは本リポジトリに含まれない。README 参照）
出力 : outputs/joint0405_survival.csv

設計判断:
- ラロキシフェン群（J04）を除外 → n=2,608
- W = 1 if ALLOCATION == 'J05テリパラチド群'
- T_VF = VF_PERSON_YEAR（生値をそのまま渡す）
- D_VF = (VF_COUNT > 0)（全観察期間でのイベント有無）
- horizon打ち切りは causal_survival_forest(horizon=2.0) 引数に委譲
- person-year の負値・ゼロを除外
"""

import pandas as pd
import numpy as np
from pathlib import Path

# code_release/ をルートとする自己完結構成。
# 原データ（data/joint0405_raw.xlsx）はプライバシー保護のため本リポジトリに含まれない。
# 入手方法は README.md の "Data availability" を参照。
PROJECT_DIR = Path(__file__).parent.parent          # = code_release/
DATA_FILE   = PROJECT_DIR / "data" / "joint0405_raw.xlsx"
OUT_DIR     = PROJECT_DIR / "outputs"
OUT_FILE    = OUT_DIR / "joint0405_survival.csv"

COVARIATES_INTEGRATED = [
    "TRIAL_ID", "AGE", "bmi", "MENOPAUSE_AGE",
    "TUG_0", "ONE_LEG_STANDING_0",
    "TRACP_5B_0", "EGFR_0",
    "LDL_0", "HDL_0", "TG_0", "TOTAL_CHOLESTEROL_0",
    "PREVALENT_VF_COUNT", "PREVALENT_VF_MAX_GR", "PREVALENT_PROXIMAL_FEMOR_FX",
    "COMORBIDITY_DM", "COMORBIDITY_HL", "COMORBIDITY_HT",
    "PRIOR_BP_TREATMENT",
]


def load_raw():
    print(f"Loading: {DATA_FILE}")
    df = pd.read_excel(DATA_FILE, sheet_name=0)
    print(f"  Raw shape: {df.shape}")
    return df


def exclude_raloxifene(df):
    """ラロキシフェン群を除外"""
    before = len(df)
    df = df[df["ALLOCATION"] != "J04ラロキシフェン群"].copy()
    after = len(df)
    print(f"  Raloxifene excluded: {before} → {after} (removed {before-after})")
    return df


def make_treatment(df):
    """W: テリパラチド=1, BP=0"""
    df["W"] = (df["ALLOCATION"] == "J05テリパラチド群").astype(int)
    n_pth = df["W"].sum()
    n_bp  = (df["W"] == 0).sum()
    print(f"  W=1(PTH)={n_pth}, W=0(BP)={n_bp}")
    assert n_pth == 489, f"PTH count mismatch: {n_pth}"
    assert len(df) == 2608, f"n mismatch: {len(df)}"
    return df


def make_site_code(df):
    df["SITE_CODE"] = df["ID"].astype(str).str.split("-").str[0]
    print(f"  Sites: {df['SITE_CODE'].nunique()}")
    return df


def clean_person_year(df, col):
    """負値・ゼロを NaN に置換（除外フラグとして使用）"""
    n_bad = (df[col] <= 0).sum()
    if n_bad > 0:
        print(f"  {col}: {n_bad} 件の負値/ゼロ → NaN")
        df[col] = df[col].where(df[col] > 0, np.nan)
    return df


def make_survival_outcomes(df):
    """
    E1(VF), E3(複合VF+NVF), E4(大腿骨) の time/event を生成
    horizon打ち切りは causal_survival_forest(horizon=2.0) に委譲するため、
    生の person-year と全観察期間でのイベント有無をそのまま渡す。
    """
    df = clean_person_year(df, "VF_PERSON_YEAR")
    df = clean_person_year(df, "NVF_PERSON_YEAR")

    # ---- E1: 椎体骨折 ----
    vfpy_valid = df["VF_PERSON_YEAR"].notna()
    df["T_VF"] = np.where(vfpy_valid, df["VF_PERSON_YEAR"], np.nan)
    df["D_VF"] = np.where(vfpy_valid,
                           (df["VF_COUNT"] > 0).astype(float),
                           np.nan)
    n_event_vf = (df["D_VF"] == 1).sum()
    n_valid_vf = vfpy_valid.sum()
    print(f"  E1(VF): valid={n_valid_vf}, total events={n_event_vf}")

    # ---- E3: 複合骨折（VF or NVF）----
    # NVF_COUNT に NA=953 あり → VF発生例は NVF NA に関わらずイベントあり
    nvf_valid = df["NVF_COUNT"].notna()
    e3_event_raw = (df["VF_COUNT"] > 0) | ((nvf_valid) & (df["NVF_COUNT"] > 0))

    # time: VF_PERSON_YEAR と NVF_PERSON_YEAR の min（より早い方）
    df["T_E3_raw"] = df[["VF_PERSON_YEAR", "NVF_PERSON_YEAR"]].min(axis=1)
    e3py_valid = df["T_E3_raw"].notna()
    df["T_E3"] = np.where(e3py_valid, df["T_E3_raw"], np.nan)
    df["D_E3"] = np.where(e3py_valid,
                           e3_event_raw.astype(float),
                           np.nan)
    n_event_e3 = (df["D_E3"] == 1).sum()
    print(f"  E3(VF+NVF): valid={e3py_valid.sum()}, total events={n_event_e3}")

    # ---- E4: 大腿骨骨折 ----
    femor_valid = df["NVF_FEMOR_COUNT"].notna() & df["NVF_PERSON_YEAR"].notna()
    df["T_E4"] = np.where(femor_valid, df["NVF_PERSON_YEAR"], np.nan)
    df["D_E4"] = np.where(femor_valid,
                           (df["NVF_FEMOR_COUNT"] > 0).astype(float),
                           np.nan)
    n_event_e4 = (df["D_E4"] == 1).sum()
    print(f"  E4(大腿骨): valid={femor_valid.sum()}, total events={n_event_e4}")

    # binary アウトカム（coherence check 用）
    df["VF_BINARY"]  = (df["VF_COUNT"] > 0).astype(int)
    df["NVF_BINARY"] = df["NVF_COUNT"].gt(0).astype(int)   # NA → 0
    df["E3_BINARY"]  = ((df["VF_COUNT"] > 0) | df["NVF_COUNT"].gt(0)).astype(int)
    df["E4_BINARY"]  = (df["NVF_FEMOR_COUNT"] > 0).astype(int)

    return df


def export(df):
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    cols_out = (
        ["ID", "W", "TRIAL_ID", "SITE_CODE", "ALLOCATION"]
        + ["VF_PERSON_YEAR", "NVF_PERSON_YEAR"]
        + ["T_VF", "D_VF", "T_E3", "D_E3", "T_E4", "D_E4"]
        + ["VF_BINARY", "NVF_BINARY", "E3_BINARY", "E4_BINARY"]
        + ["VF_COUNT", "NVF_COUNT", "NVF_FEMOR_COUNT"]
        + [c for c in COVARIATES_INTEGRATED if c != "TRIAL_ID"]
    )
    cols_out = [c for c in cols_out if c in df.columns]
    df[cols_out].to_csv(OUT_FILE, index=False)
    print(f"\n  Saved → {OUT_FILE}  ({len(df)} rows × {len(cols_out)} cols)")


def summarise(df):
    print("\n" + "=" * 55)
    print("DATA SUMMARY — 統合 n=2,608")
    print("=" * 55)
    print(f"  TRIAL_ID: J04={( df['TRIAL_ID']==4).sum()}, J05={(df['TRIAL_ID']==5).sum()}")
    print(f"  W=0(BP)={(df['W']==0).sum()}, W=1(PTH)={(df['W']==1).sum()}")
    print(f"\n  Outcomes (horizon打ち切りはCSF引数に委譲):")
    for col, label in [("D_VF","E1 VF"),("D_E3","E3 複合"),("D_E4","E4 大腿骨")]:
        n_ev = (df[col]==1).sum()
        n_va = df[col].notna().sum()
        print(f"    {label}: {n_ev}/{n_va} events ({n_ev/n_va*100:.1f}%)")
    print(f"\n  Covariate missing (%):")
    for c in COVARIATES_INTEGRATED:
        if c in df.columns:
            miss = df[c].isna().mean() * 100
            if miss > 0:
                print(f"    {c}: {miss:.1f}%")


def main():
    print("=" * 55)
    print("JOINT-04/05 Survival Preprocessing")
    print("=" * 55 + "\n")

    df = load_raw()
    df = exclude_raloxifene(df)
    df = make_treatment(df)
    df = make_site_code(df)
    df = make_survival_outcomes(df)
    export(df)
    summarise(df)
    print("\nDone.")


if __name__ == "__main__":
    main()
