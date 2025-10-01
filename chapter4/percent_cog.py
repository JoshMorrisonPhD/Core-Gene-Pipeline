#!/usr/bin/env python3


# This code is designed to run stats on COG ouputs. It needs to be within the cog output folder to run.
# The input for this file is the cog_count.tsv file in COG outputs
# The output of this is a .tsv file showing the statistical significance of COG counts.
# This code creates plots in percentage for cog categories, and runs stats as described in section 3.2.8 of my thesis 

import sys
import argparse
from pathlib import Path

import altair as alt
import pandas as pd
import numpy as np

# optional deps (graceful fallback)
try:
    from scipy.stats import binomtest
    HAVE_SCIPY = True
except Exception:
    HAVE_SCIPY = False

try:
    from statsmodels.stats.proportion import proportion_confint
    from statsmodels.stats.multitest import multipletests
    HAVE_STATSMODELS = True
except Exception:
    HAVE_STATSMODELS = False

# COGclassifier plotting (for pretty barchart)
from cogclassifier.plot import plot_cog_count_barchart

alt.renderers.enable("png")

# helpers
def infer_columns(df: pd.DataFrame):
    cand_cat = ["LETTER","category","cog","COG","cog_category","name","label","Class","class"]
    cat_col = next((c for c in cand_cat if c in df.columns), None)
    if cat_col is None:
        cat_col = next(c for c in df.columns if not pd.api.types.is_numeric_dtype(df[c]))
    cand_cnt = ["count","COUNT","n","N","freq","frequency","total"]
    cnt_col = next((c for c in cand_cnt if c in df.columns), None)
    if cnt_col is None:
        cnt_col = next(c for c in df.columns if pd.api.types.is_numeric_dtype(df[c]))
    if cat_col is None or cnt_col is None:
        raise ValueError(f"Could not infer category/count columns from: {list(df.columns)}")
    return cat_col, cnt_col

def format_num(x):
    if pd.isnull(x): return ""
    return f"{x:.3g}"

def format_pval(x):
    if pd.isnull(x): return ""
    if x < 1e-3: return "<0.001"
    return f"{x:.3g}"

def resolve_input_path(inp: Path) -> Path:
    """Accept a TSV file or a directory made by COGclassifier; if dir, autodetect a suitable TSV."""
    if inp.is_file():
        return inp
    if not inp.is_dir():
        raise FileNotFoundError(f"Input path not found: {inp}")

    candidates = [
        "cog_count.tsv",
        "cog_counts.tsv",
        "counts.tsv",
        "classification.tsv",
        "cog_classify.tsv",  # some runs
    ]
    for name in candidates:
        p = inp / name
        if p.exists():
            return p

    for p in sorted(inp.glob("*.tsv")):
        try:
            df = pd.read_csv(p, sep="\t", nrows=200)
            if df.shape[1] >= 2 and any(pd.api.types.is_numeric_dtype(df[c]) for c in df.columns):
                return p
        except Exception:
            pass

    raise FileNotFoundError(f"No suitable TSV found inside directory: {inp}")

# main
def main():
    ap = argparse.ArgumentParser(
        description="COGclassifier figure + stats table. Input can be a TSV or a COGclassifier output directory."
    )
    ap.add_argument("-i", "--input", default="cog_count.tsv",
                    help="Input TSV OR a COGclassifier output directory (default: cog_count.tsv)")
    ap.add_argument("-o", "--outbase", default=None,
                    help="Output base name for figure/table (default: inferred from input)")
    ap.add_argument("--alpha", type=float, default=0.05,
                    help="FDR threshold for significance (default: 0.05)")
    ap.add_argument("--table", default=None,
                    help="Output TSV for the rich table (default: <outbase>_percent_significance.tsv)")
    args = ap.parse_args()

    here = Path(".").resolve()
    given = Path(args.input)
    in_path = resolve_input_path(given)

    # derive outputs
    if args.outbase is None:
        outbase = in_path.with_suffix("").name if in_path.suffix else in_path.name
    else:
        outbase = args.outbase
    out_html = here / f"{outbase}.html"
    out_png  = here / f"{outbase}.png"
    out_tsv  = here / (args.table if args.table else f"{outbase}_percent_significance.tsv")

    # load counts
    df_all = pd.read_csv(in_path, sep="\t")
    cat_col, cnt_col = infer_columns(df_all)
    df_all[cnt_col] = pd.to_numeric(df_all[cnt_col], errors="coerce").fillna(0).astype(int)

    # ---- zero-count safe path: no plots, minimal table, exit cleanly ----
    if df_all[cnt_col].sum() == 0:
        table = df_all[[cat_col, cnt_col]].rename(columns={cat_col:"category", cnt_col:"count"}).copy()
        for col in ["percent","pct_ci_low","pct_ci_high","expected_pct",
                    "enrichment_ratio","log2_enrichment","p_value","q_value","significant"]:
            table[col] = ""
        table.to_csv(out_tsv, sep="\t", index=False)
        print(f"[info] Zero COG hits in {in_path.name}; skipped plotting and wrote {out_tsv.name}")
        sys.exit(0)

    # figure (bar chart)
    chart = plot_cog_count_barchart(str(in_path), percent_style=True, sort=True)
    chart.save(out_html)
    chart.save(out_png)
    print(f"Saved chart:\n - {out_png}\n - {out_html}")

    # stats on nonzero categories only (if deps available)
    mask_nz = df_all[cnt_col] > 0
    df_nz = df_all[mask_nz].copy()

    # if missing deps, write counts-only table (with blanks for stats)
    if not (HAVE_SCIPY and HAVE_STATSMODELS) or len(df_nz) == 0:
        table = df_all[[cat_col, cnt_col]].rename(columns={cat_col:"category", cnt_col:"count"}).copy()
        for col in ["percent","pct_ci_low","pct_ci_high","expected_pct",
                    "enrichment_ratio","log2_enrichment","p_value","q_value","significant"]:
            table[col] = ""
        table = table.sort_values("count", ascending=False).reset_index(drop=True)
        table.to_csv(out_tsv, sep="\t", index=False)
        msg = "SciPy/statsmodels not available" if not (HAVE_SCIPY and HAVE_STATSMODELS) else "No non-zero categories"
        print(f"[info] {msg}; wrote counts-only table {out_tsv.name}")
        sys.exit(0)

    # compute stats
    total = int(df_nz[cnt_col].sum())
    K = len(df_nz)
    p0 = 1.0 / K

    pct = (df_nz[cnt_col] / total).to_numpy()
    ci_low, ci_high = proportion_confint(df_nz[cnt_col].to_numpy(), total, alpha=0.05, method="wilson")

    expected = np.full(K, p0)
    with np.errstate(divide="ignore", invalid="ignore"):
        enrich = pct / expected
        log2_enrich = np.log2(enrich)
        enrich[~np.isfinite(enrich)] = np.nan
        log2_enrich[~np.isfinite(log2_enrich)] = np.nan

    pvals = np.array([binomtest(int(x), total, p=p0, alternative="two-sided").pvalue for x in df_nz[cnt_col]])
    _, qvals, _, _ = multipletests(pvals, alpha=args.alpha, method="fdr_bh")
    sig = qvals <= args.alpha

    stats_nz = pd.DataFrame({
        "category": df_nz[cat_col].astype(str).to_numpy(),
        "count": df_nz[cnt_col].astype(int).to_numpy(),
        "percent": pct * 100.0,
        "pct_ci_low": ci_low * 100.0,
        "pct_ci_high": ci_high * 100.0,
        "expected_pct": expected * 100.0,
        "enrichment_ratio": enrich,
        "log2_enrichment": log2_enrich,
        "p_value": pvals,
        "q_value": qvals,
        "significant": sig
    })

    # format numbers
    for col in ["percent","pct_ci_low","pct_ci_high","expected_pct","enrichment_ratio","log2_enrichment"]:
        stats_nz[col] = stats_nz[col].apply(format_num)
    stats_nz["p_value"] = stats_nz["p_value"].apply(format_pval)
    stats_nz["q_value"] = stats_nz["q_value"].apply(format_pval)

    # merge back in original order, keeping zero rows with blank stats
    stats_map = stats_nz.set_index("category").to_dict(orient="index")
    rows = []
    for _, r in df_all.iterrows():
        cat = str(r[cat_col]); cnt = int(r[cnt_col])
        if cnt > 0 and cat in stats_map:
            rows.append({"category": cat, **stats_map[cat]})
        else:
            rows.append({
                "category": cat, "count": cnt,
                "percent": "", "pct_ci_low": "", "pct_ci_high": "",
                "expected_pct": "", "enrichment_ratio": "", "log2_enrichment": "",
                "p_value": "", "q_value": "", "significant": ""
            })

    table = pd.DataFrame(rows, columns=[
        "category","count","percent","pct_ci_low","pct_ci_high",
        "expected_pct","enrichment_ratio","log2_enrichment","p_value","q_value","significant"
    ]).sort_values("count", ascending=False).reset_index(drop=True)

    table.to_csv(out_tsv, sep="\t", index=False)
    print("Wrote stats table:", out_tsv.name)

if __name__ == "__main__":
    main()
