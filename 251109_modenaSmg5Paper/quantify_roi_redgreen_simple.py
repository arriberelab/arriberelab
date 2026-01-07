#!/usr/bin/env python3
# Simplified ROI quantification using *fuzzy/averaged* values for BOTH channels (RFP and GFP)
#
# Usage:
#   python quantify_roi_redgreen_simple.py C1_RED.tif C2_GREEN.tif outPrefix \
#       --roi x1 y1 x2 y2 [--r 5] [--cut 5000]
#
import argparse
from pathlib import Path
import numpy as np
import tifffile
import matplotlib.pyplot as plt
from scipy.stats import linregress
import matplotlib.patches as patches

def fuzzy_avg(row: int, col: int, arr: np.ndarray, r: int) -> float:
    """Average of arr values in a disk of radius r centered at (row, col)."""
    vals = []
    H, W = arr.shape
    r2 = r * r
    for rr in range(row - r, row + r + 1):
        if rr < 0 or rr >= H:
            continue
        dy2 = (rr - row) * (rr - row)
        for cc in range(col - r, col + r + 1):
            if cc < 0 or cc >= W:
                continue
            dx2 = (cc - col) * (cc - col)
            if dx2 + dy2 <= r2:
                vals.append(arr[rr, cc])
    return float(np.average(vals)) if vals else 0.0

def rfp_summary(r_vals: np.ndarray) -> dict:
    """Compute rich summary metrics for RFP (fuzzy) values."""
    n = int(r_vals.size)
    if n == 0:
        return dict(
            red_mean=np.nan, red_median=np.nan, red_std=np.nan, red_mad=np.nan,
            red_min=np.nan, red_max=np.nan, red_q10=np.nan, red_q25=np.nan,
            red_q75=np.nan, red_q90=np.nan, red_sum=np.nan, red_cv=np.nan
        )
    med = float(np.median(r_vals))
    mad = float(np.median(np.abs(r_vals - med)))
    mean = float(np.mean(r_vals))
    std = float(np.std(r_vals, ddof=1)) if n > 1 else float('nan')
    cv = (std / mean) if (n > 1 and mean != 0.0) else float('nan')
    return dict(
        red_mean=mean,
        red_median=med,
        red_std=std,
        red_mad=mad,
        red_min=float(np.min(r_vals)),
        red_max=float(np.max(r_vals)),
        red_q10=float(np.percentile(r_vals, 10)),
        red_q25=float(np.percentile(r_vals, 25)),
        red_q75=float(np.percentile(r_vals, 75)),
        red_q90=float(np.percentile(r_vals, 90)),
        red_sum=float(np.sum(r_vals)),
        red_cv=float(cv)
    )

def linreg(g_vals: np.ndarray, r_vals: np.ndarray):
    """Linear regression RFP ~ GFP, return dict of slope/intercept/r2 and n."""
    n = int(r_vals.size)
    if n < 2:
        return dict(slope=np.nan, intercept=np.nan, r2=np.nan, n=n)
    slope, intercept, r, p, se = linregress(g_vals, r_vals)
    return dict(slope=float(slope), intercept=float(intercept), r2=float(r**2), n=n)

def save_overlay(bg_vals: np.ndarray, roi_box, path: str, title: str):
    """Save an overlay image with ROI rectangle (white)."""
    x1, y1, x2, y2 = roi_box
    fig, ax = plt.subplots()
    im = ax.imshow(bg_vals, origin="upper")
    plt.colorbar(im, label="Fuzzy RFP (ROI, GFP-pass)")
    rect = patches.Rectangle((x1, y1), x2 - x1, y2 - y1,
                             linewidth=1.5, edgecolor="white", facecolor="none")
    ax.add_patch(rect)
    ax.set_title(title)
    plt.savefig(path, dpi=150, bbox_inches="tight")
    plt.close()

def save_scatter(g_vals: np.ndarray, r_vals: np.ndarray, fit: dict, path: str):
    """Scatter of GFP_fuzzy vs RFP_fuzzy with regression line and R^2 annotation."""
    fig, ax = plt.subplots()
    # Downsample for plotting if huge
    N = g_vals.size
    if N > 50000:
        idx = np.random.choice(N, 50000, replace=False)
        gx, ry = g_vals[idx], r_vals[idx]
    else:
        gx, ry = g_vals, r_vals
    ax.plot(gx, ry, '.', markersize=1)
    if not np.isnan(fit['slope']):
        xline = np.linspace(np.min(gx), np.max(gx), 100)
        yline = fit['slope'] * xline + fit['intercept']
        ax.plot(xline, yline, linewidth=1)
        ax.set_title(f"RFP_fuzzy ~ GFP_fuzzy   n={fit['n']}   R^2={fit['r2']:.4f}")
    else:
        ax.set_title(f"RFP_fuzzy ~ GFP_fuzzy   n={fit['n']}   R^2=NA")
    ax.set_xlabel('GFP (fuzzy)')
    ax.set_ylabel('RFP (fuzzy)')
    plt.tight_layout()
    plt.savefig(path, dpi=200)
    plt.close()

def parse_args():
    ap = argparse.ArgumentParser(description="Simplified ROI quant using *fuzzy* (RFP,GFP) values + overlay + scatter + rich metrics")
    ap.add_argument("c1_red", help="Red channel TIFF (c1)")
    ap.add_argument("c2_green", help="Green channel TIFF (c2)")
    ap.add_argument("outPrefix", help="Prefix for outputs (directories created as needed)")
    ap.add_argument("--roi", nargs=4, type=int, metavar=("x1","y1","x2","y2"),
                    required=True, help="Rectangle in image coords (x=cols, y=rows)")
    ap.add_argument("--r", type=int, default=5, help="Fuzzy radius (default 5)")
    ap.add_argument("--cut", type=float, default=5000, help="GFP cutoff on *fuzzy* GFP (default 5000)")
    return ap.parse_args()

def main():
    args = parse_args()
    x1, y1, x2, y2 = args.roi
    if not (x2 > x1 and y2 > y1):
        raise SystemExit("ROI must satisfy x1<x2 and y1<y2")

    # Load channels (raw)
    red_raw = tifffile.imread(args.c1_red)
    green_raw = tifffile.imread(args.c2_green)
    if red_raw.shape != green_raw.shape:
        raise SystemExit("Red and Green images must have identical dimensions")
    H, W = red_raw.shape

    # Clip ROI to bounds
    x1 = max(0, x1); y1 = max(0, y1)
    x2 = min(W, x2); y2 = min(H, y2)

    r = int(args.r)
    cut = float(args.cut)

    # Build ROI mask
    roi_mask = np.zeros((H, W), dtype=bool)
    roi_mask[y1:y2, x1:x2] = True

    # Compute fuzzy maps for BOTH channels within ROI region (skip borders like original)
    fuzzy_g = np.zeros((H, W), dtype=float)
    fuzzy_r = np.zeros((H, W), dtype=float)
    for row in range(r, H - r):
        for col in range(r, W - r):
            if not roi_mask[row, col]:
                continue
            fuzzy_g[row, col] = fuzzy_avg(row, col, green_raw, r)
            fuzzy_r[row, col] = fuzzy_avg(row, col, red_raw, r)

    # GFP-pass mask based on fuzzy GFP
    passed = np.zeros((H, W), dtype=bool)
    passed[y1:y2, x1:x2] = (fuzzy_g[y1:y2, x1:x2] >= cut)

    # Collect *fuzzy* values at passed pixels
    rows, cols = np.where(passed)
    g_all = fuzzy_g[rows, cols]
    r_all = fuzzy_r[rows, cols]

    # Summaries
    fit = linreg(g_all, r_all)
    rsum = rfp_summary(r_all)

    outPrefix = Path(args.outPrefix)
    outPrefix.parent.mkdir(parents=True, exist_ok=True)

    # Save ROI metrics TSV (append)
    tsv = outPrefix.with_suffix("").as_posix() + ".roi_metrics.tsv"
    header = [
        "x1","y1","x2","y2","radius","gfp_cut",
        "n_passed","slope","intercept","r2",
        "red_mean","red_median","red_std","red_mad","red_min","red_max",
        "red_q10","red_q25","red_q75","red_q90","red_sum","red_cv",
        "values_type"
    ]
    values = [
        x1, y1, x2, y2, r, cut,
        fit["n"], fit["slope"], fit["intercept"], fit["r2"],
        rsum["red_mean"], rsum["red_median"], rsum["red_std"], rsum["red_mad"],
        rsum["red_min"], rsum["red_max"], rsum["red_q10"], rsum["red_q25"],
        rsum["red_q75"], rsum["red_q90"], rsum["red_sum"], rsum["red_cv"],
        "fuzzy"
    ]
    write_header = not Path(tsv).exists()
    with open(tsv, "a") as f:
        if write_header:
            f.write("\t".join(header) + "\n")
        f.write("\t".join(str(v) for v in values) + "\n")

    # Overlay based on *fuzzy* RFP
    just_reds = np.zeros_like(red_raw, dtype=float)
    if r_all.size:
        just_reds[rows, cols] = r_all
    save_overlay(just_reds, (x1, y1, x2, y2),
                 str(outPrefix) + ".roi_overlay.png",
                 "ROI overlay — *fuzzy* RFP (GFP-pass)")

    # Scatter plot with regression
    save_scatter(g_all, r_all, fit, str(outPrefix) + ".roi_scatter.png")

    # Console summary
    print(f"[ROI x1={x1} y1={y1} x2={x2} y2={y2}] radius={r} gfp_cut={cut}")
    print(f"  Passed pixels: n={fit['n']}  (fuzzy values)")
    print(f"  r2={fit['r2']}  slope={fit['slope']}  intercept={fit['intercept']}")
    print(f"  RFP: median={rsum['red_median']}, mean={rsum['red_mean']}, std={rsum['red_std']}, MAD={rsum['red_mad']}, CV={rsum['red_cv']}")

if __name__ == "__main__":
    main()
