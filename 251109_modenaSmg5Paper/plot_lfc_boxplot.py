#!/usr/bin/env python3
import sys
import pandas as pd
import matplotlib.pyplot as plt

def load_gene_list(path):
    """Load one gene per line."""
    with open(path, 'r') as f:
        return {line.strip() for line in f if line.strip()}

def load_deseq_file(path):
    """Load DESeq2-like CSV file with gene IDs in column 0."""
    df = pd.read_csv(path)
    df = df.rename(columns={df.columns[0]: "gene"})
    return df.set_index("gene")

def drop_outliers_iqr(s, k=1.5):
    """
    Drop outliers using Tukey's rule (k * IQR).
    Returns a filtered Series (outliers removed).
    """
    s = pd.to_numeric(s, errors="coerce").dropna()
    if s.empty:
        return s
    q1 = s.quantile(0.25)
    q3 = s.quantile(0.75)
    iqr = q3 - q1
    if iqr == 0:
        return s  # nothing to filter meaningfully
    lo = q1 - k * iqr
    hi = q3 + k * iqr
    return s[(s >= lo) & (s <= hi)]

def main():
    if len(sys.argv) != 6:
        print("Usage: python plot_common_lfc_boxplot.py common_genes.txt file1.csv file2.csv file3.csv output.svg")
        sys.exit(1)

    common_file, f1, f2, f3, out_svg = sys.argv[1:]

    # Load data
    common_genes = load_gene_list(common_file)
    df1 = load_deseq_file(f1)
    df2 = load_deseq_file(f2)
    df3 = load_deseq_file(f3)

    # Extract log2FC for only the common genes and flip to positive
    idx1 = df1.index.intersection(common_genes)
    idx2 = df2.index.intersection(common_genes)
    idx3 = df3.index.intersection(common_genes)

    lfc1 = pd.to_numeric(df1.loc[idx1, "log2FoldChange"], errors="coerce").abs()
    lfc2 = pd.to_numeric(df2.loc[idx2, "log2FoldChange"], errors="coerce").abs()
    lfc3 = pd.to_numeric(df3.loc[idx3, "log2FoldChange"], errors="coerce").abs()

    # Build dataframe (aligned by gene) and drop genes missing in any file
    plot_df = pd.DataFrame({"File1": lfc1, "File2": lfc2, "File3": lfc3}).dropna()
    n_before = plot_df.shape[0]

    # Drop outliers *from the data* (per column) before plotting
    f1_vals = drop_outliers_iqr(plot_df["File1"])
    f2_vals = drop_outliers_iqr(plot_df["File2"])
    f3_vals = drop_outliers_iqr(plot_df["File3"])

    med1 = f1_vals.median() if len(f1_vals) else float("nan")
    med2 = f2_vals.median() if len(f2_vals) else float("nan")
    med3 = f3_vals.median() if len(f3_vals) else float("nan")

    print("Median absolute log2FoldChange (after outlier removal):")
    print("  File 1: %.4f" % med1)
    print("  File 2: %.4f" % med2)
    print("  File 3: %.4f" % med3)

    # Plot
    plt.figure(figsize=(8, 6))

    positions = [1.0, 1.8, 2.6]

    plt.boxplot(
        [f1_vals.values, f2_vals.values, f3_vals.values],
        positions=positions,
        widths=0.5,
        labels=["File 1", "File 2", "File 3"],
        showfliers=False
    )

    plt.ylabel("absolute log2FoldChange")
    plt.title("Absolute log2FoldChange for Common WBGenes (outliers dropped)")

    plt.ylim(bottom=0)

    ymax = max(f1_vals.max(), f2_vals.max(), f3_vals.max())
    for y in range(0, int(ymax) + 1):
        plt.axhline(y, linestyle=":", linewidth=0.8, alpha=0.5)

    plt.tight_layout()
    plt.savefig(out_svg, format="svg")

    print(f"Saved plot to: {out_svg}")
    print(f"Used {n_before} common genes (complete across all files).")
    print(f"After outlier dropping (1.5×IQR): File1={len(f1_vals)}, File2={len(f2_vals)}, File3={len(f3_vals)}")

if __name__ == "__main__":
    main()
