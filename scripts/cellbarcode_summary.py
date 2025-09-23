import pandas as pd
import argparse

parser = argparse.ArgumentParser(description="Summarize cell barcode statistics from STARsolo output.")
parser.add_argument("--count-table", required=True, help="Path to the matrix.mtx file.")
parser.add_argument("--output", required=True, help="Path to the output summary file.")
args = parser.parse_args()

# Load the count matrix
count_table = pd.read_csv(args.count_table, sep="\t", header=0, index_col=0)

total_reads = count_table.sum(axis=0)
detected_genes = (count_table > 0).sum(axis=0)
avg_per_gene = count_table.mean(axis=0)
median_per_gene = count_table.median(axis=0)

# Combine into one dataframe
metrics = pd.DataFrame({
    "Total_reads": total_reads,
    "Detected_genes": detected_genes,
    "Avg_reads_per_gene": avg_per_gene,
    "Median_reads_per_gene": median_per_gene
})
metrics = metrics.loc[~(metrics==0).all(axis=1)]
metrics.to_csv(args.output, sep="\t", index_label="Cell_Barcode")
