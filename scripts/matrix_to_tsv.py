#!/usr/bin/env python3
from scipy.io import mmread
import pandas as pd
import argparse

parser = argparse.ArgumentParser(description="Convert matrix to TSV")
parser.add_argument("--matrix", required=True, help="Path to matrix.mtx file")
parser.add_argument("--features", required=True, help="Path to features.tsv file")
parser.add_argument("--barcodes", required=True, help="Path to barcodes.tsv file")
parser.add_argument("--output", required=True, help="Path to output TSV file")
args = parser.parse_args()

matrix = mmread(args.matrix).tocsc()
features = pd.read_csv(args.features, header=None, sep="\t")
barcodes = pd.read_csv(args.barcodes, header=None, sep="\t", names=["barcode", "alias"])
print(barcodes.head())

# Use alias if non-empty, otherwise keep barcode
columns = barcodes["alias"].where(
    barcodes["alias"].notna() & (barcodes["alias"] != ""), barcodes["barcode"]
)

pd.DataFrame(
    matrix.toarray(),
    index=features[0],
    columns=columns,
).to_csv(args.output, sep="\t", index=True, header=True)
