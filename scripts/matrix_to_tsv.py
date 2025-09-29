#!/usr/bin/env python3
from scipy.io import mmread
import pandas as pd
import argparse

parser = argparse.ArgumentParser(description="Convert matrix to TSV")
parser.add_argument("--matrix", required=True, help="Path to matrix.mtx file")
parser.add_argument("--features", required=True, help="Path to features.tsv file")
parser.add_argument("--barcodes", required=True, help="Path to barcodes.tsv file")
parser.add_argument("--aliases", required=True, help="Path to aliases.tsv file")
parser.add_argument("--output", required=True, help="Path to output TSV file")
args = parser.parse_args()

matrix = mmread(args.matrix).tocsc()
features = pd.read_csv(args.features, header=None, sep="\t")
barcodes = pd.read_csv(args.barcodes, header=None, sep="\t", names=["barcode"])
aliases = pd.read_csv(args.aliases, header=None, sep="\t", names=["barcode", "alias"])
aliases["alias"] = aliases["alias"].replace("", pd.NA).fillna(aliases["barcode"])
print(aliases.head())

matrix = pd.DataFrame(
    matrix.toarray(),
    index=features[0],
    columns=barcodes["barcode"],
)

for col in matrix.columns:
    match = aliases.loc[aliases["barcode"] == col, "alias"]
    if not match.empty:
        matrix = matrix.rename(columns={col: match.iloc[0]})

matrix.index.name = "gene"
matrix.to_csv(args.output, sep="\t", index=True, header=True)
