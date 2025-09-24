import pysam
import argparse
import pandas as pd
from collections import Counter
from multiprocessing import Pool, cpu_count


def count_bam_region(bam_path, reference):
    """Count mapped/unmapped reads per CB for one reference (chromosome)."""
    mapped, unmapped = Counter(), Counter()
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for read in bam.fetch(reference):
            if not read.has_tag("CB"):
                continue
            cb = read.get_tag("CB", with_value_type=False)
            if read.is_unmapped:
                unmapped[cb] += 1
            else:
                mapped[cb] += 1
    return mapped, unmapped


def merge_counts(counts_list):
    """Merge a list of (mapped, unmapped) counters."""
    mapped_total, unmapped_total = Counter(), Counter()
    for mapped, unmapped in counts_list:
        mapped_total.update(mapped)
        unmapped_total.update(unmapped)
    return mapped_total, unmapped_total


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Summarize cell barcode statistics from STARsolo output."
    )
    parser.add_argument("--bam", required=True, help="Path to the BAM file.")
    parser.add_argument(
        "--cellbarcodes", required=True, help="Path to the cell barcode file"
    )
    parser.add_argument(
        "--output", required=True, help="Path to the output summary file."
    )
    parser.add_argument(
        "--threads", type=int, default=cpu_count(), help="Number of worker processes"
    )
    args = parser.parse_args()

    # load barcode map
    barcode_map = pd.read_csv(
        args.cellbarcodes, sep="\t", header=None, names=["barcode", "alias"]
    )

    # get references (chromosomes/contigs)
    with pysam.AlignmentFile(args.bam, "rb") as bam:
        references = list(bam.references)

    # parallel count
    with Pool(processes=args.threads) as pool:
        results = pool.starmap(count_bam_region, [(args.bam, r) for r in references])

    mapped_counts, unmapped_counts = merge_counts(results)

    # build output
    all_barcodes = (
        set(barcode_map["barcode"])
        | set(mapped_counts.keys())
        | set(unmapped_counts.keys())
    )

    rows = []
    for bc in all_barcodes:
        # find alias in barcode_map, if any
        alias_row = barcode_map[barcode_map["barcode"] == bc]
        if not alias_row.empty:  # Barcode is in cell_barcode file
            alias = (
                alias_row.iloc[0]["alias"]
                if pd.notna(alias_row.iloc[0]["alias"])
                else bc
            )
        else:
            alias = f"unknown:{bc}"
        rows.append(
            {
                "sample": alias,
                "mapped": mapped_counts.get(bc, 0),
                "unmapped": unmapped_counts.get(bc, 0),
            }
        )

    pd.DataFrame(rows).to_csv(args.output, sep="\t", index=False)
