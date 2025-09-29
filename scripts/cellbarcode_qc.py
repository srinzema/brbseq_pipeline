import argparse
import pandas as pd
import gzip
from collections import Counter
from multiprocessing import Pool, cpu_count

def hamming_distance(s1, s2):
    return sum(c1 != c2 for c1, c2 in zip(s1, s2))

def find_closest_barcode(barcodes, query):
    distances = [hamming_distance(barcode, query) for barcode in barcodes]
    min_dist = min(distances)
    if distances.count(min_dist) > 1:
        return "ambiguous"
    return barcodes[distances.index(min_dist)]

def process_chunk(reads, barcodes):
    bc_len = len(barcodes[0])
    counts = Counter()
    for read_seq in reads:
        read_bc = read_seq[:bc_len]
        closest = find_closest_barcode(barcodes, read_bc)
        counts[closest] += 1
    return counts

def read_fastq_chunks(filename, chunk_size=100000):
    open_func = gzip.open if filename.endswith(".gz") else open
    with open_func(filename, "rt") as f:
        chunk = []
        for i, line in enumerate(f):
            if i % 4 == 1:  # sequence line
                chunk.append(line.strip())
                if len(chunk) >= chunk_size:
                    yield chunk
                    chunk = []
        if chunk:
            yield chunk

parser = argparse.ArgumentParser()
parser.add_argument("--barcode-file", help="CSV/TSV file with 'barcode' and 'alias' or no header")
parser.add_argument("--fastq", help="FASTQ file (.gz or plain) to scan for barcodes")
parser.add_argument("--outfile", default="barcode_counts.csv", help="Output CSV file for barcode counts")
parser.add_argument("--threads", type=int, default=cpu_count(), help="Number of threads to use")
args = parser.parse_args()

# Load barcode file
try:
    df = pd.read_csv(args.barcode_file, sep=None, engine='python')
    if "barcode" not in df.columns or "alias" not in df.columns:
        raise ValueError
except:
    df = pd.read_csv(args.barcode_file, sep=None, engine='python', header=None)
    df.columns = ["barcode", "alias"]

barcodes = df["barcode"].tolist()
counts = Counter()

# Parallel processing
pool = Pool(args.threads)
results = []

for chunk in read_fastq_chunks(args.fastq):
    results.append(pool.apply_async(process_chunk, args=(chunk, barcodes)))

pool.close()
pool.join()

# Combine results
for r in results:
    counts.update(r.get())

# Print counts per barcode with alias
barcode_to_alias = dict(zip(df["barcode"], df["alias"]))
rows = []

for bc, count in counts.items():
    alias = barcode_to_alias.get(bc, bc)  # use bc itself for ambiguous
    rows.append([bc, alias, count])

df_counts = pd.DataFrame(rows, columns=["barcode", "alias", "count"])
df_counts.to_csv(args.outfile, index=False)