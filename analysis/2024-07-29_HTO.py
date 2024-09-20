# %%
import pandas as pd

# %%
import barqs
import fastq

# %%
barcode_size = 16
umi_size = 12

# %%
raw_features = [
    "TGGTATCAGTCTCGT",
    "ACGTAGTGCACTAGA",
    "ATTGCCTCAGGGTTT",
    "AGAGCGACAGAATAC",
    "TGAGAGAAGGGACTA",
    "GCATCTCTCAAGCAA",
    "CACCCAATAGTAAGG",
    "TTAACCACTGTGTAC",
    "TCGTTATAACATCTG",
    "GCAGTTAAAGTCGGG",
    "AGAATAGGAGCCGCA",
    "ATACACAGTTTACCG",
    "CATCTAGACCTAACT",
    "GAGATACTGCCATAG",
]

# %%
features = [
    (seq, seq) for seq in raw_features
]

# %%
feature_map = {seq: name for name, seq in features}

# %%
read1_file = "/Users/dcollins/workspace/data/jyun/HTO_R1.fastq.gz"
read1 = fastq.load(read1_file)

# %%
read2_file = "/Users/dcollins/workspace/data/jyun/HTO_R2.fastq.gz"
read2 = fastq.load(read2_file)

# %%
identifiers = (
    barqs.extract(read, umi_size=umi_size, barcode_size=barcode_size)
    for read in read1
)

# %%
reads = (
    barqs.tag(read, barcode, umi, trim = False)
    for read, (barcode, umi) in zip(read2, identifiers)
)

# %%
trimmed_reads = (
    barqs.trim_by_index(
        read, 
        feature_lookup=feature_map,
        region=(0, 15)
    ) 
    for read in reads
)

# %%
trimmed_reads = (
    barqs.trim_by_regex(
        read, 
        feature_lookup=feature_map,
        tolerance=1,
    )
    for read in trimmed_reads
)

# %%
observed_features = (barqs.filter_duplicates(trimmed_reads))

# %%
counts = barqs.quantify(observed_features, features)

# %%
counts_df = pd.DataFrame(counts).fillna(0)

# %%
counts_df.head()
# %%
