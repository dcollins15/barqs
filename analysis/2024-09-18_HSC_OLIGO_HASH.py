import pandas as pd

import barqs
import fastq

barcode_size = 16
umi_size = 12

feature_file = "/Users/dcollins/workspace/data/hto_demux/20240918/HTO_barcodes.xlsx"
feature_df = pd.read_excel(feature_file)
feature_df = feature_df.dropna()
feature_df = feature_df[feature_df.Oligo == "HSC_OLIGO_HASH"]

features = list(zip(feature_df.Condition, feature_df.Barcode))

feature_map = {seq: name for name, seq in features}

read1_file = "/Users/dcollins/workspace/data/hto_demux/20240918/HSC_OLIGO_HASH_S7_R1_001.fastq.gz"
read1 = fastq.load(read1_file)

read2_file = "/Users/dcollins/workspace/data/hto_demux/20240918/HSC_OLIGO_HASH_S7_R2_001.fastq.gz"
read2 = fastq.load(read2_file)

identifiers = (
    barqs.extract(read, umi_size=umi_size, barcode_size=barcode_size)
    for read in read1
)

reads = (
    barqs.tag(read, barcode, umi, trim = False)
    for read, (barcode, umi) in zip(read2, identifiers)
)

trimmed_reads = (
    barqs.trim_by_index(
        read, 
        feature_lookup=feature_map,
        region=(30, 45)
    ) 
    for read in reads
)

observed_features = (barqs.filter_duplicates(trimmed_reads))

counts = barqs.quantify(observed_features, features)

counts_df = pd.DataFrame(counts).fillna(0)

counts_df.head()


