import barqs
import fastq


path_to_fastq = "/Users/dcollins/workspace/barqs/tests/test_data/small.fastq"
raw_reads = list(fastq.load(path_to_fastq))

umi_size = 12
barcode_size = 16

identifiers = (
    barqs.extract(read, umi_size=umi_size, barcode_size=barcode_size)
    for read in raw_reads
)

reads = list(
    barqs.tag(read, barcode, umi)
    for read, (barcode, umi) in zip(raw_reads, identifiers)
)

features = [
    ("DAVE", "TGGAACTGATTA"),
    ("ELIA", "TTGCGCTGATTA"),
]

feature_lookup = {seq: name for name, seq in features}

trimmed_reads = list(
    barqs.trim_by_index(
        read, 
        feature_lookup=feature_lookup,
        region=(0, 12)
    ) 
    for read in reads
)

observed_features = list(barqs.filter_duplicates(trimmed_reads))

counts = barqs.quantify(observed_features, features)
print(counts)