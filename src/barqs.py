import regex
from collections import defaultdict
from typing import Mapping, Optional

import fasta
import fastq


def extract(
    read1: fastq.FASTQRecord,
    barcode_length: int,
    umi_length: int,
    read2: Optional[fastq.FASTQRecord] = None,
) -> fastq.FASTQRecord:

    header, seq, quality_scores = read1

    barcode = seq[:barcode_length]
    umi = seq[barcode_length : barcode_length + umi_length]

    if read2:
        header, seq, quality_scores = read2
    else:
        seq = seq[barcode_length + umi_length :]
        quality_scores = quality_scores[barcode_length + umi_length :]

    header = f"{header} {barcode}:{umi}"

    return (
        header,
        seq,
        quality_scores,
    )


def trim_by_regex(
    read: fastq.FASTQRecord, 
    features: fasta.FASTAish,
    tolerance = 3,
) -> fastq.FASTQRecord:
    
    header, seq, quality_scores = read

    fuzzy_pattern = f"{{e<={tolerance}}}"
    feature_matches = {
        feature_name: match
        for feature_name, feature_seq in features
        if (match := regex.search(f"{feature_seq}{fuzzy_pattern}", seq))
    }

    if any(feature_matches):
        match_name, top_match = min(
            feature_matches.items(), key=lambda x: sum(x[1].fuzzy_counts)
        )

        header = f"{header} {match_name}"
        seq = top_match.group()
        quality_scores = quality_scores[top_match.start(): top_match.end()]

    return (
        header,
        seq,
        quality_scores,
    )


def trim_by_index(
    read: fastq.FASTQRecord, 
    features_lookup: Mapping[str, str],
    region: tuple[str, str]
) -> fastq.FASTQRecord:

    header, seq, quality_scores = read

    start, end = region
    seq = seq[start: end]
    quality_scores = quality_scores[start: end]

    feature_name = features_lookup.get(seq)
    if feature_name:
        header = f"{header} {feature_name}"

    return (
        header,
        seq,
        quality_scores,
    )


def filter_duplicates(reads: fastq.FASTQish) -> fasta.FASTAStream:
    observed = set()
    for header, seq, _ in reads:
        key = get_barcodes(header)

        if key not in observed:
            yield (header, seq)

            observed.add(key)


def quantify(
    seqs: fasta.FASTAish, features: fasta.FASTAish
) -> dict[str, dict[str, int]]:
    
    umis_by_barcode_feature = defaultdict(lambda: defaultdict(set))
    for header, _ in seqs:
        barcode, umi = get_barcodes(header)

        feature_matches = [
            feature_name for feature_name, _ in features if feature_name in header
        ]

        if feature_matches:
            assert len(feature_matches) == 1
            umis_by_barcode_feature[barcode][feature_matches[0]].add(umi)

    return {
        barcode: {feature: len(umis) for feature, umis in features.items()}
        for barcode, features in umis_by_barcode_feature.items()
    }


def get_barcodes(header: str) -> tuple[str, str]:
    pattern = rf" ([{fastq.IUPAC_DNA}]+):([{fastq.IUPAC_DNA}]+)( |$)"

    pattern_match = regex.search(pattern, header)

    if pattern_match:
        barcode = pattern_match.group(1)
        umi = pattern_match.group(2)
    else:
        raise ValueError

    return (barcode, umi)
