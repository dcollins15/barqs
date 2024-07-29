import regex
from collections import defaultdict
from typing import Callable, Hashable, Mapping

import fasta
import fastq


def extract(
    read: fastq.FASTQRecord,
    barcode_size: int,
    umi_size: int,
) -> tuple[str, str]:

    _, seq, _ = read

    barcode = seq[:barcode_size]
    umi = seq[barcode_size : barcode_size + umi_size]

    return (barcode, umi)


def tag(
    read: fastq.FASTQRecord,
    barcode: str,
    umi: str,
    trim = True,
) -> fastq.FASTQRecord:

    header, seq, quality_scores = read
    
    if trim:
        seq = seq[len(barcode) + len(umi) :]
        quality_scores = quality_scores[len(barcode) + len(umi) :]

    header = f"{header} {barcode}:{umi}"

    return (
        header,
        seq,
        quality_scores,
    )


def trim_by_regex(
    read: fastq.FASTQRecord, 
    feature_lookup: Mapping[str, str],
    tolerance = 3,
) -> fastq.FASTQRecord:
    
    header, seq, quality_scores = read

    match_lookup = {
        match: feature
        for feature in feature_lookup.keys()
        if (match := regex.search(f"{feature}{{e<={tolerance}}}", seq))
    }
    if any(match_lookup):
        match = min(match_lookup.keys(), key=lambda x: x.fuzzy_counts)
        feature = match_lookup[match]
        feature_name = feature_lookup[feature]

        header = f"{header} {feature_name}"
        seq = match.group()
        quality_scores = quality_scores[match.start(): match.end()]

    return (
        header,
        seq,
        quality_scores,
    )


def trim_by_index(
    read: fastq.FASTQRecord, 
    feature_lookup: Mapping[str, str],
    region: tuple[str, str]
) -> fastq.FASTQRecord:

    header, seq, quality_scores = read

    start, end = region
    seq = seq[start: end]
    quality_scores = quality_scores[start: end]

    feature_name = feature_lookup.get(seq)
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

    return observed


def quantify(
    seqs: fasta.FASTAish, 
    features: fasta.FASTAish,
) -> dict[str, dict[str, int]]:
    
    umis_by_barcode_feature = defaultdict(lambda: defaultdict(set))
    for header, _ in seqs:
        barcode, umi = get_barcodes(header)

        feature_matches = [
            feature_name 
            # checking for the presence of the feature name isn't robust enough
            # will create spurious matches if feature names are subsets of eachother
            for feature_name, _ in features if feature_name in header
        ]

        if feature_matches:
            assert len(feature_matches) == 1
            umis_by_barcode_feature[barcode][feature_matches[0]].add(umi)

    return {
        barcode: {feature: len(umis) for feature, umis in features.items()}
        for barcode, features in umis_by_barcode_feature.items()
    }


def get_barcodes(header: str) -> tuple[str, str]:
    pattern = r" ([ATGCN]+):([ATCGN]+)( |$)"

    pattern_match = regex.search(pattern, header)

    if pattern_match:
        barcode = pattern_match.group(1)
        umi = pattern_match.group(2)
    else:
        raise ValueError(header)

    return (barcode, umi)
