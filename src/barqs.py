from collections import defaultdict
from typing import Iterable, Iterator, Mapping

import regex

from fastq import FASTQish, FASTQRecord

FASTARecord = tuple[str, str]
FASTAStream = Iterator[FASTARecord]
FASTALike = Iterable[FASTARecord]
FASTAish = dict[str, str] | Iterable[str] | FASTALike


def extract(
    read: FASTQRecord,
    barcode_size: int,
    umi_size: int,
) -> tuple[str, str]:
    """Extracts a cell barcode and UMI from the 3' end of the FASTQ record.

    Args:
        read (FASTQRecord):
            A tuple containing the FASTQ record (header, sequence, quality scores).
        barcode_size (int):
            The length of the cell barcode sequence.
        umi_size (int):
            The length of the UMI sequence.

    Returns:
        tuple[str, str]:
            The extracted cell barcode and UMI.
    """

    _, seq, _ = read

    barcode = seq[:barcode_size]
    umi = seq[barcode_size : barcode_size + umi_size]

    return (barcode, umi)


def tag(
    read: FASTQRecord,
    barcode: str,
    umi: str,
    trim=True,
) -> FASTQRecord:
    """Append barcode and UMI to the FASTQ record header and optionally trim
    them from the sequence.

    Args:
        read (FASTQRecord):
            A tuple containing the FASTQ record (header, sequence, quality scores).
        barcode (str):
            The cell barcode sequence.
        umi (str):
            The UMI sequence.
        trim (bool, optional):
            Whether to remove the barcode and UMI from the sequence and
            quality scores (assumes they appear on the 3' end of the read).
            Defaults to True.

    Returns:
        FASTQRecord:
            The updated FASTQ record.
    """

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


def trim_by_index(
    read: FASTQRecord, feature_lookup: Mapping[str, str], region: tuple[int, int]
) -> FASTQRecord:
    """Trim the FASTQ record to the specified region and append a feature name
    to the header if there is a sequence match.

    Args:
        read (FASTQRecord):
            A tuple containing the FASTQ record (header, sequence, quality scores).
        feature_lookup (Mapping[str, str]):
            A dictionary mapping sequences to feature names.
        region (tuple[int, int]):
            [blank].

    Returns:
        FASTQRecord:
            The updated FASTQ record.
    """

    header, seq, quality_scores = read

    start, end = region
    seq = seq[start:end]
    quality_scores = quality_scores[start:end]

    feature_name = feature_lookup.get(seq)
    if feature_name:
        header = f"{header} {feature_name}"

    return (
        header,
        seq,
        quality_scores,
    )


def trim_by_regex(
    read: FASTQRecord,
    feature_lookup: Mapping[str, str],
    tolerance=3,
) -> FASTQRecord:
    """Perform a fuzzy match between the FASTQ sequence and [blank]—if there
    is a match, trim the sequence down to the corresponding region and append
    the feature name to the header.

    Args:
        read (FASTQRecord):
            A tuple containing the FASTQ record (header, sequence, quality scores).
        feature_lookup (Mapping[str, str]):
            A dictionary mapping feature sequences to their names.
        tolerance (int, optional):
            The maximum number of mismatches allowed. Defaults to 3.

    Returns:
        FASTQRecord:
            The updated FASTQ record.
    """
    header, seq, quality_scores = read

    match_lookup = {
        match: feature
        for feature in feature_lookup.keys()
        if (match := regex.search(f"(?:{feature}){{s<={tolerance}}}", seq))
    }
    if any(match_lookup):
        match = min(match_lookup.keys(), key=lambda x: x.fuzzy_counts)
        feature = match_lookup[match]
        feature_name = feature_lookup[feature]

        header = f"{header} {feature_name}"
        seq = match.group()
        quality_scores = quality_scores[match.start() : match.end()]

    return (
        header,
        seq,
        quality_scores,
    )


def filter_duplicates(reads: FASTQish) -> FASTAStream:
    """Remove duplicate sequences based on barcode information.

    Args:
        reads (FASTQish):
            An iterable of FASTQ records.

    Returns:
        FASTAStream:
            An iterator of unique FASTA records.
    """

    observed = set()
    for header, seq, _ in reads:
        key = get_barcodes(header)

        if key not in observed:
            yield (header, seq)

            observed.add(key)


def quantify(
    seqs: FASTAish,
    features: FASTAish,
) -> dict[str, dict[str, int]]:
    """Count unique UMIs per cell barcode, feature pair.

    Args:
        seqs (FASTAish):
            An iterable of FASTA records.
        features (FASTAish):
            An iterable of known feature sequences.

    Returns:
        dict[str, dict[str, int]]:
            A nested dictionary with barcode-feature pairs and UMI counts.
    """

    umis_by_barcode_feature = defaultdict(lambda: defaultdict(set))
    for header, _ in seqs:
        barcode, umi = get_barcodes(header)

        feature_matches = [
            feature_name
            # Checking for the presence of the feature name isn't robust enough
            # since it will create spurious matches if feature names are
            # subsets of each other.
            for feature_name, _ in features
            if feature_name in header
        ]

        if feature_matches:
            assert len(feature_matches) == 1
            umis_by_barcode_feature[barcode][feature_matches[0]].add(umi)

    return {
        barcode: {feature: len(umis) for feature, umis in features.items()}
        for barcode, features in umis_by_barcode_feature.items()
    }


def get_barcodes(header: str) -> tuple[str, str]:
    """Pull the cell barcode and UMI from the header of a FASTQ record.

    Args:
        header (str):
            The header string from a FASTQ record.

    Returns:
        tuple[str, str]:
            The extracted barcode and UMI.

    Raises:
        ValueError: If the header does not contain a valid barcode-UMI pair.
    """

    pattern = r" ([ATGCN]+):([ATCGN]+)( |$)"

    pattern_match = regex.search(pattern, header)

    if pattern_match:
        barcode = pattern_match.group(1)
        umi = pattern_match.group(2)
    else:
        raise ValueError(header)

    return (barcode, umi)
