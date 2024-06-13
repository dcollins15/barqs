import re
from collections import defaultdict

from cutadapt.adapters import AnywhereAdapter

import fasta
import fastq


def extract(
    record: fastq.FASTQRecord,
    barcode_length: tuple[int, int],
    umi_length: tuple[int, int],
) -> fastq.FASTQRecord:

    header, seq, quality_scores = record

    barcode = seq[:barcode_length]
    umi = seq[barcode_length : barcode_length + umi_length]

    seq = seq[barcode_length + umi_length :]
    quality_scores = quality_scores[barcode_length + umi_length :]

    header = f"{header} {barcode}:{umi}"

    return (
        header,
        seq,
        quality_scores,
    )


def trim(record: fastq.FASTQRecord, features: fasta.FASTAish) -> fastq.FASTQRecord:

    header, seq, quality_scores = record

    feature_matches = {
        feature_name: AnywhereAdapter(feature).match_to(seq)
        for feature_name, feature in features
    }

    if any(feature_matches.values()):
        match_name, top_match = max(
            feature_matches.items(), key=lambda x: x[1].score if x[1] else 0
        )

        header = f"{header} {match_name}"

        feature_start = top_match.rstart
        feature_stop = top_match.rstop

        seq = seq[feature_start:feature_stop]

        quality_scores = quality_scores[feature_start:feature_stop]

    return (
        header,
        seq,
        quality_scores,
    )


def filter_duplicates(reads: fastq.FASTQish) -> fasta.FASTAStream:
    pattern = rf" ([{fastq.IUPAC_DNA}]+):([{fastq.IUPAC_DNA}]+)( |$)"

    observed = set()
    for header, seq, _ in reads:
        pattern_match = re.search(pattern, header)

        if pattern_match:
            barcode = pattern_match.group(1)
            umi = pattern_match.group(2)
        else:
            raise ValueError()

        key = (barcode, umi)
        if key not in observed:
            yield (header, seq)

            observed.add(key)


def quantify(
    seqs: fasta.FASTAish, features: fasta.FASTAish
) -> dict[str, dict[str, int]]:
    id_pattern = rf" ([{fastq.IUPAC_DNA}]+):([{fastq.IUPAC_DNA}]+)( |$)"

    umis_by_barcode_feature = defaultdict(lambda: defaultdict(set))
    for header, _ in seqs:
        id_match = re.search(id_pattern, header)

        if id_match:
            barcode = id_match.group(1)
            umi = id_match.group(2)
        else:
            raise ValueError()

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
