import re
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

    if any(feature_matches):
        match_name, top_match = max(
            feature_matches.items(), key=lambda x: x[1].score if x else 0
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
    pattern = fr" ([{fastq.IUPAC_DNA}]+):([{fastq.IUPAC_DNA}]+)( |$)"

    observed = set()
    for header, seq, _ in reads:
        pattern_match = re.serach(pattern, header)

        if pattern_match:
            barcode = pattern_match.group(1)
            umi = pattern_match.group(2)
        else:
            raise ValueError()
        
        key = tuple(barcode, umi)
        if key not in observed:
            yield (header, seq)

            observed.add(key)
