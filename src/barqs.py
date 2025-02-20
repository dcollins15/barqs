from typing import Mapping

from fastq import FASTQRecord


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
