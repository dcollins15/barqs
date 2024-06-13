from cutadapt.adapters import
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
