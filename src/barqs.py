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
