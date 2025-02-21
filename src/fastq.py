from typing import Iterable, Iterator, Mapping

from _io_utils import FileOrPathish, maybe_open

FASTQRecord = tuple[str, str, str]
FASTQStream = Iterator[FASTQRecord]
FASTQLike = Iterable[FASTQRecord]
FASTQish = Mapping[str, tuple[str, str]] | FASTQLike


def load(fp: FileOrPathish) -> FASTQStream:
    """Deserialize `fp` to an iterator of FASTQ records.

    Args:
        fp (FileOrPathish):
            A path or file-like object.

    Returns:
        FASTQStream:
            An iterator of (header, sequence, quality score) tuples.
    """

    with maybe_open(fp) as fastq_file:
        if fastq_file.read(0) != "":
            raise IOError("FASTQ files must be opened in text mode")

        fastq_by_line = enumerate(fastq_file)

        try:
            line_number, next_header = next(fastq_by_line)
        except StopIteration:
            raise EOFError("File is empty")

        while True:
            header = next_header
            if not header.startswith("@"):
                raise ValueError(f"Invalid header at line {line_number + 1}")

            header = header[1:].rstrip()

            seq_list = []
            for line_number, next_line in fastq_by_line:
                if next_line.startswith("+"):
                    break
                seq_list.append(next_line.rstrip())
            else:
                raise EOFError(f"Missing quality header: @{header}")

            seq = "".join(seq_list)

            quality_header = next_line[1:].rstrip()
            if quality_header and quality_header != header:
                raise EOFError(f"Header mismatch: {header} != {quality_header}")

            seq_length = len(seq)
            quality_list = []
            for line_number, next_line in fastq_by_line:
                if (
                    next_line.startswith("@")
                    and sum(len(quality_line) for quality_line in quality_list)
                    >= seq_length
                ):
                    break
                quality_list.append(next_line.rstrip())
            else:
                if quality_list:
                    next_line = None
                else:
                    raise EOFError(f"Empty quality scores: @{header}")

            quality_scores = "".join(quality_list)
            if seq_length != (quality_length := len(quality_scores)):
                raise EOFError(
                    f"Lengths of sequence and quality scores differs "
                    f"({seq_length} != {quality_length}): @{header}"
                )

            yield (header, seq, quality_scores)

            if next_line is None:
                break
            else:
                next_header = next_line


def dump(data: FASTQish, fp: FileOrPathish):
    """Serialize `data` as a FASTQ formatted stream to `fp`.

    Args:
        data (FASTQish):
            An iterable of FASTQ records.
        fp:
            A path or file-like object.
    """

    with maybe_open(fp, "wt") as fastq_file:
        for header, seq, quality_scores in data:
            fastq_file.write(
                f"@{header}\n" + f"{seq}\n" + "+\n" + f"{quality_scores}\n"
            )
