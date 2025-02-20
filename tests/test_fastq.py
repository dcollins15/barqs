import gzip
import os
import pathlib

import fastq

working_directory = os.path.dirname(__file__)
path_to_fastq = os.path.join(working_directory, "test_data", "small.fastq")


def test_fastq_unzip(tmp_path: pathlib.Path):
    """Test that reading a .fastq with `fastq.load` and then compressing it to
    a .fastq.gz with `fastq.dump` leaves the contents otherwise unchanged.
    """

    reads = fastq.load(path_to_fastq)

    path_to_zip = tmp_path / "small.fastq.gz"
    fastq.dump(reads, path_to_zip)

    assert open(path_to_fastq).read() == gzip.open(path_to_zip, mode="rt").read()


def test_maybe_open(tmp_path: pathlib.Path):
    """Test that IO streams can be safely passed to `fastq.load` and `fastq.dump`."""

    with open(path_to_fastq) as file_in:
        reads = list(fastq.load(file_in))

    path_to_zip = tmp_path / "small.fastq.gz"
    with gzip.open(path_to_zip, "wt") as file_out:
        fastq.dump(reads, file_out)

    assert open(path_to_fastq).read() == gzip.open(path_to_zip, mode="rt").read()
