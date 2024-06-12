import os
import gzip
import pathlib
import pytest

import fastq

working_directory = os.path.dirname(__file__)
path_to_fastq = os.path.join(working_directory, "test_data", "small.fastq")


def test_fastq_unzip(tmp_path: pathlib.Path):
    """Test that loading .fastq.gz and then dumping it to an uncompressed
    .fastq file leaves the contents otherwise unchanged.
    """

    reads = fastq.load(path_to_fastq)
    fastq_out = tmp_path / "result.fastq.gz"
    fastq.dump(reads, fastq_out)

    assert open(path_to_fastq).read() == gzip.open(fastq_out, mode="rt").read()
