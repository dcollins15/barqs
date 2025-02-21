import random

import pytest

import barqs


@pytest.fixture
def barcodes() -> list[tuple[str, str]]:
    """Return a Fasta-like list of (name, sequence) tuples to use as barcodes
    for testing.
    """

    return [
        ("DAVE", "GATGCGGTGGAA"),
        ("ELIA", "GAACTGATTGCG"),
    ]


@pytest.fixture
def features() -> list[tuple[str, str]]:
    """Return a FASTA like list of (name, sequence) tuples to use as features
    for testing.
    """

    return [
        ("MAGIC", "ATGGCGGGCATTTGC"),
        ("WATER", "TGGGCGACCGAACGC"),
        ("NIGHT", "AACATTGGCCATACC"),
    ]


@pytest.fixture
def reads(
    barcodes: list[tuple[str, str]], features: list[tuple[str, str]]
) -> list[tuple[str, str, str]]:
    """Returns a FASTQ-like list of (header, sequence, quality score) tuples
    to use for testing.
    """

    # Generate three random UMI sequences 8 bases long.
    umis = ["".join(random.choices("ATGC", k=8)) for _ in range(3)]
    # Define a
    constant = "CTGATTAAAGAAAGC"

    barcode_dict = dict(barcodes)
    features_lookup = dict(features)

    # Introduce additional features to simulate simple sequencing errors.
    # Three nucleotide edits from "MAGIC".
    features_lookup["MANIC"] = "ATGGCGAATATTTGC"
    # Two nucleotide edits from "WATER".
    features_lookup["WAGER"] = "TGGGCGGGCGAACGC"
    # Two nucleotide edits from "LIGHT".
    features_lookup["LIGHT"] = "CTCATTGGCCATACC"

    seqs = [
        # "DAVE" sees "MAGIC" in a single UMI.
        barcode_dict["DAVE"] + umis[0] + constant + features_lookup["MAGIC"],
        barcode_dict["DAVE"] + umis[0] + constant + features_lookup["MAGIC"],
        # "DAVE" sees "MANIC" in a second UMI.
        barcode_dict["DAVE"] + umis[1] + constant + features_lookup["MANIC"],
        # "DAVE" sees "WATER" in a single UMI.
        barcode_dict["DAVE"] + umis[2] + constant + features_lookup["WATER"],
        # "ELIA" sees "WATER" in two UMIs.
        barcode_dict["ELIA"] + umis[0] + constant + features_lookup["WATER"],
        # The first time "ELIA" sees this UMI, the feature is "WAGER".
        barcode_dict["ELIA"] + umis[1] + constant + features_lookup["WAGER"],
        # The second time "ELIA" sees this UMI, the feature is "WATER".
        barcode_dict["ELIA"] + umis[1] + constant + features_lookup["WATER"],
        # "ELIA" sees "NIGHT" in a single UMI.
        barcode_dict["ELIA"] + umis[1] + constant + features_lookup["NIGHT"],
        # "ELIA" sees "LIGHT" in a second UMI.
        barcode_dict["ELIA"] + umis[2] + constant + features_lookup["LIGHT"],
    ]

    # Return a list of FASTQ (header, sequence, quality score) tuples.
    return [
        (
            f"READ{i}",
            seq,
            # Generate a random quality score for each read.
            "".join(
                random.choices("!#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHI", k=len(seq))
            ),
        )
        for i, seq in enumerate(seqs)
    ]


def test_barqs(reads, barcodes, features):
    """Test that the methods from `barqs` work as expected."""

    barcode_dict = dict(barcodes)

    # Extract barcode and UMI identifiers from FASTQ sequence and append them
    # onto the header.
    identifiers = (barqs.extract(read, barcode_size=12, umi_size=8) for read in reads)
    reads = (
        barqs.tag(read, barcode, umi)
        for read, (barcode, umi) in zip(reads, identifiers)
    )
    trimmed_reads = list(barqs.trim(read, start=15, end=30) for read in reads)

    # Quantify feature occurrences using the default mismatch tolerance.
    counts = barqs.quantify(trimmed_reads, features)
    assert counts[barcode_dict["DAVE"]]["MAGIC"] == 2
    assert counts[barcode_dict["DAVE"]]["WATER"] == 1
    assert "NIGHT" not in counts[barcode_dict["DAVE"]]
    assert "MAGIC" not in counts[barcode_dict["ELIA"]]
    assert counts[barcode_dict["ELIA"]]["WATER"] == 2
    assert counts[barcode_dict["ELIA"]]["NIGHT"] == 1

    # Quantify feature occurrences allowing up to 2 mismatches.
    counts = barqs.quantify(trimmed_reads, features, tolerance=2)
    assert counts[barcode_dict["DAVE"]]["MAGIC"] == 1
    assert counts[barcode_dict["DAVE"]]["WATER"] == 1
    assert "NIGHT" not in counts[barcode_dict["DAVE"]]
    assert "MAGIC" not in counts[barcode_dict["ELIA"]]
    assert counts[barcode_dict["ELIA"]]["WATER"] == 2
    assert counts[barcode_dict["ELIA"]]["NIGHT"] == 1

    # Quantify feature occurrences with strict matching (zero mismatches allowed).
    counts = barqs.quantify(trimmed_reads, features, tolerance=0)
    assert counts[barcode_dict["DAVE"]]["MAGIC"] == 1
    assert counts[barcode_dict["DAVE"]]["WATER"] == 1
    assert "NIGHT" not in counts[barcode_dict["DAVE"]]
    assert "MAGIC" not in counts[barcode_dict["ELIA"]]
    assert counts[barcode_dict["ELIA"]]["WATER"] == 2
    assert "NIGHT" not in counts[barcode_dict["ELIA"]]
