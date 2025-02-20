import random
from functools import partial

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

    # Extract the cell barcode and UMI from each FASTQ sequence.
    identifiers = (barqs.extract(read, barcode_size=12, umi_size=8) for read in reads)
    # Strip the cell barcode and UMI from the FASTQ record and append them to
    # the header.
    reads = list(
        barqs.tag(read, barcode, umi)
        for read, (barcode, umi) in zip(reads, identifiers)
    )

    # Create a dictionary mapping feature sequences to their names.
    feature_lookup = {seq: name for name, seq in features}

    # Use `partial` to setup equivalent calls to `barqs.trim_by_index` and
    # `barqs.trim_by_regex`.
    ways_to_trim = [
        # Specify the region to look for exact sequence matches apriori.
        partial(barqs.trim_by_index, feature_lookup=feature_lookup, region=(15, 30)),
        # Require exact sequence matches, search using `regex`.
        partial(barqs.trim_by_regex, feature_lookup=feature_lookup, tolerance=0),
    ]
    for trim in ways_to_trim:
        # For any reads containing a feature of interest, trim it down to the
        # matching subsequence and append the feature's name to the header.
        trimmed_reads = (trim(read) for read in reads)
        # Keep the first instance of each cell barcode, UMI pair and ignore
        # any subsequent reads.
        observed_features = barqs.filter_duplicates(trimmed_reads)
        # Count the number of UMIs per unique cell barcode, feature pair.
        counts = barqs.quantify(observed_features, features)

        # Check the expected counts.
        barcode_dict = dict(barcodes)
        assert counts[barcode_dict["DAVE"]]["MAGIC"] == 1
        assert counts[barcode_dict["DAVE"]]["WATER"] == 1
        assert "NIGHT" not in counts[barcode_dict["DAVE"]]
        assert "MAGIC" not in counts[barcode_dict["ELIA"]]
        assert counts[barcode_dict["ELIA"]]["WATER"] == 1
        assert "NIGHT" not in counts[barcode_dict["ELIA"]]

    # Run `barqs.trim_by_regex` using the default `tolerance=3`
    trimmed_reads = (
        barqs.trim_by_regex(read, feature_lookup=feature_lookup) for read in reads
    )
    observed_features = barqs.filter_duplicates(trimmed_reads)
    counts = barqs.quantify(observed_features, features)

    # Check that the expected values changed.
    barcode_dict = dict(barcodes)
    # "MANIC" matched to "MAGIC"
    assert counts[barcode_dict["DAVE"]]["MAGIC"] == 2
    # "WAGER" matched to "WATER"
    assert counts[barcode_dict["ELIA"]]["WATER"] == 2
    # LIGHT matched to "LIGHT"
    assert counts[barcode_dict["ELIA"]]["NIGHT"] == 1

    # Run `barqs.trim_by_regex` using `tolerance=2`
    trimmed_reads = (
        barqs.trim_by_regex(read, feature_lookup=feature_lookup, tolerance=2)
        for read in reads
    )
    observed_features = barqs.filter_duplicates(trimmed_reads)
    counts = barqs.quantify(observed_features, features)

    # Check that the expected values changed.
    barcode_dict = dict(barcodes)
    # "MANIC" no longer matches "MAGIC"
    assert counts[barcode_dict["DAVE"]]["MAGIC"] == 1
    # The other two errors should still be within tolerance.
    assert counts[barcode_dict["ELIA"]]["WATER"] == 2
    assert counts[barcode_dict["ELIA"]]["NIGHT"] == 1
