import random
import pytest

import barqs
from fastq import IUPAC_DNA


@pytest.fixture
def barcodes() -> list[tuple[str, str]]:
    return [
        ("DAVE", "GATGCGGTGGAA"),
        ("ELIA", "GAACTGATTGCG"),
    ]


@pytest.fixture
def features() -> list[tuple[str, str]]:
    return [
        ("MAGIC", "ATGGCGGGCATTTGC"),
        ("WATER", "TGGGCGACCGAACGC"),
        ("NIGHT", "AACATTGGCCATACC"),
    ]


@pytest.fixture
def fake_reads(
    barcodes: list[tuple[str, str]], features: list[tuple[str, str]]
) -> list[tuple[str, str, str]]:

    umis = ["".join(random.choices(IUPAC_DNA, k=8)) for _ in range(2)]
    constant = "CTGATTAAAGAAAGC"

    barcode_dict = dict(barcodes)
    features_dict = dict(features)

    seqs = [
        barcode_dict["DAVE"] + umis[0] + constant + features_dict["MAGIC"],
        barcode_dict["DAVE"] + umis[0] + constant + features_dict["MAGIC"],
        barcode_dict["DAVE"] + umis[1] + constant + features_dict["WATER"],
        barcode_dict["ELIA"] + umis[0] + constant + features_dict["WATER"],
        barcode_dict["ELIA"] + umis[1] + constant + features_dict["WATER"],
        barcode_dict["ELIA"] + umis[1] + constant + features_dict["NIGHT"],
    ]

    return [
        (
            f"READ{i}",
            seq,
            "".join(
                random.choices("!#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHI", k=len(seq))
            ),
        )
        for i, seq in enumerate(seqs)
    ]


def test_basic_workflow(fake_reads, barcodes, features):
    reads = (
        barqs.extract(read, barcode_length=12, umi_length=8) 
        for read in fake_reads
    )
    trimmed_reads = (barqs.trim(read, features=features) for read in reads)
    observed_features = barqs.filter_duplicates(trimmed_reads)
    counts = barqs.quantify(observed_features, features)

    barcode_dict = dict(barcodes)
    assert counts[barcode_dict["DAVE"]]["MAGIC"] == 1
    assert counts[barcode_dict["DAVE"]]["WATER"] == 1
    assert counts[barcode_dict["ELIA"]]["WATER"] == 2
