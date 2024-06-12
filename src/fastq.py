import os

from typing import AnyStr, IO, Iterable, Iterator, Union


FileOrPathish = Union[AnyStr, os.PathLike, IO]

FASTQRecord = tuple[str, str, str]
FASTQStream = Iterator[FASTQRecord]
FASTQLike = Iterable[FASTQRecord]
FASTQish = Union[dict[str, tuple[str, str]], FASTQLike]


def load(fp: FileOrPathish) -> FASTQStream:
    """Deserialize `fp` to an iterator of key, sequence, quality score triples.

    Args:
        fp (FileOrPathish): 
            A path or file-like object.

    Returns:
        FASTQStream: 
            An iterator of `FASTQRecord`s.
    """


def dump(data: FASTQish, fp: FileOrPathish):
    """Serialize `data` as a FASTQ formatted stream to `fp`.

    Args:
        data (FASTQish): 
            An iterable of `FASTQRecord`s.
        fp: 
            A path or file-like object.
    """
