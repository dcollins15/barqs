from typing import Iterable, Iterator, Union

from _io_utils import FileOrPathish


FASTARecord = tuple[str, str]
FASTAStream = Iterator[FASTARecord]
FASTALike = Iterable[FASTARecord]
FASTAish = Union[
    dict[str, str],
    Iterable[str],
    FASTALike,
]


def load(fp: FileOrPathish) -> FASTAStream:
    """Deserialize `fp` to an iterator of key, sequence pairs.

    ...

    Args:
        fp (FileOrPathish): 
            A path or file-like object.

    Returns:
        FASTAStream: 
            An iterator of `FASTQRecord`s.
    """


def dump(data: FASTAish, fp: FileOrPathish):
    """Serialize `data` as a FASTA formatted stream to `fp`.

    Args:
        data (FASTQish): 
            An iterable of `FASTARecord`s.
        fp: 
            A path or file-like object.
    """
