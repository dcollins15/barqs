from __future__ import annotations

import contextlib
import io
import os
import typing

import fsspec

FileOrPathish = typing.IO | os.PathLike | str


@contextlib.contextmanager
def maybe_open(
    fp: FileOrPathish,
    mode: str = "r",
    compression: str = "infer",
    encoding: str = "utf8",
    **kwargs,
) -> typing.Generator[typing.IO]:
    """Given a path, or file-like object, return a file-like object.

    Args:
        fp (FileOrPathish):
            A path or a file-like object.
        mode (str):
            "rb", "wt", etc.
        compression (str, optional):
            If given, open file using compression codec. Can either be a
            compression name (a key in ``fsspec.compression.compr``) or “infer”
            to guess the compression from the filename suffix.
            Defaults to "infer".
        encoding (str):
            For text mode only.
        **kwargs:
            Additional arguments passed to ``fsspec.open``.

    Returns:
        IO: A file-like object.
    """

    # If `fp` has already been opened, hand it back.
    if isinstance(fp, io.IOBase):
        yield typing.cast(typing.IO, fp)
    else:
        # Otherwise, open the file and return its contents.
        with fsspec.open_files(
            [fp], mode=mode, encoding=encoding, compression=compression, **kwargs
        ) as (file,):
            yield file
