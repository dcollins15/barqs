from __future__ import annotations

import contextlib
import os
import fsspec
from typing import AnyStr, IO, Optional, Protocol, Union, runtime_checkable


FileOrPathish = Union[AnyStr, os.PathLike, IO[str]]


@runtime_checkable
class IOBaseProtocol(Protocol[AnyStr]):
    """Defines a minimal interface for file-like objects.

    This protocol is useful for detecting file-like instances that do not
    exclicitly inheret from the `io.IOBase`, like those returned by `fsspec`.
    """

    def close(self) -> None: ...

    def flush(self) -> None: ...

    def isatty(self) -> bool: ...

    def readable(self) -> bool: ...

    def readline(self, limit: int = -1) -> AnyStr: ...

    def readlines(self, hint: int = -1) -> list[AnyStr]: ...

    def seek(self, offset: int, whence: int = 0) -> int: ...

    def seekable(self) -> bool: ...

    def tell(self) -> int: ...

    def truncate(self, size: Optional[int] = None) -> int: ...

    def writable(self) -> bool: ...

    def writelines(self, lines: list[AnyStr]) -> None: ...

    def read(self, size: int = -1) -> AnyStr: ...

    def write(self, s: AnyStr) -> int: ...

    def __enter__(self) -> IOBaseProtocol[AnyStr]: ...

    def __exit__(self, exc_type, exc_value, traceback) -> Optional[bool]: ...


@contextlib.contextmanager
def maybe_open(
    fp: FileOrPathish,
    mode: str = "r",
    compression: str = "infer",
    encoding: str = "utf8",
    **kwargs,
):
    """Given a path, or file-like object, return a file-like object.

    Args:
        fp (FileOrPathish):
            A path or a file-like object.
        mode (str):
            "rb", "wt", etc.
        compression (str | None):
            If given, open file using compression codec. Can either be a
            compression name (a key in fsspec.compression.compr) or “infer”
            to guess the compression from the filename suffix.
        encoding (str):
            For text mode only.
        **kwargs:
            Additional arguments passed to `fsspec.open`.

    Returns:
        A file-like object.
    """

    if isinstance(fp, IOBaseProtocol):
        yield fp

    with fsspec.open(
        fp, mode=mode, encoding=encoding, compression=compression, **kwargs
    ) as file:
        yield file


class EOFError(Exception): ...
