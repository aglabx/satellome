"""Reader for the UCSC .2bit sequence format.

Why a reader of our own: the alternatives (`py2bit`, `twoBitToFa`) are a C
extension and an external binary respectively, and both would become one more
thing that is "not installed" on a cluster node. The format is small and fully
specified, so reading it directly costs less than depending on either.

Layout (all integers unsigned, little-endian unless the signature says
otherwise; see http://genome.ucsc.edu/FAQ/FAQformat.html#format7):

    header      signature uint32 = 0x1A412743, version, sequenceCount, reserved
    index       per sequence: nameSize uint8, name, offset uint32
    record      dnaSize, nBlockCount, nBlockStarts[], nBlockSizes[],
                maskBlockCount, maskBlockStarts[], maskBlockSizes[],
                reserved, packedDna

DNA is packed four bases per byte, most significant bits first, with
T=0, C=1, A=2, G=3. Runs of N are not stored in the packed stream at all: they
occupy positions (as whatever bases happen to be encoded there) and are listed
separately as N blocks, which is why decoding *must* overwrite them afterwards.
Mask blocks mark soft-masked (lowercase) regions.
"""

import struct
from pathlib import Path
from typing import Dict, Iterator, List, Optional, Tuple

SIGNATURE = 0x1A412743
SIGNATURE_SWAPPED = 0x4327411A
# T=0, C=1, A=2, G=3 - the order is part of the format, not a choice.
BASES = "TCAG"
DEFAULT_CHUNK = 1 << 22  # bases per decoded chunk (~4 MB of output)
DEFAULT_LINE_WIDTH = 60


class TwoBitError(Exception):
    """Malformed or unsupported .2bit file."""


def _decode_table() -> List[str]:
    """byte value -> its four bases, precomputed (256 entries)."""
    table = []
    for value in range(256):
        table.append(
            BASES[(value >> 6) & 3]
            + BASES[(value >> 4) & 3]
            + BASES[(value >> 2) & 3]
            + BASES[value & 3]
        )
    return table


_BYTE_TO_BASES = _decode_table()


def is_twobit(path) -> bool:
    """True if the file starts with a .2bit signature, whatever its name."""
    try:
        with open(str(path), "rb") as handle:
            raw = handle.read(4)
    except OSError:
        return False
    if len(raw) < 4:
        return False
    value = struct.unpack("<I", raw)[0]
    return value in (SIGNATURE, SIGNATURE_SWAPPED)


class TwoBitFile:
    """Random-access reader over a .2bit file.

    Sequences are decoded in chunks, so a 250 Mb chromosome never has to exist
    as one Python string unless the caller asks for it.
    """

    def __init__(self, path):
        self.path = str(path)
        self._handle = open(self.path, "rb")
        try:
            self._read_header()
            self._read_index()
        except Exception:
            self._handle.close()
            raise

    # -- context manager -------------------------------------------------
    def __enter__(self):
        return self

    def __exit__(self, *exc_info):
        self.close()
        return False

    def close(self) -> None:
        if not self._handle.closed:
            self._handle.close()

    # -- header / index --------------------------------------------------
    def _read_header(self) -> None:
        raw = self._handle.read(16)
        # Check the signature before the length: handing someone "truncated
        # header" for a file that is simply a FASTA sends them looking for
        # corruption that is not there.
        if len(raw) >= 4:
            signature = struct.unpack("<I", raw[:4])[0]
            if signature == SIGNATURE:
                self.byte_order = "<"
            elif signature == SIGNATURE_SWAPPED:
                self.byte_order = ">"
            else:
                raise TwoBitError(
                    f"{self.path}: not a .2bit file (signature 0x{signature:08X})"
                )
        if len(raw) < 16:
            raise TwoBitError(f"{self.path}: truncated header ({len(raw)} bytes)")

        _, self.version, self.sequence_count, _ = struct.unpack(
            self.byte_order + "IIII", raw
        )
        if self.version not in (0, 1):
            raise TwoBitError(
                f"{self.path}: unsupported .2bit version {self.version}"
            )
        # Version 1 uses 64-bit offsets. Refuse it explicitly rather than
        # reading 32 bits of a 64-bit field and seeking to a nonsense offset.
        self._offset_format = self.byte_order + ("Q" if self.version == 1 else "I")
        self._offset_size = 8 if self.version == 1 else 4

    def _read_index(self) -> None:
        self._offsets: Dict[str, int] = {}
        self.names: List[str] = []
        for _ in range(self.sequence_count):
            size_raw = self._handle.read(1)
            if not size_raw:
                raise TwoBitError(f"{self.path}: index ends early")
            name_size = size_raw[0]
            name = self._handle.read(name_size).decode("utf-8", "replace")
            offset_raw = self._handle.read(self._offset_size)
            if len(offset_raw) < self._offset_size:
                raise TwoBitError(f"{self.path}: index ends early at '{name}'")
            (offset,) = struct.unpack(self._offset_format, offset_raw)
            if name in self._offsets:
                raise TwoBitError(f"{self.path}: duplicate sequence name '{name}'")
            self._offsets[name] = offset
            self.names.append(name)

    def _read_uint32s(self, count: int) -> Tuple[int, ...]:
        raw = self._handle.read(4 * count)
        if len(raw) < 4 * count:
            raise TwoBitError(f"{self.path}: truncated block list")
        return struct.unpack(f"{self.byte_order}{count}I", raw)

    def _read_record_header(self, name: str):
        if name not in self._offsets:
            raise KeyError(f"{self.path}: no sequence named '{name}'")
        self._handle.seek(self._offsets[name])

        raw = self._handle.read(4)
        if len(raw) < 4:
            raise TwoBitError(f"{self.path}: truncated record for '{name}'")
        (dna_size,) = struct.unpack(self.byte_order + "I", raw)

        (n_count,) = struct.unpack(self.byte_order + "I", self._handle.read(4))
        n_starts = self._read_uint32s(n_count)
        n_sizes = self._read_uint32s(n_count)

        (mask_count,) = struct.unpack(self.byte_order + "I", self._handle.read(4))
        mask_starts = self._read_uint32s(mask_count)
        mask_sizes = self._read_uint32s(mask_count)

        self._handle.read(4)  # reserved
        dna_offset = self._handle.tell()
        return dna_size, list(zip(n_starts, n_sizes)), list(zip(mask_starts, mask_sizes)), dna_offset

    # -- public API ------------------------------------------------------
    def sequence_sizes(self) -> Dict[str, int]:
        """name -> length in bases, without decoding any sequence."""
        return {name: self._read_record_header(name)[0] for name in self.names}

    def iter_chunks(
        self, name: str, chunk_size: int = DEFAULT_CHUNK, soft_mask: bool = True
    ) -> Iterator[str]:
        """Yield the sequence in chunks of at most ``chunk_size`` bases."""
        dna_size, n_blocks, mask_blocks, dna_offset = self._read_record_header(name)
        if dna_size == 0:
            return

        position = 0
        while position < dna_size:
            end = min(position + chunk_size, dna_size)
            yield self._decode_range(
                dna_offset, position, end, n_blocks,
                mask_blocks if soft_mask else [],
            )
            position = end

    def read(self, name: str, soft_mask: bool = True) -> str:
        """The whole sequence as one string."""
        return "".join(self.iter_chunks(name, soft_mask=soft_mask))

    def _decode_range(self, dna_offset, start, end, n_blocks, mask_blocks) -> str:
        """Decode bases [start, end), applying N blocks and soft masking."""
        first_byte = start // 4
        last_byte = (end + 3) // 4
        self._handle.seek(dna_offset + first_byte)
        packed = self._handle.read(last_byte - first_byte)
        if len(packed) < last_byte - first_byte:
            raise TwoBitError(
                f"{self.path}: truncated DNA (wanted {last_byte - first_byte} "
                f"bytes at {dna_offset + first_byte})"
            )

        decoded = "".join(_BYTE_TO_BASES[byte] for byte in packed)
        # Trim the padding introduced by rounding to byte boundaries.
        lead = start - first_byte * 4
        bases = list(decoded[lead:lead + (end - start)])

        # N blocks overwrite whatever the packed stream held at those offsets;
        # the format stores no bases for them, so skipping this step yields
        # plausible-looking but wrong sequence.
        for block_start, block_size in n_blocks:
            overlap_start = max(block_start, start)
            overlap_end = min(block_start + block_size, end)
            for i in range(overlap_start, overlap_end):
                bases[i - start] = "N"

        for block_start, block_size in mask_blocks:
            overlap_start = max(block_start, start)
            overlap_end = min(block_start + block_size, end)
            for i in range(overlap_start, overlap_end):
                bases[i - start] = bases[i - start].lower()

        return "".join(bases)


def sc_iter_twobit(file_name, soft_mask: bool = True) -> Iterator[Tuple[str, str]]:
    """Iterate a .2bit file as ``(header, sequence)``, like sc_iter_fasta_brute.

    The header carries the leading '>' so callers can treat both readers the
    same way.
    """
    with TwoBitFile(file_name) as handle:
        for name in handle.names:
            yield f">{name}", handle.read(name, soft_mask=soft_mask)


def twobit_to_fasta(
    file_name,
    output_file,
    line_width: int = DEFAULT_LINE_WIDTH,
    soft_mask: bool = True,
) -> int:
    """Write a .2bit out as FASTA. Returns the number of sequences written.

    Streams chunk by chunk so a large genome is never held in memory, and
    publishes atomically: a reader must never see a half-converted genome and
    treat it as the input.
    """
    from satellome.core_functions.tools.atomic_io import atomic_output

    count = 0
    with TwoBitFile(file_name) as source, atomic_output(str(output_file), "w") as out:
        for name in source.names:
            out.write(f">{name}\n")
            column = 0
            for chunk in source.iter_chunks(name, soft_mask=soft_mask):
                start = 0
                while start < len(chunk):
                    take = line_width - column
                    piece = chunk[start:start + take]
                    out.write(piece)
                    column += len(piece)
                    start += len(piece)
                    if column == line_width:
                        out.write("\n")
                        column = 0
            if column:
                out.write("\n")
            count += 1
    return count
