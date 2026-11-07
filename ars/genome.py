"""Genome FASTA parsing and landmark handling."""

import gzip
import re
from pathlib import Path

ROMAN = ["I", "II", "III", "IV", "V", "VI", "VII", "VIII", "IX", "X",
         "XI", "XII", "XIII", "XIV", "XV", "XVI"]


def read_fasta(path: Path) -> dict[str, str]:
    """Read a (optionally gzipped) FASTA and return {chromosome: sequence}.

    NCBI headers contain ``chromosome I`` etc.; these are mapped to "1".."16".
    The mitochondrial contig is named "MT".
    """
    path = Path(path)
    opener = gzip.open if path.suffix == ".gz" else open
    seqs: dict[str, str] = {}
    name = None
    chunks: list[str] = []
    with opener(path, "rt") as handle:
        for line in handle:
            if line.startswith(">"):
                if name is not None:
                    seqs[name] = "".join(chunks).upper()
                header = line[1:].strip()
                if "mitochondrion" in header.lower():
                    name = "MT"
                else:
                    match = re.search(r"chromosome ([IVX]+)", header)
                    name = str(ROMAN.index(match.group(1)) + 1) if match else header.split()[0]
                chunks = []
            else:
                chunks.append(line.strip())
    if name is not None:
        seqs[name] = "".join(chunks).upper()
    return seqs


def read_centromeres(path: Path) -> dict[str, int]:
    """Parse SGD_features.tab and return {chromosome: centromere midpoint}.

    The file has no header. Column indices (0-based):
    0 SGDID, 1 feature type, 2 qualifier, 3 feature name, 4 gene name,
    5 alias, 6 parent, 7 secondary SGDID, 8 chromosome, 9 start, 10 stop,
    11 strand, 12-13 dates, 14 description.
    """
    centromeres: dict[str, int] = {}
    with open(path, encoding="utf-8") as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 11 or fields[1] != "centromere":
                continue
            chrom = fields[8]
            start, stop = int(fields[9]), int(fields[10])
            centromeres[chrom] = (start + stop) // 2
    return centromeres


def landmark_positions(chrom: str, length: int, centromere: int | None) -> list[float]:
    """Return positions used for the distance feature.

    Telomeres are approximated by the chromosome ends (0 and length), which is
    the standard coarse proxy for the terminal repeat arrays.
    """
    positions = [0.0, float(length)]
    if centromere is not None:
        positions.append(float(centromere))
    return positions


def distance_to_nearest_landmark(mid: float, chrom: str, length: int,
                                 centromere: int | None) -> float:
    return min(abs(mid - p) for p in landmark_positions(chrom, length, centromere))
