"""Build the labelled window dataset: positives + two negative sets.

Negative set A ("random"): windows drawn uniformly across the nuclear genome.
Negative set B ("distance_matched"): one window per positive, placed at a
matched distance from the nearest centromere / telomere. This removes the
genomic-location signal that dominates the random-negative comparison.

Both negative sets exclude windows within `exclusion` bp of any annotated
origin, and windows with >2% ambiguous bases.
"""

import random
from pathlib import Path

from . import genome


def load_origins(path: Path) -> list[dict]:
    """Read the OriDB TSV and return confirmed origins.

    Columns: chr, start, end, name, othernames, status.
    """
    origins = []
    with open(path, encoding="utf-8") as handle:
        header = next(handle).rstrip("\n").split("\t")
        if header[0].lstrip("\ufeff") != "chr":
            raise ValueError(f"Unexpected OriDB header: {header}")
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 6 or fields[5] != "Confirmed":
                continue
            chrom, start, end = fields[0], int(fields[1]), int(fields[2])
            if not chrom.isdigit():
                continue
            origins.append({
                "chr": chrom,
                "start": start,
                "end": end,
                "name": fields[3] or f"ARS_{chrom}_{start}",
                "mid": (start + end) // 2,
            })
    return origins


def positive_windows(origins: list[dict], seqs: dict[str, str],
                     centromeres: dict[str, int], window: int) -> list[dict]:
    """Centre a window on each origin midpoint, merging close origins."""
    by_chr: dict[str, list[dict]] = {}
    for origin in origins:
        by_chr.setdefault(origin["chr"], []).append(origin)

    records = []
    for chrom, items in by_chr.items():
        length = len(seqs[chrom])
        kept_mids: list[int] = []
        for origin in sorted(items, key=lambda x: x["mid"]):
            mid = origin["mid"]
            if mid - window // 2 < 0 or mid + window // 2 > length:
                continue
            if kept_mids and mid - kept_mids[-1] < window:
                continue  # overlapping window: keep the first origin only
            kept_mids.append(mid)
            start = mid - window // 2
            records.append({
                "window_id": f"pos_{chrom}_{mid}",
                "chr": chrom,
                "start": start,
                "end": start + window,
                "mid": mid,
                "label": 1,
                "negative_type": "positive",
                "origin_name": origin["name"],
                "sequence": seqs[chrom][start:start + window],
                "dist_bp": genome.distance_to_nearest_landmark(
                    mid, chrom, length, centromeres.get(chrom)),
            })
    return records


def _ambiguous_fraction(seq: str) -> float:
    return sum(c not in "ACGT" for c in seq) / len(seq)


def _too_close(mid: int, chrom: str, positive_mids: dict[str, list[int]],
               min_dist: int) -> bool:
    return any(abs(mid - p) < min_dist for p in positive_mids.get(chrom, ()))


def random_negatives(seqs: dict[str, str], positives: list[dict],
                     centromeres: dict[str, int], window: int, seed: int,
                     exclusion: int) -> list[dict]:
    rng = random.Random(seed)
    positive_mids: dict[str, list[int]] = {}
    for rec in positives:
        positive_mids.setdefault(rec["chr"], []).append(rec["mid"])

    chroms = [c for c in seqs if c.isdigit()]
    weights = [len(seqs[c]) for c in chroms]
    n_needed = len(positives)
    records = []
    attempts = 0
    while len(records) < n_needed and attempts < 200 * n_needed:
        attempts += 1
        chrom = rng.choices(chroms, weights=weights, k=1)[0]
        length = len(seqs[chrom])
        start = rng.randrange(0, length - window)
        mid = start + window // 2
        seq = seqs[chrom][start:start + window]
        if _ambiguous_fraction(seq) > 0.02:
            continue
        if _too_close(mid, chrom, positive_mids, exclusion):
            continue
        records.append({
            "window_id": f"negrand_{chrom}_{mid}",
            "chr": chrom,
            "start": start,
            "end": start + window,
            "mid": mid,
            "label": 0,
            "negative_type": "random",
            "origin_name": "",
            "sequence": seq,
            "dist_bp": genome.distance_to_nearest_landmark(
                mid, chrom, length, centromeres.get(chrom)),
        })
    if len(records) < n_needed:
        raise RuntimeError(
            f"Only {len(records)} random negatives generated ({n_needed} needed)")
    return records


def distance_matched_negatives(seqs: dict[str, str], positives: list[dict],
                               centromeres: dict[str, int], window: int, seed: int,
                               exclusion: int) -> list[dict]:
    """Place one negative per positive at a matched landmark distance."""
    rng = random.Random(seed + 1)
    positive_mids: dict[str, list[int]] = {}
    for rec in positives:
        positive_mids.setdefault(rec["chr"], []).append(rec["mid"])

    used_starts: set[tuple[str, int]] = set()
    records = []
    for pos in positives:
        chrom, length = pos["chr"], len(seqs[pos["chr"]])
        target = pos["dist_bp"]
        centromere = centromeres.get(chrom)
        landmarks = genome.landmark_positions(chrom, length, centromere)

        chosen = None
        # Build candidates at the target distance from each landmark, then
        # verify the *actual* distance to the nearest landmark (another
        # landmark may be closer). Tolerances widen only if needed.
        for tol in (0.02, 0.05, 0.10, 0.25, 0.50):
            candidates = []
            for landmark in landmarks:
                for sign in (-1, 1):
                    for step in (0.0, -1.0, 1.0, -0.5, 0.5):
                        candidates.append(landmark + sign * target * (1 + step * tol))
            rng.shuffle(candidates)
            for cand in candidates:
                start = int(round(cand - window / 2))
                if start < 0 or start + window > length:
                    continue
                if (chrom, start) in used_starts:
                    continue
                mid = start + window // 2
                seq = seqs[chrom][start:start + window]
                if _ambiguous_fraction(seq) > 0.02:
                    continue
                if _too_close(mid, chrom, positive_mids, exclusion):
                    continue
                actual = genome.distance_to_nearest_landmark(
                    mid, chrom, length, centromere)
                if abs(actual - target) <= max(window / 2, tol * target):
                    chosen = start
                    break
            if chosen is not None:
                break

        if chosen is None:
            # Fallback: closest safe window among a random candidate pool.
            best = None
            for _ in range(20000):
                start = rng.randrange(0, length - window)
                mid = start + window // 2
                seq = seqs[chrom][start:start + window]
                if (chrom, start) in used_starts or _ambiguous_fraction(seq) > 0.02:
                    continue
                if _too_close(mid, chrom, positive_mids, exclusion):
                    continue
                dist = genome.distance_to_nearest_landmark(mid, chrom, length, centromere)
                score = abs(dist - target)
                if best is None or score < best[0]:
                    best = (score, start)
            if best is None:
                raise RuntimeError(f"Could not place a matched negative on chr {chrom}")
            chosen = best[1]

        used_starts.add((chrom, chosen))
        mid = chosen + window // 2
        records.append({
            "window_id": f"negmatch_{chrom}_{mid}",
            "chr": chrom,
            "start": chosen,
            "end": chosen + window,
            "mid": mid,
            "label": 0,
            "negative_type": "distance_matched",
            "origin_name": "",
            "sequence": seqs[chrom][chosen:chosen + window],
            "dist_bp": genome.distance_to_nearest_landmark(
                mid, chrom, length, centromere),
        })
    return records
