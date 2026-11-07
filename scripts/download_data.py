"""Download all public inputs used by the pipeline.

Sources
-------
* Reference genome: NCBI RefSeq GCF_000146045.2 (S. cerevisiae S288C, R64/sacCer3)
* ARS coordinates:  OriDB (S. cerevisiae), status == "Confirmed"
* Centromeres:      SGD_features.tab (SGD)

Raw files land in data/raw/ and are git-ignored. A manifest with URLs, byte
sizes, SHA-256 digests and the download date is written next to them.
"""

import argparse
import gzip
import hashlib
import json
import re
import sys
from datetime import date
from pathlib import Path
from urllib.parse import urlencode

import requests

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from ars import config  # noqa: E402

HEADERS = {"User-Agent": "yeast-ars-prediction/0.1 (academic use)"}
TIMEOUT = 60

SOURCES = [
    {
        "name": "genome_fasta",
        "url": config.GENOME_URL,
        "path": config.GENOME_FILE,
        "description": "S288C R64 reference genome, RefSeq GCF_000146045.2",
    },
    {
        "name": "assembly_report",
        "url": config.ASSEMBLY_REPORT_URL,
        "path": config.ASSEMBLY_REPORT_FILE,
        "description": "RefSeq assembly report (sequence-name mapping)",
    },
    {
        "name": "sgd_features",
        "url": config.SGD_FEATURES_URL,
        "path": config.SGD_FEATURES_FILE,
        "description": "SGD chromosomal features (centromere coordinates)",
    },
]


def sha256_of(path: Path) -> str:
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def download(url: str, dest: Path, force: bool = False) -> None:
    if dest.exists() and not force:
        print(f"  exists, skipping: {dest.name}")
        return
    dest.parent.mkdir(parents=True, exist_ok=True)
    last_error = None
    for attempt in range(1, 4):
        try:
            with requests.get(url, headers=HEADERS, timeout=TIMEOUT,
                              stream=True) as response:
                response.raise_for_status()
                tmp = dest.with_suffix(dest.suffix + ".part")
                with open(tmp, "wb") as handle:
                    for chunk in response.iter_content(1 << 16):
                        handle.write(chunk)
                tmp.replace(dest)
            print(f"  downloaded: {dest.name} ({dest.stat().st_size:,} bytes)")
            return
        except Exception as exc:  # noqa: BLE001 - report and retry
            last_error = exc
            print(f"  attempt {attempt} failed for {url}: {exc}", file=sys.stderr)
    raise RuntimeError(f"Could not download {url}: {last_error}")


def download_oridb(dest: Path, force: bool = False) -> None:
    """Fetch the OriDB origin table and unwrap the HTML <textarea> payload."""
    if dest.exists() and not force:
        print(f"  exists, skipping: {dest.name}")
        return
    response = requests.get(config.ORIDB_BASE, params=config.ORIDB_PARAMS,
                            headers=HEADERS, timeout=90)
    response.raise_for_status()
    match = re.search(r'<textarea id="output"[^>]*>(.*?)</textarea>',
                      response.text, re.S)
    if not match:
        raise RuntimeError(
            "OriDB did not return the expected table payload; the download "
            "endpoint may have changed. See docs/data-card.md for fallbacks.")
    payload = match.group(1).strip() + "\n"
    dest.parent.mkdir(parents=True, exist_ok=True)
    dest.write_text(payload, encoding="utf-8")
    print(f"  downloaded: {dest.name} ({dest.stat().st_size:,} bytes)")


def validate_genome(path: Path) -> None:
    with gzip.open(path, "rt") as handle:
        seen = set()
        for line in handle:
            if line.startswith(">"):
                seen.add(line.split()[0])
            if len(seen) >= 17:
                break
    if not seen:
        raise RuntimeError("Genome FASTA looks empty")


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--force", action="store_true",
                        help="re-download files even if they already exist")
    args = parser.parse_args(argv)

    config.RAW_DIR.mkdir(parents=True, exist_ok=True)
    print("Downloading reference inputs")

    manifest = {"downloaded_on": date.today().isoformat(), "files": []}
    for source in SOURCES:
        print(f"- {source['name']}")
        download(source["url"], source["path"], force=args.force)
        manifest["files"].append({
            "name": source["name"],
            "description": source["description"],
            "url": source["url"],
            "path": str(source["path"].relative_to(config.ROOT)),
            "bytes": source["path"].stat().st_size,
            "sha256": sha256_of(source["path"]),
        })

    print("- oridb_origins")
    download_oridb(config.ORIDB_FILE, force=args.force)
    manifest["files"].append({
        "name": "oridb_origins",
        "description": "OriDB confirmed/likely/dubious origin coordinates",
        "url": f"{config.ORIDB_BASE}?{urlencode(config.ORIDB_PARAMS)}",
        "path": str(config.ORIDB_FILE.relative_to(config.ROOT)),
        "bytes": config.ORIDB_FILE.stat().st_size,
        "sha256": sha256_of(config.ORIDB_FILE),
    })

    validate_genome(config.GENOME_FILE)
    manifest_path = config.RAW_DIR / "download_manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    print(f"Manifest written to {manifest_path}")
    print("Raw data is git-ignored; see docs/data-card.md for provenance.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
