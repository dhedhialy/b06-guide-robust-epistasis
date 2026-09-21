"""Merge/audit provenance (plan Sec 17): hashes, versions, seeds, commit."""

from __future__ import annotations

import hashlib
import json
import platform
import subprocess
import sys
from pathlib import Path

SOURCES = {
    "raw/MOESM3_design.xlsx": {
        "url": "https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fs41467-025-67256-9/MediaObjects/41467_2025_67256_MOESM3_ESM.xlsx",
        "role": "Burgold pilot library: design + counts (Supplementary Data 1)",
    },
    "raw/MOESM6_results.xlsx": {
        "url": "https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fs41467-025-67256-9/MediaObjects/41467_2025_67256_MOESM6_ESM.xlsx",
        "role": "Figure source data",
    },
}


def _sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def _commit() -> str:
    try:
        return subprocess.check_output(
            ["git", "rev-parse", "HEAD"], stderr=subprocess.DEVNULL, text=True
        ).strip()
    except Exception:
        return "uncommitted"


def provenance(repo: Path) -> dict:
    repo = Path(repo)
    report = {
        "method": "B06 v0.1.0",
        "git_commit": _commit(),
        "python": platform.python_version(),
        "platform": platform.platform(),
        "package_versions": {
            name: __import__(name).__version__
            for name in ("pandas", "numpy", "openpyxl")
        },
        "seed": 0,
        "estimand": "gene-level interaction via matched dual-guide factorial contrast on FC_500x (log2 FC D14/D3), HT-29 pilot",
        "sources": {k: {**v, "sha256": _sha256(repo / k)} for k, v in SOURCES.items()},
        "artifacts": {},
    }
    for p in sorted((repo / "results").rglob("*")):
        if p.is_file():
            rel = p.relative_to(repo).as_posix()
            report["artifacts"][rel] = {"sha256": _sha256(p), "bytes": p.stat().st_size}
    return report


def main(argv=None) -> int:
    import argparse

    ap = argparse.ArgumentParser()
    ap.add_argument("--repo", default=".")
    ap.add_argument("--out", default="results/provenance.json")
    args = ap.parse_args(argv)
    Path(args.out).write_text(json.dumps(provenance(args.repo), indent=2) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())