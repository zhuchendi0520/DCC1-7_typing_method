#!/usr/bin/env python3
"""Assign an M. abscessus isolate to DCC1-DCC7 using clade-defining SNPs."""

from __future__ import annotations

import argparse
import csv
from collections import defaultdict
from pathlib import Path


EXPECTED_DCCS = {
    "ABS": {"DCC1", "DCC2", "DCC4", "DCC5"},
    "MAS": {"DCC3", "DCC6", "DCC7"},
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("markers", type=Path, help="Subspecies-specific marker TSV")
    parser.add_argument("typesnp", type=Path, help="Position/Ref/Var sample TSV")
    parser.add_argument(
        "--subspecies",
        required=True,
        choices=sorted(EXPECTED_DCCS),
        help="ABS for subsp. abscessus or MAS for subsp. massiliense",
    )
    parser.add_argument(
        "--threshold",
        type=float,
        default=95.0,
        help="Minimum marker match percentage for assignment (default: 95)",
    )
    parser.add_argument("--sample", help="Sample name (default: inferred from input)")
    return parser.parse_args()


def load_markers(path: Path, allowed_dccs: set[str]) -> dict[str, set[tuple[int, str]]]:
    markers: dict[str, set[tuple[int, str]]] = defaultdict(set)
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"Lineage", "Position", "Derived_Call"}
        if not reader.fieldnames or not required.issubset(reader.fieldnames):
            raise ValueError(f"Marker file must contain: {', '.join(sorted(required))}")
        for row in reader:
            lineage = row["Lineage"].strip()
            if lineage not in allowed_dccs:
                raise ValueError(
                    f"Unexpected lineage {lineage!r} for the selected subspecies"
                )
            markers[lineage].add((int(row["Position"]), row["Derived_Call"].upper()))

    missing = allowed_dccs.difference(markers)
    if missing:
        raise ValueError(f"Marker file is missing DCCs: {', '.join(sorted(missing))}")
    return markers


def load_sample_calls(path: Path) -> set[tuple[int, str]]:
    calls: set[tuple[int, str]] = set()
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"Position", "Var"}
        if not reader.fieldnames or not required.issubset(reader.fieldnames):
            raise ValueError(f"Typing file must contain: {', '.join(sorted(required))}")
        for row in reader:
            calls.add((int(row["Position"]), row["Var"].upper()))
    return calls


def infer_sample_name(path: Path) -> str:
    name = path.name
    for suffix in (".cns.typesnp.tsv", ".typesnp.tsv", ".typesnp"):
        if name.endswith(suffix):
            return name[: -len(suffix)]
    return path.stem


def main() -> None:
    args = parse_args()
    allowed_dccs = EXPECTED_DCCS[args.subspecies]
    markers = load_markers(args.markers, allowed_dccs)
    calls = load_sample_calls(args.typesnp)

    scores = []
    for lineage, marker_set in markers.items():
        matched = len(calls.intersection(marker_set))
        total = len(marker_set)
        percentage = 100.0 * matched / total
        scores.append((lineage, percentage, matched, total))

    best_percentage = max(score[1] for score in scores)
    best_matches = sorted(score for score in scores if score[1] == best_percentage)
    sample = args.sample or infer_sample_name(args.typesnp)

    print("Sample\tSubspecies\tBest_match\tMatch_percentage\tMatched_markers\tTotal_markers\tAssignment")
    for lineage, percentage, matched, total in best_matches:
        assignment = lineage if percentage >= args.threshold else "Non_DCC"
        print(
            f"{sample}\t{args.subspecies}\t{lineage}\t{percentage:.2f}"
            f"\t{matched}\t{total}\t{assignment}"
        )


if __name__ == "__main__":
    main()
