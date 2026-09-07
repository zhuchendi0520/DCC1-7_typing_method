#!/usr/bin/env python3
"""Convert a VarScan mpileup2cns table to the SNP format used for DCC typing."""

from __future__ import annotations

import argparse
from pathlib import Path


VALID_BASES = {"A", "C", "G", "T"}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Convert a VarScan CNS file to Position/Ref/Var typing format."
    )
    parser.add_argument("cns", type=Path, help="VarScan mpileup2cns output file")
    parser.add_argument(
        "-o", "--output", type=Path, help="Output TSV (default: <input>.typesnp.tsv)"
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    output = args.output or Path(f"{args.cns}.typesnp.tsv")

    with args.cns.open() as source, output.open("w") as destination:
        destination.write("Position\tRef\tVar\n")
        for line_number, raw_line in enumerate(source, start=1):
            line = raw_line.strip()
            if not line:
                continue

            fields = line.split()
            if fields[0].lower() in {"chrom", "chromosome", "position"}:
                continue
            if len(fields) < 4:
                raise ValueError(
                    f"Malformed CNS row at line {line_number}: expected at least 4 columns"
                )

            try:
                position = int(fields[1])
            except ValueError as error:
                raise ValueError(
                    f"Invalid genomic position at line {line_number}: {fields[1]}"
                ) from error

            reference = fields[2].upper()
            consensus = fields[3].upper()
            if consensus == ".":
                consensus = reference

            if reference not in VALID_BASES:
                continue
            if consensus not in VALID_BASES:
                continue

            destination.write(f"{position}\t{reference}\t{consensus}\n")

    print(output)


if __name__ == "__main__":
    main()
