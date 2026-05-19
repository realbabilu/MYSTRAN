#!/usr/bin/env python3
from __future__ import annotations

import argparse
import re
from pathlib import Path

SOL112_RE = re.compile(r"^\s*SOL\s+112\s*$", re.IGNORECASE)
SUBCASE_RE = re.compile(r"^\s*SUBCASE\b", re.IGNORECASE)
HAS_FREQ_RE = re.compile(r"^\s*FREQ\s*=", re.IGNORECASE)
BEGIN_BULK_RE = re.compile(r"^\s*BEGIN\s+BULK\s*$", re.IGNORECASE)
FREQ1_RE = re.compile(r"^\s*FREQ1\s*,\s*(\d+)\s*,", re.IGNORECASE)
TITLE_SOL112_RE = re.compile(r"(\bSOL\s*)112(\b)", re.IGNORECASE)
ACCEL_REQ_RE = re.compile(r"^\s*ACCELERATION\s*\(", re.IGNORECASE)


def convert(src: Path, dst: Path, freq_sid: int | None, clean_sol111: bool) -> None:
    lines = src.read_text(encoding="utf-8", errors="ignore").splitlines()

    out = []
    in_case = True
    in_subcase = False
    subcase_has_freq = False
    first_subcase_written = False
    detected_freq_sid = None

    # pre-scan bulk for FREQ1 sid if not provided
    if freq_sid is None:
        for l in lines:
            m = FREQ1_RE.match(l)
            if m:
                detected_freq_sid = int(m.group(1))
                break
    use_freq_sid = freq_sid if freq_sid is not None else (detected_freq_sid if detected_freq_sid is not None else 40)

    for i, l in enumerate(lines):
        raw = l
        s = raw.strip()

        if SOL112_RE.match(raw):
            out.append("SOL 111")
            continue
        if TITLE_SOL112_RE.search(raw):
            out.append(TITLE_SOL112_RE.sub(r"\g<1>111\2", raw))
            continue

        if BEGIN_BULK_RE.match(raw):
            # close previous subcase with auto FREQ if missing
            if in_subcase and not subcase_has_freq:
                out.append(f"  FREQ = {use_freq_sid}")
            in_case = False
            in_subcase = False
            out.append(raw)
            # add marker in bulk for compatibility traceability
            out.append("$ ! --- response_spectrum_mystran_add begin --- !")
            out.append("$ SOL112->SOL111 compatibility conversion applied")
            out.append("$ ! --- response_spectrum_mystran_add end --- !")
            continue

        if in_case:
            if SUBCASE_RE.match(raw):
                # finalize previous subcase
                if in_subcase and not subcase_has_freq:
                    out.append(f"  FREQ = {use_freq_sid}")
                in_subcase = True
                subcase_has_freq = False
                first_subcase_written = True
                out.append(raw)
                continue

            if in_subcase and HAS_FREQ_RE.match(raw):
                subcase_has_freq = True
                out.append(raw)
                continue

            if clean_sol111 and ACCEL_REQ_RE.match(raw):
                out.append("$ " + raw + "  $ cleaned by sol112_to_sol111_compat")
                continue

            out.append(raw)
            continue

        out.append(raw)

    # EOF finalize
    if in_subcase and not subcase_has_freq:
        out.append(f"  FREQ = {use_freq_sid}")

    dst.write_text("\n".join(out) + "\n", encoding="ascii")


def main() -> None:
    ap = argparse.ArgumentParser(description="Compatibility translator: SOL112-style deck -> SOL111 deck")
    ap.add_argument("input", type=Path, help="source .bdf/.dat")
    ap.add_argument("output", type=Path, help="target .bdf/.dat")
    ap.add_argument("--freq-sid", type=int, default=None, help="force FREQ sid in SUBCASE if absent")
    ap.add_argument(
        "--clean-sol111",
        action="store_true",
        help="comment unsupported/verbose SOL111 case-control requests (currently ACCELERATION)",
    )
    args = ap.parse_args()

    convert(args.input, args.output, args.freq_sid, args.clean_sol111)
    print(args.output)


if __name__ == "__main__":
    main()
