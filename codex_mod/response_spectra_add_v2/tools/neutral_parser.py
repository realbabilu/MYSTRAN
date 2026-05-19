#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import re
from pathlib import Path

SET_ID_RE = re.compile(r"^\s*(\d+)\s*,\s*$")
VEC_HDR_RE = re.compile(r"^\s*(\d+)\s*,\s*(\d+)\s*,\s*1\s*,\s*$")


def parse_neu(path: Path) -> dict:
    lines = path.read_text(encoding="utf-8", errors="ignore").splitlines()
    in_450 = False
    in_451 = False
    sets = []
    vectors = []

    i = 0
    while i < len(lines):
        s = lines[i].strip()

        if s == "450" and i >= 2 and lines[i - 1].strip() == "-1" and lines[i - 2].strip() == "-1":
            in_450 = True
            in_451 = False
            i += 1
            continue
        if s == "451" and i >= 2 and lines[i - 1].strip() == "-1" and lines[i - 2].strip() == "-1":
            in_450 = False
            in_451 = True
            i += 1
            continue

        if in_450:
            m = SET_ID_RE.match(s)
            if m and i + 1 < len(lines):
                sets.append({"set_id": int(m.group(1)), "name": lines[i + 1].strip()})
            i += 1
            continue

        if in_451:
            m = VEC_HDR_RE.match(s)
            if m:
                vectors.append({
                    "set_id": int(m.group(1)),
                    "vector_id": int(m.group(2)),
                    "title": lines[i + 1].strip() if i + 1 < len(lines) else "",
                })
            i += 1
            continue

        i += 1

    return {
        "file": str(path),
        "set_count": len(sets),
        "vector_count": len(vectors),
        "sets": sets,
        "vectors": vectors,
    }


def main() -> None:
    ap = argparse.ArgumentParser(description="Parse FEMAP neutral set/vector headers")
    ap.add_argument("neu", type=Path)
    ap.add_argument("--out", type=Path, default=None)
    args = ap.parse_args()

    result = parse_neu(args.neu)
    text = json.dumps(result, indent=2)
    if args.out:
        args.out.write_text(text + "\n", encoding="utf-8")
    else:
        print(text)


if __name__ == "__main__":
    main()
