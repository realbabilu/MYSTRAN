#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import re
from pathlib import Path

HDR_RE = re.compile(r"^\s*(\d+)\s*,\s*(\d+)\s*,\s*1\s*,\s*$")
NODE_RE = re.compile(r"^\s*(\d+)\s*,\s*([+-]?[0-9.]+(?:E[+-]?\d+)?)\s*,\s*$", re.I)


def parse_neu(path: Path) -> dict[tuple[int, int], dict[int, float]]:
    lines = path.read_text(encoding="utf-8", errors="ignore").splitlines()
    out: dict[tuple[int, int], dict[int, float]] = {}
    i = 0
    while i < len(lines):
        m = HDR_RE.match(lines[i].strip())
        if not m:
            i += 1
            continue
        sid = int(m.group(1))
        vid = int(m.group(2))
        i += 1
        vals: dict[int, float] = {}
        while i < len(lines):
            s = lines[i].strip()
            if s.startswith("-1,"):
                i += 1
                break
            nm = NODE_RE.match(s)
            if nm:
                nid = int(nm.group(1))
                if nid > 0:
                    vals[nid] = float(nm.group(2))
            i += 1
        out[(sid, vid)] = vals
    return out


def main() -> None:
    ap = argparse.ArgumentParser(description="Compare NEU set/vector nodal values")
    ap.add_argument("left", type=Path, help="left NEU")
    ap.add_argument("right", type=Path, help="right NEU")
    ap.add_argument("--abs-tol", type=float, default=0.0)
    ap.add_argument("--rel-tol", type=float, default=0.0)
    ap.add_argument("--json-out", type=Path, default=None)
    args = ap.parse_args()

    left = parse_neu(args.left)
    right = parse_neu(args.right)

    common = sorted(set(left.keys()) & set(right.keys()))
    left_only = sorted(set(left.keys()) - set(right.keys()))
    right_only = sorted(set(right.keys()) - set(left.keys()))

    max_abs = 0.0
    max_rel = 0.0
    worst = None

    for key in common:
        a = left[key]
        b = right[key]
        nids = set(a.keys()) | set(b.keys())
        for nid in nids:
            va = a.get(nid, 0.0)
            vb = b.get(nid, 0.0)
            d = abs(va - vb)
            r = d / (abs(va) + 1.0e-30)
            if d > max_abs:
                max_abs = d
                worst = {"set": key[0], "vec": key[1], "nid": nid, "left": va, "right": vb, "abs": d, "rel": r}
            if r > max_rel:
                max_rel = r

    status = "same" if (max_abs <= args.abs_tol and max_rel <= args.rel_tol and not left_only and not right_only) else "different"
    result = {
        "left": str(args.left),
        "right": str(args.right),
        "common_vectors": len(common),
        "left_only_vectors": len(left_only),
        "right_only_vectors": len(right_only),
        "max_abs_diff": max_abs,
        "max_rel_diff": max_rel,
        "worst": worst,
        "abs_tol": args.abs_tol,
        "rel_tol": args.rel_tol,
        "status": status,
    }

    print(json.dumps(result, indent=2))
    if args.json_out is not None:
        args.json_out.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
