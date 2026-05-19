#!/usr/bin/env python3
from __future__ import annotations

import argparse
import math
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Tuple

HDR_RE = re.compile(r"^\s*(\d+)\s*,\s*(\d+)\s*,\s*1\s*,\s*$")
NODE_RE = re.compile(r"^\s*(\d+)\s*,\s*([+-]?[0-9.]+(?:E[+-]?\d+)?)\s*,\s*$")


@dataclass
class VecRec:
    set_id: int
    vec_id: int
    lines: List[str]


def lines_of(path: Path) -> List[str]:
    return path.read_text(encoding="utf-8", errors="ignore").splitlines()


def first_450_start(lines: List[str]) -> int:
    for i, l in enumerate(lines):
        if l.strip() == "450" and i >= 2 and lines[i - 1].strip() == "-1" and lines[i - 2].strip() == "-1":
            return i - 2
    return -1


def geom_prefix(geom_neu: Path) -> List[str]:
    lines = lines_of(geom_neu)
    i = first_450_start(lines)
    return lines if i < 0 else lines[:i]


def parse_vecs(neu: Path, set_id: int) -> List[VecRec]:
    lines = lines_of(neu)
    out: List[VecRec] = []
    i = 0
    while i < len(lines):
        m = HDR_RE.match(lines[i].strip())
        if not m:
            i += 1
            continue
        sid = int(m.group(1))
        vid = int(m.group(2))
        rec = [lines[i]]
        i += 1
        while i < len(lines):
            rec.append(lines[i])
            if lines[i].strip().startswith("-1,"):
                i += 1
                break
            i += 1
        if sid == set_id:
            out.append(VecRec(sid, vid, rec))
    return out


def parse_freqs(f06: Path, n: int = 4) -> List[float]:
    txt = f06.read_text(encoding="utf-8", errors="ignore")
    vals: List[float] = []
    for m in re.finditer(r"^\s*\d+\s+\d+\s+[0-9.E+-]+\s+[0-9.E+-]+\s+([0-9.E+-]+)", txt, re.M):
        vals.append(float(m.group(1)))
        if len(vals) >= n:
            break
    while len(vals) < n:
        vals.append(float(len(vals) + 1))
    return vals


def rename_title(s: str) -> str:
    k = s.strip()
    mp = {
        "RSS translation": "Total Translation",
        "T1  translation": "T1 Translation",
        "T2  translation": "T2 Translation",
        "T3  translation": "T3 Translation",
        "RSS rotation": "Total Rotation",
        "R1  rotation": "R1 Rotation",
        "R2  rotation": "R2 Rotation",
        "R3  rotation": "R3 Rotation",
    }
    return mp.get(k, k)


def node_map(rec: VecRec) -> Dict[int, float]:
    m: Dict[int, float] = {}
    for l in rec.lines:
        mm = NODE_RE.match(l.strip())
        if mm:
            nid = int(mm.group(1))
            if nid > 0:
                m[nid] = float(mm.group(2))
    return m


def replace_ref_ids(line: str) -> str:
    parts = [p.strip() for p in line.split(",")]
    out = []
    for p in parts:
        if not p:
            continue
        try:
            iv = int(p)
            if 10001 <= iv <= 10008:
                out.append(str(iv - 10000))
            else:
                out.append(str(iv))
        except ValueError:
            out.append(p)
    return ",".join(out) + ","


def normalize_set(vecs: List[VecRec], target_set: int) -> List[List[str]]:
    out: List[List[str]] = []
    for i, rec in enumerate(vecs, start=1):
        lines = rec.lines[:]
        lines[0] = f"{target_set},{i},1,"
        if len(lines) > 1:
            lines[1] = rename_title(lines[1])
        if len(lines) > 4:
            lines[3] = replace_ref_ids(lines[3])
            lines[4] = replace_ref_ids(lines[4])
        out.append(lines)
    return out


def combine(rx: List[VecRec], ry: List[VecRec], set_id: int, name: str, a: float, b: float) -> Tuple[List[str], List[List[str]]]:
    if len(rx) != len(ry):
        raise RuntimeError("RSX/RSY vector count mismatch")
    set_block = [
        f"{set_id},",
        name,
        "36,12,",
        "0.,",
        "4,",
        "From: MYSTRAN",
        "Date : 2026-05-10",
        "<NULL>",
        "SEISMIC ANALYSIS",
    ]
    out_vecs: List[List[str]] = []
    for i, (vx, vy) in enumerate(zip(rx, ry), start=1):
        lx = vx.lines[:]
        ly = vy.lines[:]
        nx = node_map(vx)
        ny = node_map(vy)
        rec = lx[:]
        rec[0] = f"{set_id},{i},1,"
        rec[1] = rename_title(rec[1])
        rec[3] = replace_ref_ids(rec[3])
        rec[4] = replace_ref_ids(rec[4])
        vals = []
        for j, l in enumerate(rec):
            m = NODE_RE.match(l.strip())
            if not m:
                continue
            nid = int(m.group(1))
            if nid <= 0:
                continue
            v = a * nx.get(nid, 0.0) + b * ny.get(nid, 0.0)
            rec[j] = f"{nid},{v:.7E},"
            vals.append(v)
        if vals:
            vmin, vmax = min(vals), max(vals)
            vabs = max(abs(vmin), abs(vmax))
            rec[2] = f"{vmin:.7E},{vmax:.7E},{vabs:.7E},"
        out_vecs.append(rec)
    return set_block, out_vecs


def flatten(records: List[List[str]]) -> List[str]:
    o: List[str] = []
    for r in records:
        o.extend(r)
    return o


def parse_combo_spec(spec: str, sid: int) -> Tuple[int, str, float, float]:
    # format: NAME:a:b
    parts = spec.split(":")
    if len(parts) != 3:
        raise ValueError(f"invalid combo spec '{spec}'. expected NAME:a:b")
    name = parts[0].strip()
    if not name:
        raise ValueError(f"invalid combo spec '{spec}'. empty NAME")
    try:
        a = float(parts[1].strip())
        b = float(parts[2].strip())
    except ValueError as exc:
        raise ValueError(f"invalid combo spec '{spec}'. a/b must be numeric") from exc
    return sid, name, a, b


def main() -> None:
    ap = argparse.ArgumentParser(description="Build FEMAP RS NEU (geometry+modal+RS+combos)")
    ap.add_argument("--geometry", required=True, type=Path)
    ap.add_argument("--modal-neu", required=True, type=Path)
    ap.add_argument("--modal-f06", required=True, type=Path)
    ap.add_argument("--rsx-neu", required=True, type=Path)
    ap.add_argument("--rsy-neu", required=True, type=Path)
    ap.add_argument("--out", required=True, type=Path)
    ap.add_argument(
        "--no-default-combos",
        action="store_true",
        help="disable built-in combo set generation",
    )
    ap.add_argument(
        "--combo",
        action="append",
        default=[],
        help="custom combo in form NAME:a:b (repeatable), computed as a*RSX + b*RSY",
    )
    args = ap.parse_args()

    prefix = geom_prefix(args.geometry)
    freqs = parse_freqs(args.modal_f06)

    m1 = parse_vecs(args.modal_neu, 1)
    m2 = parse_vecs(args.modal_neu, 2)
    m3 = parse_vecs(args.modal_neu, 3)
    m4 = parse_vecs(args.modal_neu, 4)
    rsx = parse_vecs(args.rsx_neu, 1)
    rsy = parse_vecs(args.rsy_neu, 1)

    sets_450: List[str] = []
    base_defs = [
        (1, f"Mode 1, {freqs[0]:.6g} Hz", "MODAL ANALYSIS"),
        (2, f"Mode 2, {freqs[1]:.6g} Hz", "MODAL ANALYSIS"),
        (3, f"Mode 3, {freqs[2]:.6g} Hz", "MODAL ANALYSIS"),
        (4, f"Mode 4, {freqs[3]:.6g} Hz", "MODAL ANALYSIS"),
        (5, "RS X Direction", "SEISMIC ANALYSIS"),
        (6, "RS Y Direction", "SEISMIC ANALYSIS"),
    ]
    for sid, name, atype in base_defs:
        sets_450.extend([f"{sid},", name, "36,12,", "0.,", "4,", "From: MYSTRAN", "Date : 2026-05-10", "<NULL>", atype])

    combo_specs: List[Tuple[int, str, float, float]] = []
    next_sid = 7
    if not args.no_default_combos:
        combo_specs += [
            (next_sid, "RSX+RSY Combo", 1.0, 1.0),
            (next_sid + 1, "RSX-RSY Combo", 1.0, -1.0),
            (next_sid + 2, "-RSX-RSY Combo", -1.0, -1.0),
            (next_sid + 3, "-RSX+RSY Combo", -1.0, 1.0),
            (next_sid + 4, "1.0RSX+0.3RSY", 1.0, 0.3),
            (next_sid + 5, "0.3RSX+1.0RSY", 0.3, 1.0),
            (next_sid + 6, "1.0RSX-0.3RSY", 1.0, -0.3),
            (next_sid + 7, "-0.3RSX+1.0RSY", -0.3, 1.0),
        ]
        next_sid += 8
    for spec in args.combo:
        combo_specs.append(parse_combo_spec(spec, next_sid))
        next_sid += 1

    vec_blocks: List[List[str]] = []
    vec_blocks += normalize_set(m1, 1)
    vec_blocks += normalize_set(m2, 2)
    vec_blocks += normalize_set(m3, 3)
    vec_blocks += normalize_set(m4, 4)
    vec_blocks += normalize_set(rsx, 5)
    vec_blocks += normalize_set(rsy, 6)

    for sid, name, a, b in combo_specs:
        s450, v = combine(rsx, rsy, sid, name, a, b)
        sets_450.extend(s450)
        vec_blocks += v

    out_lines: List[str] = []
    out_lines += prefix
    out_lines += ["   -1", "   -1", "   450"]
    out_lines += sets_450
    out_lines += ["   -1", "   -1", "   451"]
    out_lines += flatten(vec_blocks)
    out_lines += ["   -1"]

    args.out.write_text("\n".join(out_lines) + "\n", encoding="utf-8")
    print(args.out)


if __name__ == "__main__":
    main()
