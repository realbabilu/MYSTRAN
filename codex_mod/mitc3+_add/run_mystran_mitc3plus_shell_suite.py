from __future__ import annotations

from pathlib import Path
import argparse
import re
import sys
from typing import Iterable


ROOT = Path(r"D:\mystran2")
RUN_DIR = ROOT / "MYSTRANSolver-18.0.0" / "run_debug" / "dkmq" / "shell_static_mct"
MYSTRAN_EXE = ROOT / "MYSTRANSolver-18.0.0" / "Binaries" / "mystran.exe"
OUT_MD = ROOT / "codex_mod" / "mitc3+_add" / "mitc3plus_mystran_shell_suite.md"


BENCHMARKS = [
    {
        "label": "Static-22",
        "template": "static-22_mitc4.bdf",
        "out_bdf": "static-22_{tag}.bdf",
        "out_f06": "static-22_{tag}.F06",
        "theory": -3.09e-01,
        "extract": ("node_dof", 49, 3),
    },
    {
        "label": "Static-23",
        "template": "static-23_mitc4_case1.bdf",
        "out_bdf": "static-23_{tag}_case1.bdf",
        "out_f06": "static-23_{tag}_case1.F06",
        "theory": 4.5197e-04,
        "extract": ("max_abs_dof", 2),
    },
    {
        "label": "Static-24",
        "template": "static-24_mitc4_case1.bdf",
        "out_bdf": "static-24_{tag}_case1.bdf",
        "out_f06": "static-24_{tag}_case1.F06",
        "theory": 9.4e-02,
        "extract": ("node_abs_dof", 1, 1),
    },
    {
        "label": "Static-33 In-plane",
        "template": "static-33_mitc4_in_plane_shear.bdf",
        "out_bdf": "static-33_{tag}_in_plane_shear.bdf",
        "out_f06": "static-33_{tag}_in_plane_shear.F06",
        "theory": -5.424e-03,
        "extract": ("node_dof", 123, 3),
    },
    {
        "label": "Static-33 Out-of-plane",
        "template": "static-33_mitc4_out_of_plane_shear.bdf",
        "out_bdf": "static-33_{tag}_out_of_plane_shear.bdf",
        "out_f06": "static-33_{tag}_out_of_plane_shear.F06",
        "theory": -1.754e-03,
        "extract": ("node_dof", 123, 2),
    },
]


FLOAT_RE = r"(?:[+-]?\d+\.\d+(?:E[+-]\d+)?|[+-]?\d+\.(?:0)?)"
F06_ROW_RE = re.compile(
    rf"^\s*(\d+)\s+\d+\s+"
    rf"({FLOAT_RE})\s+({FLOAT_RE})\s+({FLOAT_RE})\s+"
    rf"({FLOAT_RE})\s+({FLOAT_RE})\s+({FLOAT_RE})\s*$"
)


def split_quad_to_tria(n1: int, n2: int, n3: int, n4: int) -> list[tuple[int, int, int]]:
    return [(n1, n2, n3), (n1, n3, n4)]


def convert_bdf(template_path: Path, out_path: Path, tria3typ: str, title_tag: str) -> None:
    lines = template_path.read_text().splitlines()
    out_lines: list[str] = []
    next_tri_eid = 1

    for line in lines:
        stripped = line.strip()
        if stripped.startswith("TITLE ="):
            out_lines.append(line.replace("MITC4", title_tag))
            continue
        if stripped.startswith("ID "):
            out_lines.append(line.replace("MITC4", title_tag))
            continue
        if stripped.startswith("PARAM,QUAD4TYP,"):
            continue
        if stripped.startswith("PARAM,K6ROT,"):
            out_lines.append(line)
            out_lines.append(f"PARAM,TRIA3TYP,{tria3typ}")
            continue
        if stripped.startswith("CQUAD4,"):
            fields = [part.strip() for part in line.split(",")]
            pid = int(fields[2])
            n1, n2, n3, n4 = map(int, fields[3:7])
            for tri_nodes in split_quad_to_tria(n1, n2, n3, n4):
                out_lines.append(
                    f"CTRIA3,{next_tri_eid},{pid},{tri_nodes[0]},{tri_nodes[1]},{tri_nodes[2]}"
                )
                next_tri_eid += 1
            continue
        out_lines.append(line)

    out_path.write_text("\n".join(out_lines) + "\n")


def parse_displacements(f06_path: Path) -> dict[int, tuple[float, float, float, float, float, float]]:
    data: dict[int, tuple[float, float, float, float, float, float]] = {}
    for line in f06_path.read_text(errors="ignore").splitlines():
        match = F06_ROW_RE.match(line)
        if not match:
            continue
        grid = int(match.group(1))
        vals = tuple(float(match.group(i)) for i in range(2, 8))
        data[grid] = vals
    if not data:
        raise RuntimeError(f"No displacement rows parsed from {f06_path}")
    return data


def extract_metric(
    disp: dict[int, tuple[float, float, float, float, float, float]],
    extract_spec: tuple,
) -> float:
    kind = extract_spec[0]
    if kind == "node_dof":
        _, node_id, dof = extract_spec
        return disp[node_id][dof - 1]
    if kind == "node_abs_dof":
        _, node_id, dof = extract_spec
        return abs(disp[node_id][dof - 1])
    if kind == "max_abs_dof":
        _, dof = extract_spec
        idx = dof - 1
        return max(abs(row[idx]) for row in disp.values())
    raise ValueError(f"Unknown extract spec: {extract_spec}")


def format_rows(rows: Iterable[dict[str, float | str]], family_label: str) -> str:
    lines = [
        f"| Load Case | Theory | MYSTRAN {family_label} | Error % |",
        "|---|---:|---:|---:|",
    ]
    for row in rows:
        lines.append(
            f"| `{row['label']}` | {row['theory']:.8e} | {row['result']:.8e} | {row['err']:.3f} |"
        )
    return "\n".join(lines)


def generate_all(tria3typ: str, tag: str, title_tag: str) -> None:
    for case in BENCHMARKS:
        template_path = RUN_DIR / case["template"]
        out_bdf_path = RUN_DIR / case["out_bdf"].format(tag=tag)
        convert_bdf(template_path, out_bdf_path, tria3typ=tria3typ, title_tag=title_tag)
        print(f"generated {out_bdf_path.name}")


def parse_all(tag: str) -> list[dict[str, float | str]]:
    rows: list[dict[str, float | str]] = []
    for case in BENCHMARKS:
        out_f06_path = RUN_DIR / case["out_f06"].format(tag=tag)
        disp = parse_displacements(out_f06_path)
        result = extract_metric(disp, case["extract"])
        theory = float(case["theory"])
        err = 100.0 * abs(result - theory) / abs(theory)
        rows.append(
            {
                "label": case["label"],
                "theory": theory,
                "result": result,
                "err": err,
            }
        )
        print(f"{case['label']}: {result:.8e} err={err:.3f}%")
    return rows


def write_summary(rows: list[dict[str, float | str]], label: str, out_md: Path) -> None:
    md = (
        f"# {label} MYSTRAN Shell Suite\n\n"
        "Generated from legacy `static-*mitc4*.bdf` shell decks by splitting each `CQUAD4` into two `CTRIA3`\n"
        f"and injecting `PARAM,TRIA3TYP,{label}`.\n\n"
        + format_rows(rows, label)
        + "\n"
    )
    out_md.write_text(md)
    print(f"\nWrote summary to {out_md}")


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--generate-only", action="store_true")
    parser.add_argument("--parse-only", action="store_true")
    parser.add_argument("--tria3typ", default="MITC3+")
    parser.add_argument("--tag", default="mitc3p")
    parser.add_argument("--title-tag", default="MITC3+")
    parser.add_argument("--out-md", default=str(OUT_MD))
    args = parser.parse_args(argv)
    out_md = Path(args.out_md)

    if args.generate_only:
        generate_all(args.tria3typ, args.tag, args.title_tag)
        return 0

    if args.parse_only:
        rows = parse_all(args.tag)
        write_summary(rows, args.tria3typ, out_md)
        return 0

    generate_all(args.tria3typ, args.tag, args.title_tag)
    rows = parse_all(args.tag)
    write_summary(rows, args.tria3typ, out_md)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
