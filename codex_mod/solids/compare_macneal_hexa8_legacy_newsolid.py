from __future__ import annotations

import re
import subprocess
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
RUN_DIR = ROOT / "run_debug" / "solid_newsolid_compare"
EXE = ROOT / "Binaries" / "mystran.exe"
OUT_MD = ROOT / "codex_mod" / "solids" / "macneal_hexa8_legacy_newsolid.md"

CASES = [
    ("inplane_z", 3, -0.005424),
    ("outplane_y", 2, -0.001754),
]


def run_mystran(deck: Path) -> Path:
    result = subprocess.run([str(EXE), str(deck)], cwd=ROOT, text=True, capture_output=True)
    if result.returncode != 0:
        print(result.stdout[-4000:])
        print(result.stderr[-4000:])
        result.check_returncode()
    return deck.with_suffix(".F06")


def make_legacy_deck(newsolid_deck: Path) -> Path:
    legacy = newsolid_deck.with_name(newsolid_deck.name.replace("_newsolid.dat", "_legacy.dat"))
    lines = []
    for line in newsolid_deck.read_text(encoding="ascii").splitlines():
        if line.strip().upper() == "PARAM,SOLIDTYP,NEWSOLID":
            continue
        lines.append(line.replace("NEWSOLID", "LEGACY"))
    legacy.write_text("\n".join(lines) + "\n", encoding="ascii")
    return legacy


def force_nodes(deck: Path) -> list[int]:
    nodes: list[int] = []
    for line in deck.read_text(encoding="ascii").splitlines():
        if line.upper().startswith("FORCE,"):
            fields = [x.strip() for x in line.split(",")]
            nodes.append(int(fields[2]))
    return sorted(set(nodes))


def parse_avg_disp(f06: Path, nodes: list[int], comp: int) -> float:
    values: dict[int, tuple[float, float, float]] = {}
    pat = re.compile(r"^\s*(\d+)\s+0\s+([-+0-9.E]+)\s+([-+0-9.E]+)\s+([-+0-9.E]+)")
    for line in f06.read_text(errors="ignore").splitlines():
        m = pat.match(line)
        if m:
            nid = int(m.group(1))
            if nid in nodes:
                values[nid] = (float(m.group(2)), float(m.group(3)), float(m.group(4)))
    missing = sorted(set(nodes) - set(values))
    if missing:
        raise RuntimeError(f"{f06.name}: missing displacement rows for {missing}")
    return sum(values[nid][comp - 1] for nid in nodes) / len(nodes)


def main() -> None:
    lines = [
        "# MacNeal Static-37 CHEXA8 Legacy vs NEWSOLID",
        "",
        "Same generated CHEXA8 decks, compared with and without `PARAM,SOLIDTYP,NEWSOLID`.",
        "",
        "| mesh | case | reference | legacy avg | legacy/ref | NEWSOLID avg | NEWSOLID/ref | NEWSOLID/legacy |",
        "|---|---|---:|---:|---:|---:|---:|---:|",
    ]
    for nx in (4, 8, 12):
        for case, comp, ref in CASES:
            newsolid = RUN_DIR / f"macneal_eas9_n{nx}_{case}_newsolid.dat"
            legacy = make_legacy_deck(newsolid)
            nodes = force_nodes(newsolid)
            legacy_avg = parse_avg_disp(run_mystran(legacy), nodes, comp)
            new_avg = parse_avg_disp(run_mystran(newsolid), nodes, comp)
            lines.append(
                f"| `{nx}x2x1` | `{case}` | `{ref:.9e}` | `{legacy_avg:.9e}` | `{legacy_avg / ref:.3f}` | "
                f"`{new_avg:.9e}` | `{new_avg / ref:.3f}` | `{new_avg / legacy_avg:.3f}` |"
            )
    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(OUT_MD.read_text(encoding="utf-8"))


if __name__ == "__main__":
    main()
