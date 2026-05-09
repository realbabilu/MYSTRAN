from __future__ import annotations

import math
import re
import subprocess
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
RUN_DIR = ROOT / "run_debug" / "solid_newsolid_compare"
EXE = ROOT / "Binaries" / "mystran.exe"
OUT_MD = ROOT / "codex_mod" / "solids" / "macneal_cpenta6_legacy_newsolid.md"

CASES = [("inplane_z", 3, -0.005424), ("outplane_y", 2, -0.001754)]


def bdf_float(value: float) -> str:
    text = f"{value:.8g}"
    if "." not in text and "e" not in text.lower():
        text += ".0"
    return text


def card(name: str, fields: list[object]) -> str:
    return ",".join([name, *[str(f) for f in fields]])


def macneal_point(xi: float, eta: float, zeta: float) -> tuple[float, float, float]:
    length, width, thickness = 12.0, 1.1, 0.32
    x = length * xi
    s = width * (eta - 0.5)
    q = thickness * (zeta - 0.5)
    theta = 0.5 * math.pi * xi
    y = s * math.cos(theta) - q * math.sin(theta)
    z = s * math.sin(theta) + q * math.cos(theta)
    return x, y, z


def run_mystran(deck: Path) -> Path:
    result = subprocess.run([str(EXE), str(deck)], cwd=ROOT, text=True, capture_output=True)
    if result.returncode != 0:
        print(result.stdout[-4000:])
        print(result.stderr[-4000:])
        result.check_returncode()
    return deck.with_suffix(".F06")


def write_deck(nx: int, case: str, newsolid: bool) -> tuple[Path, list[int]]:
    ny, nz = 2, 1
    node: dict[tuple[int, int, int], int] = {}
    coords: dict[int, tuple[float, float, float]] = {}
    nid = 1
    for i in range(nx + 1):
        for j in range(ny + 1):
            for k in range(nz + 1):
                node[(i, j, k)] = nid
                coords[nid] = macneal_point(i / nx, j / ny, k / nz)
                nid += 1

    elems: list[list[int]] = []
    for i in range(nx):
        for j in range(ny):
            for k in range(nz):
                c000 = node[(i, j, k)]
                c100 = node[(i + 1, j, k)]
                c110 = node[(i + 1, j + 1, k)]
                c010 = node[(i, j + 1, k)]
                c001 = node[(i, j, k + 1)]
                c101 = node[(i + 1, j, k + 1)]
                c111 = node[(i + 1, j + 1, k + 1)]
                c011 = node[(i, j + 1, k + 1)]
                elems.append([c000, c100, c110, c001, c101, c111])
                elems.append([c000, c110, c010, c001, c111, c011])

    tag = "newsolid" if newsolid else "legacy"
    deck = RUN_DIR / f"macneal_cpenta6_n{nx}_{case}_{tag}.dat"
    lines = [
        "SOL 101",
        "CEND",
        f"TITLE = MACNEAL STATIC37 CPENTA6 {tag.upper()} {case} NX{nx}",
        "SUBCASE 1",
        "  SPC = 1",
        "  LOAD = 1",
        "  DISPLACEMENT = ALL",
        "BEGIN BULK",
        "PARAM,POST,-1",
    ]
    if newsolid:
        lines.append("PARAM,SOLIDTYP,NEWSOLID")
    for gid in sorted(coords):
        x, y, z = coords[gid]
        lines.append(card("GRID", [gid, "", bdf_float(x), bdf_float(y), bdf_float(z)]))
    for eid, conn in enumerate(elems, start=1):
        lines.append(card("CPENTA", [eid, 1, *conn]))
    lines.append("PSOLID,1,1,,3,,FULL")
    lines.append("MAT1,1,2.9+7,,0.22,1.0")
    for j in range(ny + 1):
        for k in range(nz + 1):
            gid = node[(0, j, k)]
            lines.append(card("SPC", [1, gid, 13, "0."]))
            if j == ny // 2:
                lines.append(card("SPC", [1, gid, 2, "0."]))
    free_face = [node[(nx, j, k)] for j in range(ny + 1) for k in range(nz + 1)]
    direction = (0.0, 0.0, -1.0) if case == "inplane_z" else (0.0, -1.0, 0.0)
    load = 1.0 / len(free_face)
    for gid in free_face:
        lines.append(card("FORCE", [1, gid, "", f"{load:.8f}", *[bdf_float(v) for v in direction]]))
    lines.append("ENDDATA")
    deck.write_text("\n".join(lines) + "\n", encoding="ascii")
    return deck, free_face


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
        "# MacNeal Static-37 CPENTA6 Legacy vs NEWSOLID",
        "",
        "Same CPENTA6 wedge mesh; legacy means standard CPENTA6, NEWSOLID means CPENTA6_EAS9.",
        "",
        "| mesh | case | reference | legacy avg | legacy/ref | NEWSOLID avg | NEWSOLID/ref | NEWSOLID/legacy |",
        "|---|---|---:|---:|---:|---:|---:|---:|",
    ]
    for nx in (4, 8, 12):
        for case, comp, ref in CASES:
            legacy_deck, nodes = write_deck(nx, case, newsolid=False)
            new_deck, _ = write_deck(nx, case, newsolid=True)
            legacy_avg = parse_avg_disp(run_mystran(legacy_deck), nodes, comp)
            new_avg = parse_avg_disp(run_mystran(new_deck), nodes, comp)
            lines.append(
                f"| `{nx}x2x1 blocks -> 2 penta/block` | `{case}` | `{ref:.9e}` | `{legacy_avg:.9e}` | `{legacy_avg / ref:.3f}` | "
                f"`{new_avg:.9e}` | `{new_avg / ref:.3f}` | `{new_avg / legacy_avg:.3f}` |"
            )
    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(OUT_MD.read_text(encoding="utf-8"))


if __name__ == "__main__":
    main()
