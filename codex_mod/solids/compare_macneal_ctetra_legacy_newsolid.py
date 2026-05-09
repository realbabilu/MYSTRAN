from __future__ import annotations

import math
import re
import subprocess
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
RUN_DIR = ROOT / "run_debug" / "solid_newsolid_compare"
EXE = ROOT / "Binaries" / "mystran.exe"
OUT_MD = ROOT / "codex_mod" / "solids" / "macneal_ctetra_legacy_newsolid.md"

CASES = [("inplane_z", 3, -0.005424), ("outplane_y", 2, -0.001754)]


def bdf_float(value: float) -> str:
    text = f"{value:.8g}"
    if "." not in text and "e" not in text.lower():
        text += ".0"
    return text


def card(name: str, fields: list[object]) -> str:
    return ",".join([name, *[str(f) for f in fields]])


def cont(fields: list[object]) -> str:
    return ",".join(["", *[str(f) for f in fields]])


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


def orient_tet(conn: list[int], coords: dict[int, tuple[float, float, float]]) -> list[int]:
    a, b, c, d = conn
    ax, ay, az = coords[a]
    bx, by, bz = coords[b]
    cx, cy, cz = coords[c]
    dx, dy, dz = coords[d]
    m = (
        (bx - ax, by - ay, bz - az),
        (cx - ax, cy - ay, cz - az),
        (dx - ax, dy - ay, dz - az),
    )
    det = (
        m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
        - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
        + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0])
    )
    if det < 0.0:
        return [a, c, b, d]
    return conn


def add_mid(
    a: int,
    b: int,
    coords: dict[int, tuple[float, float, float]],
    edge_mid: dict[tuple[int, int], int],
    next_nid: list[int],
) -> int:
    key = tuple(sorted((a, b)))
    if key not in edge_mid:
        ax, ay, az = coords[a]
        bx, by, bz = coords[b]
        nid = next_nid[0]
        next_nid[0] += 1
        coords[nid] = ((ax + bx) * 0.5, (ay + by) * 0.5, (az + bz) * 0.5)
        edge_mid[key] = nid
    return edge_mid[key]


def write_deck(kind: str, nx: int, case: str, newsolid: bool) -> tuple[Path, list[int]]:
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

    tets4: list[list[int]] = []
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
                tets4.extend(
                    [
                        [c000, c100, c010, c001],
                        [c100, c110, c010, c111],
                        [c100, c010, c001, c111],
                        [c100, c101, c001, c111],
                        [c010, c001, c011, c111],
                        [c100, c010, c111, c001],
                    ]
                )
    tets4 = [orient_tet(conn, coords) for conn in tets4]

    elems: list[list[int]] = []
    if kind == "ctetra4":
        elems = tets4
    elif kind == "ctetra10":
        edge_mid: dict[tuple[int, int], int] = {}
        next_nid = [nid]
        for a, b, c, d in tets4:
            elems.append(
                [
                    a,
                    b,
                    c,
                    d,
                    add_mid(a, b, coords, edge_mid, next_nid),
                    add_mid(b, c, coords, edge_mid, next_nid),
                    add_mid(a, c, coords, edge_mid, next_nid),
                    add_mid(a, d, coords, edge_mid, next_nid),
                    add_mid(b, d, coords, edge_mid, next_nid),
                    add_mid(c, d, coords, edge_mid, next_nid),
                ]
            )
    else:
        raise ValueError(kind)

    tag = "newsolid" if newsolid else "legacy"
    deck = RUN_DIR / f"macneal_{kind}_n{nx}_{case}_{tag}.dat"
    lines = [
        "SOL 101",
        "CEND",
        f"TITLE = MACNEAL STATIC37 {kind.upper()} {tag.upper()} {case} NX{nx}",
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
        if kind == "ctetra4":
            lines.append(card("CTETRA", [eid, 1, *conn]))
        else:
            label = f"+T{eid}"
            lines.append(card("CTETRA", [eid, 1, *conn[:6], label]))
            lines.append(",".join([label, *[str(x) for x in conn[6:10]]]))
    lines.append("PSOLID,1,1")
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
        "# MacNeal Static-37 CTETRA Legacy vs NEWSOLID",
        "",
        "Same tetra meshes; CTETRA4 NEWSOLID is smooth alpha 0.9, CTETRA10 NEWSOLID is the baseline quadratic guard.",
        "",
        "| type | mesh | case | reference | legacy avg | legacy/ref | NEWSOLID avg | NEWSOLID/ref | NEWSOLID/legacy |",
        "|---|---|---|---:|---:|---:|---:|---:|---:|",
    ]
    for kind in ("ctetra4", "ctetra10"):
        for nx in (4, 8, 12):
            if kind == "ctetra10" and nx > 8:
                continue
            for case, comp, ref in CASES:
                legacy_deck, nodes = write_deck(kind, nx, case, newsolid=False)
                new_deck, _ = write_deck(kind, nx, case, newsolid=True)
                legacy_avg = parse_avg_disp(run_mystran(legacy_deck), nodes, comp)
                new_avg = parse_avg_disp(run_mystran(new_deck), nodes, comp)
                lines.append(
                    f"| `{kind.upper()}` | `{nx}x2x1 blocks -> 6 tet/block` | `{case}` | `{ref:.9e}` | "
                    f"`{legacy_avg:.9e}` | `{legacy_avg / ref:.3f}` | `{new_avg:.9e}` | `{new_avg / ref:.3f}` | `{new_avg / legacy_avg:.3f}` |"
                )
    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(OUT_MD.read_text(encoding="utf-8"))


if __name__ == "__main__":
    main()
