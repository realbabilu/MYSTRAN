from __future__ import annotations

import math
import re
import subprocess
import sys
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[2]
WORKSPACE = ROOT.parents[0]
SOLID_PY = WORKSPACE / "codex_mod" / "solid-python"
RUN_DIR = ROOT / "run_debug" / "solid_newsolid_compare"
EXE = ROOT / "Binaries" / "mystran.exe"
OUT_MD = ROOT / "codex_mod" / "solids" / "macneal_python_vs_mystran.md"

if str(SOLID_PY) not in sys.path:
    sys.path.insert(0, str(SOLID_PY))

from core import Model  # noqa: E402
from solid3d_chexa8_eas import CHEXA8EAS  # noqa: E402
from solid3d_chexa8_eas9_frozen import CHEXA8EAS9Frozen  # noqa: E402
from solid3d_cpenta6 import CPENTA6  # noqa: E402
from solid3d_cpenta6_eas import CPENTA6EAS  # noqa: E402
from solid3d_cpyra5 import CPYRA5  # noqa: E402
from solid3d_cpyra5_eas import CPYRA5EAS54  # noqa: E402
from solid3d_ctetra4 import CTETRA4  # noqa: E402
from solid3d_ctetra4_smooth import assemble_blended_smooth_stiffness  # noqa: E402
from solid3d_ctetra10 import CTETRA10  # noqa: E402


CASES = [("inplane_z", 2, -0.005424), ("outplane_y", 1, -0.001754)]
E = 2.9e7
NU = 0.22


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


def orient_tet(conn: list[int], coords: dict[int, tuple[float, float, float]]) -> list[int]:
    a, b, c, d = conn
    ax, ay, az = coords[a]
    bx, by, bz = coords[b]
    cx, cy, cz = coords[c]
    dx, dy, dz = coords[d]
    m = ((bx - ax, by - ay, bz - az), (cx - ax, cy - ay, cz - az), (dx - ax, dy - ay, dz - az))
    det = (
        m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
        - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
        + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0])
    )
    return [a, c, b, d] if det < 0.0 else conn


def add_mid(a: int, b: int, coords: dict[int, tuple[float, float, float]], mids: dict[tuple[int, int], int], next_nid: list[int]) -> int:
    key = tuple(sorted((a, b)))
    if key not in mids:
        ax, ay, az = coords[a]
        bx, by, bz = coords[b]
        nid = next_nid[0]
        next_nid[0] += 1
        coords[nid] = ((ax + bx) * 0.5, (ay + by) * 0.5, (az + bz) * 0.5)
        mids[key] = nid
    return mids[key]


def base_grid(nx: int, ny: int = 2, nz: int = 1):
    node: dict[tuple[int, int, int], int] = {}
    coords: dict[int, tuple[float, float, float]] = {}
    nid = 1
    for i in range(nx + 1):
        for j in range(ny + 1):
            for k in range(nz + 1):
                node[(i, j, k)] = nid
                coords[nid] = macneal_point(i / nx, j / ny, k / nz)
                nid += 1
    return node, coords, nid


def apply_bcs_loads(model: Model, node: dict, nx: int, case: str, ny: int = 2, nz: int = 1) -> list[int]:
    for j in range(ny + 1):
        for k in range(nz + 1):
            nid = node[(0, j, k)]
            model.fix_node(nid, [0, 2])
            if j == ny // 2:
                model.fix_node(nid, [1])
    face = [node[(nx, j, k)] for j in range(ny + 1) for k in range(nz + 1)]
    dof = 2 if case == "inplane_z" else 1
    for nid in face:
        model.add_load(nid, dof, -1.0 / len(face))
    return face


def solve_python(model: Model, face: list[int], dof: int) -> float:
    u = model.solve_static()
    return sum(float(u[model.nodes[nid].dofs[dof]]) for nid in face) / len(face)


def parse_avg_disp(f06: Path, nodes: list[int], comp1: int) -> float:
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
    return sum(values[nid][comp1 - 1] for nid in nodes) / len(nodes)


def add_nodes(model: Model, coords: dict[int, tuple[float, float, float]]) -> None:
    for nid in sorted(coords):
        model.add_node(nid, *coords[nid])


def hexa_conn(node: dict, nx: int, ny: int = 2, nz: int = 1):
    for i in range(nx):
        for j in range(ny):
            for k in range(nz):
                yield [
                    node[(i, j, k)],
                    node[(i + 1, j, k)],
                    node[(i + 1, j + 1, k)],
                    node[(i, j + 1, k)],
                    node[(i, j, k + 1)],
                    node[(i + 1, j, k + 1)],
                    node[(i + 1, j + 1, k + 1)],
                    node[(i, j + 1, k + 1)],
                ]


def penta_conn(node: dict, nx: int, ny: int = 2, nz: int = 1):
    for h in hexa_conn(node, nx, ny, nz):
        c000, c100, c110, c010, c001, c101, c111, c011 = h
        yield [c000, c100, c110, c001, c101, c111]
        yield [c000, c110, c010, c001, c111, c011]


def pyra_conn(node: dict, coords: dict[int, tuple[float, float, float]], next_nid: list[int], nx: int, ny: int = 2, nz: int = 1):
    for h in hexa_conn(node, nx, ny, nz):
        center = next_nid[0]
        next_nid[0] += 1
        coords[center] = tuple(sum(coords[n][a] for n in h) / 8.0 for a in range(3))
        c000, c100, c110, c010, c001, c101, c111, c011 = h
        yield [c000, c010, c011, c001, center]
        yield [c100, c101, c111, c110, center]
        yield [c000, c001, c101, c100, center]
        yield [c010, c110, c111, c011, center]
        yield [c000, c100, c110, c010, center]
        yield [c001, c011, c111, c101, center]


def tet4_conn(node: dict, coords: dict[int, tuple[float, float, float]], nx: int, ny: int = 2, nz: int = 1):
    for h in hexa_conn(node, nx, ny, nz):
        c000, c100, c110, c010, c001, c101, c111, c011 = h
        raw = [
            [c000, c100, c010, c001],
            [c100, c110, c010, c111],
            [c100, c010, c001, c111],
            [c100, c101, c001, c111],
            [c010, c001, c011, c111],
            [c100, c010, c111, c001],
        ]
        for conn in raw:
            yield orient_tet(conn, coords)


def tet10_conn(tets4: list[list[int]], coords: dict[int, tuple[float, float, float]], next_nid: list[int]):
    mids: dict[tuple[int, int], int] = {}
    for a, b, c, d in tets4:
        yield [
            a,
            b,
            c,
            d,
            add_mid(a, b, coords, mids, next_nid),
            add_mid(b, c, coords, mids, next_nid),
            add_mid(a, c, coords, mids, next_nid),
            add_mid(a, d, coords, mids, next_nid),
            add_mid(b, d, coords, mids, next_nid),
            add_mid(c, d, coords, mids, next_nid),
        ]


def python_model(kind: str, branch: str, nx: int, case: str):
    node, coords, nid = base_grid(nx)
    next_nid = [nid]
    conns: list[list[int]]
    cls: type
    if kind == "chexa8":
        conns = list(hexa_conn(node, nx))
        cls = CHEXA8EAS9Frozen if branch == "newsolid" else (lambda eid, nodes, E, nu: CHEXA8EAS(eid, nodes, E, nu, use_eas=False))  # type: ignore
    elif kind == "cpenta6":
        conns = list(penta_conn(node, nx))
        cls = CPENTA6EAS if branch == "newsolid" else CPENTA6
    elif kind == "cpyra5":
        conns = list(pyra_conn(node, coords, next_nid, nx))
        cls = CPYRA5EAS54 if branch == "newsolid" else CPYRA5
    elif kind == "ctetra4":
        conns = list(tet4_conn(node, coords, nx))
        cls = CTETRA4
    elif kind == "ctetra10":
        tets = list(tet4_conn(node, coords, nx))
        conns = list(tet10_conn(tets, coords, next_nid))
        cls = CTETRA10
    else:
        raise ValueError(kind)
    model = Model(ndof_per_node=3)
    add_nodes(model, coords)
    for eid, conn in enumerate(conns, start=1):
        model.add_element(cls(eid, [model.nodes[n] for n in conn], E, NU))
    face = apply_bcs_loads(model, node, nx, case)
    if kind == "ctetra4" and branch == "newsolid":
        model.build()
        node_index = {nid2: i for i, nid2 in enumerate(model.nodes)}
        model.K = assemble_blended_smooth_stiffness(model.K, node_index, model.elements, E, NU, alpha=0.9)
    return model, face


def mystran_deck(kind: str, branch: str, nx: int, case: str) -> Path:
    if kind == "chexa8":
        return RUN_DIR / f"macneal_eas9_n{nx}_{case}_{branch}.dat"
    if kind == "cpenta6":
        return RUN_DIR / f"macneal_cpenta6_n{nx}_{case}_{branch}.dat"
    if kind == "cpyra5":
        return RUN_DIR / f"macneal_cpyra5_n{nx}_{case}_{branch}.dat"
    if kind in {"ctetra4", "ctetra10"}:
        return RUN_DIR / f"macneal_{kind}_n{nx}_{case}_{branch}.dat"
    raise ValueError(kind)


def main() -> None:
    lines = [
        "# MacNeal Python vs MYSTRAN",
        "",
        "Same generated MacNeal Static-37 meshes/decks. Python values are computed from the local `codex_mod/solid-python` element classes.",
        "",
        "| type | branch | mesh | case | Python avg | MYSTRAN avg | delta | ref ratio Python | ref ratio MYSTRAN |",
        "|---|---|---|---|---:|---:|---:|---:|---:|",
    ]
    specs = [
        ("chexa8", ["legacy", "newsolid"], [4, 8, 12]),
        ("cpenta6", ["legacy", "newsolid"], [4, 8, 12]),
        ("cpyra5", ["legacy", "newsolid"], [4, 8, 12]),
        ("ctetra4", ["legacy", "newsolid"], [4, 8, 12]),
        ("ctetra10", ["legacy", "newsolid"], [4, 8]),
    ]
    for kind, branches, meshes in specs:
        for branch in branches:
            for nx in meshes:
                for case, dof, ref in CASES:
                    model, face = python_model(kind, branch, nx, case)
                    py = solve_python(model, face, dof)
                    comp1 = dof + 1
                    my = parse_avg_disp(run_mystran(mystran_deck(kind, branch, nx, case)), face, comp1)
                    lines.append(
                        f"| `{kind}` | `{branch}` | `{nx}x2x1` | `{case}` | `{py:.9e}` | `{my:.9e}` | `{my - py:+.3e}` | "
                        f"`{py / ref:.3f}` | `{my / ref:.3f}` |"
                    )
    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(OUT_MD.read_text(encoding="utf-8"))


if __name__ == "__main__":
    main()
