from __future__ import annotations

import re
import subprocess
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[3]
SOLID_PY = ROOT / "codex_mod" / "solid-python"
MYSTRAN = ROOT / "MYSTRANSolver-18.0.0"
RUN_DIR = MYSTRAN / "run_debug" / "solid_newsolid_compare"
EXE = MYSTRAN / "Binaries" / "mystran.exe"

if str(SOLID_PY) not in sys.path:
    sys.path.insert(0, str(SOLID_PY))

from core import Model  # noqa: E402
from solid3d_ctetra4_validation import add_tet_grid  # noqa: E402
from solid3d_ctetra4_smooth import assemble_blended_smooth_stiffness  # noqa: E402


def bdf_float(value: float) -> str:
    text = f"{value:.8g}"
    if "." not in text and "e" not in text.lower():
        text += ".0"
    return text


def card(name: str, fields: list[object]) -> str:
    return ",".join([name, *[str(f) for f in fields]])


def build_case(nx: int, smooth: bool) -> tuple[Model, list[int]]:
    model = Model(ndof_per_node=3)
    e, nu = 1.0e7, 0.3
    node = add_tet_grid(model, nx, 2, 2, e, nu)
    for j in range(3):
        for k in range(3):
            model.fix_node(node[(0, j, k)], [0, 1, 2])
    face = [node[(nx, j, k)] for j in range(3) for k in range(3)]
    for nid in face:
        model.add_load(nid, 2, -1.0 / len(face))
    if smooth:
        model.build()
        node_index = {nid: i for i, nid in enumerate(model.nodes)}
        model.K = assemble_blended_smooth_stiffness(model.K, node_index, model.elements, e, nu, alpha=0.9)
    return model, sorted(face)


def write_deck(nx: int, model: Model, face: list[int]) -> Path:
    RUN_DIR.mkdir(parents=True, exist_ok=True)
    deck = RUN_DIR / f"ctetra4_n{nx}_legacy_target.dat"
    lines = [
        "SOL 101",
        "CEND",
        f"TITLE = CTETRA4 LEGACY TARGET NX{nx}",
        "SUBCASE 1",
        "  SPC = 1",
        "  LOAD = 1",
        "  DISPLACEMENT = ALL",
        "BEGIN BULK",
        "PARAM,POST,-1",
        "PARAM,SOLIDTYP,NEWSOLID",
    ]
    for nid in sorted(model.nodes):
        node = model.nodes[nid]
        lines.append(card("GRID", [nid, "", bdf_float(node.x), bdf_float(node.y), bdf_float(node.z)]))
    for elem in model.elements:
        ids = [node.nid for node in elem.nodes]
        lines.append(card("CTETRA", [elem.eid, 1, *ids]))
    lines.append("PSOLID,1,1")
    lines.append("MAT1,1,1.0+7,,0.30,1.0")
    fixed = [nid for nid, node in model.nodes.items() if abs(node.x) < 1.0e-12]
    for nid in sorted(fixed):
        lines.append(card("SPC1", [1, 123, nid]))
    load = -1000.0 / len(face)
    for nid in face:
        lines.append(card("FORCE", [1, nid, "", f"{load:.8f}", "0.", "0.", "1."]))
    lines.append("ENDDATA")
    deck.write_text("\n".join(lines) + "\n", encoding="ascii")
    return deck


def run_mystran(deck: Path) -> Path:
    result = subprocess.run([str(EXE), str(deck)], cwd=MYSTRAN, text=True, capture_output=True)
    if result.returncode != 0:
        print(result.stdout[-4000:])
        print(result.stderr[-4000:])
        result.check_returncode()
    return deck.with_suffix(".F06")


def parse_avg_t3(f06: Path, face: list[int]) -> float:
    values: dict[int, float] = {}
    pat = re.compile(r"^\s*(\d+)\s+0\s+([-+0-9.E]+)\s+([-+0-9.E]+)\s+([-+0-9.E]+)")
    for line in f06.read_text(errors="ignore").splitlines():
        match = pat.match(line)
        if match:
            nid = int(match.group(1))
            if nid in face:
                values[nid] = float(match.group(4))
    missing = sorted(set(face) - set(values))
    if missing:
        raise RuntimeError(f"{f06.name}: missing displacement rows for {missing[:8]}")
    return sum(values[nid] for nid in face) / len(face)


def python_tip(model: Model, face: list[int]) -> float:
    u = model.solve_static()
    return sum(float(u[model.nodes[nid].dofs[2]]) for nid in face) / len(face)


def main() -> None:
    print("| mesh | Python standard T3 x1000 | Python smooth a09 T3 x1000 | MYSTRAN NEWSOLID T3 | MYSTRAN/standard | smooth/standard |")
    print("|---|---:|---:|---:|---:|---:|")
    for nx in (4, 8, 12):
        standard, face = build_case(nx, smooth=False)
        smooth, _ = build_case(nx, smooth=True)
        py_standard = python_tip(standard, face) * 1000.0
        py_smooth = python_tip(smooth, face) * 1000.0
        deck = write_deck(nx, standard, face)
        mystran = parse_avg_t3(run_mystran(deck), face)
        print(
            f"| n{nx} | {py_standard:.9e} | {py_smooth:.9e} | {mystran:.9e} | "
            f"{mystran / py_standard:.6f} | {py_smooth / py_standard:.6f} |"
        )


if __name__ == "__main__":
    main()
