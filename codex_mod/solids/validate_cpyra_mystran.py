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
from solid3d_cpyra5_eas import CPYRA5EAS54  # noqa: E402
from solid3d_cpyra5_validation import add_pyra_grid  # noqa: E402
from solid3d_cpyra14_liu_validation import add_liu_grid  # noqa: E402


def bdf_float(value: float) -> str:
    text = f"{value:.8g}"
    if "." not in text and "e" not in text.lower():
        text += ".0"
    return text


def card(name: str, fields: list[object]) -> str:
    return ",".join([name, *[str(f) for f in fields]])


def cont(fields: list[object]) -> str:
    return ",".join(["", *[str(f) for f in fields]])


def build_case(kind: str, nx: int) -> tuple[Model, list[int]]:
    model = Model(ndof_per_node=3)
    e, nu = 1.0e7, 0.3
    if kind == "cpyra5":
        add_pyra_grid(model, nx, 2, 2, e, nu)
        model.elements = [CPYRA5EAS54(elem.eid, elem.nodes, e, nu) for elem in model.elements]
        model._built = False
    elif kind == "cpyra14":
        add_liu_grid(model, nx, 2, 2, e, nu)
    else:
        raise ValueError(kind)

    face = []
    for nid, node in model.nodes.items():
        if abs(node.x) < 1.0e-12:
            model.fix_node(nid, [0, 1, 2])
        if abs(node.x - 10.0) < 1.0e-12:
            face.append(nid)
    for nid in face:
        model.add_load(nid, 2, -1.0 / len(face))
    return model, sorted(face)


def element_cards(kind: str, elem) -> list[str]:
    ids = [node.nid for node in elem.nodes]
    if kind == "cpyra5":
        return [card("CPYRA", [elem.eid, 1, *ids])]
    if kind == "cpyra14":
        return [
            card("CPYRA", [elem.eid, 1, *ids[:5]]),
            cont(ids[5:13]),
            cont(ids[13:14]),
        ]
    raise ValueError(kind)


def write_deck(kind: str, nx: int, model: Model, face: list[int]) -> Path:
    RUN_DIR.mkdir(parents=True, exist_ok=True)
    deck = RUN_DIR / f"{kind}_n{nx}_newsolid.dat"
    lines = [
        "SOL 101",
        "CEND",
        f"TITLE = {kind.upper()} NEWSOLID NX{nx}",
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
        lines.extend(element_cards(kind, elem))
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
    print("| case | Python T3 x1000 | MYSTRAN T3 | delta |")
    print("|---|---:|---:|---:|")
    for kind, nx in (("cpyra5", 4), ("cpyra14", 2)):
        model, face = build_case(kind, nx)
        py = python_tip(model, face) * 1000.0
        deck = write_deck(kind, nx, model, face)
        f06 = run_mystran(deck)
        my = parse_avg_t3(f06, face)
        print(f"| {kind} n{nx} | {py:.9e} | {my:.9e} | {my - py:+.3e} |")


if __name__ == "__main__":
    main()
