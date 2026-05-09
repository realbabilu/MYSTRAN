from __future__ import annotations

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
        add_pyra_grid(model, nx, 1, 1, e, nu, rho=1.0)
    elif kind == "cpyra14":
        add_liu_grid(model, nx, 1, 1, e, nu, rho=1.0)
    else:
        raise ValueError(kind)
    right = [nid for nid, node in model.nodes.items() if abs(node.x - 10.0) < 1.0e-12]
    return model, sorted(right)


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


def write_deck(kind: str, nx: int, model: Model, right: list[int]) -> Path:
    RUN_DIR.mkdir(parents=True, exist_ok=True)
    deck = RUN_DIR / f"{kind}_n{nx}_buckling_smoke.dat"
    lines = [
        "ID CPYRA BUCKLING",
        "SOL 5",
        "CEND",
        f"TITLE = {kind.upper()} NEWSOLID BUCKLING SMOKE",
        "DISP = ALL",
        "SUBCASE 1",
        "  LOAD = 10",
        "SUBCASE 2",
        "  METHOD = 20",
        "BEGIN BULK",
        "PARAM,POST,-1",
        "PARAM,SOLIDTYP,NEWSOLID",
        "EIGRL   20                      2                       10.     MAX     +E1",
        "+E1             DGB",
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
    load = -1000.0 / len(right)
    for nid in right:
        lines.append(card("FORCE", [10, nid, "", f"{load:.8f}", "1.", "0.", "0."]))
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


def has_eigen_output(f06: Path) -> bool:
    text = f06.read_text(errors="ignore")
    return "R E A L   E I G E N V A L U E S" in text or "E I G E N V A L U E" in text


def main() -> None:
    print("| case | run | eigen table |")
    print("|---|---:|---:|")
    for kind, nx in (("cpyra5", 2), ("cpyra14", 1)):
        model, right = build_case(kind, nx)
        deck = write_deck(kind, nx, model, right)
        f06 = run_mystran(deck)
        print(f"| {kind} n{nx} | ok | {str(has_eigen_output(f06)).lower()} |")


if __name__ == "__main__":
    main()
