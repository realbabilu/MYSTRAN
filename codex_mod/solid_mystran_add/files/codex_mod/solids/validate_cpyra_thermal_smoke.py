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


def build_case(kind: str) -> Model:
    model = Model(ndof_per_node=3)
    if kind == "cpyra5":
        add_pyra_grid(model, 1, 1, 1, 1.0e7, 0.3)
    elif kind == "cpyra14":
        add_liu_grid(model, 1, 1, 1, 1.0e7, 0.3)
    else:
        raise ValueError(kind)
    return model


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


def write_deck(kind: str, model: Model) -> Path:
    RUN_DIR.mkdir(parents=True, exist_ok=True)
    deck = RUN_DIR / f"{kind}_thermal_smoke.dat"
    lines = [
        "SOL 101",
        "CEND",
        f"TITLE = {kind.upper()} NEWSOLID THERMAL SMOKE",
        "SUBCASE 1",
        "  SPC = 1",
        "  TEMP = 100",
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
    lines.append("MAT1,1,1.0+7,,0.30,1.0,1.0-5,0.0")
    fixed = [nid for nid, node in model.nodes.items() if abs(node.x) < 1.0e-12]
    for nid in sorted(fixed):
        lines.append(card("SPC1", [1, 123, nid]))
    lines.append("TEMPD,100,25.0")
    lines.append("ENDDATA")
    deck.write_text("\n".join(lines) + "\n", encoding="ascii")
    return deck


def run_mystran(deck: Path) -> Path:
    result = subprocess.run([str(EXE), str(deck)], cwd=MYSTRAN, text=True, capture_output=True)
    if result.returncode != 0:
        print(result.stdout[-4000:])
        print(result.stderr[-4000:])
        result.check_returncode()
    return deck.with_suffix(".ERR")


def has_fatal(err: Path) -> bool:
    return "*FATAL" in err.read_text(errors="ignore")


def main() -> None:
    print("| case | run | fatal |")
    print("|---|---:|---:|")
    for kind in ("cpyra5", "cpyra14"):
        model = build_case(kind)
        err = run_mystran(write_deck(kind, model))
        print(f"| {kind} | ok | {str(has_fatal(err)).lower()} |")


if __name__ == "__main__":
    main()
