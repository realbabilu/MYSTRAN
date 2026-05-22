from __future__ import annotations

import math
import re
import subprocess
from pathlib import Path


WORKDIR = Path(r"D:\user\doc\New project")
OUTDIR = WORKDIR / "benchmark_1_020_modal"
MYSTRAN_EXE = Path(r"D:\mystran2\MYSTRANSolver-18.0.0\Binaries\mystran.exe")

E = 3000.0
NU = 0.2
G = E / (2.0 * (1.0 + NU))
A_BIG = 100000.0
J_TORSION = 1.0
I_BIG = 2000.0
I_SMALL = 1000.0

JOINTS = {
    1: (-120.0, 0.0, 0.0),
    2: (-120.0, 0.0, 120.0),
    3: (-120.0, 0.0, 240.0),
    4: (120.0, 0.0, 0.0),
    5: (120.0, 0.0, 120.0),
    6: (120.0, 0.0, 240.0),
    7: (0.0, 0.0, 120.0),
    8: (0.0, 0.0, 240.0),
}

ELEMENTS = {
    1: (1, 2, "BIG"),
    2: (2, 3, "SMALL"),
    3: (4, 5, "BIG"),
    4: (5, 6, "SMALL"),
    5: (2, 7, "BIG"),
    6: (7, 5, "BIG"),
    7: (3, 8, "SMALL"),
    8: (8, 6, "SMALL"),
}

MASSES = {
    7: 1.0364,  # 2m
    8: 0.5182,  # m
}

EIG_ROW_RE = re.compile(
    r"^\s*(\d+)\s+(\d+)\s+([\-+0-9.EedD]+)\s+([\-+0-9.EedD]+)\s+([\-+0-9.EedD]+)\s+([\-+0-9.EedD]+)\s+([\-+0-9.EedD]+)\s*$"
)


def fmt(value: float) -> str:
    if abs(value) < 1.0e-12:
        return "0."
    text = f"{value:.12g}"
    if "e" in text:
        mantissa, exp = text.split("e", 1)
        if "." not in mantissa:
            mantissa += ".0"
        return f"{mantissa}E{exp}"
    if "." not in text:
        return text + "."
    return text


def line(*fields: object) -> str:
    return ",".join(str(field) for field in fields)


def pbeam_block(pid: int, inertia: float) -> list[str]:
    return [
        line("PBEAM", pid, 1, fmt(A_BIG), fmt(inertia), fmt(inertia), fmt(J_TORSION), fmt(0.0), fmt(0.0)),
        line("", fmt(0.0), fmt(0.0), fmt(0.0), fmt(0.0), fmt(0.0), fmt(0.0), fmt(0.0), fmt(0.0)),
        line("", "YES", fmt(1.0), fmt(A_BIG), fmt(inertia), fmt(inertia), fmt(J_TORSION), fmt(0.0), fmt(0.0)),
        line("", fmt(0.0), fmt(0.0), fmt(0.0), fmt(0.0), fmt(0.0), fmt(0.0), fmt(0.0), fmt(0.0)),
        line("", fmt(0.0), fmt(0.0)),
    ]


def build_dat() -> str:
    lines = [
        "ID,PROBLEM_1_020_MODAL",
        "SOL 103",
        "CEND",
        "TITLE = Problem 1-020 modal validation",
        "ECHO = NONE",
        "DISPLACEMENT(PRINT) = ALL",
        "MEFFMASS = ALL",
        "MPFACTOR = ALL",
        "OEF1 = ALL",
        "SUBCASE 1",
        "  METHOD = 1",
        "  SPC = 1",
        "BEGIN BULK",
        "PARAM,POST,-1",
        "EIGRL,1,,,2",
        ",DENSE,,64",
        line("MAT1", 1, fmt(E), fmt(G), fmt(NU), fmt(0.0)),
    ]
    lines.extend(pbeam_block(100, I_BIG))
    lines.extend(pbeam_block(200, I_SMALL))
    for jid, (x, y, z) in JOINTS.items():
        lines.append(line("GRID", jid, "", fmt(x), fmt(y), fmt(z)))
    for eid, (ga, gb, family) in ELEMENTS.items():
        pid = 100 if family == "BIG" else 200
        # Keep a consistent local frame; with in-plane DOFs only, equal Iyy/Izz
        # removes axis-mapping sensitivity.
        orient = (1.0, 0.0, 0.0) if JOINTS[ga][0] == JOINTS[gb][0] else (0.0, 1.0, 0.0)
        lines.append(line("CBEAM", eid, pid, ga, gb, fmt(orient[0]), fmt(orient[1]), fmt(orient[2])))
    lines.append(line("SPC1", 1, 123456, 1, 4))
    # Restrict inactive SAP DOFs everywhere else: U2, R1, R3
    lines.append(line("SPC1", 1, 246, 2, 3, 5, 6, 7, 8))
    for eid, (jid, mass) in enumerate(MASSES.items(), start=1):
        lines.append(line("CMASS2", 800 + eid, fmt(mass), jid, 1, 0, 0))
    lines.append("ENDDATA")
    return "\n".join(lines) + "\n"


def parse_cycles(f06_path: Path) -> list[float]:
    cycles: list[float] = []
    for raw in f06_path.read_text(encoding="utf-8", errors="ignore").splitlines():
        match = EIG_ROW_RE.match(raw)
        if match:
            cycles.append(float(match.group(5).replace("D", "E").replace("d", "e")))
    if len(cycles) < 2:
        raise RuntimeError(f"Could not extract two modal cycles from {f06_path}")
    return cycles[:2]


def main() -> None:
    OUTDIR.mkdir(exist_ok=True)
    dat_path = OUTDIR / "problem_1_020_modal.dat"
    dat_path.write_text(build_dat(), encoding="utf-8")
    subprocess.run([str(MYSTRAN_EXE), str(dat_path)], cwd=str(OUTDIR), check=False)
    cycles = parse_cycles(dat_path.with_suffix(".F06"))
    periods = [1.0 / c for c in cycles]
    for i, (cycle, period) in enumerate(zip(cycles, periods), start=1):
        print(f"mode {i}: f = {cycle:.6f} Hz, T = {period:.6f} s")
    print(f"dat: {dat_path}")
    print(f"f06: {dat_path.with_suffix('.F06')}")


if __name__ == "__main__":
    main()
