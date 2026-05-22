from __future__ import annotations

import math
import re
from pathlib import Path

import benchmark_1_020_modal as modal


F06_PATH = modal.OUTDIR / "problem_1_020_modal.F06"

SPECTRUM = [
    (0.03, 0.5),
    (0.125, 1.355),
    (0.5868, 1.355),
    (0.66, 1.355),
    (1.562, 0.576),
    (4.0, 0.2185),
    (10.0, 0.037),
]
G_ACCEL = 386.098

REFERENCE_PERIODS = [1.5620, 0.5868]
REFERENCE_UX = {2: 7.576, 3: 18.84}
REFERENCE_M = {
    (1, 1): 12636.0,
    (1, 2): 6793.0,
    (2, 2): 6023.0,
    (2, 3): 5222.0,
    (5, 2): 9810.0,
    (6, 5): 9810.0,
    (7, 3): 5222.0,
    (8, 6): 5222.0,
}

EIG_ROW_RE = re.compile(
    r"^\s*(\d+)\s+(\d+)\s+([\-+0-9.EedD]+)\s+([\-+0-9.EedD]+)\s+([\-+0-9.EedD]+)\s+([\-+0-9.EedD]+)\s+([\-+0-9.EedD]+)\s*$"
)
GRID_ROW_RE = re.compile(
    r"^\s*(\d+)\s+(\d+)\s+([\-+0-9.EedD]+|0\.0)\s+([\-+0-9.EedD]+|0\.0)\s+([\-+0-9.EedD]+|0\.0)\s+([\-+0-9.EedD]+|0\.0)\s+([\-+0-9.EedD]+|0\.0)\s+([\-+0-9.EedD]+|0\.0)\s*$"
)
def fnum(text: str) -> float:
    return float(text.replace("D", "E").replace("d", "e"))


def spectrum_accel(period: float) -> float:
    if period <= SPECTRUM[0][0]:
        return SPECTRUM[0][1] * G_ACCEL
    for (t0, a0), (t1, a1) in zip(SPECTRUM, SPECTRUM[1:]):
        if period <= t1 + 1.0e-12:
            if abs(period - t0) < 1.0e-12:
                return a0 * G_ACCEL
            if abs(period - t1) < 1.0e-12:
                return a1 * G_ACCEL
            frac = (period - t0) / (t1 - t0)
            return (a0 + frac * (a1 - a0)) * G_ACCEL
    return SPECTRUM[-1][1] * G_ACCEL


def parse_f06():
    text = F06_PATH.read_text(encoding="utf-8", errors="ignore").splitlines()
    cycles: list[float] = []
    eigvecs: dict[int, dict[int, tuple[float, float, float]]] = {}
    mpf_t1: dict[int, float] = {}
    current_mode: int | None = None
    in_mpf = False
    for raw in text:
        match = EIG_ROW_RE.match(raw)
        if match:
            cycles.append(fnum(match.group(5)))
            continue
        if "OUTPUT FOR EIGENVECTOR" in raw:
            current_mode = int(raw.split()[-1])
            eigvecs[current_mode] = {}
            in_mpf = False
            continue
        if "P A R T I C I P A T I O N" in raw:
            in_mpf = True
            current_mode = None
            continue
        if in_mpf:
            if "E F F E C T I V E" in raw:
                in_mpf = False
                continue
            parts = raw.split()
            if len(parts) == 8 and parts[0].isdigit():
                mode = int(parts[0])
                mpf_t1[mode] = fnum(parts[2])
            continue
        if current_mode is not None:
            match = GRID_ROW_RE.match(raw)
            if match:
                grid = int(match.group(1))
                ux = fnum(match.group(3))
                uz = fnum(match.group(5))
                ry = fnum(match.group(7))
                eigvecs[current_mode][grid] = (ux, uz, ry)
    return cycles[:2], eigvecs, mpf_t1


def local_stiffness(E: float, A: float, I: float, L: float) -> list[list[float]]:
    ea = E * A / L
    e12 = 12.0 * E * I / (L**3)
    e6 = 6.0 * E * I / (L**2)
    e4 = 4.0 * E * I / L
    e2 = 2.0 * E * I / L
    return [
        [ea, 0.0, 0.0, -ea, 0.0, 0.0],
        [0.0, e12, e6, 0.0, -e12, e6],
        [0.0, e6, e4, 0.0, -e6, e2],
        [-ea, 0.0, 0.0, ea, 0.0, 0.0],
        [0.0, -e12, -e6, 0.0, e12, -e6],
        [0.0, e6, e2, 0.0, -e6, e4],
    ]


def mat_vec(mat: list[list[float]], vec: list[float]) -> list[float]:
    return [sum(a * b for a, b in zip(row, vec)) for row in mat]


def modal_scale(mode: int, period: float, cycle: float, mpf: dict[int, float]) -> float:
    omega = 2.0 * math.pi * cycle
    sd = spectrum_accel(period) / (omega * omega)
    return mpf[mode] * sd


def element_end_moment(
    eid: int,
    mode: int,
    eigvecs: dict[int, dict[int, tuple[float, float, float]]],
    scale: float,
) -> tuple[float, float]:
    ga, gb, family = modal.ELEMENTS[eid]
    xa, _, za = modal.JOINTS[ga]
    xb, _, zb = modal.JOINTS[gb]
    dx = xb - xa
    dz = zb - za
    L = math.hypot(dx, dz)
    c = dx / L
    s = dz / L
    theta_sign = -1.0 if abs(dx) < 1.0e-12 else 1.0
    da = eigvecs[mode][ga]
    db = eigvecs[mode][gb]
    u_local = [
        c * da[0] + s * da[1],
        -s * da[0] + c * da[1],
        theta_sign * da[2],
        c * db[0] + s * db[1],
        -s * db[0] + c * db[1],
        theta_sign * db[2],
    ]
    u_local = [scale * value for value in u_local]
    inertia = modal.I_BIG if family == "BIG" else modal.I_SMALL
    force = mat_vec(local_stiffness(modal.E, modal.A_BIG, inertia, L), u_local)
    return force[2], force[5]


def main() -> None:
    cycles, eigvecs, mpf = parse_f06()
    periods = [1.0 / cycle for cycle in cycles]
    print("modal periods")
    for i, (period, ref) in enumerate(zip(periods, REFERENCE_PERIODS), start=1):
        print(f"  mode {i}: T = {period:.6f} s vs ref {ref:.6f} (dT = {period - ref:+.6e})")

    print("response-spectrum Ux (SRSS)")
    for grid in (2, 3):
        total = 0.0
        for mode, (period, cycle) in enumerate(zip(periods, cycles), start=1):
            scale = modal_scale(mode, period, cycle, mpf)
            total += (eigvecs[mode][grid][0] * scale) ** 2
        value = math.sqrt(total)
        ref = REFERENCE_UX[grid]
        print(f"  joint {grid}: {value:.6f} in vs ref {ref:.3f} (d = {value - ref:+.6e})")

    print("response-spectrum M33 end moments (SRSS)")
    for (eid, joint), ref in REFERENCE_M.items():
        total = 0.0
        ga, gb, _ = modal.ELEMENTS[eid]
        for mode, (period, cycle) in enumerate(zip(periods, cycles), start=1):
            scale = modal_scale(mode, period, cycle, mpf)
            mi, mj = element_end_moment(eid, mode, eigvecs, scale)
            value = mi if joint == ga else mj
            total += value * value
        moment = math.sqrt(total)
        print(f"  elem {eid} @ joint {joint}: {abs(moment):.3f} k-in vs ref {ref:.1f} (d = {abs(moment) - ref:+.6e})")


if __name__ == "__main__":
    main()
