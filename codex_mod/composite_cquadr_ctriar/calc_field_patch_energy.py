import csv
import re
import sys
from pathlib import Path


DOFS = ("T1", "T2", "T3", "R1", "R2", "R3")


def read_enforced(path: Path):
    data = {}
    with path.open(newline="") as f:
        reader = csv.DictReader(f, skipinitialspace=True)
        for row in reader:
            nid = int(row["nid"])
            data[nid] = {
                "T1": float(row["ux"]),
                "T2": float(row["uy"]),
                "T3": float(row["uz"]),
                "R1": float(row["rx"]),
                "R2": float(row["ry"]),
                "R3": float(row["rz"]),
            }
    return data


GRID_RE = re.compile(r"FORCE BALANCE FOR GRID POINT\s+(\d+)\s+IN GLOBAL COORD SYSTEM", re.I)
SPC_RE = re.compile(
    r"SPC FORCE\s+([+-]?\d+\.\d+E[+-]\d+)\s+([+-]?\d+\.\d+E[+-]\d+)\s+([+-]?\d+\.\d+E[+-]\d+)\s+"
    r"([+-]?\d+\.\d+E[+-]\d+)\s+([+-]?\d+\.\d+E[+-]\d+)\s+([+-]?\d+\.\d+E[+-]\d+)",
    re.I,
)


def read_spc_forces(path: Path):
    reactions = {}
    current_grid = None
    with path.open(errors="ignore") as f:
        for line in f:
            grid_match = GRID_RE.search(line)
            if grid_match:
                current_grid = int(grid_match.group(1))
                continue
            if current_grid is None:
                continue
            spc_match = SPC_RE.search(line)
            if spc_match:
                reactions[current_grid] = {
                    dof: float(val) for dof, val in zip(DOFS, spc_match.groups())
                }
                current_grid = None
    return reactions


def compute_energy(enforced, reactions):
    contributions = []
    for nid, u in enforced.items():
        r = reactions.get(nid)
        if r is None:
            continue
        work = sum(u[dof] * r[dof] for dof in DOFS)
        contributions.append((nid, work))
    total = 0.5 * sum(work for _, work in contributions)
    return total, contributions


def main():
    if len(sys.argv) != 3:
        print("usage: py -3 calc_field_patch_energy.py <enf.csv> <f06>")
        raise SystemExit(2)
    enf_path = Path(sys.argv[1])
    f06_path = Path(sys.argv[2])
    enforced = read_enforced(enf_path)
    reactions = read_spc_forces(f06_path)
    total, contributions = compute_energy(enforced, reactions)
    print(f"enforced_nodes={len(enforced)}")
    print(f"reaction_nodes={len(reactions)}")
    print(f"energy_half_uTf={total:.12e}")
    top = sorted(contributions, key=lambda item: abs(item[1]), reverse=True)[:8]
    print("top_node_work:")
    for nid, work in top:
        print(f"  nid={nid:4d} work={work:.12e}")


if __name__ == "__main__":
    main()
