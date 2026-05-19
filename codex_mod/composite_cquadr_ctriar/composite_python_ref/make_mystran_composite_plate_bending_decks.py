from pathlib import Path


ROOT = Path(r"D:\fortran\mystran3\codex_mod\composite_plate_bending_mystran")
ROOT.mkdir(parents=True, exist_ok=True)


E1 = 135.0e9
E2 = 10.0e9
NU12 = 0.28
G12 = 5.0e9
G1Z = 5.0e9
G2Z = 3.8e9
RHO = 0.0

A = 1.0
H = 0.01
PLY_T = H / 4.0
PRESSURE = 1000.0


def nid_ij(i: int, j: int, n: int) -> int:
    return j * (n + 1) + i + 1


def fmt(x: float) -> str:
    s = f"{x:.10g}"
    if "." not in s and "E" not in s and "e" not in s:
        s += ".0"
    return s


def grid_lines(n: int) -> list[str]:
    lines = []
    for j in range(n + 1):
        y = A * j / n
        for i in range(n + 1):
            x = A * i / n
            nid = nid_ij(i, j, n)
            lines.append(f"GRID,{nid},,{fmt(x)},{fmt(y)},0.0")
    return lines


def cquadr_lines(n: int) -> list[str]:
    lines = []
    eid = 1
    for j in range(n):
        for i in range(n):
            n1 = nid_ij(i, j, n)
            n2 = nid_ij(i + 1, j, n)
            n3 = nid_ij(i + 1, j + 1, n)
            n4 = nid_ij(i, j + 1, n)
            lines.append(f"CQUADR,{eid},1,{n1},{n2},{n3},{n4}")
            eid += 1
    return lines


def ctriar_lines(n: int) -> list[str]:
    lines = []
    eid = 1
    for j in range(n):
        for i in range(n):
            n1 = nid_ij(i, j, n)
            n2 = nid_ij(i + 1, j, n)
            n3 = nid_ij(i + 1, j + 1, n)
            n4 = nid_ij(i, j + 1, n)
            lines.append(f"CTRIAR,{eid},1,{n1},{n2},{n3}")
            eid += 1
            lines.append(f"CTRIAR,{eid},1,{n1},{n3},{n4}")
            eid += 1
    return lines


def spc_lines(n: int) -> list[str]:
    lines = []
    # w = 0 on all edges
    for j in range(n + 1):
        for i in range(n + 1):
            if i == 0 or i == n or j == 0 or j == n:
                nid = nid_ij(i, j, n)
                lines.append(f"SPC1,1,3,{nid}")
    # corner (0,0): ux=uy=0
    lines.append("SPC1,1,12,1")
    # corner (a,0): uy=0
    lines.append(f"SPC1,1,2,{nid_ij(n, 0, n)}")
    return lines


def nodal_force_map(n: int, tri: bool) -> dict[int, float]:
    loads = {nid_ij(i, j, n): 0.0 for j in range(n + 1) for i in range(n + 1)}
    dx = A / n
    dy = A / n
    if tri:
        elem_area = 0.5 * dx * dy
        share = -PRESSURE * elem_area / 3.0
        for j in range(n):
            for i in range(n):
                n1 = nid_ij(i, j, n)
                n2 = nid_ij(i + 1, j, n)
                n3 = nid_ij(i + 1, j + 1, n)
                n4 = nid_ij(i, j + 1, n)
                for nid in (n1, n2, n3):
                    loads[nid] += share
                for nid in (n1, n3, n4):
                    loads[nid] += share
    else:
        elem_area = dx * dy
        share = -PRESSURE * elem_area / 4.0
        for j in range(n):
            for i in range(n):
                n1 = nid_ij(i, j, n)
                n2 = nid_ij(i + 1, j, n)
                n3 = nid_ij(i + 1, j + 1, n)
                n4 = nid_ij(i, j + 1, n)
                for nid in (n1, n2, n3, n4):
                    loads[nid] += share
    return loads


def force_lines(n: int, tri: bool) -> list[str]:
    lines = []
    for nid, fz in nodal_force_map(n, tri).items():
        if abs(fz) > 0.0:
            lines.append(f"FORCE,10,{nid},0,{fz:.6E},0.0,0.0,1.0")
    return lines


def header(title: str) -> list[str]:
    return [
        f"ID {title}",
        "SOL 1",
        "CEND",
        f"TITLE = {title}",
        "ECHO = NONE",
        "DISPLACEMENT(PRINT) = ALL",
        "SPCFORCE(PRINT) = ALL",
        "GPFORCE(PRINT) = ALL",
        "OLOAD(PRINT) = ALL",
        "ELDATA(0,PRINT) = ALL",
        "SUBCASE 1",
        f"  LABEL = {title}",
        "  SPC = 1",
        "  LOAD = 10",
        "BEGIN BULK",
        "PARAM,POST,-1",
        "PARAM,AUTOSPC,N",
        "PARAM,GRDPNT,0",
        "PARAM,K6ROT,1.0E+00",
        "PARAM,WTMASS,1.0",
        f"MAT8,1,{E1:.6E},{E2:.6E},{NU12:.6E},{G12:.6E},{G1Z:.6E},{G2Z:.6E},{RHO:.6E}",
        "PCOMP,1,,0.0",
        f",1,{PLY_T:.6E},0.0",
        f",1,{PLY_T:.6E},90.0",
        f",1,{PLY_T:.6E},90.0",
        f",1,{PLY_T:.6E},0.0",
    ]


def write_deck(kind: str, n: int) -> None:
    tri = kind == "ctriar"
    title = f"COMP PLATE BENDING {kind.upper()} PCOMP N{n}"
    lines = header(title)
    lines.extend(grid_lines(n))
    if tri:
        lines.extend(ctriar_lines(n))
    else:
        lines.extend(cquadr_lines(n))
    lines.extend(force_lines(n, tri))
    lines.extend(spc_lines(n))
    lines.append("ENDDATA")
    out = ROOT / f"comp_plate_bending_{kind}_pcomp_n{n}.dat"
    out.write_text("\n".join(lines) + "\n", encoding="ascii")


def main() -> None:
    for kind in ("cquadr", "ctriar"):
        for n in (4, 8):
            write_deck(kind, n)


if __name__ == "__main__":
    main()
