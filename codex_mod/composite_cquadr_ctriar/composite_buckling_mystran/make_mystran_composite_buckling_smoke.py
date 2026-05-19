from __future__ import annotations

from pathlib import Path


OUTDIR = Path(r"D:\fortran\mystran3\codex_mod\composite_buckling_mystran")
OUTDIR.mkdir(parents=True, exist_ok=True)

A = 1.0
H = 0.01
N = 4

# Python composite_buckling_v5_operator anchors for isotropic / oneply_iso
NXX = -7.692307692308e3
NYY = -2.307692307692e3


def nid(i: int, j: int) -> int:
    return j * (N + 1) + i + 1


def edge_lumped(total_line_force_mag: float) -> list[float]:
    dx = A / N
    vals = [0.0] * (N + 1)
    for i in range(N):
        seg = total_line_force_mag * dx
        vals[i] += 0.5 * seg
        vals[i + 1] += 0.5 * seg
    return vals


def shell_common_header(title: str) -> list[str]:
    return [
        f"ID {title}",
        "SOL 5",
        "CEND",
        f"TITLE = {title}",
        "ECHO = NONE",
        "DISP(PRINT) = ALL",`r`n        "SPCFORCE(PRINT) = ALL",`r`n        "GPFORCE(PRINT) = ALL",`r`n        "OLOAD(PRINT) = ALL",`r`n        "STRESS(PRINT) = ALL",`r`n        "STRAIN(PRINT) = ALL",
        "SUBCASE 1",
        "  SPC = 1",
        "  LOAD = 10",
        "SUBCASE 2",
        "  METHOD = 20",
        "BEGIN BULK",
        "PARAM,POST,-1",
        "PARAM,SOLLIB,BANDED",
        "PARAM,AUTOSPC,N",
        "PARAM,GRDPNT,0",
        "EIGRL,20,,,3,,,,MAX",
    ]


def isotropic_property_block(pid: int = 1, mid: int = 1) -> list[str]:
    return [
        "MAT1,1,7.000000E+10,,3.000000E-01",
        f"PSHELL,{pid},{mid},{H:.6E},{mid},1.0,{mid},8.333333E-01",
    ]


def oneply_iso_property_block(pid: int = 1, mid: int = 1) -> list[str]:
    g = 7.0e10 / (2.0 * (1.0 + 0.3))
    return [
        f"MAT8,{mid},7.000000E+10,7.000000E+10,3.000000E-01,{g:.6E},{g:.6E},{g:.6E},0.0",
        f"PCOMP,{pid},,0.0",
        f",{mid},{H:.6E},0.0",
    ]


def grid_block() -> list[str]:
    lines: list[str] = []
    for j in range(N + 1):
        y = A * j / N
        for i in range(N + 1):
            x = A * i / N
            lines.append(f"GRID,{nid(i,j)},,{x:.6E},{y:.6E},0.0")
    return lines


def cquadr_block(pid: int = 1) -> list[str]:
    lines: list[str] = []
    eid = 1
    for j in range(N):
        for i in range(N):
            n1 = nid(i, j)
            n2 = nid(i + 1, j)
            n3 = nid(i + 1, j + 1)
            n4 = nid(i, j + 1)
            lines.append(f"CQUADR,{eid},{pid},{n1},{n2},{n3},{n4}")
            eid += 1
    return lines


def ctriar_block(pid: int = 1) -> list[str]:
    lines: list[str] = []
    eid = 1
    for j in range(N):
        for i in range(N):
            n1 = nid(i, j)
            n2 = nid(i + 1, j)
            n3 = nid(i + 1, j + 1)
            n4 = nid(i, j + 1)
            lines.append(f"CTRIAR,{eid},{pid},{n1},{n2},{n3}")
            eid += 1
            lines.append(f"CTRIAR,{eid},{pid},{n1},{n3},{n4}")
            eid += 1
    return lines


def spc_block() -> list[str]:
    edge_nodes = []
    seen = set()
    for j in range(N + 1):
        for i in range(N + 1):
            if i == 0 or i == N or j == 0 or j == N:
                g = nid(i, j)
                if g not in seen:
                    edge_nodes.append(g)
                    seen.add(g)
    lines = []
    chunk = []
    for g in edge_nodes:
        chunk.append(str(g))
        if len(chunk) == 5:
            lines.append("SPC1,1,3," + ",".join(chunk))
            chunk = []
    if chunk:
        lines.append("SPC1,1,3," + ",".join(chunk))
    lines.append(f"SPC1,1,12,{nid(0,0)}")
    lines.append(f"SPC1,1,2,{nid(N,0)}")
    return lines


def load_block(load_sid: int = 10) -> list[str]:
    fx = edge_lumped(abs(NXX))
    fy = edge_lumped(abs(NYY))
    lines: list[str] = []
    lid = 1
    # left/right edges: x compression
    for j in range(N + 1):
        gl = nid(0, j)
        gr = nid(N, j)
        if fx[j] != 0.0:
            lines.append(f"FORCE,{load_sid},{gl},0,{fx[j]:.6E},1.0,0.0,0.0")
            lines.append(f"FORCE,{load_sid},{gr},0,{fx[j]:.6E},-1.0,0.0,0.0")
        lid += 2
    # bottom/top edges: y compression
    for i in range(N + 1):
        gb = nid(i, 0)
        gt = nid(i, N)
        if fy[i] != 0.0:
            lines.append(f"FORCE,{load_sid},{gb},0,{fy[i]:.6E},0.0,1.0,0.0")
            lines.append(f"FORCE,{load_sid},{gt},0,{fy[i]:.6E},0.0,-1.0,0.0")
        lid += 2
    return lines


def write_case(name: str, prop_lines: list[str], elem_lines: list[str]) -> None:
    lines = shell_common_header(name)
    lines.extend(prop_lines)
    lines.extend(grid_block())
    lines.extend(elem_lines)
    lines.extend(spc_block())
    lines.extend(load_block())
    lines.append("ENDDATA")
    (OUTDIR / f"{name}.dat").write_text("\n".join(lines) + "\n", encoding="ascii")



def sym_crossply_property_block(pid: int = 1, mid: int = 1) -> list[str]:
    return [
        f"MAT8,{mid},1.350000E+11,1.000000E+10,2.800000E-01,5.000000E+09,5.000000E+09,3.800000E+09,0.0",
        f"PCOMP,{pid},,0.0",
        f",{mid},2.500000E-03,0.0",
        f",{mid},2.500000E-03,90.0",
        f",{mid},2.500000E-03,90.0",
        f",{mid},2.500000E-03,0.0",
    ]


def load_block_sym_crossply(load_sid: int = 10) -> list[str]:
    fx = edge_lumped(7.292349645348e3)
    fy = edge_lumped(2.816355725100e2)
    lines: list[str] = []
    for j in range(N + 1):
        gl = nid(0, j)
        gr = nid(N, j)
        if fx[j] != 0.0:
            lines.append(f"FORCE,{load_sid},{gl},0,{fx[j]:.6E},1.0,0.0,0.0")
            lines.append(f"FORCE,{load_sid},{gr},0,{fx[j]:.6E},-1.0,0.0,0.0")
    for i in range(N + 1):
        gb = nid(i, 0)
        gt = nid(i, N)
        if fy[i] != 0.0:
            lines.append(f"FORCE,{load_sid},{gb},0,{fy[i]:.6E},0.0,1.0,0.0")
            lines.append(f"FORCE,{load_sid},{gt},0,{fy[i]:.6E},0.0,-1.0,0.0")
    return lines


def write_case_with_loads(name: str, prop_lines: list[str], elem_lines: list[str], load_lines: list[str]) -> None:
    lines = shell_common_header(name)
    lines.extend(prop_lines)
    lines.extend(grid_block())
    lines.extend(elem_lines)
    lines.extend(spc_block())
    lines.extend(load_lines)
    lines.append("ENDDATA")
    (OUTDIR / f"{name}.dat").write_text("\n".join(lines) + "\n", encoding="ascii")
def main() -> None:
    write_case("comp_buckling_cquadr_iso_n4", isotropic_property_block(), cquadr_block())
    write_case("comp_buckling_cquadr_oneply_iso_n4", oneply_iso_property_block(), cquadr_block())
    write_case("comp_buckling_ctriar_iso_n4", isotropic_property_block(), ctriar_block())
    write_case("comp_buckling_ctriar_oneply_iso_n4", oneply_iso_property_block(), ctriar_block())
    write_case_with_loads("comp_buckling_cquadr_pcomp_sym_n4", sym_crossply_property_block(), cquadr_block(), load_block_sym_crossply())
    write_case_with_loads("comp_buckling_ctriar_pcomp_sym_n4", sym_crossply_property_block(), ctriar_block(), load_block_sym_crossply())


if __name__ == "__main__":
    main()




