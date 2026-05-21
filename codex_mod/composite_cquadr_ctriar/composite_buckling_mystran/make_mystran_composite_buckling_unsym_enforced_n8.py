from __future__ import annotations

from pathlib import Path


OUTDIR = Path(
    r"D:\fortran\mystran3\MYSTRANSolver-18.0.0\codex_mod\composite_cquadr_ctriar\composite_buckling_mystran"
)
OUTDIR.mkdir(parents=True, exist_ok=True)

A = 1.0
H = 0.01
N = 8
EPS_X = -1.0e-5


def nid(i: int, j: int) -> int:
    return j * (N + 1) + i + 1


def shell_common_header(title: str, enf_name: str) -> list[str]:
    return [
        f"ID {title}",
        "SOL 5",
        "CEND",
        f"TITLE = {title}",
        "ECHO = NONE",
        f"ENFORCED = {enf_name}",
        "DISP(PRINT) = ALL",
        "SPCFORCE(PRINT) = ALL",
        "GPFORCE(PRINT) = ALL",
        "OLOAD(PRINT) = ALL",
        "STRESS(PRINT) = ALL",
        "STRAIN(PRINT) = ALL",
        "SUBCASE 1",
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


def unsym_crossply_property_block(pid: int = 1, mid: int = 1) -> list[str]:
    return [
        f"MAT8,{mid},1.350000E+11,1.000000E+10,2.800000E-01,5.000000E+09,5.000000E+09,3.800000E+09,0.0",
        f"PCOMP,{pid},,0.0",
        f",{mid},5.000000E-03,0.0",
        f",{mid},5.000000E-03,90.0",
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


def enforced_lines() -> list[str]:
    lines = ["UNSYM prescribed eps_x state"]
    for j in range(N + 1):
        y = A * j / N
        for i in range(N + 1):
            x = A * i / N
            g = nid(i, j)
            ux = EPS_X * x
            uy = 0.0 * y
            uz = 0.0
            rx = 0.0
            ry = 0.0
            rz = 0.0
            lines.append(
                f"{g}, {ux:.9E}, {uy:.9E}, {uz:.9E}, {rx:.9E}, {ry:.9E}, {rz:.9E}"
            )
    return lines


def write_case(name: str, elem_lines: list[str]) -> None:
    enf_name = f"{name}.enf"
    lines = shell_common_header(name, enf_name)
    lines.extend(unsym_crossply_property_block())
    lines.extend(grid_block())
    lines.extend(elem_lines)
    lines.append("FORCE,10,1,0,1.000000E-12,1.0,0.0,0.0")
    lines.append("ENDDATA")
    (OUTDIR / f"{name}.dat").write_text("\n".join(lines) + "\n", encoding="ascii")
    (OUTDIR / enf_name).write_text("\n".join(enforced_lines()) + "\n", encoding="ascii")


def main() -> None:
    write_case("comp_buckling_cquadr_pcomp_unsym_n8_enforced", cquadr_block())
    write_case("comp_buckling_ctriar_pcomp_unsym_n8_enforced", ctriar_block())


if __name__ == "__main__":
    main()
