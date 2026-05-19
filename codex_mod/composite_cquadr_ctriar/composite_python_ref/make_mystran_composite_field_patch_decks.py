from pathlib import Path


ROOT = Path(
    r"D:\fortran\mystran3\MYSTRANSolver-18.0.0\codex_mod\composite_cquadr_ctriar\composite_field_patch_mystran"
)
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

EPS_COUPLED = (1.0e-4, -0.5e-4, 0.25e-4)
KAP_COUPLED = (2.0e-3, -1.0e-3, 0.5e-3)
EPS_MEMBRANE = EPS_COUPLED
KAP_MEMBRANE = (0.0, 0.0, 0.0)
EPS_BENDING = (0.0, 0.0, 0.0)
KAP_BENDING = KAP_COUPLED


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


def displacement_field(x: float, y: float, eps, kap):
    ex, ey, gxy = eps
    kx, ky, kxy = kap

    ux = ex * x + 0.5 * gxy * y
    uy = ey * y + 0.5 * gxy * x

    w = -0.5 * kx * x * x - 0.5 * ky * y * y - 0.5 * kxy * x * y
    dwdx = -kx * x - 0.5 * kxy * y
    dwdy = -ky * y - 0.5 * kxy * x

    rx = dwdy
    ry = -dwdx
    rz = 0.0

    return ux, uy, w, rx, ry, rz


def enf_lines(n: int, eps, kap) -> list[str]:
    lines = ["nid, ux, uy, uz, rx, ry, rz"]
    for j in range(n + 1):
        y = A * j / n
        for i in range(n + 1):
            x = A * i / n
            nid = nid_ij(i, j, n)
            vals = displacement_field(x, y, eps, kap)
            txt = ",".join([str(nid)] + [f"{v:.12e}" for v in vals])
            lines.append(txt)
    return lines


def header(title: str, enf_name: str) -> list[str]:
    return [
        f"ID {title}",
        "SOL 1",
        "CEND",
        f"TITLE = {title}",
        f"ENFORCED = {enf_name}",
        "ECHO = NONE",
        "DISPLACEMENT(PRINT) = ALL",
        "SPCFORCE(PRINT) = ALL",
        "GPFORCE(PRINT) = ALL",
        "OLOAD(PRINT) = ALL",
        "STRESS(PRINT) = ALL",
        "STRAIN(PRINT) = ALL",
        "ELDATA(0,PRINT) = ALL",
        "SUBCASE 1",
        f"  LABEL = {title}",
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
        f",1,{PLY_T:.6E},0.0",
        f",1,{PLY_T:.6E},90.0",
        "ENDDATA",
    ]


def laminate_lines(laminate: str) -> list[str]:
    if laminate == "sym_crossply":
        return [
            "PCOMP,1,,0.0",
            f",1,{PLY_T:.6E},0.0",
            f",1,{PLY_T:.6E},90.0",
            f",1,{PLY_T:.6E},90.0",
            f",1,{PLY_T:.6E},0.0",
        ]
    if laminate == "unsym_crossply":
        t = H / 2.0
        return [
            "PCOMP,1,,0.0",
            f",1,{t:.6E},0.0",
            f",1,{t:.6E},90.0",
        ]
    raise ValueError(laminate)


def case_eps_kap(case: str):
    if case == "coupled":
        return EPS_COUPLED, KAP_COUPLED
    if case == "membrane":
        return EPS_MEMBRANE, KAP_MEMBRANE
    if case == "bending":
        return EPS_BENDING, KAP_BENDING
    raise ValueError(case)


def write_deck(kind: str, laminate: str, case: str, n: int) -> None:
    title = f"COMP FIELD PATCH {kind.upper()} {laminate.upper()} {case.upper()} N{n}"
    stem = f"comp_field_patch_{kind}_{laminate}_{case}_n{n}"
    enf_name = f"{stem}.enf"
    eps, kap = case_eps_kap(case)

    lines = [
        f"ID {title}",
        "SOL 1",
        "CEND",
        f"TITLE = {title}",
        f"ENFORCED = {enf_name}",
        "ECHO = NONE",
        "DISPLACEMENT(PRINT) = ALL",
        "SPCFORCE(PRINT) = ALL",
        "GPFORCE(PRINT) = ALL",
        "OLOAD(PRINT) = ALL",
        "STRESS(PRINT) = ALL",
        "STRAIN(PRINT) = ALL",
        "ELDATA(0,PRINT) = ALL",
        "SUBCASE 1",
        f"  LABEL = {title}",
        "BEGIN BULK",
        "PARAM,POST,-1",
        "PARAM,AUTOSPC,N",
        "PARAM,GRDPNT,0",
        "PARAM,K6ROT,1.0E+00",
        "PARAM,WTMASS,1.0",
        f"MAT8,1,{E1:.6E},{E2:.6E},{NU12:.6E},{G12:.6E},{G1Z:.6E},{G2Z:.6E},{RHO:.6E}",
    ]
    lines.extend(laminate_lines(laminate))
    lines.extend(grid_lines(n))
    if kind == "cquadr":
        lines.extend(cquadr_lines(n))
    else:
        lines.extend(ctriar_lines(n))
    lines.append("ENDDATA")

    (ROOT / f"{stem}.dat").write_text("\n".join(lines) + "\n", encoding="ascii")
    (ROOT / enf_name).write_text("\n".join(enf_lines(n, eps, kap)) + "\n", encoding="ascii")


def main() -> None:
    for kind in ("cquadr", "ctriar"):
        for laminate in ("sym_crossply", "unsym_crossply"):
            for case in ("membrane", "bending", "coupled"):
                write_deck(kind, laminate, case, 4)


if __name__ == "__main__":
    main()
