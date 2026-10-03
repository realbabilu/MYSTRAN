"""
test_result_2_001_t6_disp.py
=============================
Compare T6 displacements: Python (testpython.txt) vs MYSTRAN (F06 + OP2).

Only T6 elements are compared (MITC6, MacNeal_MH6T, SimoT6).
Q8 elements are skipped.

Sources:
  - testpython.txt  : Python reference
  - duel3a_mh6t_1.F06 : MYSTRAN MH6T membrane displacement (SUBCASE 1)
  - duel3a_mitc6_1.F06 : MYSTRAN MITC6 membrane displacement (SUBCASE 1)
  - duel3a_simo_1.F06 : MYSTRAN Simo T6 membrane displacement (SUBCASE 1)
  - duel3a_mh6t_1.OP2 : MYSTRAN MH6T bending displacement (SUBCASE 2)
  - duel3a_mitc6_1.OP2 : MYSTRAN MITC6 bending displacement (SUBCASE 2)
  - duel3a_simo_1.OP2 : MYSTRAN Simo T6 bending displacement (SUBCASE 2)

Output: result_disp_t6_2_001.txt
"""

import sys, os, re, struct

try:
    sys.stdout.reconfigure(encoding='utf-8')
except AttributeError:
    pass

# ═══════════════════════════════════════════════════════════════════
# PATHS
# ═══════════════════════════════════════════════════════════════════
PYTHON_TXT = r"C:\PROJECTAI\18a\python\testpython.txt"
TEST_DIR   = r"C:\PROJECTAI\18a\test"
OUT_TXT    = r"C:\PROJECTAI\18a\python\result_disp_t6_2_001.txt"

F06_FILES = {
    'MH6T':  os.path.join(TEST_DIR, 'duel3a_mh6t_1.F06'),
    'MITC6': os.path.join(TEST_DIR, 'duel3a_mitc6_1.F06'),
    'Simo':  os.path.join(TEST_DIR, 'duel3a_simo_1.F06'),
}

OP2_FILES = {
    'MH6T':  os.path.join(TEST_DIR, 'duel3a_mh6t_1.OP2'),
    'MITC6': os.path.join(TEST_DIR, 'duel3a_mitc6_1.OP2'),
    'Simo':  os.path.join(TEST_DIR, 'duel3a_simo_1.OP2'),
}

# ═══════════════════════════════════════════════════════════════════
# PARSE PYTHON testpython.txt
# ═══════════════════════════════════════════════════════════════════
def parse_python_txt(path):
    """Extract T6 GRID NODE DISPLACEMENTS blocks from testpython.txt."""
    with open(path, 'r', encoding='utf-8') as f:
        lines = f.readlines()

    results = {}  # {name: {mode: {grid: [6 values]}}}

    i = 0
    while i < len(lines):
        line = lines[i]
        # Match: "      MITC6 GRID NODE DISPLACEMENTS (membrane):"
        m = re.match(r'\s+(\S+)\s+GRID NODE DISPLACEMENTS\s+\((\w+)\):', line)
        if m:
            name = m.group(1)
            mode = m.group(2)
            if name not in results:
                results[name] = {}
            results[name][mode] = {}

            # Parse displacement rows until blank line or next section
            i += 1
            while i < len(lines):
                row = lines[i].strip()
                if not row or row.startswith('(') or row.startswith('[') or row.startswith('='):
                    break
                # Match: "        301    1.0000   2.0000   0.000000E+00   0.000000E+00 ..."
                parts = row.split()
                if len(parts) >= 9:
                    try:
                        grid = int(parts[0])
                        vals = [float(x) for x in parts[3:9]]
                        results[name][mode][grid] = vals
                    except ValueError:
                        pass
                i += 1
        else:
            i += 1

    return results

# ═══════════════════════════════════════════════════════════════════
# PARSE F06 DISPLACEMENT (SUBCASE 1 = membrane)
# ═══════════════════════════════════════════════════════════════════
def parse_f06_displacement(path):
    """Extract displacement table from F06 (SUBCASE 1 = membrane)."""
    with open(path, 'r', encoding='utf-8', errors='replace') as f:
        lines = f.readlines()

    disp = {}
    in_disp = False
    for i, line in enumerate(lines):
        if 'D I S P L A C E M E N T S' in line:
            in_disp = True
            continue
        if in_disp:
            # Match: "            301        0  0.0           0.0           0.0           0.0           0.0          -3.545151E-17"
            m = re.match(r'\s+(\d+)\s+0\s+([-\d.E+]+)\s+([-\d.E+]+)\s+([-\d.E+]+)\s+([-\d.E+]+)\s+([-\d.E+]+)\s+([-\d.E+]+)', line)
            if m:
                grid = int(m.group(1))
                vals = [float(m.group(j)) for j in range(2, 8)]
                disp[grid] = vals
            elif 'MAX*' in line or 'MIN*' in line or 'ABS*' in line:
                break
    return disp

# ═══════════════════════════════════════════════════════════════════
# PARSE OP2 DISPLACEMENT (SUBCASE 2 = bending)
# ═══════════════════════════════════════════════════════════════════
def parse_op2_displacement(path):
    """Extract displacement from OP2 binary (SUBCASE 2 = bending)."""
    with open(path, 'rb') as f:
        data = f.read()

    disp = {}
    pos = 0
    while pos < len(data) - 12:
        # OP2 record: 4 bytes length, 4 bytes name_len, 4 bytes name
        rec_len = struct.unpack('<i', data[pos:pos+4])[0]
        if rec_len <= 0 or pos + rec_len > len(data):
            pos += 1
            continue
        name_len = struct.unpack('<i', data[pos+4:pos+8])[0]
        name = data[pos+8:pos+8+name_len].decode('ascii', errors='replace')

        if name == 'OUGV1':
            # Displacement table
            data_start = pos + 8 + name_len
            data_end = pos + rec_len
            p = data_start
            while p < data_end - 4:
                grid = struct.unpack('<i', data[p:p+4])[0]
                if grid == 0:
                    p += 4
                    continue
                # 6 DOF values (4 bytes each)
                vals = struct.unpack('<6f', data[p+4:p+28])
                disp[grid] = list(vals)
                p += 28
            break
        pos += rec_len
    return disp

# ═══════════════════════════════════════════════════════════════════
# COMPARISON
# ═══════════════════════════════════════════════════════════════════
def compare_displacements(py_disp, my_disp, label_py, label_my, mode):
    """Compare two displacement dicts and return report lines."""
    lines = []
    lines.append(f"  {label_py} vs {label_my} ({mode}):")
    lines.append(f"    {'GRID':>5}  {'DOF':>4}  {'Python':>14}  {'MYSTRAN':>14}  {'Diff':>12}  {'%Err':>8}")
    lines.append(f"    {'-'*5}  {'-'*4}  {'-'*14}  {'-'*14}  {'-'*12}  {'-'*8}")

    all_grids = sorted(set(py_disp.keys()) | set(my_disp.keys()))
    max_err = 0.0
    n_match = 0
    n_total = 0

    for grid in all_grids:
        py_vals = py_disp.get(grid, [0.0]*6)
        my_vals = my_disp.get(grid, [0.0]*6)
        for dof in range(6):
            pv = py_vals[dof]
            mv = my_vals[dof]
            diff = abs(pv - mv)
            denom = max(abs(pv), abs(mv), 1e-30)
            pct = diff / denom * 100
            if pct > max_err:
                max_err = pct
            n_total += 1
            if pct < 0.01:
                n_match += 1
            lines.append(f"    {grid:>5}  {dof:>4}  {pv:14.6E}  {mv:14.6E}  {diff:12.4E}  {pct:8.4f}%")

    lines.append(f"    {'-'*70}")
    lines.append(f"    Summary: {n_match}/{n_total} DOFs match (<0.01% err), max err = {max_err:.4f}%")
    if max_err < 0.01:
        lines.append(f"    RESULT: MATCH")
    else:
        lines.append(f"    RESULT: MISMATCH")
    lines.append("")
    return lines, max_err

# ═══════════════════════════════════════════════════════════════════
# MAIN
# ═══════════════════════════════════════════════════════════════════
def main():
    out = []
    out.append("=" * 72)
    out.append("  T6 Displacement Comparison: Python vs MYSTRAN")
    out.append("  Problem: SAP2000 2-001 (MacNeal-Harder patch test)")
    out.append("=" * 72)
    out.append("")

    # Parse Python results
    py_results = parse_python_txt(PYTHON_TXT)
    out.append(f"Parsed Python results from: {PYTHON_TXT}")
    out.append(f"  T6 elements found: {list(py_results.keys())}")
    out.append("")

    # Parse F06 (membrane) and OP2 (bending)
    f06_disp = {}
    op2_disp = {}
    for name, path in F06_FILES.items():
        if os.path.exists(path):
            f06_disp[name] = parse_f06_displacement(path)
            out.append(f"Parsed F06 {name}: {len(f06_disp[name])} grids")
    for name, path in OP2_FILES.items():
        if os.path.exists(path):
            op2_disp[name] = parse_op2_displacement(path)
            out.append(f"Parsed OP2 {name}: {len(op2_disp[name])} grids")
    out.append("")

    # Compare membrane (Python vs F06) - ONLY T6 elements
    out.append("=" * 72)
    out.append("  MEMBRANE DISPLACEMENT COMPARISON (SUBCASE 1)")
    out.append("=" * 72)
    out.append("")

    membrane_results = {}
    for py_name in py_results:
        # Only compare T6 elements
        if 'T6' not in py_name and 'tri' not in py_name.lower():
            continue
        if 'membrane' not in py_results[py_name]:
            continue
        py_mem = py_results[py_name]['membrane']

        # Map Python name to MYSTRAN name
        my_name = None
        if 'MITC6' in py_name:
            my_name = 'MITC6'
        elif 'MH6T' in py_name:
            my_name = 'MH6T'
        elif 'Simo' in py_name or 'SIMO' in py_name:
            my_name = 'Simo'

        if my_name and my_name in f06_disp:
            lines, max_err = compare_displacements(py_mem, f06_disp[my_name], py_name, f"F06-{my_name}", "membrane")
            out.extend(lines)
            membrane_results[py_name] = max_err

    # Compare bending (Python vs OP2) - ONLY T6 elements
    out.append("=" * 72)
    out.append("  BENDING DISPLACEMENT COMPARISON (SUBCASE 2)")
    out.append("=" * 72)
    out.append("")

    bending_results = {}
    for py_name in py_results:
        # Only compare T6 elements
        if 'T6' not in py_name and 'tri' not in py_name.lower():
            continue
        if 'bending' not in py_results[py_name]:
            continue
        py_bend = py_results[py_name]['bending']

        my_name = None
        if 'MITC6' in py_name:
            my_name = 'MITC6'
        elif 'MH6T' in py_name:
            my_name = 'MH6T'
        elif 'Simo' in py_name or 'SIMO' in py_name:
            my_name = 'Simo'

        if my_name and my_name in op2_disp:
            lines, max_err = compare_displacements(py_bend, op2_disp[my_name], py_name, f"OP2-{my_name}", "bending")
            out.extend(lines)
            bending_results[py_name] = max_err

    # Summary
    out.append("=" * 72)
    out.append("  SUMMARY")
    out.append("=" * 72)
    out.append("")
    out.append("  Membrane (Python vs F06):")
    for name, err in membrane_results.items():
        status = "MATCH" if err < 0.01 else "MISMATCH"
        out.append(f"    {name:20s}: max err = {err:.4f}%  [{status}]")
    out.append("")
    out.append("  Bending (Python vs OP2):")
    for name, err in bending_results.items():
        status = "MATCH" if err < 0.01 else "MISMATCH"
        out.append(f"    {name:20s}: max err = {err:.4f}%  [{status}]")
    out.append("")
    out.append("=" * 72)

    # Write output
    report = "\n".join(out)
    with open(OUT_TXT, 'w', encoding='utf-8') as f:
        f.write(report)

    print(report)
    print(f"\nReport saved to: {OUT_TXT}")

if __name__ == '__main__':
    main()
