"""
compare_f06_vs_python.py
========================
Compare MYSTRAN F06 local stress (center element) vs Python v4c output.

Usage:
    python test_patch2001_q8_t6_v4c.py > hasilpython.txt
    python compare_f06_vs_python.py

Reads:
    hasilpython.txt          — Python v4c output (redirected)
    duel3a_mh6t_1.F06        — MYSTRAN F06 for MH6T
    duel3a_mitc6_1.F06       — MYSTRAN F06 for MITC6
    duel3a_simo_1.F06        — MYSTRAN F06 for Simo

F06 local stress format:
    E L E M E N T   S T R E S S E S   I N   L O C A L   E L E M E N T   C O O R D I N A T E   S Y S T E M
    Element  Location  Stress  Stresses In Reported Coord System
       51    CENTER  -5.00000E-04  3.21872E+02  9.16097E+02 -7.42781E-02 ...
                        5.00000E-04  3.21872E+02  9.16097E+02 -7.42781E-02 ...

Python v4c center stress format:
       EID            s11            s22            s12  |          sXX            sYY            sXY
        51     1.0575E+03     1.6092E+03    -2.8966E-02  |     1.3333E+03     1.3333E+03     4.0000E+02

Comparison:
    F06 Normal-X  vs Python s11 (LOCAL)
    F06 Normal-Y  vs Python s22 (LOCAL)
    F06 Shear-XY  vs Python s12 (LOCAL)
"""

import re, sys, os

TEST_DIR = r"C:\PROJECTAI\18a\test"
PYTHON_OUT = os.path.join(TEST_DIR, "hasilpython.txt")

F06_FILES = {
    'MH6T':  os.path.join(TEST_DIR, "duel3a_mh6t_1.F06"),
    'MITC6': os.path.join(TEST_DIR, "duel3a_mitc6_1.F06"),
    'Simo':  os.path.join(TEST_DIR, "duel3a_simo_1.F06"),
}

# ═══════════════════════════════════════════════════════════════════
# PARSE F06 — extract center-element local stress
# ═══════════════════════════════════════════════════════════════════
def parse_f06_center_stress(f06_path):
    """Parse F06 file, return dict: eid -> (s11, s22, s12) from CENTER row."""
    results = {}
    with open(f06_path, 'r') as f:
        lines = f.readlines()

    in_local = False
    for i, line in enumerate(lines):
        if 'E L E M E N T   S T R E S S E S   I N   L O C A L' in line:
            in_local = True
            continue
        if in_local:
            # Look for CENTER rows: "       51    CENTER     -5.00000E-04  ..."
            m = re.match(
                r'\s*(\d+)\s+CENTER\s+'
                r'[-+]?\d+\.?\d*E[-+]?\d+\s+'
                r'([-+]?\d+\.?\d*E[-+]?\d+)\s+'
                r'([-+]?\d+\.?\d*E[-+]?\d+)\s+'
                r'([-+]?\d+\.?\d*E[-+]?\d+)',
                line
            )
            if m:
                eid = int(m.group(1))
                s11 = float(m.group(2))
                s22 = float(m.group(3))
                s12 = float(m.group(4))
                # Only store first CENTER row (top fiber) — bottom is same for membrane
                if eid not in results:
                    results[eid] = (s11, s22, s12)
            # Stop at next major section
            if 'S T R E S S E S   A T   G R I D' in line:
                break
    return results

# ═══════════════════════════════════════════════════════════════════
# PARSE PYTHON OUTPUT — extract center stress per formulation
# ═══════════════════════════════════════════════════════════════════
def parse_python_center_stress(python_out_path):
    """Parse hasilpython.txt, return dict: name -> {eid: (s11, s22, s12)}."""
    with open(python_out_path, 'r') as f:
        content = f.read()

    results = {}
    # Split by formulation blocks
    # Pattern: "      {name} CENTER stress (membrane):"
    # More robust: find the header line, then read data rows until blank line
    lines = content.split('\n')
    i = 0
    while i < len(lines):
        line = lines[i]
        m = re.match(r'\s+(\S+)\s+CENTER stress \(membrane\):', line)
        if m:
            name = m.group(1)
            eid_stress = {}
            # Skip header lines (EID, ----, ====)
            j = i + 1
            while j < len(lines):
                data_line = lines[j]
                # Stop at blank line or next section
                if not data_line.strip() or 'CENTER stress' in data_line or 'GRID NODE' in data_line:
                    break
                parts = data_line.split()
                if len(parts) >= 4:
                    try:
                        eid = int(parts[0])
                        s11 = float(parts[1])
                        s22 = float(parts[2])
                        s12 = float(parts[3])
                        eid_stress[eid] = (s11, s22, s12)
                    except ValueError:
                        pass
                j += 1
            results[name] = eid_stress
            i = j
        else:
            i += 1
    return results

# ═══════════════════════════════════════════════════════════════════
# COMPARE
# ═══════════════════════════════════════════════════════════════════
def compare(f06_data, py_data, label):
    """Compare F06 vs Python center stress. Print table."""
    print(f"\n{'='*80}")
    print(f"  {label}")
    print(f"{'='*80}")
    print(f"  {'EID':>4}  {'F06 s11':>12} {'Py s11':>12} {'err%':>7}  "
          f"{'F06 s22':>12} {'Py s22':>12} {'err%':>7}  "
          f"{'F06 s12':>12} {'Py s12':>12} {'err%':>7}")
    print(f"  {'-'*4}  {'-'*12} {'-'*12} {'-'*7}  "
          f"{'-'*12} {'-'*12} {'-'*7}  "
          f"{'-'*12} {'-'*12} {'-'*7}")

    all_match = True
    for eid in sorted(f06_data.keys()):
        f06 = f06_data[eid]
        py = py_data.get(eid)
        if py is None:
            print(f"  {eid:>4}  MISSING in Python")
            all_match = False
            continue

        errs = []
        for k in range(3):
            denom = max(abs(f06[k]), 1e-14)
            errs.append(abs(f06[k] - py[k]) / denom * 100.0)

        match = all(e < 1.0 for e in errs)
        if not match:
            all_match = False

        print(f"  {eid:>4}  {f06[0]:12.4E} {py[0]:12.4E} {errs[0]:7.2f}  "
              f"{f06[1]:12.4E} {py[1]:12.4E} {errs[1]:7.2f}  "
              f"{f06[2]:12.4E} {py[2]:12.4E} {errs[2]:7.2f}  "
              f"{'✓' if match else '✗'}")

    print(f"\n  Overall: {'✓ ALL MATCH (<1%)' if all_match else '✗ MISMATCH DETECTED'}")
    return all_match

# ═══════════════════════════════════════════════════════════════════
# MAIN
# ═══════════════════════════════════════════════════════════════════
if __name__ == '__main__':
    print("="*80)
    print("  MYSTRAN F06 vs Python v4c — Center Element Local Stress Comparison")
    print("="*80)

    # Check files exist
    if not os.path.exists(PYTHON_OUT):
        print(f"\n  ERROR: {PYTHON_OUT} not found.")
        print("  Run: python test_patch2001_q8_t6_v4c.py > hasilpython.txt")
        sys.exit(1)

    for label, path in F06_FILES.items():
        if not os.path.exists(path):
            print(f"\n  ERROR: {path} not found.")
            sys.exit(1)

    # Parse Python output
    py_results = parse_python_center_stress(PYTHON_OUT)
    print(f"\n  Python formulations found: {list(py_results.keys())}")

    # Parse F06 files and compare
    all_ok = True
    for label, f06_path in F06_FILES.items():
        f06_data = parse_f06_center_stress(f06_path)
        print(f"\n  F06 {label}: {len(f06_data)} elements parsed")

        # Find matching Python formulation
        py_key = None
        for k in py_results:
            if label.lower() in k.lower():
                py_key = k
                break
        if py_key is None:
            print(f"  WARNING: No Python formulation matching '{label}' found")
            continue

        py_data = py_results[py_key]
        ok = compare(f06_data, py_data, f"F06 {label} vs Python {py_key}")
        if not ok:
            all_ok = False

    print(f"\n{'='*80}")
    print(f"  FINAL: {'✓ ALL MATCH' if all_ok else '✗ MISMATCH — see above'}")
    print(f"{'='*80}")
