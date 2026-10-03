"""
compare_disp_v4c.py
===================
Compare grid-node displacements: Python v4c vs MYSTRAN OP2.
Usage: python compare_disp_v4c.py
"""
import sys, numpy as np
sys.path.insert(0, '/mnt/user-data/uploads')
sys.path.insert(0, '/md/user-data/outputs')
sys.path.insert(0, 'C:/PROJECTAI/18a/python')

from pyNastran.op2.op2 import read_op2

# Import test module
import test_patch2001_q8_t6_v4c as t6

def get_mystran_disp(op2_file):
    """Extract displacements from MYSTRAN OP2 file.
    Returns dict: {subcase_id: {node_id: [6 DOF values]}}
    """
    model = read_op2(op2_file, build_dataframe=False, debug=False)
    result = {}
    for sc_id, disp in model.displacements.items():
        node_ids = disp.node_gridtype[:, 0]
        data = disp.data[0]  # first time step
        result[sc_id] = {}
        for i, nid in enumerate(node_ids):
            result[sc_id][int(nid)] = data[i]
    return result

def get_python_disp(elem_cls, node_xy, conn_dict, bnd_nodes, x0, y0, mode):
    """Run Python v4c and extract grid-node displacements.
    Returns dict: {node_id: [6 DOF values]}
    """
    fem = t6.Model(ndof_per_node=6, autospc=True, autospc_tol=1e-10)
    for nid, (xn, yn) in node_xy.items():
        fem.add_node(nid, xn, yn, 0.0)
    elems = {}
    for eid, conn in conn_dict.items():
        e = elem_cls(eid, [fem.nodes[n] for n in conn], t6.E_, t6.NU_, t6.H_)
        fem.add_element(e)
        elems[eid] = e

    for nid in bnd_nodes:
        xn, yn = node_xy[nid]
        if mode == 'membrane':
            ux, uy = t6.mem_disp(xn, yn, x0, y0)
            fem.add_bc(nid, 0, ux); fem.add_bc(nid, 1, uy)
            for d in [2, 3, 4, 5]: fem.add_bc(nid, d, 0.0)
        else:
            uz, rx, ry = t6.bend_disp(xn, yn, x0, y0)
            fem.add_bc(nid, 0, 0.0); fem.add_bc(nid, 1, 0.0)
            fem.add_bc(nid, 2, uz)
            fem.add_bc(nid, 3, rx); fem.add_bc(nid, 4, ry)
            fem.add_bc(nid, 5, 0.0)
    fem.solve_static()

    result = {}
    u = fem.u
    for nid in node_xy:
        dofs = fem.nodes[nid].dofs
        result[nid] = np.array([u[dofs[i]] for i in range(6)])
    return result

def compare_displacements(py_disp, my_disp, label, mode):
    """Compare two displacement dicts and print differences."""
    print(f"\n{'='*72}")
    print(f"  {label} — {mode.upper()}")
    print(f"{'='*72}")
    print(f"  {'GRID':>5}  {'DOF':>4}  {'Python v4c':>14}  {'MYSTRAN OP2':>14}  {'Diff':>12}  {'%Err':>8}")
    print(f"  {'-'*5}  {'-'*4}  {'-'*14}  {'-'*14}  {'-'*12}  {'-'*8}")

    max_err = 0.0
    max_err_info = None
    n_compared = 0

    all_nodes = sorted(set(py_disp.keys()) | set(my_disp.keys()))
    for nid in all_nodes:
        if nid not in py_disp or nid not in my_disp:
            print(f"  {nid:>5}  MISSING in one source")
            continue
        py_vals = py_disp[nid]
        my_vals = my_disp[nid]
        for dof in range(6):
            pv = py_vals[dof]
            mv = my_vals[dof]
            diff = abs(pv - mv)
            denom = max(abs(pv), abs(mv), 1e-14)
            pct = diff / denom * 100
            if pct > max_err:
                max_err = pct
                max_err_info = (nid, dof, pv, mv, diff, pct)
            if diff > 1e-10:
                print(f"  {nid:>5}  {dof:>4}  {pv:14.6E}  {mv:14.6E}  {diff:12.4E}  {pct:7.2f}%")
            n_compared += 1

    print(f"\n  Total DOFs compared: {n_compared}")
    if max_err_info:
        nid, dof, pv, mv, diff, pct = max_err_info
        print(f"  Max error: {pct:.4f}% at Node {nid}, DOF {dof}")
        print(f"    Python:  {pv:.6E}")
        print(f"    MYSTRAN: {mv:.6E}")
    else:
        print(f"  Max error: 0.0000% — EXACT MATCH")
    return max_err

if __name__ == '__main__':
    print("="*72)
    print("  Displacement Comparison: Python v4c vs MYSTRAN OP2")
    print("="*72)

    # MYSTRAN OP2 files
    op2_files = {
        'MH6T':  'C:/PROJECTAI/18a/test/duel3a_mh6t_1.OP2',
        'MITC6': 'C:/PROJECTAI/18a/test/duel3a_mitc6_1.OP2',
        'SimoT6': 'C:/PROJECTAI/18a/test/duel3a_simo_1.OP2',
    }

    # Python element classes
    from MITC6_Tri_v4 import MITC6_Tri_v4
    from MacNeal_MH6T_Tri_v3 import MacNeal_MH6T_Tri_v3
    try:
        from Simo1993_Tri6_ShellElement_v2 import Simo1993_Tri6_ShellElement_v2
        HAS_SIMO = True
    except ImportError:
        HAS_SIMO = False

    py_classes = {
        'MH6T':  MacNeal_MH6T_Tri_v3,
        'MITC6': MITC6_Tri_v4,
    }
    if HAS_SIMO:
        py_classes['SimoT6'] = Simo1993_Tri6_ShellElement_v2

    for label, op2_file in op2_files.items():
        if label not in py_classes:
            print(f"\nSkipping {label} — Python class not available")
            continue

        print(f"\n{'#'*72}")
        print(f"  {label}")
        print(f"{'#'*72}")

        # MYSTRAN displacements
        my_disp = get_mystran_disp(op2_file)

        # Python displacements
        cls = py_classes[label]
        for mode in ['membrane', 'bending']:
            py_disp = get_python_disp(
                cls, t6.NODE_T6, t6.CONN_T6, t6.BND_T6,
                t6.X0_T6, t6.Y0_T6, mode
            )

            # Compare per subcase
            for sc_id in sorted(my_disp.keys()):
                my_sc = my_disp[sc_id]
                # Determine which subcase this is
                my_mode = 'membrane' if sc_id == 1 else 'bending'
                if my_mode != mode:
                    continue
                compare_displacements(py_disp, my_sc, label, mode)

    print(f"\n{'='*72}")
    print("  Comparison complete")
    print(f"{'='*72}")
