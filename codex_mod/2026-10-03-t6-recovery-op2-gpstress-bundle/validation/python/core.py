"""
core.py — FEM Engine dengan AutoSPC
=====================================
Node, DOF, Boundary Conditions, Assembler, Solver, AutoSPC.

AutoSPC (Automatic Single Point Constraint):
  Mendeteksi dan mengunci DOF yang tidak terkonstrain (zero-stiffness rows)
  di K global sebelum solve. Penting untuk model shell 3D di mana:
    - Elemen flat memberikan zero stiffness untuk in-plane rotation (drilling)
    - Node bebas di udara tanpa elemen yang mengisi semua DOF
    - Struktur yang secara geometri singular (mechanism)
  AutoSPC mencegah LinAlgError singular dan memberikan warning yang informatif.
"""

import numpy as np
from dataclasses import dataclass
from typing import Optional


# ═══════════════════════════════════════════════════════════════
#  NODE & DATA CLASSES
# ═══════════════════════════════════════════════════════════════

class Node:
    def __init__(self, nid: int, x: float, y: float, z: float):
        self.nid = nid
        self.x   = float(x)
        self.y   = float(y)
        self.z   = float(z)
        self.dofs: list = []

    @property
    def coords(self) -> np.ndarray:
        return np.array([self.x, self.y, self.z])

    def __repr__(self):
        return f"Node({self.nid}: [{self.x:.4g}, {self.y:.4g}, {self.z:.4g}])"


@dataclass
class BC:
    dof  : int
    value: float = 0.0

@dataclass
class NodalLoad:
    dof  : int
    value: float

@dataclass
class _PendingBC:
    node_id  : int
    dof_local: int
    value    : float = 0.0

@dataclass
class _PendingLoad:
    node_id  : int
    dof_local: int
    value    : float


# ═══════════════════════════════════════════════════════════════
#  ELEMENT BASE CLASS
# ═══════════════════════════════════════════════════════════════

class Element:
    eid  : int
    nodes: list

    @property
    def ndof_per_node(self) -> int:
        raise NotImplementedError

    def k_local(self) -> np.ndarray:
        raise NotImplementedError

    def m_local(self) -> Optional[np.ndarray]:
        return None

    def T_matrix(self) -> np.ndarray:
        n = len(self.nodes) * self.ndof_per_node
        return np.eye(n)

    def k_global(self) -> np.ndarray:
        T  = self.T_matrix()
        Kl = self.k_local()
        return T.T @ Kl @ T

    def m_global(self) -> Optional[np.ndarray]:
        Ml = self.m_local()
        if Ml is None:
            return None
        T = self.T_matrix()
        return T.T @ Ml @ T

    def global_dof_indices(self) -> list:
        dofs = []
        for node in self.nodes:
            dofs.extend(node.dofs[:self.ndof_per_node])
        return dofs

    # ── Stress/moment recovery (center + corners, local -> global) ──
    #
    # Template method: this is the SAME for every element class (loop
    # over natural-coordinate sample points, compute local strain via
    # Bm/Bb/shear, apply isotropic constitutive law, rotate to global).
    # What genuinely differs between element classes is WHICH local frame
    # each class's Bm/Bb actually used at a given point -- verified in
    # this codebase to vary: DKMQ24_ShellElement_RHR and
    # DKMQ24_MystranCQUADR_ShellElement_RHR recompute a frame fresh at
    # every Gauss point (_local_basis_at_point), Simo1993_ShellElement_v1p6
    # uses ONE frame fixed at the centroid for the whole element, and
    # DKMQ24_EAS4_ShellElement_RHR is conditional on its own flatness
    # check. Reporting stress in the WRONG frame for a given class is
    # exactly the class of bug that broke the EAS patch test earlier in
    # this codebase's history -- so this is factored out as an explicit,
    # overridable hook (bm_frame) rather than assumed.
    #
    # Subclasses only need to override:
    #   - bm_frame(r1, r2): if their Bm/Bb use something other than
    #     _local_basis_at_point(r1, r2) (e.g. a frame fixed at centroid).
    #   - _shear_B(r1, r2): if their shear B-matrix method isn't named
    #     _compute_Bs (e.g. Simo1993's _compute_Bs_ANS).
    # k_local() etc. do NOT need to be touched.

    CORNER_POINTS = {
        "center": (0.0, 0.0),
        "corner1": (-1.0, -1.0),
        "corner2": (1.0, -1.0),
        "corner3": (1.0, 1.0),
        "corner4": (-1.0, 1.0),
    }

    def bm_frame(self, r1: float, r2: float):
        """Local (e1, e2, e3) frame that this element's Bm/Bb actually
        used at natural coordinates (r1, r2). Default: delegate to
        _local_basis_at_point (correct for DKMQ24_ShellElement_RHR and
        DKMQ24_MystranCQUADR_ShellElement_RHR). Override in any subclass
        whose Bm/Bb use a different frame (fixed-centroid, conditional,
        etc.) -- do not assume the default is correct without checking
        how that class's _compute_Bm is actually built."""
        return self._local_basis_at_point(r1, r2)

    def _shear_B(self, r1: float, r2: float) -> Optional[np.ndarray]:
        """Uniform accessor for the shear B-matrix, regardless of what a
        given subclass happens to name its own method. Returns None if
        the element has no shear term at all (caller should skip shear
        recovery in that case)."""
        if hasattr(self, "_compute_Bs"):
            return self._compute_Bs(r1, r2)
        if hasattr(self, "_compute_Bs_ANS"):
            return self._compute_Bs_ANS(r1, r2)
        return None

    def _isotropic_constitutive(self):
        """Standard isotropic membrane/bending/shear constitutive
        matrices from E, nu, h -- identical formulas across every element
        class in this codebase (DKMQ24 family and Simo1993 alike), so
        provided once here rather than duplicated per class."""
        h_avg = float(np.mean(self.h)) if hasattr(self.h, "__len__") else float(self.h)
        factor = self.E / (1.0 - self.nu ** 2)
        C_mb = np.array([
            [factor,           factor * self.nu, 0.0],
            [factor * self.nu, factor,            0.0],
            [0.0,              0.0,               factor * (1.0 - self.nu) / 2.0],
        ])
        C_b = C_mb * (h_avg ** 2 / 12.0)  # bending constitutive = membrane * h^2/12 for isotropic
        k_shear = getattr(self, "k_shear", 5.0 / 6.0)
        G = getattr(self, "G", self.E / (2.0 * (1.0 + self.nu)))
        C_s = k_shear * G * np.eye(2)
        return C_mb, C_b, C_s, h_avg

    def _element_u(self, u_global: np.ndarray) -> np.ndarray:
        """Extract this element's local DOF vector from the full global
        displacement vector, via T_matrix (identity for every class in
        this codebase currently, but kept general)."""
        idx = self.global_dof_indices()
        u_elem = u_global[idx]
        T = self.T_matrix()
        return T @ u_elem

    def recover_output(self, u_global: np.ndarray, points: Optional[dict] = None) -> dict:
        """Recover membrane stress, bending (fiber) stress, membrane
        force, bending moment, and shear force at center + 4 corners
        (NASTRAN CQUAD4/CQUADR "CORNER" output convention), both in this
        element's own local frame AND rotated to global.

        NOTE on the constitutive matrices used here: every element class
        in this codebase builds C_mb = E/(1-nu^2)*[[1,nu,0],...] WITHOUT
        a thickness factor -- the h_avg factor lives in k_local's volume
        weight (dV_m = detJ*w*h_avg), not in C_mb. That means C_mb @
        strain is already STRESS (matches this codebase's own patch-test
        reference values directly, e.g. sxx=1333.33 for the SAP2000
        2-001 case), not a force resultant. An earlier version of this
        method conflated the two (labeled C_mb @ strain as "Nxx" and used
        C_mb*h^2/12 for "moment", which is dimensionally wrong -- true
        bending rigidity is D_b = C_mb * h^3/12). Fixed here: STRESS and
        FORCE/MOMENT resultants are both returned, computed separately
        and correctly.

        Returns:
            { point_name: {
                "stress_local":  {"Sxx","Syy","Sxy","Mx_top","My_top",
                                   "Mxy_top" (bending fiber stress at
                                   z=+h/2), "Qx","Qy" (shear stress)},
                "stress_global": {same keys, rotated to global},
                "force_local":   {"Nxx","Nyy","Nxy" (membrane force/length),
                                   "Mxx","Myy","Mxy" (moment/length),
                                   "Qx","Qy" (shear force/length)},
                "force_global":  {same keys, rotated to global},
              }, ... }
        """
        if points is None:
            points = self.CORNER_POINTS

        u = self._element_u(u_global)
        C_mb, _, C_s, h_avg = self._isotropic_constitutive()
        D_b = C_mb * (h_avg ** 3 / 12.0)   # true bending rigidity

        out = {}
        for name, (r1, r2) in points.items():
            Bm = self._compute_Bm(r1, r2)
            Bb = self._compute_Bb(r1, r2)
            Bs = self._shear_B(r1, r2)

            eps_m = Bm @ u
            kappa = Bb @ u
            gamma = Bs @ u if Bs is not None else np.zeros(2)

            # --- stress (membrane stress + bending fiber stress at z=+h/2) ---
            S_mem = C_mb @ eps_m
            S_bend = C_mb @ kappa * (h_avg / 2.0)
            Q_stress = C_s @ gamma

            # --- force / moment resultants ---
            N_force = C_mb @ eps_m * h_avg
            M_force = D_b @ kappa
            Q_force = C_s @ gamma * h_avg

            e1, e2, e3 = self.bm_frame(r1, r2)

            def _rotate_shear_to_global(qx, qy, e1, e2, e3):
                """Shear is a rank-1 (vector) quantity in the tangent
                plane, not a rank-2 tensor -- rotate it as a 3-vector via
                Q @ [qx, qy, 0], Q = [e1 e2 e3] (columns), then return the
                full 3D global vector. The through-thickness (e3)
                component of a properly-defined transverse shear is zero
                by construction (Qx,Qy live entirely in the tangent
                plane), so this recovers the correct global Qx,Qy,Qz
                exactly, unlike the previous placeholder which just
                passed the local numbers through unchanged and was wrong
                for any warped or non-axis-aligned element."""
                Q = np.column_stack([e1, e2, e3])
                return Q @ np.array([qx, qy, 0.0])

            def _pack_stress(mem, bend, shear):
                g_mem = _rotate_inplane_tensor_to_global(mem[0], mem[1], mem[2], e1, e2)
                g_bend = _rotate_inplane_tensor_to_global(bend[0], bend[1], bend[2], e1, e2)
                g_shear = _rotate_shear_to_global(shear[0], shear[1], e1, e2, e3)
                local = {
                    "Sxx": mem[0], "Syy": mem[1], "Sxy": mem[2],
                    "Sxx_bend_top": bend[0], "Syy_bend_top": bend[1], "Sxy_bend_top": bend[2],
                    "Qx": shear[0], "Qy": shear[1],
                }
                glob = {
                    "Sxx": g_mem[0, 0], "Syy": g_mem[1, 1], "Sxy": g_mem[0, 1],
                    "Sxx_bend_top": g_bend[0, 0], "Syy_bend_top": g_bend[1, 1], "Sxy_bend_top": g_bend[0, 1],
                    "Qx": g_shear[0], "Qy": g_shear[1], "Qz": g_shear[2],
                }
                return local, glob

            def _pack_force(mem, moment, shear):
                g_mem = _rotate_inplane_tensor_to_global(mem[0], mem[1], mem[2], e1, e2)
                g_moment = _rotate_inplane_tensor_to_global(moment[0], moment[1], moment[2], e1, e2)
                g_shear = _rotate_shear_to_global(shear[0], shear[1], e1, e2, e3)
                local = {
                    "Nxx": mem[0], "Nyy": mem[1], "Nxy": mem[2],
                    "Mxx": moment[0], "Myy": moment[1], "Mxy": moment[2],
                    "Qx": shear[0], "Qy": shear[1],
                }
                glob = {
                    "Nxx": g_mem[0, 0], "Nyy": g_mem[1, 1], "Nxy": g_mem[0, 1],
                    "Mxx": g_moment[0, 0], "Myy": g_moment[1, 1], "Mxy": g_moment[0, 1],
                    "Qx": g_shear[0], "Qy": g_shear[1], "Qz": g_shear[2],
                }
                return local, glob

            stress_local, stress_global = _pack_stress(S_mem, S_bend, Q_stress)
            force_local, force_global = _pack_force(N_force, M_force, Q_force)

            out[name] = {
                "stress_local": stress_local, "stress_global": stress_global,
                "force_local": force_local, "force_global": force_global,
            }
        return out


def _rotate_inplane_tensor_to_global(sxx: float, syy: float, sxy: float,
                                      e1: np.ndarray, e2: np.ndarray) -> np.ndarray:
    """Embed a 2D in-plane symmetric tensor (given in local frame e1,e2,
    plane-stress convention: zero through-thickness components) as a full
    3x3 GLOBAL tensor, via the standard orthonormal-basis-change formula
    Global = Q @ Local @ Q.T, where Q's columns are the local basis
    vectors expressed in global coordinates. Q is completed with e3 =
    e1 x e2 (unused rows/columns of Local are zero, so e3's exact value
    only needs to make Q orthonormal, not carry any physical meaning)."""
    e1 = np.asarray(e1, dtype=float)
    e2 = np.asarray(e2, dtype=float)
    e3 = np.cross(e1, e2)
    Q = np.column_stack([e1, e2, e3])
    local = np.array([
        [sxx, sxy, 0.0],
        [sxy, syy, 0.0],
        [0.0, 0.0, 0.0],
    ])
    return Q @ local @ Q.T


# ═══════════════════════════════════════════════════════════════
#  MODEL (CONTAINER UTAMA)
# ═══════════════════════════════════════════════════════════════

class Model:
    """
    FEM Model dengan AutoSPC.

    Parameters
    ----------
    ndof_per_node : int
        DOF per node (default 6 untuk shell 3D umum, pakai 5 untuk DKMQ20).
    autospc : bool
        Aktifkan AutoSPC (default True). Mengunci DOF dengan stiffness nol
        otomatis sebelum solve.
    autospc_tol : float
        Threshold relatif untuk deteksi zero-stiffness DOF.
        DOF dengan K[i,i] < autospc_tol * max(K_diag) dianggap zero.
        Default: 1e-8.
    autospc_verbose : bool
        Tampilkan informasi DOF yang dikunci AutoSPC (default False).
    """

    def __init__(self,
                 ndof_per_node : int   = 6,
                 autospc        : bool  = True,
                 autospc_tol    : float = 1e-8,
                 autospc_verbose: bool  = False):

        self.ndof_per_node  = ndof_per_node
        self.autospc        = autospc
        self.autospc_tol    = float(autospc_tol)
        self.autospc_verbose= autospc_verbose

        self.nodes    : dict         = {}
        self.elements : list         = []
        self.bcs      : list         = []
        self.loads    : list         = []
        self._nodal_masses: dict     = {}

        self._built  = False
        self.ndof    = 0
        self.K       : Optional[np.ndarray] = None
        self.M       : Optional[np.ndarray] = None
        self.F       : Optional[np.ndarray] = None
        self.u       : Optional[np.ndarray] = None

        # DOF yang dikunci oleh AutoSPC (dicatat setelah build)
        self._autospc_dofs : list = []

    # ── Node / Element management ──────────────────────────────

    def add_node(self, nid: int, x: float, y: float, z: float) -> Node:
        node = Node(nid, x, y, z)
        self.nodes[nid] = node
        self._built = False
        return node

    def add_element(self, elem: Element) -> None:
        self.elements.append(elem)
        self._built = False

    def add_bc(self, node_id: int, dof_local: int, value: float = 0.0) -> None:
        self.bcs.append(_PendingBC(node_id, dof_local, value))

    def add_load(self, node_id: int, dof_local: int, value: float) -> None:
        self.loads.append(_PendingLoad(node_id, dof_local, value))

    def fix_node(self, node_id: int, dofs: list = None) -> None:
        if dofs is None:
            dofs = list(range(self.ndof_per_node))
        for d in dofs:
            self.add_bc(node_id, d, 0.0)

    def add_nodal_mass(self, node_id: int, masses: list) -> None:
        if node_id not in self._nodal_masses:
            self._nodal_masses[node_id] = np.zeros(self.ndof_per_node)
        self._nodal_masses[node_id] += np.array(masses)

    # ── Build (assembly) ─────────────────────────────────────

    def build(self) -> None:
        # 1. Assign DOF indices
        counter = 0
        for nid in sorted(self.nodes.keys()):
            node = self.nodes[nid]
            node.dofs = list(range(counter, counter + self.ndof_per_node))
            counter  += self.ndof_per_node
        self.ndof = counter

        # 2. Resolve BCs and loads
        self._resolved_bcs   = []
        self._resolved_loads = []
        for pbc in self.bcs:
            g = self.nodes[pbc.node_id].dofs[pbc.dof_local]
            self._resolved_bcs.append(BC(g, pbc.value))
        for pl in self.loads:
            g = self.nodes[pl.node_id].dofs[pl.dof_local]
            self._resolved_loads.append(NodalLoad(g, pl.value))

        # 3. Assemble K and M
        self.K = np.zeros((self.ndof, self.ndof))
        self.M = np.zeros((self.ndof, self.ndof))

        for elem in self.elements:
            Ke  = elem.k_global()
            Me  = elem.m_global()
            idx = elem.global_dof_indices()
            for i, gi in enumerate(idx):
                for j, gj in enumerate(idx):
                    self.K[gi, gj] += Ke[i, j]
                    if Me is not None:
                        self.M[gi, gj] += Me[i, j]

        # Nodal masses
        for nid, m_vals in self._nodal_masses.items():
            node = self.nodes[nid]
            for i, m in enumerate(m_vals):
                self.M[node.dofs[i], node.dofs[i]] += m

        # 4. Load vector
        self.F = np.zeros(self.ndof)
        for load in self._resolved_loads:
            self.F[load.dof] += load.value

        # 5. AutoSPC: detect zero-stiffness DOFs
        self._autospc_dofs = []
        if self.autospc:
            self._run_autospc()

        self._built = True

    def _run_autospc(self) -> None:
        """
        Deteksi DOF dengan stiffness diagonal nol (atau sangat kecil)
        yang belum di-constrain oleh user.

        Penyebab umum di model shell 3D:
          • DOF rotasi drilling (θ_z) pada elemen flat — tidak ada rigiditas
          • DOF pada node yang tidak terhubung ke elemen apapun
          • DOF in-plane yang tidak terkonstrain pada pelat datar (translasi
            dalam bidang tanpa beban/constraint in-plane)

        AutoSPC menambahkan constraint nol (BC=0) pada DOF ini.

        FIX (found via prob_2_005 thickness sweep -- see diagnosis notes):
        the threshold used to be a SINGLE value derived from the largest
        diagonal entry in the WHOLE matrix (K_max = max(diag(K))). For a
        shell element, membrane stiffness scales as ~h while bending
        stiffness scales as ~h^3 (e.g. DKMQ24_ShellElement_RHR.k_local's
        dV_b = detJ*w*h^3/12 vs dV_m = detJ*w*h). Comparing a bending DOF's
        diagonal against a GLOBAL max dominated by membrane stiffness means
        bending_diag/K_max ~ h^2 -- which crosses ANY fixed relative
        threshold once the plate is thin enough, regardless of element
        formulation. This silently auto-constrained genuinely valid
        (just naturally small) bending rotation DOFs to zero at h=1e-4 in
        prob_2_005, with NO warning since autospc_verbose defaults to
        False -- the model was quietly wrong with no error raised.

        Fix: compare each DOF only against the max diagonal AMONG ITS OWN
        DOF TYPE (i.e. same local index within a node: 0=Ux,1=Uy,2=Uz,
        3=Rx,4=Ry,5=Rz for a 6-dof node), not against the whole matrix.
        This keeps AutoSPC's original purpose (catching genuinely
        unconnected/singular DOFs of a given type) while no longer
        conflating physically different stiffness families that have
        different power-law scaling with thickness.
        """
        if self.ndof == 0:
            return

        # DOF yang sudah di-constrain user
        user_fixed = {bc.dof for bc in self._resolved_bcs}

        K_diag = np.abs(np.diag(self.K))
        n = self.ndof_per_node

        # Per-DOF-type max, instead of one global max. Assumes uniform
        # ndof_per_node and contiguous per-node DOF blocks (true given how
        # build() assigns node.dofs), so global DOF index d has local type
        # d % n.
        group_max = np.zeros(n)
        for loc in range(n):
            vals = K_diag[loc::n]
            group_max[loc] = vals.max() if vals.size else 0.0

        if group_max.max() < 1e-30:
            return   # Model kosong / semua nol

        spc_dofs = []
        for dof in range(self.ndof):
            if dof in user_fixed:
                continue
            loc = dof % n
            gmax = group_max[loc]
            if gmax < 1e-30:
                # An entire DOF-type with ZERO stiffness at every node in
                # the whole model (e.g. an element whose drilling rotation
                # genuinely carries no stiffness anywhere, by design, and
                # relies entirely on AutoSPC -- MacNealQ8_1992_native is
                # exactly this case) is the case that MOST needs locking,
                # not one to skip. The original `continue` here silently
                # left every DOF of that type completely unconstrained,
                # producing a singular K (LinAlgError, or NaN wherever the
                # caller catches that exception and reports it that way)
                # for any element relying on AutoSPC rather than an
                # artificial drilling penalty. Lock the whole group here
                # instead of skipping it.
                spc_dofs.append(dof)
                continue
            if K_diag[dof] < self.autospc_tol * gmax:
                spc_dofs.append(dof)

        self._autospc_dofs = spc_dofs

        if spc_dofs:
            # Tambahkan ke resolved BCs
            for dof in spc_dofs:
                self._resolved_bcs.append(BC(dof, 0.0))

            # Always print a one-line warning, even if autospc_verbose is
            # False -- silent auto-constraining is exactly what caused the
            # prob_2_005 bug to go undetected. Full per-DOF detail still
            # gated behind autospc_verbose.
            print(f"[AutoSPC] WARNING: locked {len(spc_dofs)} zero-stiffness "
                  f"DOF(s) (autospc_tol={self.autospc_tol:.1e}, per-DOF-type "
                  f"threshold). Set autospc_verbose=True for details, or "
                  f"autospc=False to disable and see the real singular-"
                  f"matrix error instead.")

            if self.autospc_verbose:
                print(f"[AutoSPC] Mengunci {len(spc_dofs)} DOF zero-stiffness:")
                for dof in spc_dofs:
                    # Identifikasi node dan DOF lokal
                    for nid in sorted(self.nodes.keys()):
                        nd = self.nodes[nid]
                        if dof in nd.dofs:
                            loc = nd.dofs.index(dof)
                            dof_names = ['U','V','W','θ₁','θ₂','θ_z']
                            name = dof_names[loc] if loc < len(dof_names) else f'DOF{loc}'
                            print(f"  Node {nid} ({name}) — K_diag={np.diag(self.K)[dof]:.3e}")
                            break

    @property
    def autospc_count(self) -> int:
        """Jumlah DOF yang dikunci oleh AutoSPC."""
        return len(self._autospc_dofs)

    # ── Solve Static ─────────────────────────────────────────

    def solve_static(self) -> np.ndarray:
        """
        Solve K·u = F (linear static).
        AutoSPC diterapkan otomatis di build() sebelum solve.
        Returns u (ndof,) — displacement vector.
        """
        if not self._built:
            self.build()

        K = self.K.copy()
        F = self.F.copy()

        fixed_dofs = [bc.dof   for bc in self._resolved_bcs]
        fixed_vals = [bc.value for bc in self._resolved_bcs]
        free_dofs  = [d for d in range(self.ndof) if d not in set(fixed_dofs)]

        # Modifikasi RHS untuk prescribed displacements
        for dof, val in zip(fixed_dofs, fixed_vals):
            if val != 0.0:
                F[free_dofs] -= K[np.ix_(free_dofs, [dof])].flatten() * val

        K_ff = K[np.ix_(free_dofs, free_dofs)]
        F_f  = F[free_dofs]

        u = np.zeros(self.ndof)
        for dof, val in zip(fixed_dofs, fixed_vals):
            u[dof] = val

        try:
            u[free_dofs] = np.linalg.solve(K_ff, F_f)
        except np.linalg.LinAlgError as e:
            # Tambahkan info diagnostik
            cond = np.linalg.cond(K_ff)
            raise np.linalg.LinAlgError(
                f"K singular (kondisi={cond:.2e}). "
                f"Cek BC, autospc={self.autospc}, "
                f"autospc_dofs={len(self._autospc_dofs)}. "
                f"Original: {e}"
            ) from e

        self.u = u
        return u

    # ── Solve Eigen ──────────────────────────────────────────

    def solve_eigen(self, n_modes: int = 6) -> tuple:
        """
        Eigenvalue analysis: K·φ = λ·M·φ
        Returns (freqs_Hz, mode_shapes).
        """
        from scipy.linalg import eigh

        if not self._built:
            self.build()

        fixed_dofs = [bc.dof for bc in self._resolved_bcs]
        free_dofs  = [d for d in range(self.ndof) if d not in set(fixed_dofs)]

        K_ff = self.K[np.ix_(free_dofs, free_dofs)]
        M_ff = self.M[np.ix_(free_dofs, free_dofs)]

        # Regularisasi M untuk stabilitas numerik
        M_ff += np.eye(len(free_dofs)) * 1e-12

        n_modes = min(n_modes, len(free_dofs))
        eigenvalues, eigenvectors = eigh(
            K_ff, M_ff, subset_by_index=[0, n_modes-1]
        )

        eigenvalues = np.maximum(eigenvalues, 0.0)
        omega = np.sqrt(eigenvalues)
        freqs = omega / (2 * np.pi)

        modes_global = np.zeros((self.ndof, n_modes))
        for i in range(n_modes):
            modes_global[free_dofs, i] = eigenvectors[:, i]

        self.freqs = freqs
        self.modes = modes_global
        return freqs, modes_global

    # ── Utilities ────────────────────────────────────────────

    def get_displacement(self, node_id: int) -> np.ndarray:
        """Displacement vector untuk node_id (panjang ndof_per_node)."""
        node = self.nodes[node_id]
        return self.u[node.dofs[:self.ndof_per_node]]

    def get_reaction(self, node_id: int) -> np.ndarray:
        """
        Reaction force di node (untuk DOF yang di-constrain).
        R = K·u - F  di DOF terkonstrain.
        """
        node = self.nodes[node_id]
        dofs = node.dofs[:self.ndof_per_node]
        R = np.zeros(self.ndof_per_node)
        for i, d in enumerate(dofs):
            R[i] = (self.K[d, :] @ self.u) - self.F[d]
        return R

    def summary(self) -> str:
        """Ringkasan model."""
        n_elem  = len(self.elements)
        n_node  = len(self.nodes)
        n_bc    = len([b for b in self._resolved_bcs
                       if b.dof not in self._autospc_dofs]) if self._built else len(self.bcs)
        n_spc   = len(self._autospc_dofs)
        n_load  = len(self.loads)
        lines   = [
            f"Model Summary",
            f"  Nodes    : {n_node}",
            f"  Elements : {n_elem}",
            f"  DOF total: {self.ndof}",
            f"  BC (user): {n_bc}",
            f"  AutoSPC  : {n_spc} DOF",
            f"  Loads    : {n_load}",
        ]
        return "\n".join(lines)
