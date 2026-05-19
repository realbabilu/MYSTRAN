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
        """
        if self.ndof == 0:
            return

        # DOF yang sudah di-constrain user
        user_fixed = {bc.dof for bc in self._resolved_bcs}

        # Deteksi: diagonal K < tol * max_diagonal
        K_diag = np.abs(np.diag(self.K))
        K_max  = K_diag.max()

        if K_max < 1e-30:
            return   # Model kosong / semua nol

        threshold = self.autospc_tol * K_max
        spc_dofs  = []

        for dof in range(self.ndof):
            if dof in user_fixed:
                continue
            if K_diag[dof] < threshold:
                spc_dofs.append(dof)

        self._autospc_dofs = spc_dofs

        if spc_dofs:
            # Tambahkan ke resolved BCs
            for dof in spc_dofs:
                self._resolved_bcs.append(BC(dof, 0.0))

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
