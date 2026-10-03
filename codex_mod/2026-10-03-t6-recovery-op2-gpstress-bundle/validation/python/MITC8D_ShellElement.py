import numpy as np
from MITC8_ShellElement_v1p2 import MITC8_ShellElement_v1p2

class MITC8D_ShellElement(MITC8_ShellElement_v1p2):
    """
    MITC8 with pure Wilson drilling penalty (only theta_z is penalized).
    """




    def __init__(self, eid, nodes, E, nu, h,
                 kshear=5.0/6.0, beta_drill=1e-4):
        # beta_drill is now the coefficient for the pure penalty
        super().__init__(eid, nodes, E, nu, h,
                         kshear=kshear, beta_drill=beta_drill)
        self.nodal_normals = None  # opsional

    def _compute_Bdrill(self, xi, eta):
        N, _ = self._shape_q8(xi, eta)
        Bd = np.zeros((1, 48))
        for i in range(8):
            col = 6*i
            # Only theta_z (index col+5) is active
            Bd[0, col+5] = N[i]
        return Bd

    def k_local(self):
        # You can either:
        # 1. Call the parent's k_local and then adjust the penalty, or
        # 2. Reimplement the whole k_local with the new Bdrill.
        # Option 1 is simpler but parent's k_local already includes the old penalty.
        # So better to reimplement k_local without the old Bdrill.
        # Here we copy the parent's k_local but replace the drilling term.
        h  = self.h
        Cm = self._C_membrane()
        Cb = Cm
        Cs = self._C_shear()
        alpha_drill = self.beta_drill * self.G   # still uses beta_drill

        K = np.zeros((48, 48))
        gp = np.array([-np.sqrt(3.0/5.0), 0.0, np.sqrt(3.0/5.0)])
        gw = np.array([5.0/9.0, 8.0/9.0, 5.0/9.0])

        for i, xi in enumerate(gp):
            for j, eta in enumerate(gp):
                J   = self._jacobian(xi, eta)
                detJ = np.linalg.det(J)
                if detJ < 1e-14:
                    continue
                w = gw[i] * gw[j] * detJ

                Bm = self._compute_Bm(xi, eta)
                Bb = self._compute_Bb_direct(xi, eta)
                Bs = self._compute_Bs(xi, eta)
                Bd = self._compute_Bdrill(xi, eta)   # now the pure one

                K += h     * w * (Bm.T @ Cm @ Bm)
                K += (h**3/12.0) * w * (Bb.T @ Cb @ Bb)
                K += h     * w * (Bs.T @ Cs @ Bs)
                K += h     * w * alpha_drill * (Bd.T @ Bd)

        return K

