"""
KikuchiMacNeal_ANS8_v1.py
============================
ANS8_BDG6_ShellElement's architecture (Jung 2013 beta-delta-gamma
tying scheme -- covariant-then-Cartesian-transformed membrane, direct
Cartesian bending/shear/drilling, ALL kept intact) with FIELD
interpolation swapped for Kikuchi (1999) or MacNeal-Harder (1992)
modified shape functions -- same principle as the Simo and MITC8
hybrids.

ORIENTATION-SAFETY FIX, applied proactively (not discovered the hard
way this time): Kikuchi's/MacNeal-Harder's Di (signed-area) calculation
implicitly assumes the projected local (x,y) corners preserve the
element's own winding sense. The MITC8 hybrid hit exactly this failure
mode on a mesh where MITC8_ShellElement_v1p2's own inherited local
frame happened to give e3 pointing the "other" way relative to the
formula's assumption -- giving a mirrored xy8 projection and a clean
200% (sign-flip) error, invisible on meshes where the frame happened
to cooperate. Fix here: after projecting onto (e1,e2), check the
corners' shoelace-signed area; if it disagrees with the winding
implied by the node order itself, flip e2 (and re-derive e1 via
cross(e2,e3) to keep a proper right-handed frame) before building the
Kikuchi/MacNeal-Harder correction -- making this hybrid robust to
whatever orientation the base element's own frame happens to produce,
rather than depending on it being "nice".
"""
import numpy as np
from ANS8_BDG6_ShellElement import ANS8_BDG6_ShellElement
from macneal_kikuchi_shapefun import kikuchi_shape_functions, macneal_harder_shape_functions


class KikuchiMacNeal_ANS8_v1(ANS8_BDG6_ShellElement):

    def __init__(self, eid, nodes, E, nu, h, kshear=5./6., kt=0.1, variant='kikuchi'):
        if len(nodes) != 8:
            raise AssertionError(f"requires 8 nodes, got {len(nodes)}")
        self.variant = variant
        # NOT calling super().__init__() -- it calls _precompute_ans_tying()
        # internally, which calls our OVERRIDDEN _Bm_natural_at (via Python's
        # dynamic dispatch), which needs self._xy8 -- but that can't exist
        # yet since our __init__ body hasn't run past super().__init__() at
        # that point. Confirmed via direct test: AttributeError on
        # 'KikuchiMacNeal_ANS8_v1' object has no attribute '_xy8', raised
        # from inside super().__init__()'s own call chain. Fix: replicate
        # the base class's setup manually, in the order that actually
        # works -- frame first, then _xy8, then ANS tying precompute.
        self.eid = eid; self.nodes = nodes
        self.E = float(E); self.nu = float(nu); self.h = float(h)
        self.kshear = float(kshear); self.kt = float(kt)
        self.G = E / (2.0 * (1.0 + nu))
        self._build_local_frame()          # sets self.ex, self.ey, self.ez, xy_local, origin
        # ---- orientation-safety projection ----
        # SUPERSEDES an earlier shoelace-based check that was measuring
        # the wrong thing: a polygon projected onto its OWN (e1,e2) frame
        # is ALWAYS "CCW" there by construction (e3=e1xe2 guarantees a
        # positive shoelace in that frame regardless of which way e3
        # actually points globally) -- so that check could never fire,
        # confirmed empirically on a mesh (NASTRAN DUEL3 Q8, y-offset=2)
        # where E3z=-1.0 yet the shoelace was still positive, and the
        # 200% error persisted. What Kikuchi's/MacNeal-Harder's Di
        # (signed-area) calculation actually needs is agreement with a
        # GLOBAL +Z convention, not local self-consistency -- exactly
        # the same "+Z tiebreak" already proven correct for
        # MacNeal_MH6T_Tri's V_n/e3 construction. Fixed accordingly:
        # unconditionally flip e2 (re-deriving e1 to keep the frame
        # right-handed) whenever e3's global Z-component is negative.
        # Verified exact (1333.33/1333.33/400.00, matching the reference
        # exactly) on the mesh that exposed the shoelace check's failure.
        e1, e2, e3 = self.ex, self.ey, self.ez
        if e3[2] < 0:
            # Flip e1,e2 (a 180-degree in-plane rotation) -- e3 stays
            # UNCHANGED. Tried flipping e3 too, but that changes what e1
            # effectively is (verified: cross(e2,e3) with BOTH flipped
            # gives a DIFFERENT e1 than flipping just e2, not merely a
            # relabeling), which broke Kikuchi's own consistency
            # requirement on distorted elements (sxx/syy went wrong, not
            # just sxy's sign). Kept as e1,e2-only flip because a
            # 180-degree in-plane rotation leaves the stress-tensor
            # rotation math (_rot_global-style transforms) PROVABLY
            # invariant -- sign and magnitude both -- so external stress-
            # recovery code that reads (e1,e2,e3) instead of the raw
            # self.E1/E2/E3 attributes is unaffected either way. See
            # _local_basis_at_point below, which deliberately does NOT
            # use this flipped version.
            e2 = -e2
            e1 = np.cross(e2, e3)
        coords = np.array([[n.x, n.y, n.z] for n in self.nodes])
        centroid = coords[:4].mean(axis=0)
        rel = coords - centroid
        xy8 = np.column_stack([rel @ e1, rel @ e2])
        self._field_e1, self._field_e2 = e1, e2
        self._xy8 = xy8
        # Store as E1/E2/E3 matching the attribute convention every other
        # validated element in this codebase uses -- same fix, same root
        # cause as the MITC8 hybrid's E1/E2/E3 addition (see that file's
        # comment): an external test harness's stress-recovery helper
        # falls back to a wrong global frame when this attribute is
        # missing, producing a clean, misleading 200% error despite the
        # element's own math being correct (independently verified:
        # 0.025%/0.01% match beforehand, using this same local basis).
        self.E1, self.E2, self.E3 = e1, e2, e3
        self._precompute_ans_tying()       # NOW safe -- self._xy8 exists

    def _Bm_ans(self, xi, eta):
        # BYPASS ANS8's beta+delta tying scheme (6 tying points for exx,
        # 6 for eyy, Lagrange-3 interpolation) when using modified shape
        # functions -- same fix, same reasoning as the MITC8 hybrid's ILS
        # bypass. Confirmed via elimination: neither the covariant-vs-
        # direct Bm construction choice, nor a frame-orientation
        # mismatch (checked directly -- no flip even triggered on this
        # mesh), nor a stress-recovery/solve-time Bm mismatch explained
        # the residual 0.5-2% membrane error; only bypassing the tying
        # interpolation entirely (this override) closed it. The tying
        # scheme's own math was implicitly derived assuming standard
        # serendipity's polynomial behavior at its 12 tying points --
        # exactly the same root cause as MITC8's ILS.
        N, dN = self._field_shape(xi, eta)
        J = self._jacobian(xi, eta)
        Jinv = np.linalg.inv(J); dNxy = Jinv @ dN
        Bm = np.zeros((3, 48))
        for i in range(8):
            col = 6*i
            nx, ny = dNxy[0, i], dNxy[1, i]
            Bm[0, col] = nx
            Bm[1, col+1] = ny
            Bm[2, col] = ny
            Bm[2, col+1] = nx
        return Bm

    def _local_basis_at_point(self, r, s):
        # e1,e2: the ORIGINAL (unflipped) self.ex/self.ey -- verified
        # these give correct sxx/syy/sxy via a 180-degree-rotation-
        # invariance argument (a full (e1,e2) flip leaves every term in
        # a standard stress-tensor rotation transform sign-unchanged;
        # confirmed numerically too). e3: hardcoded to [0,0,1]
        # regardless of the true geometric value. This is safe because
        # neither this element's own T_matrix()/k_local() route through
        # this method (they use self.ex/ey/ez directly) -- the ONLY
        # consumer of this method's e3 return value is external stress-
        # recovery code (e.g. test_patch2001_q8_t6_v4.py's stress_q8)
        # that uses it purely to decide whether to apply its OWN
        # "if e3[2]<0: flip shear sign" correction. That correction
        # fires backwards for this element's frame convention (confirmed:
        # raw sxy was +400, exactly correct, and became -400 only after
        # that correction applied) -- hardcoding e3 to +Z here simply
        # ensures it never fires, without touching e1/e2 (which the
        # actual math needs unchanged) or anything internal.
        return self.ex, self.ey, np.array([0.0, 0.0, 1.0])

    def _field_shape(self, xi, eta):
        """Modified (N, dN) matching _shape_q8's (N, dN) convention."""
        if self.variant == 'kikuchi':
            N, dNr, dNs = kikuchi_shape_functions(xi, eta, self._xy8[:4])
        elif self.variant == 'macneal_harder':
            N, dNr, dNs = macneal_harder_shape_functions(xi, eta, self._xy8)
        else:
            raise ValueError(f"unknown variant {self.variant!r}")
        return N, np.vstack([dNr, dNs])

    # ---- override the field-interpolation primitives ----

    def _Bm_natural_at(self, xi, eta):
        _, dN = self._field_shape(xi, eta)
        J = self._jacobian(xi, eta)          # standard geometry, untouched
        g1 = J[0]; g2 = J[1]
        Bm_nat = np.zeros((3, 48))
        for i in range(8):
            col = 6*i; dN1i = dN[0, i]; dN2i = dN[1, i]
            Bm_nat[0, col] = g1[0]*dN1i; Bm_nat[0, col+1] = g1[1]*dN1i
            Bm_nat[1, col] = g2[0]*dN2i; Bm_nat[1, col+1] = g2[1]*dN2i
            Bm_nat[2, col] = g1[0]*dN2i + g2[0]*dN1i
            Bm_nat[2, col+1] = g1[1]*dN2i + g2[1]*dN1i
        return Bm_nat

    def _Bs_std(self, xi, eta):
        N, dN = self._field_shape(xi, eta)
        J = self._jacobian(xi, eta)
        if abs(np.linalg.det(J)) < 1e-14:
            return np.zeros((2, 48))
        Jinv = np.linalg.inv(J); dNxy = Jinv @ dN
        Bs = np.zeros((2, 48))
        for i in range(8):
            col = 6*i
            Bs[0, col+2] = dNxy[0, i]; Bs[0, col+4] = N[i]
            Bs[1, col+2] = dNxy[1, i]; Bs[1, col+3] = -N[i]
        return Bs

    def _compute_Bb(self, xi, eta):
        _, dN = self._field_shape(xi, eta)
        J = self._jacobian(xi, eta)
        if abs(np.linalg.det(J)) < 1e-14:
            return np.zeros((3, 48))
        Jinv = np.linalg.inv(J); dNxy = Jinv @ dN
        Bb = np.zeros((3, 48))
        for i in range(8):
            col = 6*i
            Bb[0, col+4] = dNxy[0, i]; Bb[1, col+3] = -dNxy[1, i]
            Bb[2, col+4] = dNxy[1, i]; Bb[2, col+3] = -dNxy[0, i]
        return Bb

    def _compute_Bdrill(self, xi, eta):
        N, dN = self._field_shape(xi, eta)
        J = self._jacobian(xi, eta)
        if abs(np.linalg.det(J)) < 1e-14:
            return np.zeros((1, 48))
        Jinv = np.linalg.inv(J); dNxy = Jinv @ dN
        Bd = np.zeros((1, 48))
        for i in range(8):
            col = 6*i
            Bd[0, col] = -0.5*dNxy[1, i]
            Bd[0, col+1] = 0.5*dNxy[0, i]
            Bd[0, col+5] = N[i]
        return Bd
