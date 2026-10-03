"""
Modified 8-node shape functions (drop-in for standard Q8 serendipity's field
interpolation), implementing:

  - MacNeal & Harder (1992), "Eight nodes or nine?" -- general numeric
    construction via a ninth-node constraint solved from an 8x8 metric+
    parametric basis.
  - Kikuchi, Okabe & Fujio (1999) -- closed-form alternative, built from the
    9-node LAGRANGE shape functions (not the standalone serendipity ones),
    with nicer degenerate-to-6-node-triangle behavior.

IMPORTANT ARCHITECTURAL POINT (from both papers): the modification only
touches how *field variables* (u,v,w,theta) are interpolated from nodal
values. The element's GEOMETRY / Jacobian mapping stays on the standard
(bilinear-subparametric or standard serendipity) map. Do not use these
functions to interpolate x,y,z -- only to build Bm/Bb/Bs.

Node/parametric convention matches Simo1993_Q8_ShellElement_v1p8_standalone.py:
  corners: 1(-1,-1) 2(1,-1) 3(1,1) 4(-1,1)
  midsides: 5(0,-1)[1-2]  6(1,0)[2-3]  7(0,1)[3-4]  8(-1,0)[4-1]

Both return (N, dNr, dNs) exactly like Simo1993's _shape_functions(r,s).
"""
import numpy as np

_PARAM_RS = np.array([[-1.,-1.],[1.,-1.],[1.,1.],[-1.,1.],
                       [0.,-1.],[1.,0.],[0.,1.],[-1.,0.]])

def _N8_std(r, s):
    """Standalone 8-node serendipity (unchanged geometry mapping / baseline)."""
    N = np.array([0.25*(1-r)*(1-s)*(-1-r-s),
                  0.25*(1+r)*(1-s)*(-1+r-s),
                  0.25*(1+r)*(1+s)*(-1+r+s),
                  0.25*(1-r)*(1+s)*(-1-r+s),
                  0.5*(1-r*r)*(1-s),
                  0.5*(1+r)*(1-s*s),
                  0.5*(1-r*r)*(1+s),
                  0.5*(1-r)*(1-s*s)])
    dNr = np.array([0.25*(1-s)*(2*r+s), 0.25*(1-s)*(2*r-s),
                     0.25*(1+s)*(2*r+s), 0.25*(1+s)*(2*r-s),
                     -r*(1-s), 0.5*(1-s*s), -r*(1+s), -0.5*(1-s*s)])
    dNs = np.array([0.25*(1-r)*(r+2*s), 0.25*(1+r)*(-r+2*s),
                     0.25*(1+r)*(r+2*s), 0.25*(1-r)*(-r+2*s),
                     -0.5*(1-r*r), -s*(1+r), 0.5*(1-r*r), -s*(1-r)])
    return N, dNr, dNs

def _N9(r, s):
    return (1-r*r)*(1-s*s)

def _dN9(r, s):
    return -2*r*(1-s*s), -2*s*(1-r*r)

def _N8_lagrange(r, s):
    """The 8 boundary shape functions of the 9-node LAGRANGE element
    (biquadratic tensor product, node 9 omitted). These sum to (1 - N9),
    NOT to 1 on their own -- that's the point: Kikuchi/MacNeal reconstitute
    a valid 8-dof interpolant by folding N9's contribution back in via a
    constraint, instead of just dropping it (which is what produces plain
    serendipity)."""
    l_m1 = lambda x: 0.5*x*(x-1.0)   # Lagrange node at x=-1
    l_0  = lambda x: 1.0-x*x         # Lagrange node at x=0
    l_p1 = lambda x: 0.5*x*(x+1.0)   # Lagrange node at x=+1
    dl_m1 = lambda x: x - 0.5
    dl_0  = lambda x: -2.0*x
    dl_p1 = lambda x: x + 0.5

    Lr = {-1: l_m1(r), 0: l_0(r), 1: l_p1(r)}
    Ls = {-1: l_m1(s), 0: l_0(s), 1: l_p1(s)}
    dLr = {-1: dl_m1(r), 0: dl_0(r), 1: dl_p1(r)}
    dLs = {-1: dl_m1(s), 0: dl_0(s), 1: dl_p1(s)}

    N = np.empty(8); dNr = np.empty(8); dNs = np.empty(8)
    for k, (r0, s0) in enumerate(_PARAM_RS):
        r0i, s0i = int(r0), int(s0)
        N[k]   = Lr[r0i]  * Ls[s0i]
        dNr[k] = dLr[r0i] * Ls[s0i]
        dNs[k] = Lr[r0i]  * dLs[s0i]
    return N, dNr, dNs


def macneal_harder_shape_functions(r, s, xy8):
    """
    MacNeal & Harder (1992) modified shape functions -- eqs (9)-(15).
    xy8: (8,2) array -- LOCAL PLANAR (x,y) coords of the 8 nodes, i.e. nodes'
         positions projected onto the element's tangent plane (dot with
         E1,E2 as in Simo1993's local basis). Used ONLY to build the
         constraint T_i; geometry mapping is untouched elsewhere.
    Returns N, dNr, dNs (each length-8) for interpolating FIELD variables.

    Per the paper (eq. 16): the (x,y) origin used to build the constraint
    MUST be placed at the physical location of the eliminated 9th node,
    x9 = sum_i N_i^(8)(0,0) x_i, else eq. (14)'s "evaluate X_m at the
    origin" shortcut is invalid. We re-center xy8 internally so the caller
    doesn't have to worry about this.
    """
    xy8 = np.asarray(xy8, dtype=float)
    N8_0, _, _ = _N8_std(0.0, 0.0)          # -0.25 corners, 0.5 mids
    x9y9 = N8_0 @ xy8                        # eq. (16)
    xy8 = xy8 - x9y9                         # re-center so node-9 sits at (0,0)

    Xim = np.empty((8, 8))
    for i in range(8):
        x, y = xy8[i]
        xi, eta = _PARAM_RS[i]
        Xim[i, :] = [1.0, x, y, x*x, x*y, y*y, xi*xi*eta, xi*eta*eta]
    Ami = np.linalg.inv(Xim)         # [A_mi], eq. (13)
    T = Ami[0, :]                    # T_i = A_{1i}, eq. (15)

    N8_0, _, _ = _N8_std(0.0, 0.0)   # N_i^(8)(0,0): corners=-0.25, mids=0.5
    corr = T - N8_0

    N8, dN8r, dN8s = _N8_std(r, s)
    N9 = _N9(r, s)
    dN9r, dN9s = _dN9(r, s)

    N   = N8   + N9   * corr
    dNr = dN8r + dN9r * corr
    dNs = dN8s + dN9s * corr
    return N, dNr, dNs


def kikuchi_shape_functions(r, s, xy4_corners):
    """
    Kikuchi, Okabe & Fujio (1999) closed-form modified shape functions,
    eqs (35)-(37), built on top of the 9-node LAGRANGE basis (see
    _N8_lagrange), not standalone serendipity.
    xy4_corners: (4,2) array -- LOCAL PLANAR (x,y) coords of the 4 CORNER
                 nodes only (same projection convention as above).
    Returns N, dNr, dNs (length-8) for interpolating FIELD variables.
    """
    xy = np.asarray(xy4_corners, dtype=float)

    def D(i0):
        j0, m0 = (i0+1) % 4, (i0+3) % 4
        xi_, yi_ = xy[i0]; xj_, yj_ = xy[j0]; xm_, ym_ = xy[m0]
        return (xj_-xi_)*(ym_-yi_) - (xm_-xi_)*(yj_-yi_)

    Ds = np.array([D(i0) for i0 in range(4)])

    NL, dNLr, dNLs = _N8_lagrange(r, s)
    N9 = _N9(r, s)
    dN9r, dN9s = _dN9(r, s)

    coeff = np.zeros(8)   # ci for corners i=1..4, ci4 for midsides i+4
    for i0 in range(4):
        k0 = (i0+2) % 4
        m0 = (i0+3) % 4
        Di, Dk, Dm = Ds[i0], Ds[k0], Ds[m0]
        denom = (Di + Dk)
        ci  = -0.25 + (Di - Dk) / (8.0*denom)
        ci4 =  0.5  + (Dm - Di) / (4.0*denom)
        coeff[i0]   = ci
        coeff[i0+4] = ci4

    N   = NL    + N9   * coeff
    dNr = dNLr  + dN9r * coeff
    dNs = dNLs  + dN9s * coeff
    return N, dNr, dNs
