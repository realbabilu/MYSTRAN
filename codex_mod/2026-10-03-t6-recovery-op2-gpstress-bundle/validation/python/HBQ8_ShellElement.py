"""
HBQ8_ShellElement.py  v2 — BUG FIX
=====================================
Fix: hapus invalid drilling-membrane coupling dari _compute_Bm().
Sebelumnya Bm punya terms Bm[0/1/2, col+5] yang mencoupling εxx/εyy/γxy
dengan DOF drilling θz — ini menyebabkan:
  • rank k_local = 41 (harusnya 42) → 1 spurious zero mode extra
  • respon ~200M% salah di problem 2-002 (LC1/LC2/LC5)
Fix: drilling hanya lewat Bdrill penalty (sudah benar di k_local),
     Bm hanya berisi membran murni tanpa coupling ke θz.

Ref: Darılmaz & Kumbasar, Computers & Structures 84 (2006) 1990-2000.
"""

import numpy as np
from core import Element


class HBQ8_ShellElement(Element):
    ndof_per_node = 6

    def __init__(self, eid, nodes, E, nu, h,
                 kshear=5./6., beta_drill=1e-3):
        self.eid        = eid
        self.nodes      = nodes
        self.E          = float(E)
        self.nu         = float(nu)
        self.h          = float(h)
        self.kshear     = float(kshear)
        self.beta_drill = float(beta_drill)
        self.G          = E / (2.*(1.+nu))
        self._build_local_frame()

    def _build_local_frame(self):
        pts = np.array([n.coords for n in self.nodes])
        c   = pts[:4]
        e1  = 0.5*((c[1]-c[0]) + (c[2]-c[3]))
        e1 /= max(np.linalg.norm(e1), 1e-14)
        e3  = np.cross(c[2]-c[0], c[3]-c[1])
        e3 /= max(np.linalg.norm(e3), 1e-14)
        e2  = np.cross(e3, e1); e2 /= max(np.linalg.norm(e2), 1e-14)
        e1  = np.cross(e2, e3); e1 /= max(np.linalg.norm(e1), 1e-14)
        self.ex, self.ey, self.ez = e1, e2, e3
        origin = np.mean(pts[:4], axis=0)
        self.origin = origin
        self.xy_local = np.array([
            [np.dot(n.coords-origin, e1),
             np.dot(n.coords-origin, e2)]
            for n in self.nodes
        ])

    @staticmethod
    def _shape_q8(xi, eta):
        N=np.zeros(8); dN=np.zeros((2,8))
        xc=np.array([-1.,1.,1.,-1.]); ec=np.array([-1.,-1.,1.,1.])
        for i in range(4):
            xi_i,eta_i=xc[i],ec[i]
            N[i]=.25*(1+xi*xi_i)*(1+eta*eta_i)*(xi*xi_i+eta*eta_i-1)
            dN[0,i]=.25*xi_i*(1+eta*eta_i)*(2*xi*xi_i+eta*eta_i)
            dN[1,i]=.25*eta_i*(1+xi*xi_i)*(xi*xi_i+2*eta*eta_i)
        N[4]=.5*(1-xi**2)*(1-eta);   dN[0,4]=-xi*(1-eta);     dN[1,4]=-.5*(1-xi**2)
        N[5]=.5*(1+xi)*(1-eta**2);   dN[0,5]=.5*(1-eta**2);   dN[1,5]=-eta*(1+xi)
        N[6]=.5*(1-xi**2)*(1+eta);   dN[0,6]=-xi*(1+eta);     dN[1,6]=.5*(1-xi**2)
        N[7]=.5*(1-xi)*(1-eta**2);   dN[0,7]=-.5*(1-eta**2);  dN[1,7]=-eta*(1-xi)
        return N, dN

    def _jacobian(self, xi, eta):
        _, dN = self._shape_q8(xi, eta)
        return dN @ self.xy_local   # (2,2)

    def T_matrix(self):
        R3 = np.column_stack([self.ex, self.ey, self.ez])
        T  = np.zeros((48,48))
        for i in range(8):
            for k in range(2):
                s=6*i+3*k; T[s:s+3,s:s+3]=R3.T
        return T

    def _C_membrane(self):
        f=self.E/(1.-self.nu**2)
        return f*np.array([[1.,self.nu,0.],[self.nu,1.,0.],[0.,0.,(1.-self.nu)/2.]])

    def _C_shear(self): return self.kshear*self.G*np.eye(2)

    def _compute_Bm(self, xi, eta):
        """
        Membran B-matrix (3×48): εxx, εyy, γxy — TANPA coupling drilling.
        Bug lama: Bm[0/1/2, col+5] ada coupling θz ke strain membran.
        Fix: Bm hanya dari DOF translasi u,v (col+0, col+1).
        """
        N, dN = self._shape_q8(xi, eta)
        J = self._jacobian(xi, eta)
        if abs(np.linalg.det(J))<1e-14: return np.zeros((3,48))
        Jinv = np.linalg.inv(J)
        dNxy = Jinv @ dN    # (2,8) → dN/dx, dN/dy
        Bm = np.zeros((3,48))
        for i in range(8):
            col=6*i
            Bm[0, col]   = dNxy[0,i]   # εxx = ∂u/∂x
            Bm[1, col+1] = dNxy[1,i]   # εyy = ∂v/∂y
            Bm[2, col]   = dNxy[1,i]   # γxy = ∂u/∂y + ∂v/∂x
            Bm[2, col+1] = dNxy[0,i]
        return Bm

    def _compute_Bb(self, xi, eta):
        """Bending B-matrix (3×48): κxx, κyy, κxy."""
        _, dN = self._shape_q8(xi, eta)
        J = self._jacobian(xi, eta)
        if abs(np.linalg.det(J))<1e-14: return np.zeros((3,48))
        Jinv = np.linalg.inv(J)
        dNxy = Jinv @ dN
        Bb = np.zeros((3,48))
        for i in range(8):
            col=6*i
            Bb[0, col+4] =  dNxy[0,i]   # κxx = ∂θy/∂x
            Bb[1, col+3] = -dNxy[1,i]   # κyy = -∂θx/∂y
            Bb[2, col+4] =  dNxy[1,i]   # κxy = ∂θy/∂y - ∂θx/∂x
            Bb[2, col+3] = -dNxy[0,i]
        return Bb

    def _compute_Bs(self, xi, eta):
        """Transverse shear B-matrix (2×48): γxz=∂w/∂x+θy, γyz=∂w/∂y-θx."""
        N, dN = self._shape_q8(xi, eta)
        J = self._jacobian(xi, eta)
        if abs(np.linalg.det(J))<1e-14: return np.zeros((2,48))
        Jinv = np.linalg.inv(J)
        dNxy = Jinv @ dN
        Bs = np.zeros((2,48))
        for i in range(8):
            col=6*i
            Bs[0, col+2] = dNxy[0,i]; Bs[0, col+4] =  N[i]
            Bs[1, col+2] = dNxy[1,i]; Bs[1, col+3] = -N[i]
        return Bs

    def _compute_Bdrill(self, xi, eta):
        """Drilling penalty (1×48): θz - ½(∂v/∂x - ∂u/∂y) = 0."""
        N, dN = self._shape_q8(xi, eta)
        J = self._jacobian(xi, eta)
        if abs(np.linalg.det(J))<1e-14: return np.zeros((1,48))
        Jinv = np.linalg.inv(J)
        dNxy = Jinv @ dN
        Bd = np.zeros((1,48))
        for i in range(8):
            col=6*i
            Bd[0,col]   = -0.5*dNxy[1,i]   # -½∂N/∂y × u
            Bd[0,col+1] =  0.5*dNxy[0,i]   # +½∂N/∂x × v
            Bd[0,col+5] =  N[i]             #  N × θz
        return Bd

    def k_local(self):
        h=self.h; Cm=self._C_membrane(); Cs=self._C_shear()
        ad=self.beta_drill*self.G
        K=np.zeros((48,48))
        gp=np.array([-np.sqrt(3./5.),0.,np.sqrt(3./5.)])
        gw=np.array([5./9.,8./9.,5./9.])
        for i,(xi,wi) in enumerate(zip(gp,gw)):
            for j,(eta,wj) in enumerate(zip(gp,gw)):
                J=self._jacobian(xi,eta); detJ=np.linalg.det(J)
                if detJ<1e-14: continue
                w=wi*wj*detJ
                Bm=self._compute_Bm(xi,eta)
                Bb=self._compute_Bb(xi,eta)
                Bs=self._compute_Bs(xi,eta)
                Bd=self._compute_Bdrill(xi,eta)
                K += h          * w * (Bm.T @ Cm @ Bm)
                K += (h**3/12.) * w * (Bb.T @ Cm @ Bb)
                K += h          * w * (Bs.T @ Cs @ Bs)
                K += h          * w * ad * (Bd.T @ Bd)
        return K

    def m_local(self):
        rho=getattr(self,'rho',0.); 
        if rho==0.: return None
        h=self.h; M=np.zeros((48,48))
        R=np.diag([rho*h]*3+[rho*h**3/12.,rho*h**3/12.,0.])
        gp=np.array([-np.sqrt(3./5.),0.,np.sqrt(3./5.)])
        gw=np.array([5./9.,8./9.,5./9.])
        for xi,wi in zip(gp,gw):
            for eta,wj in zip(gp,gw):
                N,_=self._shape_q8(xi,eta); J=self._jacobian(xi,eta)
                detJ=np.linalg.det(J)
                if detJ<1e-14: continue
                w=wi*wj*detJ
                Nm=np.zeros((6,48))
                for ii in range(8):
                    for d in range(6): Nm[d,6*ii+d]=N[ii]
                M+=w*(Nm.T@R@Nm)
        return M

    # ── _local_basis_at_point (required by patch_2_001_quadratic.py) ─
    def _local_basis_at_point(self, r, s):
        """
        Return local orthonormal basis (e1,e2,e3) at parametric point (r,s).
        For flat elements e1/e2/e3 are constant (precomputed in _build_local_frame).
        For curved elements, basis is computed from surface tangents at (r,s).
        """
        J = self._jacobian(r, s)   # (2,2) in local coords
        # Tangent vectors in 3D (local → global via ex,ey)
        g1_3d = J[0,0]*self.ex + J[0,1]*self.ey
        g2_3d = J[1,0]*self.ex + J[1,1]*self.ey
        norm1 = np.linalg.norm(g1_3d)
        if norm1 < 1e-14:
            return self.ex, self.ey, self.ez
        e1 = g1_3d / norm1
        e3 = np.cross(g1_3d, g2_3d)
        norm3 = np.linalg.norm(e3)
        if norm3 < 1e-14:
            return self.ex, self.ey, self.ez
        e3 = e3 / norm3
        e2 = np.cross(e3, e1)
        e2 = e2 / np.linalg.norm(e2)
        return e1, e2, e3
