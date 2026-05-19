"""
laminate_utils.py

Small Classical Lamination Theory (CLT) helper for composite shell elements.

Conventions
-----------
Engineering strain/resultant order follows the shell elements in this folder:

    membrane strain / curvature: [xx, yy, xy_engineering]
    transverse shear strain     : [xz, yz]

Laminate stiffnesses:

    [N]   [ A  B ] [eps0]
    [M] = [ B  D ] [kappa]

    [Qshear] = As [gamma_xz, gamma_yz]

Units are whatever your model uses, as long as material moduli and thickness
are consistent.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable, Sequence
import math
import numpy as np


def _sym(a: np.ndarray) -> np.ndarray:
    return 0.5 * (a + a.T)


def q_plane_stress(E1: float, E2: float, nu12: float, G12: float) -> np.ndarray:
    """Orthotropic reduced stiffness Q in material 1-2 axes, plane stress."""
    E1 = float(E1)
    E2 = float(E2)
    nu12 = float(nu12)
    G12 = float(G12)
    nu21 = nu12 * E2 / E1
    den = 1.0 - nu12 * nu21
    if den <= 0.0:
        raise ValueError("invalid orthotropic constants: 1-nu12*nu21 must be positive")
    return np.array(
        [
            [E1 / den, nu12 * E2 / den, 0.0],
            [nu12 * E2 / den, E2 / den, 0.0],
            [0.0, 0.0, G12],
        ],
        dtype=float,
    )


def qbar_plane_stress(E1: float, E2: float, nu12: float, G12: float, theta_deg: float) -> np.ndarray:
    """
    Transformed orthotropic reduced stiffness Qbar for engineering shear strain.

    Formula follows standard CLT [xx, yy, xy_engineering] convention.
    """
    Q = q_plane_stress(E1, E2, nu12, G12)
    Q11, Q22, Q12, Q66 = Q[0, 0], Q[1, 1], Q[0, 1], Q[2, 2]

    th = math.radians(float(theta_deg))
    m = math.cos(th)
    n = math.sin(th)
    m2, n2 = m * m, n * n
    m4, n4 = m2 * m2, n2 * n2

    q11 = Q11 * m4 + 2.0 * (Q12 + 2.0 * Q66) * m2 * n2 + Q22 * n4
    q22 = Q11 * n4 + 2.0 * (Q12 + 2.0 * Q66) * m2 * n2 + Q22 * m4
    q12 = (Q11 + Q22 - 4.0 * Q66) * m2 * n2 + Q12 * (m4 + n4)
    q16 = (Q11 - Q12 - 2.0 * Q66) * m * m2 * n - (Q22 - Q12 - 2.0 * Q66) * m * n * n2
    q26 = (Q11 - Q12 - 2.0 * Q66) * m * n * n2 - (Q22 - Q12 - 2.0 * Q66) * m * m2 * n
    q66 = (Q11 + Q22 - 2.0 * Q12 - 2.0 * Q66) * m2 * n2 + Q66 * (m4 + n4)

    return _sym(np.array([[q11, q12, q16], [q12, q22, q26], [q16, q26, q66]], dtype=float))


def shear_qbar(G13: float, G23: float, theta_deg: float) -> np.ndarray:
    """
    Approximate transformed transverse shear stiffness for [xz, yz].

    This is enough for a first FSDT/Mindlin composite implementation and
    reduces exactly to G*I for isotropic plies.
    """
    th = math.radians(float(theta_deg))
    c = math.cos(th)
    s = math.sin(th)
    R = np.array([[c, -s], [s, c]], dtype=float)
    Qs = np.diag([float(G13), float(G23)])
    return _sym(R @ Qs @ R.T)


@dataclass(frozen=True)
class Ply:
    E1: float
    E2: float
    nu12: float
    G12: float
    G13: float
    G23: float
    theta: float
    t: float
    name: str = ""


@dataclass
class Laminate:
    """
    Laminate stiffness container.

    A, B, D use [xx, yy, xy_engineering].
    As uses [xz, yz].
    """

    A: np.ndarray
    B: np.ndarray
    D: np.ndarray
    As: np.ndarray
    h: float
    plies: list[Ply] | None = None
    kappa: float = 5.0 / 6.0
    E_ref: float = 1.0
    nu_ref: float = 0.3

    def __post_init__(self):
        self.A = _sym(np.asarray(self.A, dtype=float).reshape(3, 3))
        self.B = _sym(np.asarray(self.B, dtype=float).reshape(3, 3))
        self.D = _sym(np.asarray(self.D, dtype=float).reshape(3, 3))
        self.As = _sym(np.asarray(self.As, dtype=float).reshape(2, 2))
        self.h = float(self.h)
        if self.h <= 0.0:
            raise ValueError("laminate thickness h must be positive")

    @property
    def ABD(self) -> np.ndarray:
        return np.block([[self.A, self.B], [self.B, self.D]])

    def is_symmetric(self, tol: float = 1e-10) -> bool:
        scale = max(1.0, float(np.max(np.abs(self.D))))
        return float(np.max(np.abs(self.B))) <= tol * scale / max(self.h, 1e-30)

    @classmethod
    def isotropic(cls, E: float, nu: float, h: float, kappa: float = 5.0 / 6.0) -> "Laminate":
        E = float(E)
        nu = float(nu)
        h = float(h)
        G = E / (2.0 * (1.0 + nu))
        Q = E / (1.0 - nu * nu) * np.array(
            [[1.0, nu, 0.0], [nu, 1.0, 0.0], [0.0, 0.0, 0.5 * (1.0 - nu)]],
            dtype=float,
        )
        A = Q * h
        B = np.zeros((3, 3), dtype=float)
        D = Q * h**3 / 12.0
        As = kappa * G * h * np.eye(2)
        return cls(A=A, B=B, D=D, As=As, h=h, plies=None, kappa=kappa, E_ref=E, nu_ref=nu)

    @classmethod
    def from_plies(cls, plies: Sequence[Ply | dict], kappa: float = 5.0 / 6.0) -> "Laminate":
        pp: list[Ply] = []
        for p in plies:
            if isinstance(p, Ply):
                pp.append(p)
            else:
                d = dict(p)
                # Friendly aliases
                if "angle" in d and "theta" not in d:
                    d["theta"] = d.pop("angle")
                if "thickness" in d and "t" not in d:
                    d["t"] = d.pop("thickness")
                # If G13/G23 missing, use G12 as first approximation.
                d.setdefault("G13", d.get("G12"))
                d.setdefault("G23", d.get("G12"))
                pp.append(Ply(**d))

        if not pp:
            raise ValueError("at least one ply is required")

        h = float(sum(p.t for p in pp))
        if h <= 0.0:
            raise ValueError("total laminate thickness must be positive")

        # z coordinates bottom -> top, mid-plane at z=0
        z = [-0.5 * h]
        for p in pp:
            z.append(z[-1] + float(p.t))

        A = np.zeros((3, 3), dtype=float)
        B = np.zeros((3, 3), dtype=float)
        D = np.zeros((3, 3), dtype=float)
        As = np.zeros((2, 2), dtype=float)

        for k, p in enumerate(pp):
            z0, z1 = z[k], z[k + 1]
            Qb = qbar_plane_stress(p.E1, p.E2, p.nu12, p.G12, p.theta)
            Qs = shear_qbar(p.G13, p.G23, p.theta)

            A += Qb * (z1 - z0)
            B += 0.5 * Qb * (z1**2 - z0**2)
            D += (1.0 / 3.0) * Qb * (z1**3 - z0**3)
            As += kappa * Qs * (z1 - z0)

        # Rough reference constants for fictitious drilling scale only.
        E_ref = float(sum(p.E1 * p.t for p in pp) / h)
        nu_ref = float(sum(p.nu12 * p.t for p in pp) / h)

        return cls(A=A, B=B, D=D, As=As, h=h, plies=pp, kappa=kappa, E_ref=E_ref, nu_ref=nu_ref)


def make_symmetric_angle_laminate(E1, E2, nu12, G12, G13, G23, h_ply, angles) -> Laminate:
    """
    Convenience builder for symmetric [angles / reversed(angles)] laminate.
    """
    half = [
        Ply(E1=E1, E2=E2, nu12=nu12, G12=G12, G13=G13, G23=G23, theta=a, t=h_ply)
        for a in angles
    ]
    plies = half + list(reversed(half))
    return Laminate.from_plies(plies)


__all__ = [
    "Ply",
    "Laminate",
    "q_plane_stress",
    "qbar_plane_stress",
    "shear_qbar",
    "make_symmetric_angle_laminate",
]
