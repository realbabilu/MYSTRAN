"""
MacNealQ8_1992_native_v3_drill.py
====================================
Implementasi elemen MacNeal Quad8RM (1992) dengan modifikasi "Wilson Penalty" murni
pada derajat kebebasan drilling (theta_z).

Berbeda dengan Hughes-Brezzi yang mengikat theta_z pada regangan geser membran
(continuum spin), metode Wilson ini secara artifisial memberikan kekakuan fiktif
yang decoupled pada theta_z untuk menghindari singularitas (0 eigenvalue).
"""
import numpy as np
from MacNealQ8_1992_native_v2 import MacNealQ8_1992_native_v2


class MacNealQ8_1992_native_v3_drill(MacNealQ8_1992_native_v2):

    def __init__(self, eid, nodes, E, nu, h, rho=None, beta_drill=9e-6):
        """
        beta_drill: Konstanta penalty fiktif untuk drilling. 
        Disetel sangat kecil (1e-5) agar meminimalisir pengaruh fiktif pada hasil riil.
        """
        self.beta_drill = beta_drill
        super().__init__(eid, nodes, E, nu, h, rho=rho)

    def _compute_Bdrill(self, r, s):
        """
        B-Matrix khusus untuk drilling.
        Hanya mengaktifkan komponen N[k] pada DOF ke-6 (indeks 5, yaitu theta_z).
        """
        # Kita panggil _std_dxy dari base class hanya untuk mendapatkan Shape Function (N)
        N, _, _, _ = self._std_dxy(r, s)
        n = len(self.nodes)
        
        # Inisialisasi matriks (1 baris, n*6 kolom)
        Bdrill = np.zeros((1, n * 6))
        
        for k in range(n):
            c = 6 * k
            # HANYA penalize theta_z (DOF ke-6). 
            # Gradien translasi u dan v (seperti pada Hughes-Brezzi) diabaikan.
            Bdrill[0, c + 5] = N[k]
            
        return Bdrill

    def k_local(self):
        """
        Overriding fungsi k_local untuk menambahkan kekakuan drilling penalty
        ke dalam matriks kekakuan dasar elemen MacNeal Q8.
        """
        # Panggil matriks K dasar (Membrane + Bending + Transverse Shear) dari parent class
        K = super().k_local()
        
        # Persiapan integrasi Gauss 3x3 untuk penalty term
        gp3 = [-np.sqrt(3.0/5.0), 0.0, np.sqrt(3.0/5.0)]
        w3 = [5.0/9.0, 8.0/9.0, 5.0/9.0]
        h_avg = np.mean(self.h)
        
        # Koefisien penalty fiktif (beta * G)
        coeff = self.beta_drill * self.G
        
        # Loop integrasi untuk kekakuan drilling
        for r, wr in zip(gp3, w3):
            for s, ws in zip(gp3, w3):
                # Ambil determinan Jacobian dari base class
                _, _, _, detJ = self._std_dxy(r, s)
                w = wr * ws
                
                # Hitung B_drill (hanya theta_z aktif)
                Bd = self._compute_Bdrill(r, s)
                
                # Diferensial volume area elemen
                dV = abs(detJ) * w * h_avg
                
                # Akumulasi kekakuan: K = K + int( B_drill^T * coeff * B_drill ) dV
                # np.outer(Bd, Bd) ekivalen dengan Bd.T @ Bd karena Bd berukuran 1x(n*6)
                K += np.outer(Bd, Bd) * coeff * dV
                
        return K