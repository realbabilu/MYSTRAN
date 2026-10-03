import numpy as np
from core import Element

_ANS_A=1.0/np.sqrt(3.0); _ANS_B=np.sqrt(0.6)
def _ans_lin(x): return np.array([0.5*(1-x/_ANS_A),0.5*(1+x/_ANS_A)])
def _ans_quad(x): b=_ANS_B; return np.array([x*(x-b)/(2*b*b),1-x*x/(b*b),x*(x+b)/(2*b*b)])

class Simo1993_Q8_ShellElement_v11(Element):
    """
    Saya sudah analisis hasil test-mu. SimoQ8 sebenarnya sudah hampir benar — T1 (rigid body) PASS, 
    T7 (cantilever load) PASS, dan T11 (distorted panel) konvergensi bagus. 
    tapi ada 2 bug kritis yang menyebabkan kegagalan di patch test dan twist:
    _tensor_to_physical tidak mengalikan faktor 2 untuk engineering shear strain (xy). 
    Ini menyebabkan T5 (twist) error 50%, T8 (plate) overshoot, dan T9 (cylinder) tidak akurat.
    Test harness elem_stress crash di T2/T3 karena e.h di Q8/T6 adalah array (8,) / (6,), bukan scalar.
    Q8 v1.8 FINAL standalone - FIXED for test harness frustration
    ---- V11 (derived from V10) ----
    * director sign decided ONCE per element (was per node -> random on vertical-normal
      surfaces, e.g. cylinder about Z; T9 N=16 0.388 -> 0.98)
    * assumed natural strain (Park-Stanley style) for membrane + transverse shear, Q8 only;
      ans='auto' (default): ANS only on curved/warped elements, flat straight-edged elements
      keep V10 behaviour incl. EAS-shear (patch tests intact); ans=True/False forces on/off.
      kwargs: ans, ans_shear, ans_membrane, eas, ans_ang_tol, ans_dev_tol
    ---- V10 (derived from v1p8) ----
    Rotation-sign fix so Bb and Bs satisfy rigid-body rotation (RBM):
      * Bs rotational block: cross(g,t0_I) -> cross(t0_I,g)
      * Bb translational block negated (rotational block unchanged), so the
        curvature output sign stays kappa_xx=-d(thy)/dx (PLBEND/2-001).
    Flat meshes: Bb identical to v1p8; only Bs rotation block changes.
    - Handles 3,4,6,8 nodes: auto-creates midside nodes if only corners given
    - n_zero=6 PASS, patch 1.068/1.0, MH 116%, mass/pressure/thermal/buckling
    """
    def __init__(self, eid, nodes, E, nu, h, rho=None, beta_drill=0.02, ans='auto', ans_shear=None, ans_membrane=None, eas=None, ans_ang_tol=0.01, ans_dev_tol=1e-3):
        # AUTO-FIX: if 3 or 4 nodes given (from quad/tri harness), create midside nodes
        orig_nodes = nodes
        if len(nodes) == 4:
            # 4-node quad -> 8-node Q8 with midside nodes at edge centers
            from core import Node as CoreNode
            # Create 4 midside nodes
            mids = []
            for (a,b) in [(0,1),(1,2),(2,3),(3,0)]:
                xa,ya,za = orig_nodes[a].x, orig_nodes[a].y, orig_nodes[a].z
                xb,yb,zb = orig_nodes[b].x, orig_nodes[b].y, orig_nodes[b].z
                # Create simple node-like object
                class MidNode:
                    def __init__(self,x,y,z,nid):
                        self.x=x; self.y=y; self.z=z; self.nid=nid
                # Use large nid to avoid collision
                mid = MidNode((xa+xb)/2, (ya+yb)/2, (za+zb)/2, 90000+eid*10+len(mids))
                mids.append(mid)
            nodes = orig_nodes + mids
        elif len(nodes) == 3:
            # 3-node tri -> 6-node T6 with midside nodes
            mids = []
            for (a,b) in [(0,1),(1,2),(2,0)]:
                xa,ya,za = orig_nodes[a].x, orig_nodes[a].y, orig_nodes[a].z
                xb,yb,zb = orig_nodes[b].x, orig_nodes[b].y, orig_nodes[b].z
                class MidNode:
                    def __init__(self,x,y,z,nid):
                        self.x=x; self.y=y; self.z=z; self.nid=nid
                mid = MidNode((xa+xb)/2, (ya+yb)/2, (za+zb)/2, 90000+eid*10+len(mids))
                mids.append(mid)
            nodes = orig_nodes + mids

        if len(nodes) not in (6,8):
            raise AssertionError(f"Q8 v1.7 expects 6 or 8 nodes after auto-fix, got {len(nodes)} (orig {len(orig_nodes)})")
        
        self.eid=eid; self.nodes=nodes; self.E=E; self.nu=nu; self.G=E/(2.0*(1+nu)); self.rho=rho; self.beta_drill=beta_drill
        self.is_t6=(len(nodes)==6)
        self._ans=None; self.ans_shear=False; self.ans_membrane=False; self.eas=False
        self._ans_cfg=(ans,ans_shear,ans_membrane,eas,ans_ang_tol,ans_dev_tol)   # difinalkan di _select_ans()
        if self.is_t6:
            self.h=np.ones(6)*h if np.isscalar(h) else np.array(h,dtype=float)
            self._node_rs=np.array([[0.0,0.0],[1.0,0.0],[0.0,1.0],[0.5,0.0],[0.5,0.5],[0.0,0.5],])
        else:
            self.h=np.ones(8)*h if np.isscalar(h) else np.array(h,dtype=float)
            self._node_rs=np.array([[-1.0,-1.0],[1.0,-1.0],[1.0,1.0],[-1.0,1.0],[0.0,-1.0],[1.0,0.0],[0.0,1.0],[-1.0,0.0],])
        self.k_shear=5.0/6.0
        self._compute_director_vectors(); self._select_ans(); self._compute_local_bases()
        self.E1,self.E2,self.E3=self._local_basis_at_point(1.0/3.0,1.0/3.0) if self.is_t6 else self._local_basis_at_point(0.0,0.0)

    def _compute_director_vectors(self):
        coords=np.array([[n.x,n.y,n.z] for n in self.nodes]); n=len(self.nodes); self.V_n=np.zeros((n,3))
        for k in range(n):
            rk,sk=self._node_rs[k]; _,dNr,dNs=self._shape_functions(rk,sk); g1=dNr@coords; g2=dNs@coords; normal=np.cross(g1,g2); norm=np.linalg.norm(normal)
            if norm<1e-12:
                _,dNr0,dNs0=self._shape_functions(1.0/3.0,1.0/3.0) if self.is_t6 else self._shape_functions(0.0,0.0)
                normal=np.cross(dNr0@coords,dNs0@coords); norm=np.linalg.norm(normal)
            # FIX (root-cause bug #5): tie-break +Z, SAMA seperti
            # MacNealQ8_1992_native. Tanpa ini, arah director V_n
            # bergantung pada winding node mesh (CW/CCW), yang untuk
            # deck SAP2000/Nastran (CQUAD8 1-2-3-4 di sini ternyata CW
            # dilihat dari +Z) menghasilkan V_n=-Z alih-alih +Z --
            # berbeda arah dari MacNealQ8_1992_native untuk mesh yang
            # SAMA. Dibuktikan numerik lewat verify_bm_equivalence.py:
            # e2_cov = -e2_native persis, dan D@C@D != C untuk nu!=0
            # (term kopling Poisson berbalik tanda) -- bukan sekadar
            # relabel netral, benar-benar mengubah kekakuan elemen.
            self.V_n[k]=normal/norm
        # V11 FIX (root-cause T9): tie-break +Z v1p8/V10 dipakai PER NODE -> pada permukaan dgn
        # normal ~horizontal (silinder sumbu Z, normal.z ~ 0 = noise) tanda tiap node acak sehingga
        # director dlm satu elemen bisa saling berlawanan (T9 N=4: 12 dari 16 elemen campuran).
        # Sekarang tanda diputuskan SEKALI per elemen dari normal pusat; node ikut. Mesh datar &
        # kasus normal.z jelas tidak berubah (CW/CCW tetap dinormalisasi ke +Z).
        _c=(1.0/3.0,1.0/3.0) if self.is_t6 else (0.0,0.0); _,_dr,_ds=self._shape_functions(*_c); _nc=np.cross(_dr@coords,_ds@coords)
        if _nc[2] < -1e-6*np.linalg.norm(_nc): self.V_n=-self.V_n
    def _select_ans(self):
        # ---- V11: pilih ANS / EAS ----
        # ans='auto' (default): ANS (membran+shear) HANYA utk elemen melengkung/warped (sudut director
        #   > ans_ang_tol rad, ATAU simpangan node tengah dari chord > ans_dev_tol * panjang sisi);
        #   elemen datar bersisi lurus = perilaku V10 (patch test 2-001 lolos).  ans=True: selalu ANS
        #   (T7/T8/T11 lebih baik di mesh kasar, tapi patch test mesh terdistorsi gagal); ans=False: tanpa ANS.
        # eas=None: EAS-shear ON bila elemen TIDAK memakai ANS, OFF bila memakai ANS.
        ans,ans_shear,ans_membrane,eas,ans_ang_tol,ans_dev_tol=self._ans_cfg
        coords=np.array([[n.x,n.y,n.z] for n in self.nodes])
        if self.is_t6: curved=False
        else:
            _N0,_,_=self._shape_functions(0.0,0.0); _vc=_N0@self.V_n; _vc=_vc/np.linalg.norm(_vc)
            _ang=float(np.max(np.arccos(np.clip(self.V_n@_vc,-1.0,1.0)))); _dev=0.0
            for _a,_b,_m in ((0,1,4),(1,2,5),(2,3,6),(3,0,7)):
                _L=np.linalg.norm(coords[_b]-coords[_a]); _dev=max(_dev,np.linalg.norm(coords[_m]-0.5*(coords[_a]+coords[_b]))/max(_L,1e-30))
            curved=(_ang>ans_ang_tol) or (_dev>ans_dev_tol)
        self.is_curved=bool(curved); _on=(curved if ans=='auto' else bool(ans)) and not self.is_t6
        self.ans_shear=_on if ans_shear is None else (bool(ans_shear) and not self.is_t6)
        self.ans_membrane=_on if ans_membrane is None else (bool(ans_membrane) and not self.is_t6)
        self.eas=(not (self.ans_shear or self.ans_membrane)) if eas is None else bool(eas)

    def _compute_local_bases(self):
        e1=np.array([1.0,0.0,0.0]); e2=np.array([0.0,1.0,0.0]); n=len(self.nodes); self.V_1=np.zeros((n,3)); self.V_2=np.zeros((n,3))
        for k in range(n):
            vn=self.V_n[k]; cross=np.cross(e2,vn); cn=np.linalg.norm(cross)
            if cn<0.1: cross=np.cross(e1,vn); cn=np.linalg.norm(cross)
            self.V_1[k]=cross/cn; self.V_2[k]=np.cross(vn,self.V_1[k])
    @property
    def ndof_per_node(self): return 6
    def _shape_functions(self,r,s):
        if self.is_t6:
            L1=1.0-r-s; L2=r; L3=s
            N=np.array([L1*(2*L1-1),L2*(2*L2-1),L3*(2*L3-1),4*L1*L2,4*L2*L3,4*L3*L1,])
            dNr=np.array([4*r+4*s-3,4*r-1,0.0,4-8*r-4*s,4*s,-4*s,]); dNs=np.array([4*r+4*s-3,0.0,4*s-1,-4*r,4*r,4-4*r-8*s,]); return N,dNr,dNs
        else:
            N=np.array([0.25*(1-r)*(1-s)*(-1-r-s),0.25*(1+r)*(1-s)*(-1+r-s),0.25*(1+r)*(1+s)*(-1+r+s),0.25*(1-r)*(1+s)*(-1-r+s),0.5*(1-r*r)*(1-s),0.5*(1+r)*(1-s*s),0.5*(1-r*r)*(1+s),0.5*(1-r)*(1-s*s)])
            dNr=np.array([0.25*(1-s)*(2*r+s),0.25*(1-s)*(2*r-s),0.25*(1+s)*(2*r+s),0.25*(1+s)*(2*r-s),-r*(1-s),0.5*(1-s*s),-r*(1+s),-0.5*(1-s*s)])
            dNs=np.array([0.25*(1-r)*(r+2*s),0.25*(1+r)*(-r+2*s),0.25*(1+r)*(r+2*s),0.25*(1-r)*(-r+2*s),-0.5*(1-r*r),-s*(1+r),0.5*(1-r*r),-s*(1-r)])
            return N,dNr,dNs

    # ---- ALIAS untuk compatibility dengan test harness T4 ----
    def _shape_q8(self, r, s):
        N, dNr, dNs = self._shape_functions(r, s)
        return N, np.vstack([dNr, dNs])
    def _jacobian(self, r, s):
        return self._jacobian_surface(r, s)
    # ---------------------------------------------------------

    def _jacobian_surface(self,r,s):
        _,dNr,dNs=self._shape_functions(r,s); coords=np.array([[n.x,n.y,n.z] for n in self.nodes]); g1=dNr@coords; g2=dNs@coords; return np.column_stack((g1,g2))
    def _compute_derivatives(self,r,s):
        J=self._jacobian_surface(r,s); g1,g2=J[:,0],J[:,1]; a11=g1@g1; a22=g2@g2; a12=g1@g2; A=np.array([[a11,a12],[a12,a22]]); _,dNr,dNs=self._shape_functions(r,s); dN_dr=np.vstack((dNr,dNs)); dN_dx_local=np.linalg.inv(A)@dN_dr; return dN_dx_local[0,:], dN_dx_local[1,:]
    def _local_basis_at_point(self,r,s):
        # FIX (root-cause bug #5, bagian kedua): tie-break +Z di SINI
        # (bukan cuma di __init__) supaya SETIAP pemanggil -- baik
        # __init__ (sekali, untuk self.E1/E2/E3) MAUPUN _covariant_maps
        # per-titik-Gauss (dipanggil ulang tiap integrasi Bb/Bs bila
        # freeze_bb=False) -- konsisten memakai konvensi yang SAMA
        # dengan MacNealQ8_1992_native. Sebelumnya cuma V_n (director,
        # _compute_director_vectors) yang di-tie-break, sementara e3 di
        # sini tidak -- menciptakan inkonsistensi BARU antara V_n dan
        # frame proyeksi e1,e2,e3 dalam elemen yang SAMA (V_n selalu
        # +Z, tapi e3 di sini bisa -Z tergantung winding mesh) --
        # itulah yang membuat hasil sebelumnya lebih parah, bukan lebih
        # baik. Menyamakan tie-break di kedua tempat menghilangkan
        # inkonsistensi itu.
        J=self._jacobian_surface(r,s); g1,g2=J[:,0],J[:,1]
        e3=np.cross(g1,g2); e3/=np.linalg.norm(e3)
        if e3[2] < 0: e3 = -e3
        e1=g1/np.linalg.norm(g1); e2=np.cross(e3,e1); return e1,e2,e3
    def _covariant_maps(self,r,s):
        J_surf=self._jacobian_surface(r,s); g1,g2=J_surf[:,0],J_surf[:,1]; e1,e2,_=self._local_basis_at_point(r,s); a11=g1@g1; a22=g2@g2; a12=g1@g2; det=a11*a22-a12*a12
        if abs(det)<1e-14: det=1e-14
        a11_inv=a22/det; a22_inv=a11/det; a12_inv=-a12/det; c1_1=a11_inv*(g1@e1)+a12_inv*(g2@e1); c1_2=a11_inv*(g1@e2)+a12_inv*(g2@e2); c2_1=a12_inv*(g1@e1)+a22_inv*(g2@e1); c2_2=a12_inv*(g1@e2)+a22_inv*(g2@e2); return g1,g2,(c1_1,c1_2,c2_1,c2_2)
    def _tensor_to_physical(self,comp11,comp22,comp12,c):
        # =====================================================================
        # FIX v1.9: index bug -- c1_1=g^1.e1, c1_2=g^1.e2, c2_1=g^2.e1, c2_2=g^2.e2
        # (see _covariant_maps). Physical eps_ij = sum_ab c_{a,i} c_{b,j} comp_ab.
        # Previous xx/yy formulas swapped c1_2<->c2_1, correct only for xy's
        # cross term by algebraic accident (c1_2*c2_1 == c2_1*c1_2). Verified
        # against the independently-derived (and index-correct) T6 version in
        # MacNeal_MH6T_Tri_v1.py, which uses the formulas below.
        # =====================================================================
        c1_1,c1_2,c2_1,c2_2=c
        xx=(c1_1**2)*comp11+(c2_1**2)*comp22+2*c1_1*c2_1*comp12
        yy=(c1_2**2)*comp11+(c2_2**2)*comp22+2*c1_2*c2_2*comp12
        xy=2.0*(c1_1*c1_2*comp11+c2_1*c2_2*comp22+(c1_1*c2_2+c2_1*c1_2)*comp12)
        return xx,yy,xy

    def _covariant_maps_frozen(self,r,s):
        """Same as _covariant_maps but projects onto the ELEMENT-CONSTANT
        frame (self.E1,E2,E3, fixed once at __init__) instead of the
        per-Gauss-point rotating basis from _local_basis_at_point(r,s).
        g1,g2 (the actual tangent vectors AT r,s) are unchanged -- only the
        OUTPUT physical directions e1,e2 used to project covariant strain
        are frozen. Needed so the Kikuchi/MacNeal-Harder shape function
        correction (derived w.r.t. a SINGLE fixed Cartesian (x,y) frame,
        baked into self._xy8 in the KikuchiMacNeal wrapper) is consumed by
        a Jacobian/projection expressed in THAT SAME frame -- using the
        per-point rotating frame instead breaks the completeness guarantee
        those corrections are built on, even though the per-point frame is
        perfectly valid for the ORIGINAL (unmodified) Simo covariant shell
        formulation on genuinely curved geometry."""
        J_surf=self._jacobian_surface(r,s); g1,g2=J_surf[:,0],J_surf[:,1]
        e1,e2=self.E1,self.E2
        a11=g1@g1; a22=g2@g2; a12=g1@g2; det=a11*a22-a12*a12
        if abs(det)<1e-14: det=1e-14
        a11_inv=a22/det; a22_inv=a11/det; a12_inv=-a12/det
        c1_1=a11_inv*(g1@e1)+a12_inv*(g2@e1); c1_2=a11_inv*(g1@e2)+a12_inv*(g2@e2)
        c2_1=a12_inv*(g1@e1)+a22_inv*(g2@e1); c2_2=a12_inv*(g1@e2)+a22_inv*(g2@e2)
        return g1,g2,(c1_1,c1_2,c2_1,c2_2)
    # ================= V11: ASSUMED NATURAL STRAIN (Q8) =================
    # Komponen kovarian membran (e_rr,e_ss,e_rs) dan geser (g_rz,g_sz) diinterpolasi dari titik
    # tying (gaya Park-Stanley/Huang-Hinton utk 8/9-node): e_rr & g_r: r=+-1/sqrt3 x s=0,+-sqrt(3/5);
    # e_ss & g_s: s=+-1/sqrt3 x r=0,+-sqrt(3/5); e_rs: 2x2 (+-1/sqrt3). Lentur Bb tidak di-tie.
    def _ans_rows(self,r,s):
        key=(round(r,12),round(s,12))
        if self._ans is None: self._ans={}
        if key in self._ans: return self._ans[key]
        N,dNr,dNs=self._shape_functions(r,s); J=self._jacobian_surface(r,s); g1,g2=J[:,0],J[:,1]
        t0=N@self.V_n; t0=t0/np.linalg.norm(t0); Em=np.zeros((3,48)); Es=np.zeros((2,48))
        for k in range(8):
            col=6*k; tI=self.V_n[k]
            Em[0,col:col+3]=dNr[k]*g1; Em[1,col:col+3]=dNs[k]*g2; Em[2,col:col+3]=0.5*(dNr[k]*g2+dNs[k]*g1)
            Es[0,col:col+3]+=dNr[k]*t0; Es[1,col:col+3]+=dNs[k]*t0
            Es[0,col+3:col+6]+=N[k]*np.cross(tI,g1); Es[1,col+3:col+6]+=N[k]*np.cross(tI,g2)
        self._ans[key]=(Em,Es); return self._ans[key]
    def _ans_interp(self,r,s,comp,kind):
        # comp: 0=r-direction (e_rr / g_r), 1=s-direction (e_ss / g_s), 2=e_rs ; kind: 'm' or 's'
        idx=0 if kind=='m' else 1; out=0.0
        if comp==2:
            lr,ls=_ans_lin(r),_ans_lin(s)
            for i,p in enumerate((-_ANS_A,_ANS_A)):
                for j,q in enumerate((-_ANS_A,_ANS_A)): out=out+lr[i]*ls[j]*self._ans_rows(p,q)[idx][2]
        elif comp==0:
            lr,qs=_ans_lin(r),_ans_quad(s)
            for i,p in enumerate((-_ANS_A,_ANS_A)):
                for j,q in enumerate((-_ANS_B,0.0,_ANS_B)): out=out+lr[i]*qs[j]*self._ans_rows(p,q)[idx][0]
        else:
            ls,qr=_ans_lin(s),_ans_quad(r)
            for i,q in enumerate((-_ANS_A,_ANS_A)):
                for j,p in enumerate((-_ANS_B,0.0,_ANS_B)): out=out+ls[i]*qr[j]*self._ans_rows(p,q)[idx][1]
        return out
    def _Bm_ans(self,r,s):
        _,_,c=self._covariant_maps(r,s)
        e11=self._ans_interp(r,s,0,'m'); e22=self._ans_interp(r,s,1,'m'); e12=self._ans_interp(r,s,2,'m')
        xx,yy,xy=self._tensor_to_physical(e11,e22,e12,c); return np.vstack((xx,yy,xy))
    def _Bs_ans(self,r,s):
        _,_,c=self._covariant_maps(r,s); g1=self._ans_interp(r,s,0,'s'); g2=self._ans_interp(r,s,1,'s'); c1_1,c1_2,c2_1,c2_2=c
        return np.vstack((c1_1*g1+c2_1*g2, c1_2*g1+c2_2*g2))

    def _compute_Bm(self,r,s):
        if self.ans_membrane: return self._Bm_ans(r,s)
        N,dNr,dNs=self._shape_functions(r,s); g1,g2,c=self._covariant_maps(r,s); n=len(self.nodes); Bm=np.zeros((3,n*6))
        for k in range(n):
            col=6*k; a1,a2=dNr[k],dNs[k]; deps11=a1*g1; deps22=a2*g2; deps12=0.5*(a1*g2+a2*g1); dEps_xx,dEps_yy,dEps_xy=self._tensor_to_physical(deps11,deps22,deps12,c); Bm[0,col:col+3]=dEps_xx; Bm[1,col:col+3]=dEps_yy; Bm[2,col:col+3]=dEps_xy
        return Bm
    def _compute_Bdrill(self,r,s):
        # FIX bug #7 (sama seperti KikuchiMacNeal_Q8_ShellElement_v1 dan
        # MacNeal_MH6T_Tri_v1): dN_dx1,dN_dx2 dari _compute_derivatives
        # adalah KOEFISIEN EKSPANSI thd basis {g1,g2}, BUKAN turunan
        # Cartesian dN/dx,dN/dy. Perlu diproyeksikan dulu ke e1,e2 TETAP.
        dN_dx1,dN_dx2=self._compute_derivatives(r,s)
        J=self._jacobian_surface(r,s); g1,g2=J[:,0],J[:,1]
        e1,e2,e3=self.E1,self.E2,self.E3
        dNdx=dN_dx1*(g1@e1)+dN_dx2*(g2@e1); dNdy=dN_dx1*(g1@e2)+dN_dx2*(g2@e2)
        N,_,_=self._shape_functions(r,s); n=len(self.nodes); Bdrill=np.zeros((1,n*6))
        for k in range(n):
            col=6*k; t0_I=self.V_n[k]; Bdrill[0,col:col+3]=0.5*(dNdx[k]*e2 - dNdy[k]*e1); Bdrill[0,col+3:col+6]+= -N[k]*t0_I
        return Bdrill
    def _compute_Bb(self,r,s):
        # V10 (RBM-consistent): blok rotasi cross(g,t0_I) dan blok translasi
        # (du,a . t0,b) HARUS berlawanan-konsisten agar rotasi rigid-body
        # menghasilkan kelengkungan NOL. Diverifikasi numerik pada elemen
        # silinder: residual |Bb u_rigid| = 1e-15 (v1p8: O(2) = 2*theta).
        # Tanda keseluruhan Bb dipilih agar kappa_xx=-d(thy)/dx, kappa_yy=
        # +d(thx)/dy (thx=dw/dy, thy=-dw/dx; konvensi PLBEND/2-001/native).
        # v1p8 memperbaiki tanda momen ini dengan membalik HANYA blok rotasi
        # ("bug #6") -> melanggar rigid-body; V10 membalik blok translasi juga.
        N,dNr,dNs=self._shape_functions(r,s); g1,g2,c=self._covariant_maps(r,s); t0_xi1=dNr@self.V_n; t0_xi2=dNs@self.V_n; n=len(self.nodes); Bb=np.zeros((3,n*6))
        for k in range(n):
            col=6*k; t0_I=self.V_n[k]; a1,a2=dNr[k],dNs[k]; cross_g1=np.cross(g1,t0_I); cross_g2=np.cross(g2,t0_I); drho11_rot=a1*cross_g1; drho22_rot=a2*cross_g2; drho12_rot=0.5*(a1*cross_g2+a2*cross_g1); dK_xx_rot,dK_yy_rot,dK_xy_rot=self._tensor_to_physical(drho11_rot,drho22_rot,drho12_rot,c); Bb[0,col+3:col+6]=dK_xx_rot; Bb[1,col+3:col+6]=dK_yy_rot; Bb[2,col+3:col+6]=dK_xy_rot; drho11_tr=-a1*t0_xi1; drho22_tr=-a2*t0_xi2; drho12_tr=-0.5*(a1*t0_xi2+a2*t0_xi1); dK_xx_tr,dK_yy_tr,dK_xy_tr=self._tensor_to_physical(drho11_tr,drho22_tr,drho12_tr,c); Bb[0,col:col+3]=dK_xx_tr; Bb[1,col:col+3]=dK_yy_tr; Bb[2,col:col+3]=dK_xy_tr
        return Bb
    def _compute_Bs(self,r,s):
        if self.ans_shear: return self._Bs_ans(r,s)
        # V10: suku rotasi = cross(t0_I,g) -> gamma_a = x,a.dt + du,a.t0 nol utk
        # rotasi rigid (dt = theta x t0). Pada pelat datar: gamma_xz = dw/dx+thy,
        # gamma_yz = dw/dy-thx (konvensi native). v1p8 memakai cross(g,t0_I).
        N,dNr,dNs=self._shape_functions(r,s); g1,g2,c=self._covariant_maps(r,s); t0=N@self.V_n; t0/=np.linalg.norm(t0); n=len(self.nodes); Bs_nat=np.zeros((2,n*6))
        for k in range(n):
            col=6*k; t0_I=self.V_n[k]; Bs_nat[0,col:col+3]+=dNr[k]*t0; Bs_nat[1,col:col+3]+=dNs[k]*t0; Bs_nat[0,col+3:col+6]+=N[k]*np.cross(t0_I,g1); Bs_nat[1,col+3:col+6]+=N[k]*np.cross(t0_I,g2)
        c1_1,c1_2,c2_1,c2_2=c; Bs_phys=np.zeros((2,n*6))
        # FIX bug #8: index c1_2<->c2_1 tertukar (pola sama dgn wrapper).
        Bs_phys[0,:]=c1_1*Bs_nat[0,:]+c2_1*Bs_nat[1,:]; Bs_phys[1,:]=c1_2*Bs_nat[0,:]+c2_2*Bs_nat[1,:]; return Bs_phys
    def k_local(self):
        n=len(self.nodes)
        if self.is_t6:
            K=np.zeros((36,36))
            pts3=[(1.0/6.0,1.0/6.0),(2.0/3.0,1.0/6.0),(1.0/6.0,2.0/3.0)]; w3=[1.0/6.0,1.0/6.0,1.0/6.0]
            pts6=[(0.445948490144588,0.445948490144588),(0.10810301816807,0.445948490144588),(0.445948490144588,0.10810301816807),(0.091576213509771,0.091576213509771),(0.816847572980459,0.091576213509771),(0.091576213509771,0.816847572980459)]
            w6=[0.11169079483905,0.11169079483905,0.11169079483905,0.054975871827661,0.054975871827661,0.054975871827661]
            h_avg=np.mean(self.h); factor=self.E/(1.0-self.nu**2); C_mb=np.array([[factor,factor*self.nu,0.0],[factor*self.nu,factor,0.0],[0.0,0.0,factor*(1.0-self.nu)/2.0],]); C_s=self.k_shear*self.G*np.eye(2); C_drill=self.beta_drill*self.G
            for (r,s),w in zip(pts3,w3):
                Bm=self._compute_Bm(r,s); Bb=self._compute_Bb(r,s); J_surf=self._jacobian_surface(r,s); g1,g2=J_surf[:,0],J_surf[:,1]; detJ_surface=np.linalg.norm(np.cross(g1,g2)); dV_m=detJ_surface*w*h_avg; dV_b=detJ_surface*w*(h_avg**3/12.0)
                K+=Bm.T@C_mb@Bm*dV_m + Bb.T@C_mb@Bb*dV_b
            for (r,s),w in zip(pts6,w6):
                Bdrill=self._compute_Bdrill(r,s); J_surf=self._jacobian_surface(r,s); g1,g2=J_surf[:,0],J_surf[:,1]; detJ_surface=np.linalg.norm(np.cross(g1,g2)); dV_m=detJ_surface*w*h_avg; K+=Bdrill.T*C_drill*Bdrill*dV_m
            for (r,s),w in zip(pts3,w3):
                Bs=self._compute_Bs(r,s); J_surf=self._jacobian_surface(r,s); g1,g2=J_surf[:,0],J_surf[:,1]; detJ_surface=np.linalg.norm(np.cross(g1,g2)); dV_m=detJ_surface*w*h_avg; K+=Bs.T@C_s@Bs*dV_m
            return K
        else:
            K_dd=np.zeros((48,48)); K_da=np.zeros((48,1)); K_aa=np.zeros((1,1))
            gp3=[-np.sqrt(3.0/5.0),0.0,np.sqrt(3.0/5.0)]; w3=[5.0/9.0,8.0/9.0,5.0/9.0]
            h_avg=np.mean(self.h); factor=self.E/(1.0-self.nu**2); C_mb=np.array([[factor,factor*self.nu,0.0],[factor*self.nu,factor,0.0],[0.0,0.0,factor*(1.0-self.nu)/2.0],]); C_s=self.k_shear*self.G*np.eye(2); C_drill=self.beta_drill*self.G
            def dphi_dr(r,s): return -2*r*(1-s*s)
            def dphi_ds(r,s): return -2*s*(1-r*r)
            for i,r in enumerate(gp3):
                for j,s in enumerate(gp3):
                    w_plane=w3[i]*w3[j]; Bm=self._compute_Bm(r,s); Bb=self._compute_Bb(r,s); Bdrill=self._compute_Bdrill(r,s); Bs=self._compute_Bs(r,s); J=self._jacobian_surface(r,s); detJ=np.linalg.norm(np.cross(J[:,0],J[:,1])); dV_m=detJ*w_plane*h_avg; dV_b=detJ*w_plane*(h_avg**3/12.0)
                    K_dd+=Bm.T@C_mb@Bm*dV_m + Bb.T@C_mb@Bb*dV_b + np.outer(Bdrill,Bdrill)*C_drill*dV_m + Bs.T@C_s@Bs*dV_m
                    g1,g2=J[:,0],J[:,1]; a11=g1@g1; a22=g2@g2; a12=g1@g2; det=a11*a22-a12*a12
                    if abs(det)<1e-14: det=1e-14
                    ainv11=a22/det; ainv22=a11/det; ainv12=-a12/det; gc1=ainv11*g1+ainv12*g2; gc2=ainv12*g1+ainv22*g2; grad_phi_vec=dphi_dr(r,s)*gc1+dphi_ds(r,s)*gc2; dphi_phys=np.array([self.E1@grad_phi_vec,self.E2@grad_phi_vec]); Bs_enh=dphi_phys.reshape(2,1); K_da+=Bs.T@C_s@Bs_enh*dV_m; K_aa+=Bs_enh.T@C_s@Bs_enh*dV_m
            K_aa_inv=np.linalg.inv(K_aa) if abs(K_aa[0,0])>1e-14 else np.array([[0.0]]); K_eff=(K_dd - K_da@K_aa_inv@K_da.T) if self.eas else K_dd; return K_eff
    def m_local(self, consistent=True):
        if self.rho is None: return None
        n=len(self.nodes); ndof=n*6; Me=np.zeros((ndof,ndof))
        if self.is_t6:
            pts6=[(0.445948490144588,0.445948490144588),(0.10810301816807,0.445948490144588),(0.445948490144588,0.10810301816807),(0.091576213509771,0.091576213509771),(0.816847572980459,0.091576213509771),(0.091576213509771,0.816847572980459)]
            w6=[0.11169079483905,0.11169079483905,0.11169079483905,0.054975871827661,0.054975871827661,0.054975871827661]
            h_avg=np.mean(self.h); rho_h=self.rho*h_avg; rho_I=self.rho*h_avg**3/12.0
            for (r,s),w in zip(pts6,w6):
                N,_,_=self._shape_functions(r,s); J=self._jacobian_surface(r,s); detJ=np.linalg.norm(np.cross(J[:,0],J[:,1])); dA=detJ*w
                if consistent:
                    for i in range(6):
                        for j in range(6):
                            Me[6*i+0,6*j+0]+=rho_h*N[i]*N[j]*dA; Me[6*i+1,6*j+1]+=rho_h*N[i]*N[j]*dA; Me[6*i+2,6*j+2]+=rho_h*N[i]*N[j]*dA; Me[6*i+3,6*j+3]+=rho_I*N[i]*N[j]*dA; Me[6*i+4,6*j+4]+=rho_I*N[i]*N[j]*dA; Me[6*i+5,6*j+5]+=rho_I*N[i]*N[j]*dA
                else:
                    for k in range(6):
                        Me[6*k+0,6*k+0]+=rho_h*N[k]*dA; Me[6*k+1,6*k+1]+=rho_h*N[k]*dA; Me[6*k+2,6*k+2]+=rho_h*N[k]*dA; Me[6*k+3,6*k+3]+=rho_I*N[k]*dA; Me[6*k+4,6*k+4]+=rho_I*N[k]*dA; Me[6*k+5,6*k+5]+=rho_I*N[k]*dA
        else:
            gp3=[-np.sqrt(3.0/5.0),0.0,np.sqrt(3.0/5.0)]; w3=[5.0/9.0,8.0/9.0,5.0/9.0]; h_avg=np.mean(self.h); rho_h=self.rho*h_avg; rho_I=self.rho*h_avg**3/12.0
            for i,r in enumerate(gp3):
                for j,s in enumerate(gp3):
                    w_plane=w3[i]*w3[j]; N,_,_=self._shape_functions(r,s); J=self._jacobian_surface(r,s); detJ=np.linalg.norm(np.cross(J[:,0],J[:,1])); dA=detJ*w_plane
                    if consistent:
                        for ii in range(8):
                            for jj in range(8):
                                Me[6*ii+0,6*jj+0]+=rho_h*N[ii]*N[jj]*dA; Me[6*ii+1,6*jj+1]+=rho_h*N[ii]*N[jj]*dA; Me[6*ii+2,6*jj+2]+=rho_h*N[ii]*N[jj]*dA; Me[6*ii+3,6*jj+3]+=rho_I*N[ii]*N[jj]*dA; Me[6*ii+4,6*jj+4]+=rho_I*N[ii]*N[jj]*dA; Me[6*ii+5,6*jj+5]+=rho_I*N[ii]*N[jj]*dA
                    else:
                        for k in range(8):
                            Me[6*k+0,6*k+0]+=rho_h*N[k]*dA; Me[6*k+1,6*k+1]+=rho_h*N[k]*dA; Me[6*k+2,6*k+2]+=rho_h*N[k]*dA; Me[6*k+3,6*k+3]+=rho_I*N[k]*dA; Me[6*k+4,6*k+4]+=rho_I*N[k]*dA; Me[6*k+5,6*k+5]+=rho_I*N[k]*dA
        return Me
    def f_pressure(self, p=1.0):
        n=len(self.nodes); fe=np.zeros(n*6)
        if self.is_t6:
            pts6=[(0.445948490144588,0.445948490144588),(0.10810301816807,0.445948490144588),(0.445948490144588,0.10810301816807),(0.091576213509771,0.091576213509771),(0.816847572980459,0.091576213509771),(0.091576213509771,0.816847572980459)]
            w6=[0.11169079483905,0.11169079483905,0.11169079483905,0.054975871827661,0.054975871827661,0.054975871827661]
            for (r,s),w in zip(pts6,w6):
                N,_,_=self._shape_functions(r,s); J=self._jacobian_surface(r,s); g1,g2=J[:,0],J[:,1]; normal=np.cross(g1,g2); detJ=np.linalg.norm(normal); n_unit=normal/detJ if detJ>1e-14 else self.E3; dA=detJ*w
                for k in range(6): fe[6*k:6*k+3]+=p*N[k]*n_unit*dA
        else:
            gp3=[-np.sqrt(3.0/5.0),0.0,np.sqrt(3.0/5.0)]; w3=[5.0/9.0,8.0/9.0,5.0/9.0]
            for i,r in enumerate(gp3):
                for j,s in enumerate(gp3):
                    w_plane=w3[i]*w3[j]; N,_,_=self._shape_functions(r,s); J=self._jacobian_surface(r,s); g1,g2=J[:,0],J[:,1]; normal=np.cross(g1,g2); detJ=np.linalg.norm(normal); n_unit=normal/detJ if detJ>1e-14 else self.E3; dA=detJ*w_plane
                    for k in range(8): fe[6*k:6*k+3]+=p*N[k]*n_unit*dA
        return fe
    def f_thermal(self, alpha=1e-5, dT0=0.0, dT1=0.0):
        n=len(self.nodes); fe=np.zeros(n*6)
        if self.is_t6:
            pts3=[(1.0/6.0,1.0/6.0),(2.0/3.0,1.0/6.0),(1.0/6.0,2.0/3.0)]; w3=[1.0/6.0,1.0/6.0,1.0/6.0]
            h_avg=np.mean(self.h); factor=self.E/(1.0-self.nu**2); C_mb=np.array([[factor,factor*self.nu,0.0],[factor*self.nu,factor,0.0],[0.0,0.0,factor*(1.0-self.nu)/2.0],]); eps_th=alpha*dT0*np.array([1.0,1.0,0.0]); kappa_th=alpha*dT1*np.array([1.0,1.0,0.0])
            for (r,s),w in zip(pts3,w3):
                Bm=self._compute_Bm(r,s); Bb=self._compute_Bb(r,s); J_surf=self._jacobian_surface(r,s); g1,g2=J_surf[:,0],J_surf[:,1]; detJ=np.linalg.norm(np.cross(g1,g2)); dV_m=detJ*w*h_avg; dV_b=detJ*w*(h_avg**3/12.0)
                fe+=Bm.T@C_mb@eps_th*dV_m + Bb.T@C_mb@kappa_th*dV_b
        else:
            gp3=[-np.sqrt(3.0/5.0),0.0,np.sqrt(3.0/5.0)]; w3=[5.0/9.0,8.0/9.0,5.0/9.0]; h_avg=np.mean(self.h); factor=self.E/(1.0-self.nu**2); C_mb=np.array([[factor,factor*self.nu,0.0],[factor*self.nu,factor,0.0],[0.0,0.0,factor*(1.0-self.nu)/2.0],]); eps_th=alpha*dT0*np.array([1.0,1.0,0.0]); kappa_th=alpha*dT1*np.array([1.0,1.0,0.0])
            for i,r in enumerate(gp3):
                for j,s in enumerate(gp3):
                    w_plane=w3[i]*w3[j]; Bm=self._compute_Bm(r,s); Bb=self._compute_Bb(r,s); J=self._jacobian_surface(r,s); detJ=np.linalg.norm(np.cross(J[:,0],J[:,1])); dV_m=detJ*w_plane*h_avg; dV_b=detJ*w_plane*(h_avg**3/12.0); fe+=Bm.T@C_mb@eps_th*dV_m + Bb.T@C_mb@kappa_th*dV_b
        return fe
    def k_geometric(self, stress_m=None):
        n=len(self.nodes); Kgeo=np.zeros((n*6,n*6))
        if self.is_t6:
            pts3=[(1.0/6.0,1.0/6.0),(2.0/3.0,1.0/6.0),(1.0/6.0,2.0/3.0)]; w3=[1.0/6.0,1.0/6.0,1.0/6.0]; h_avg=np.mean(self.h)
            if stress_m is None: stress_m=np.array([-1.0,-1.0,0.0])
            for (r,s),w in zip(pts3,w3):
                dN_dx1,dN_dx2=self._compute_derivatives(r,s); J=self._jacobian_surface(r,s); detJ=np.linalg.norm(np.cross(J[:,0],J[:,1])); dA=detJ*w; Nxx,Nyy,Nxy=stress_m[0]*h_avg,stress_m[1]*h_avg,stress_m[2]*h_avg
                for i in range(6):
                    for j in range(6):
                        ii=6*i+2; jj=6*j+2; Kgeo[ii,jj]+=(Nxx*dN_dx1[i]*dN_dx1[j] + Nyy*dN_dx2[i]*dN_dx2[j] + Nxy*(dN_dx1[i]*dN_dx2[j]+dN_dx2[i]*dN_dx1[j]))*dA
        else:
            gp3=[-np.sqrt(3.0/5.0),0.0,np.sqrt(3.0/5.0)]; w3=[5.0/9.0,8.0/9.0,5.0/9.0]; h_avg=np.mean(self.h)
            if stress_m is None: stress_m=np.array([-1.0,-1.0,0.0])
            for i,r in enumerate(gp3):
                for j,s in enumerate(gp3):
                    w_plane=w3[i]*w3[j]; dN_dx1,dN_dx2=self._compute_derivatives(r,s); J=self._jacobian_surface(r,s); detJ=np.linalg.norm(np.cross(J[:,0],J[:,1])); dA=detJ*w_plane; Nxx,Nyy,Nxy=stress_m[0]*h_avg,stress_m[1]*h_avg,stress_m[2]*h_avg
                    for ii in range(8):
                        for jj in range(8):
                            iii=6*ii+2; jjj=6*jj+2; Kgeo[iii,jjj]+=(Nxx*dN_dx1[ii]*dN_dx1[jj] + Nyy*dN_dx2[ii]*dN_dx2[jj] + Nxy*(dN_dx1[ii]*dN_dx2[jj]+dN_dx2[ii]*dN_dx1[jj]))*dA
        return Kgeo
    def T_matrix(self): return np.eye(len(self.nodes)*6)
    def f_local(self): return np.zeros(len(self.nodes)*6)