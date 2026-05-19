"""
core_sparse.py -- sparse FEM engine compatible with core.py for large shell meshes.

Drop-in subset of core.py API used by scordelis.py:
Node, Element, Model, BC, NodalLoad.
Assembly uses scipy.sparse COO/CSR and static solve uses scipy.sparse.linalg.spsolve.
"""
from __future__ import annotations
import numpy as np
from dataclasses import dataclass
from typing import Optional
try:
    from scipy.sparse import coo_matrix, csr_matrix, lil_matrix, eye as speye
    from scipy.sparse.linalg import spsolve, eigsh
except Exception as exc:
    raise ImportError("core_sparse.py requires scipy: pip install scipy") from exc

class Node:
    def __init__(self, nid:int, x:float, y:float, z:float):
        self.nid=nid; self.x=float(x); self.y=float(y); self.z=float(z); self.dofs=[]
    @property
    def coords(self): return np.array([self.x,self.y,self.z], dtype=float)
    def __repr__(self): return f"Node({self.nid}: [{self.x:.4g}, {self.y:.4g}, {self.z:.4g}])"
@dataclass
class BC: dof:int; value:float=0.0
@dataclass
class NodalLoad: dof:int; value:float
@dataclass
class _PendingBC: node_id:int; dof_local:int; value:float=0.0
@dataclass
class _PendingLoad: node_id:int; dof_local:int; value:float

class Element:
    @property
    def ndof_per_node(self): raise NotImplementedError
    def k_local(self): raise NotImplementedError
    def m_local(self): return None
    def T_matrix(self): return np.eye(len(self.nodes)*self.ndof_per_node)
    def k_global(self):
        T=self.T_matrix(); return T.T @ self.k_local() @ T
    def m_global(self):
        M=self.m_local()
        if M is None: return None
        T=self.T_matrix(); return T.T @ M @ T
    def global_dof_indices(self):
        dofs=[]
        for node in self.nodes: dofs.extend(node.dofs[:self.ndof_per_node])
        return dofs

class Model:
    def __init__(self, ndof_per_node:int=6, autospc:bool=True, autospc_tol:float=1e-8, autospc_verbose:bool=False):
        self.ndof_per_node=ndof_per_node; self.autospc=autospc; self.autospc_tol=float(autospc_tol); self.autospc_verbose=autospc_verbose
        self.nodes={}; self.elements=[]; self.bcs=[]; self.loads=[]; self._nodal_masses={}
        self._built=False; self.ndof=0; self.K=None; self.M=None; self.F=None; self.u=None; self._autospc_dofs=[]
    def add_node(self,nid,x,y,z):
        node=Node(nid,x,y,z); self.nodes[nid]=node; self._built=False; return node
    def add_element(self,elem): self.elements.append(elem); self._built=False
    def add_bc(self,node_id,dof_local,value=0.0): self.bcs.append(_PendingBC(node_id,dof_local,value))
    def add_load(self,node_id,dof_local,value): self.loads.append(_PendingLoad(node_id,dof_local,value))
    def fix_node(self,node_id,dofs=None):
        if dofs is None: dofs=list(range(self.ndof_per_node))
        for d in dofs: self.add_bc(node_id,d,0.0)
    def add_nodal_mass(self,node_id,masses):
        if node_id not in self._nodal_masses: self._nodal_masses[node_id]=np.zeros(self.ndof_per_node)
        self._nodal_masses[node_id]+=np.array(masses,dtype=float)
    def build(self):
        counter=0
        for nid in sorted(self.nodes):
            nd=self.nodes[nid]; nd.dofs=list(range(counter,counter+self.ndof_per_node)); counter+=self.ndof_per_node
        self.ndof=counter
        self._resolved_bcs=[]; self._resolved_loads=[]
        for pbc in self.bcs: self._resolved_bcs.append(BC(self.nodes[pbc.node_id].dofs[pbc.dof_local],pbc.value))
        for pl in self.loads: self._resolved_loads.append(NodalLoad(self.nodes[pl.node_id].dofs[pl.dof_local],pl.value))
        rows=[]; cols=[]; data=[]
        for elem in self.elements:
            Ke=elem.k_global()
            if hasattr(elem, 'global_dof_indices'):
                idx=elem.global_dof_indices()
            else:
                idx=[]
                for nd in elem.nodes:
                    idx.extend(nd.dofs[:elem.ndof_per_node])
            for ii,gi in enumerate(idx):
                r=Ke[ii]
                for jj,gj in enumerate(idx):
                    v=r[jj]
                    if v!=0.0: rows.append(gi); cols.append(gj); data.append(float(v))
        self.K=coo_matrix((data,(rows,cols)),shape=(self.ndof,self.ndof)).tocsr(); self.K.sum_duplicates()
        self.M=csr_matrix((self.ndof,self.ndof))
        if self._nodal_masses:
            Mlil=lil_matrix((self.ndof,self.ndof))
            for nid,mv in self._nodal_masses.items():
                nd=self.nodes[nid]
                for i,m in enumerate(mv): Mlil[nd.dofs[i],nd.dofs[i]] += m
            self.M=Mlil.tocsr()
        self.F=np.zeros(self.ndof)
        for ld in self._resolved_loads: self.F[ld.dof]+=ld.value
        self._autospc_dofs=[]
        if self.autospc: self._run_autospc()
        self._built=True
    def _run_autospc(self):
        if self.ndof==0: return
        user_fixed={bc.dof for bc in self._resolved_bcs}; kd=np.abs(self.K.diagonal()); km=kd.max() if kd.size else 0.0
        if km<1e-30: return
        th=self.autospc_tol*km
        self._autospc_dofs=[d for d in range(self.ndof) if d not in user_fixed and kd[d]<th]
        for d in self._autospc_dofs: self._resolved_bcs.append(BC(d,0.0))
        if self._autospc_dofs and self.autospc_verbose: print(f"[AutoSPC] Mengunci {len(self._autospc_dofs)} DOF zero-stiffness")
    @property
    def autospc_count(self): return len(self._autospc_dofs)
    def solve_static(self):
        if not self._built: self.build()
        fixed_dofs=np.array([bc.dof for bc in self._resolved_bcs],dtype=int); fixed_vals=np.array([bc.value for bc in self._resolved_bcs],dtype=float)
        mask=np.zeros(self.ndof,dtype=bool); mask[fixed_dofs]=True; free=np.flatnonzero(~mask)
        F=self.F.copy(); nz=np.flatnonzero(np.abs(fixed_vals)>0.0)
        if nz.size:
            fd=fixed_dofs[nz]; fv=fixed_vals[nz]; F[free]-=self.K[free[:,None],fd].dot(fv)
        Kff=self.K[free[:,None],free].tocsc(); Ff=F[free]
        u=np.zeros(self.ndof); u[fixed_dofs]=fixed_vals; u[free]=spsolve(Kff,Ff); self.u=u; return u
    def solve_eigen(self,n_modes=6):
        if not self._built: self.build()
        fixed={bc.dof for bc in self._resolved_bcs}; free=np.array([d for d in range(self.ndof) if d not in fixed],dtype=int)
        Kff=self.K[free[:,None],free]; Mff=self.M[free[:,None],free]+speye(len(free),format='csr')*1e-12
        n_modes=min(n_modes,len(free)-2); vals,vecs=eigsh(Kff,k=n_modes,M=Mff,sigma=0.0,which='LM')
        vals=np.maximum(vals,0.0); freqs=np.sqrt(vals)/(2*np.pi); modes=np.zeros((self.ndof,n_modes)); modes[free,:]=vecs
        self.freqs=freqs; self.modes=modes; return freqs,modes
    def get_displacement(self,node_id): return self.u[self.nodes[node_id].dofs[:self.ndof_per_node]]
    def get_reaction(self,node_id): return (self.K.dot(self.u)-self.F)[self.nodes[node_id].dofs[:self.ndof_per_node]]
    def summary(self):
        nbc=len([b for b in self._resolved_bcs if b.dof not in self._autospc_dofs]) if self._built else len(self.bcs)
        return '\n'.join(['Model Summary',f'  Nodes    : {len(self.nodes)}',f'  Elements : {len(self.elements)}',f'  DOF total: {self.ndof}',f'  BC (user): {nbc}',f'  AutoSPC  : {len(self._autospc_dofs)} DOF',f'  Loads    : {len(self.loads)}',f'  K format : scipy.sparse CSR, nnz={self.K.nnz if self.K is not None else 0}'])
