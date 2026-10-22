"""Direct bordered fixed-point solver for an IPM inner profile.

No physical or rescaled time integration. All fields have (ny,nx) layout.
The Poisson boundary is supplied explicitly; changing it changes the box problem.
"""
from dataclasses import dataclass
import numpy as np
from scipy import sparse
from scipy.sparse.linalg import splu, LinearOperator, gmres


def fd_matrix(x, derivative=1, stencil=5):
    """Polynomial-exact local differences; one-sided at endpoints."""
    x = np.asarray(x, float)
    n = len(x)
    m = min(stencil, n)
    if derivative >= m:
        raise ValueError('stencil too small')
    rows, cols, vals = [], [], []
    import math
    for i in range(n):
        first = min(max(0, i-m//2), n-m)
        idx = np.arange(first, first+m)
        scale = np.max(abs(x[idx]-x[i]))
        t = (x[idx]-x[i])/scale
        rhs = np.zeros(m)
        rhs[derivative] = math.factorial(derivative)/scale**derivative
        w = np.linalg.solve(t[None, :]**np.arange(m)[:, None], rhs)
        rows.extend([i]*m); cols.extend(idx); vals.extend(w)
    return sparse.csr_matrix((vals, (rows, cols)), shape=(n, n))


@dataclass
class SolverOptions:
    max_iterations: int = 30
    residual_tolerance: float = 1e-8
    constraint_tolerance: float = 1e-7
    gmres_rtol: float = 1e-4
    gmres_restart: int = 60
    gmres_cycles: int = 8
    min_step: float = 2**-15
    method: str = 'newton'
    require_blowup_branch: bool = True


class ProfileProblem:
    def __init__(self, x, y, psi_boundary, stencil=5, inflow_mask=None,
                 inflow_values=None, rhs_forcing=None, green_reference=None):
        self.x, self.y = np.asarray(x), np.asarray(y)
        self.nx, self.ny = len(x), len(y)
        self.shape = (self.ny, self.nx)
        self.N = self.nx*self.ny
        self.X, self.Y = np.meshgrid(x, y)
        if abs(y[0]) > 1e-13 or min(abs(x)) > 1e-13:
            raise ValueError('grid must contain the wall origin')
        self.origin = int(np.argmin(abs(x)))
        ix, iy = sparse.eye(self.nx), sparse.eye(self.ny)
        self.Dx = sparse.kron(iy, fd_matrix(x, 1, stencil), format='csr')
        self.Dy = sparse.kron(fd_matrix(y, 1, stencil), ix, format='csr')
        # Separate elliptic order, intentionally independently checkable.
        self.Lap = sparse.kron(iy, fd_matrix(x, 2, stencil), format='csr') + sparse.kron(fd_matrix(y, 2, stencil), ix, format='csr')
        self.outer = np.zeros(self.shape, bool)
        self.outer[[0,-1],:] = True; self.outer[:,[0,-1]] = True
        self.interior = np.flatnonzero(~self.outer.ravel())
        self.boundary = np.flatnonzero(self.outer.ravel())
        self.A = (-self.Lap[self.interior][:, self.interior]).tocsc()
        self.lu = splu(self.A)
        p = np.asarray(psi_boundary, float).reshape(-1).copy()
        p[self.interior] = 0
        if np.max(abs(p[:self.nx])) > 1e-10:
            raise ValueError('psi must vanish on physical wall')
        self.psi_boundary = p
        self.boundary_rhs = (self.Lap @ p)[self.interior]
        self.green_kernel = None
        self.tail_closure = None
        self.tail_matching = np.zeros(self.N,bool)
        if green_reference is not None:
            # Exact linear response of the retained midpoint source quadrature;
            # the exterior source is held at its supplied reference value.
            bx = self.X.ravel()[self.boundary,None]
            by = self.Y.ravel()[self.boundary,None]
            sx = self.X.ravel()[None,self.interior]
            sy = self.Y.ravel()[None,self.interior]
            wx = np.empty(self.nx); wy = np.empty(self.ny)
            wx[1:-1] = (self.x[2:]-self.x[:-2])/2
            wy[1:-1] = (self.y[2:]-self.y[:-2])/2
            wx[[0,-1]] = 0; wy[[0,-1]] = 0
            weights = (wy[:,None]*wx[None,:]).ravel()[self.interior]
            direct = (bx-sx)**2+(by-sy)**2
            self.green_kernel = np.log1p(4*by*sy/direct)/(4*np.pi)*weights
            local_seed = self.green_kernel@(self.Dx@np.asarray(green_reference).ravel())[self.interior]
            self.psi_boundary[self.boundary] -= local_seed
        self.C = sparse.vstack([
            sparse.csr_matrix(([1.], ([0], [self.origin])), shape=(1,self.N)),
            self.Dx.getrow(self.origin),
            sparse.kron(sparse.csr_matrix(([1.],([0],[0])),shape=(1,self.ny)),
                        fd_matrix(x,2,max(stencil,7)).getrow(self.origin),format='csr'),
            sparse.kron(sparse.csr_matrix(([1.],([0],[0])),shape=(1,self.ny)),
                        fd_matrix(x,3,max(stencil,7)).getrow(self.origin),format='csr')
        ], format='csr')
        self.target = np.array([0.,1.,0.,-2.])
        self.constraint_scale = 1/np.maximum(np.sqrt(np.asarray(self.C.multiply(self.C).sum(axis=1))).ravel(),1.)
        self.inflow = np.zeros(self.N, bool) if inflow_mask is None else np.asarray(inflow_mask,bool).ravel()
        self.inflow_values = np.zeros(self.N) if inflow_values is None else np.asarray(inflow_values).ravel()
        if self.inflow[self.origin]:
            raise ValueError('origin cannot be an artificial inflow boundary')
        self.forcing = np.zeros(self.N) if rhs_forcing is None else np.asarray(rhs_forcing).ravel()

    def velocity(self, q, homogeneous=False, boundary_addition=None):
        source = self.Dx @ q
        p = np.zeros(self.N) if homogeneous else self.psi_boundary.copy()
        if self.green_kernel is not None:
            p[self.boundary] += self.green_kernel@source[self.interior]
        if boundary_addition is not None:
            p[self.boundary] += boundary_addition
        rhs = source[self.interior]+(self.Lap@p)[self.interior]
        p[self.interior] = self.lu.solve(rhs)
        ux = -self.Dy @ p
        uy = self.Dx @ p
        # Dirichlet psi=0 is exactly impermeable; do not impose Q at wall.
        uy[:self.nx] = 0
        return ux, uy, p

    def evaluate(self, z):
        q = z[:self.N]
        cl, cw, drift, offset = z[self.N:]
        tail = None if self.tail_closure is None else self.tail_closure.evaluate(q,z[self.N:])
        ux, uy, psi = self.velocity(q,boundary_addition=tail)
        vx = ux-ux[self.origin]+cl*self.X.ravel()+drift
        vy = uy+cl*self.Y.ravel()
        qx, qy = self.Dx@q, self.Dy@q
        raw = vx*qx+vy*qy-cw*q-offset-self.forcing
        r = raw.copy()
        r[self.inflow] = q[self.inflow]-self.inflow_values[self.inflow]
        matching = np.zeros(self.N)
        if self.tail_matching.any():
            matching = self.X.ravel()*qx+self.Y.ravel()*qy-(cw/cl)*(q+offset/cw)
            r[self.tail_matching] = matching[self.tail_matching]
        gauges = self.C@q-self.target
        residual = np.r_[r, gauges*self.constraint_scale]
        cache = dict(q=q, qx=qx, qy=qy, ux=ux, uy=uy, psi=psi, vx=vx,
                     vy=vy, cl=cl, cw=cw, drift=drift, offset=offset,
                     raw=raw, gauges=gauges, matching=matching)
        return residual, cache

    def jvp(self, dz, cache):
        h = dz[:self.N]
        dl, dw, dv, db = dz[self.N:]
        tail = None if self.tail_closure is None else self.tail_closure.jvp(cache['q'],
            np.array([cache['cl'],cache['cw'],cache['drift'],cache['offset']]),h,dz[self.N:])
        hx, hy, _ = self.velocity(h, homogeneous=True,boundary_addition=tail)
        jh = cache['vx']*(self.Dx@h)+cache['vy']*(self.Dy@h)-cache['cw']*h
        jh += (hx-hx[self.origin])*cache['qx']+hy*cache['qy']
        jh += dl*(self.X.ravel()*cache['qx']+self.Y.ravel()*cache['qy'])
        jh += -dw*cache['q']+dv*cache['qx']-db
        jh[self.inflow] = h[self.inflow]
        if self.tail_matching.any():
            alpha=cache['cw']/cache['cl']
            dalpha=(dw-alpha*dl)/cache['cl']
            offset=cache['offset']/cache['cw']
            doffset=(db-offset*dw)/cache['cw']
            robin=self.X.ravel()*(self.Dx@h)+self.Y.ravel()*(self.Dy@h)-alpha*(h+doffset)-dalpha*(cache['q']+offset)
            jh[self.tail_matching]=robin[self.tail_matching]
        return np.r_[jh, (self.C@h)*self.constraint_scale]

    def frozen_matrix(self, c, shift=0.):
        """Bordered Picard operator: omit only the nonlocal velocity response."""
        top = sparse.diags(c['vx'])@self.Dx+sparse.diags(c['vy'])@self.Dy-sparse.eye(self.N)*(c['cw']-shift)
        cols = np.column_stack([self.X.ravel()*c['qx']+self.Y.ravel()*c['qy'],
                                -c['q'], c['qx'], -np.ones(self.N)])
        if self.inflow.any():
            top = top.tolil()
            for i in np.flatnonzero(self.inflow):
                top.rows[i] = [i]; top.data[i] = [1.]
            top = top.tocsr(); cols[self.inflow,:] = 0
        if self.tail_matching.any():
            alpha=c['cw']/c['cl']
            robin=sparse.diags(self.X.ravel())@self.Dx+sparse.diags(self.Y.ravel())@self.Dy-alpha*sparse.eye(self.N)
            top=top.tolil();robin=robin.tocsr()
            for i in np.flatnonzero(self.tail_matching):
                rr=robin.getrow(i)
                top.rows[i]=rr.indices.tolist();top.data[i]=rr.data.tolist()
            top=top.tocsr()
            ii=self.tail_matching
            cols[ii,0]=(c['cw']*c['q'][ii]+c['offset'])/c['cl']**2
            cols[ii,1]=-c['q'][ii]/c['cl']
            cols[ii,2]=0
            cols[ii,3]=-1/c['cl']
        return sparse.bmat([[top, sparse.csr_matrix(cols)],
            [sparse.diags(self.constraint_scale)@self.C, sparse.csr_matrix((4,4))]],format='csc')

    def normalize_seed(self, q):
        """Enforce only the four gauges by localized polynomial corrections."""
        q = np.asarray(q).reshape(-1).copy()
        xx, yy = self.X.ravel(), self.Y.ravel()
        envelope = np.exp(-(xx/1.5)**8-(yy/1.5)**4)
        basis = np.column_stack([envelope*xx**k for k in range(4)])
        mat = self.C@basis
        q += basis@np.linalg.solve(mat,self.target-self.C@q)
        return q

    def metrics(self, z, c=None):
        if c is None: _,c = self.evaluate(z)
        core = (abs(self.X.ravel())<=2)&(self.Y.ravel()<=2)
        scale = max(np.max(abs(c['vx']*c['qx'])+abs(c['vy']*c['qy'])+abs(c['cw']*c['q'])+abs(c['offset'])),1.)
        boundary = c['vx'].reshape(self.shape), c['vy'].reshape(self.shape)
        outward = np.r_[-boundary[0][:,0],boundary[0][:,-1],boundary[1][-1,:]]
        poisson_r = (-self.Lap@c['psi']-self.Dx@c['q'])[self.interior]
        return dict(residual_inf=float(np.max(abs(c['raw']))),
            active_pde_residual_inf=float(np.max(abs(c['raw'][~(self.inflow|self.tail_matching)]))),
            tail_matching_residual_inf=float(np.max(abs(c['matching'][self.tail_matching]))) if self.tail_matching.any() else 0.,
            inflow_data_residual_inf=float(np.max(abs(c['q'][self.inflow]-self.inflow_values[self.inflow]))) if self.inflow.any() else 0.,
            residual_relative_inf=float(np.max(abs(c['raw']))/scale),
            core_residual_inf=float(np.max(abs(c['raw'][core]))),
            wall_residual_inf=float(np.max(abs(c['raw'][:self.nx]))),
            residual_rms=float(np.sqrt(np.mean(c['raw']**2))),
            gradient_residual_inf=float(max(np.max(abs(self.Dx@c['raw'])),np.max(abs(self.Dy@c['raw'])))),
            constraints=[float(t) for t in c['gauges']],
            poisson_residual_inf=float(np.max(abs(poisson_r))),
            cl=float(c['cl']),cw=float(c['cw']),drift=float(c['drift']),offset=float(c['offset']),
            alpha=float(c['cw']/c['cl']) if c['cl'] != 0 else None,
            kappa=float(c['cl']-c['cw']),offset_drift_defect=float(c['offset']-c['drift']),
            minimum_outward_velocity=float(np.min(outward)),
            inflow_boundary_count=int(np.sum(outward < -1e-8)),
            max_abs_Q=float(np.max(abs(c['q']))),max_abs_Qx=float(np.max(abs(c['qx']))))


def solve(problem, z0, options=None, callback=None):
    """z <- z - alpha B^{-1}H, with exact-Jacobian Newton or frozen Picard B."""
    opts = options or SolverOptions()
    z = np.asarray(z0,float).copy()
    history = []
    status = 'maximum_iterations'
    for it in range(opts.max_iterations+1):
        r,c = problem.evaluate(z)
        metrics = problem.metrics(z,c)
        merit = float(np.linalg.norm(r)/np.sqrt(r.size))
        entry = dict(iteration=it,merit=merit,**metrics)
        history.append(entry)
        if callback: callback(entry)
        gauge_inf = np.max(abs(c['gauges']))
        branch_ok = not opts.require_blowup_branch or c['cl']>c['cw']>0
        if metrics['active_pde_residual_inf'] < opts.residual_tolerance and metrics['inflow_data_residual_inf'] < opts.residual_tolerance and metrics['tail_matching_residual_inf'] < opts.residual_tolerance and gauge_inf < opts.constraint_tolerance and branch_ok:
            status = 'matched_tail_root' if problem.tail_matching.any() else 'discrete_root'; break
        if it == opts.max_iterations: break
        pre = None
        for shift in (0., 0.05, 0.5, 2.):
            try:
                pre = splu(problem.frozen_matrix(c,shift))
                break
            except RuntimeError: pass
        if pre is None:
            status = 'preconditioner_failed'; break
        if opts.method == 'picard':
            step = pre.solve(-r)
            info = 0
        elif opts.method == 'newton':
            J = LinearOperator((z.size,z.size),matvec=lambda h:problem.jvp(h,c))
            M = LinearOperator((z.size,z.size),matvec=pre.solve)
            linear_hist=[]
            step,info = gmres(J,-r,M=M,rtol=opts.gmres_rtol,atol=0.,
                restart=opts.gmres_restart,maxiter=opts.gmres_cycles,
                callback=linear_hist.append,callback_type='pr_norm')
            entry['krylov_iterations'] = len(linear_hist)
            entry['linear_residual_relative'] = float(np.linalg.norm(problem.jvp(step,c)+r)/max(np.linalg.norm(r),1e-300))
        else: raise ValueError('method must be newton or picard')
        entry['linear_info']=int(info)
        if not np.isfinite(step).all():
            status='nonfinite_step'; break
        alpha=1.
        accepted=False
        while alpha>=opts.min_step:
            trial=z+alpha*step
            cl,cw=trial[problem.N:problem.N+2]
            if opts.require_blowup_branch and not (cl>cw>0):
                alpha*=.5; continue
            rt,ct=problem.evaluate(trial)
            # A global least-squares decrease alone can hide a bad local defect.
            if np.isfinite(rt).all() and np.linalg.norm(rt)<(1-1e-4*alpha)*np.linalg.norm(r) and np.max(abs(rt[:problem.N])) <= (1+0.05*alpha)*max(np.max(abs(r[:problem.N])),opts.residual_tolerance):
                accepted=True; break
            alpha*=.5
        entry['step_length']=alpha
        if not accepted:
            status='line_search_failed'; break
        z=trial
    metrics=problem.metrics(z)
    # A box root is never automatically an infinite-domain profile candidate.
    return dict(z=z,history=history,status=status,metrics=metrics,
                infinite_domain_candidate=False,
                box_root=(status in ('discrete_root','matched_tail_root')),
                boundary_inflow_admissible=metrics['inflow_boundary_count']==0 or bool(problem.tail_matching.any()))
