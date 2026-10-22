"""Check the six highest-height-derivative columns of the potential reduction.

Finite-interval numerical rank and derivative checks only; these do not show
that the lower-order transport residual belongs to this six-function span.
"""
from pathlib import Path
import json
import numpy as np
from scipy.integrate import cumulative_simpson
from scipy.interpolate import CubicSpline
from streamfunction_shape_fit import mother, BETA

ROOT=Path(__file__).resolve().parents[2]
OUT=ROOT/'runs/streamfunction_reduction_20260928'
GEOMETRY=ROOT/'runs/streamfunction_geometry_20260928/shape_fit.json'


class Primitives:
    def __init__(self, parameters, nodes=4001):
        self.parameters=np.asarray(parameters)
        xi=np.linspace(-2,2,nodes); midpoint=(nodes-1)//2
        s=np.sinh(BETA*xi); c=np.cosh(BETA*xi); g=mother(xi,self.parameters)
        def primitive(f):
            v=cumulative_simpson(f,x=xi,initial=0)
            return v-v[midpoint]
        j0=primitive(g*c); j1=primitive(g*s)
        self.J0=CubicSpline(xi,j0,extrapolate=False)
        self.J1=CubicSpline(xi,j1,extrapolate=False)
        self.K0=CubicSpline(xi,primitive(j0*c),extrapolate=False)
        self.K1=CubicSpline(xi,primitive(j0*s+j1*c),extrapolate=False)
        self.K2=CubicSpline(xi,primitive(j1*s),extrapolate=False)

    def values(self,xi,k):
        s=np.sinh(BETA*xi); c=np.cosh(BETA*xi)
        G=(s+k*(c-1))/BETA; Gp=c+k*s; Gk=(c-1)/BETA
        J=self.J0(xi)+k*self.J1(xi)
        K=self.K0(xi)+k*self.K1(xi)+k*k*self.K2(xi)
        Kk=self.K1(xi)+2*k*self.K2(xi)
        return G,Gp,Gk,J,K,Kk

    def phi(self,x,y,theta,chi=0.):
        c,B,n,W,k=theta
        xi=(np.arcsinh((BETA*(x-n)/W+k)/np.sqrt(1-k*k))-np.arctanh(k))/BETA
        G,Gp,Gk,J,K,Kk=self.values(xi,k)
        return chi+y*(c*W*G+B*W*W*K)

    def columns(self,xi,y,theta):
        c,B,n,W,k=theta
        G,Gp,Gk,J,K,Kk=self.values(xi,k);H=c+B*W*J
        return np.column_stack([y*W*G,y*W*W*K,-y*H,
            y*B*W*(2*K-J*G),y*B*W*W*(Kk-J*Gk),np.ones_like(xi)])


def projection(columns,rhs,weights=None):
    """Solve projected highest derivatives while retaining the full remainder.

    `columns` multiplies (c''',B''',n''',W''',k''',chi''').
    Return coefficients and rhs-columns@coefficients, never discard the latter.
    """
    weights=np.ones(len(rhs)) if weights is None else np.asarray(weights)
    if np.any(weights<=0):
        raise ValueError('Projection weights must be positive')
    root=np.sqrt(weights); weighted=columns*root[:,None]
    scale=np.linalg.norm(weighted,axis=0)
    if np.any(scale==0):
        raise ValueError('Degenerate highest-derivative column; use the unsolved DAE')
    coeff,_,rank,singular=np.linalg.lstsq(weighted/scale,root*rhs,rcond=None)
    coeff/=scale
    return dict(coefficients=coeff,remainder=rhs-columns@coeff,rank=int(rank),singular=singular,
                columnScale=scale)


def main():
    data=json.loads(GEOMETRY.read_text());xi=np.linspace(-1.5,1.5,601)
    records=[]; checks=[]
    for family in ('student','student_skew_even'):
        fit=next(z for z in data['fits'] if z['coordinate']=='sinh' and z['family']==family)
        primitive=Primitives(fit['parameters'])
        for geo in data['geometry']:
            B,W,k,n,y=(geo[s] for s in ('B','W','k','n','yInner'))
            xiz=(np.arcsinh((BETA*(geo['s']-n)/W+k)/np.sqrt(1-k*k))-np.arctanh(k))/BETA
            c=-B*W*primitive.values(np.array(xiz),k)[3]
            theta=np.array([c,B,n,W,k]); columns=primitive.columns(xi,y,theta)
            g=mother(xi,fit['parameters']); lam=.60735924
            A=y*(lam+B*g[:,None].ravel())[:,None]*columns
            norms=np.linalg.norm(A,axis=0); singular=np.linalg.svd(A/norms,compute_uv=False)
            records.append(dict(family=family,Ynative=geo['Ynative'],
                normalizedRank=int(np.sum(singular>singular[0]*1e-10)),
                singularValues=singular.tolist(),normalizedCondition=float(singular[0]/singular[-1]),
                V2Minimum=float(np.min(y*(lam+B*g))),includesVelocityPeak=True))
        # Check fixed-x parameter derivatives on a moderate synthetic geometry,
        # independent of observed height derivatives or fitted transport.
        theta=np.array([-.17,.43,-.31,1.28,-.61]); y=.7
        x=theta[2]+theta[3]*primitive.values(xi,theta[4])[0]
        analytic=primitive.columns(xi,y,theta)
        errors=[]
        for i in range(5):
            h=1e-5*max(1.,abs(theta[i])); offset=np.eye(5)[i]*h
            derivative=(primitive.phi(x,y,theta+offset)-primitive.phi(x,y,theta-offset))/(2*h)
            errors.append(float(np.linalg.norm(derivative-analytic[:,i])/np.linalg.norm(analytic[:,i])))
        eps=1e-5
        px=(primitive.phi(x+eps,y,theta)-primitive.phi(x-eps,y,theta))/(2*eps)
        H=theta[0]+theta[1]*theta[3]*primitive.values(xi,theta[4])[3]
        phi_x_error=float(np.linalg.norm(px-y*H)/np.linalg.norm(y*H))
        # Check quadrature refinement of the actual matrix, not PDE convergence.
        refined=Primitives(fit['parameters'],8001).columns(xi,y,theta)
        refinement=float(np.linalg.norm(refined-analytic)/np.linalg.norm(refined))
        assert max(errors)<2e-6 and phi_x_error<2e-6 and refinement<2e-7
        checks.append(dict(family=family,fixedXParameterDerivativeRelativeErrors=errors,
                           phiXEqualsPsiRelativeError=phi_x_error,quadratureRefinement=refinement))
    out=dict(parameterOrder=['c','B','n','W','k','chi'],
        potential='Phi=chi+y[c W G+B W^2 K], J_xi=g G_xi, K_xi=J G_xi',
        residual='E=V dot grad(Delta Phi)-mu Delta Phi+v; highest derivative coefficient is +V2 Phi_parameter',
        highestColumns='V2*(yWG, yW^2K, -yH, yBW(2K-JG), yBW^2(K_k-JG_k), 1)',
        xi=xi.tolist(),checks=checks,rankDiagnostics=records,
        interpretation='Finite-interval rank and differentiation checks. No transport zero, invariant manifold, admissible global tail, or integration of the reduced height ODE is claimed.')
    OUT.mkdir(exist_ok=True)
    (OUT/'basis_checks.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(dict(checks=checks,conditionRange=[min(v['normalizedCondition'] for v in records),max(v['normalizedCondition'] for v in records)]),indent=2))


if __name__=='__main__':
    main()
