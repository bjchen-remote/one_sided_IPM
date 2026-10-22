"""Exact finite algebraic-relation gate for the normalized cubic benchmark.

Set w=s*y/(1+s*s)**2, h=y*y/(1+s*s)**4, Z=(X-F0)/y-k.
The verified parity pattern writes Z as a weighted series in w,h.
Weights: w=1, h=2, Z=1. No numerical tolerance or fitted data enter.
"""
import json
from pathlib import Path
import sympy as sp

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT/'runs/normal_lift_20260928'


def gate():
    data = json.loads((OUT/'cubic_exact_jets.json').read_text())
    s = sp.symbols('s',real=True)
    order = len(data['coefficients'])-1
    Z = {}
    for row in data['coefficients'][1:]:
        n = row['n']
        expr = sp.sympify(row['a_n'],locals={'s':s})
        poly = sp.Poly(sp.cancel(expr*(1+s*s)**(2*n-2)),s)
        for (i,), value in poly.terms():
            assert (n-1-i)%2==0
            Z[(i,(n-1-i)//2)] = value
    def multiply(a,b):
        out = {}
        for (i,j),x in a.items():
            for (k,l),y in b.items():
                if i+k+2*(j+l)<=order:
                    pair = (i+k,j+l)
                    out[pair] = out.get(pair,0)+x*y
        return {pair:value for pair,value in out.items() if value}
    powers = [{(0,0):sp.Integer(1)}]
    for _ in range(4):
        powers.append(multiply(powers[-1],Z))
    rows = [(i,j) for j in range(order//2+1) for i in range(order-2*j+1)]
    results = []
    for degree in range(1,5):
        monomials = [(i,j,k) for j in range(degree//2+1)
                     for i in range(degree-2*j+1)
                     for k in range(degree-2*j-i+1)]
        columns = []
        for i,j,k in monomials:
            column = {(u+i,v+j):a for (u,v),a in powers[k].items() if u+i+2*(v+j)<=order}
            columns.append([column.get(row,0) for row in rows])
        # Full rank modulo a prime is an exact full-rank certificate over Q.
        # Avoid expensive fraction growth in generic symbolic matrix.rank().
        prime = 1000003
        matrix = [[int(columns[j][i].p) * pow(int(columns[j][i].q),-1,prime)%prime
                   if isinstance(columns[j][i],sp.Rational) else int(columns[j][i])%prime
                   for j in range(len(columns))] for i in range(len(rows))]
        rank = 0
        for col in range(len(columns)):
            pivot = next((i for i in range(rank,len(rows)) if matrix[i][col]),None)
            if pivot is None:
                continue
            matrix[rank],matrix[pivot] = matrix[pivot],matrix[rank]
            scale = pow(matrix[rank][col],-1,prime)
            matrix[rank] = [v*scale%prime for v in matrix[rank]]
            for i in range(rank+1,len(rows)):
                factor = matrix[i][col]
                matrix[i] = [(x-factor*y)%prime for x,y in zip(matrix[i],matrix[rank])]
            rank += 1
        results.append(dict(weightedDegree=degree,unknowns=len(monomials),
            conditions=len(rows),modularRank=rank,prime=prime,
            fullColumnRankCertified=rank==len(monomials)))
        print(results[-1],flush=True)
    result = dict(k=data['k'],knownNormalOrder=order+1,
        variables=dict(w='s*y/(1+s^2)^2',h='y^2/(1+s^2)^4',Z='(X-F0)/y-k'),
        weights=dict(w=1,h=2,Z=1),gates=results,
        scope='No nonzero polynomial relation within the stated weighted degree for this k=1 cubic-wall lift. Does not exclude other walls, higher degrees, nonpolynomial coefficients, or differential/special-function closure.')
    (OUT/'algebraic_gate.json').write_text(json.dumps(result,indent=2)+'\n')


if __name__=='__main__':
    gate()
