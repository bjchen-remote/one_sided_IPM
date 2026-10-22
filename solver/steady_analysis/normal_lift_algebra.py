"""Truncated exact polynomial arithmetic for the full hodograph curl identity."""


def curl_coefficients(j, js, jss, k, ks, ky, s, lam, mu, order, simplify=lambda x:x):
    """E(y) modulo y**(order+1), arrays may hold symbolic or NumPy values.

    j = X_s, k = X_y. E/J**3 is the physical curl defect.
    """
    n=order+1
    def pad(a):
        return list(a[:n])+[0]*max(0,n-len(a))
    def add(*arrays):
        return [simplify(sum(pad(a)[i] for a in arrays)) for i in range(n)]
    def mul(*arrays):
        out=[1]+[0]*order
        for a in arrays:
            a=pad(a)
            out=[simplify(sum(out[z]*a[i-z] for z in range(i+1))) for i in range(n)]
        return out
    def scale(a,b):
        return [simplify(b*v) for v in pad(a)]
    j,js,jss,k,ks,ky=map(pad,(j,js,jss,k,ks,ky))
    nu=2*lam-mu
    f=[0]+[(nu*j[i]-mu*s*js[i])/(i+1) for i in range(order)]
    fs=[0]+[((nu-mu)*js[i]-mu*s*jss[i])/(i+1) for i in range(order)]
    one_k2=add([1],mul(k,k))
    return add(mul(j,one_k2,fs),
               mul(f,add(scale(mul(j,k,ks),2),scale(mul(one_k2,js),-1),scale(mul(j,j,ky),-1))),
               scale(mul(j,j,j,k),-2*(lam-mu)),
               scale(mul(j,j,k,js),2*mu*s),
               scale(mul(j,j,j,ks),-mu*s),mul(j,j))
