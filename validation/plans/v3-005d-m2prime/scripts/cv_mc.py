import numpy as np, math
rng=np.random.default_rng(12345)
def formula(x):  # x=E_F/E
    y=1/x  # E/E_F
    f=1-1.4*x
    if y<2: f+= 0.4*x*(2-y)**2.5
    return f if y>1 else 0.0
def mc(x,n=4_000_000):
    k0=np.array([0,0,1.0]); kf=math.sqrt(x)
    u=rng.random(n); d=rng.normal(size=(n,3)); d/=np.linalg.norm(d,axis=1)[:,None]
    k2=d*kf*u[:,None]**(1/3)
    P=(k0+k2)/2; q=np.linalg.norm(k0-k2,axis=1)/2
    e=rng.normal(size=(n,3)); e/=np.linalg.norm(e,axis=1)[:,None]
    k1p=P+q[:,None]*e; k2p=P-q[:,None]*e
    ok=(np.linalg.norm(k1p,axis=1)>kf)&(np.linalg.norm(k2p,axis=1)>kf)
    vrel=2*q  # |k0-k2| / |k0| with k0=1
    return float(np.mean(vrel*ok)), float(np.mean(ok)), float(np.mean(vrel))
print(" E/EF   formula   MC<vrel*Q>/v0  MC<Q>   <vrel>/v0")
for y in (1.0,1.05,1.2,1.5,1.8,2.0,2.5,3.0,4.0,6.0):
    a,b,c=mc(1/y); print(f"{y:5.2f}  {formula(1/y):.4f}   {a:.4f}        {b:.4f}  {c:.4f}")
