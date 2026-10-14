import numpy as np, math, sys
from numpy import trapz
d=np.load(sys.argv[1]); E=d['grid_e_mev']; Z=list(d['target_z'])
tO=Z.index(8); tH=Z.index(1)
from math import pi
def j1(x):
    t=np.linspace(0,pi,401)
    return trapz(np.cos(t[None,:]-x[:,None]*np.sin(t[None,:])),t,axis=1)/pi
alpha,p=0.0022,1.77  # Bragg-Kleeman water, cm, MeV
Rng=lambda e: alpha*np.maximum(e,0)**p
nO=3.343e22; nH=2*nO
mu=np.linspace(-1,1,4001)
M=d['target_mass_mev'][tO]; mp=938.272
def shape(e,Rfac=1.0,s=0.0,gauss=False):
    i=np.argmin(abs(E-e)); k=d['p_cm_mev'][tO,i]/197.327; R=d['radius_fm'][tO,i]*Rfac
    q=2*k*np.sqrt((1-mu)/2)
    if gauss: w=np.exp(-(q*R)**2/4)
    else:
        x=np.maximum(q*R,1e-9); w=(2*j1(x)/x)**2*np.exp(-(q*s)**2)
    return w/trapz(w,mu), k
def run(E0, w_mm=(1,2,3,5), **kw):
    R0=Rng(E0); zs=np.linspace(0,R0*0.98,120); res={}
    acc={'P':0,'tr':0,'cos':0}; accw={x:0 for x in w_mm}
    for j in range(len(zs)-1):
        z=0.5*(zs[j]+zs[j+1]); dz=zs[j+1]-zs[j]
        e=(  (R0-z)/alpha)**(1/p)
        if e<7: continue
        i=np.argmin(abs(E-e)); sig=d['sigma_barn'][tO,i]*1e-24
        if sig<=0: continue
        wgt,k=shape(e,**kw)
        # CM->lab for elastic p+A: proton lab angle and energy
        g=mp/M
        thl=np.arctan2(np.sqrt(1-mu**2), mu+g)
        # recoil kinetic energy (nonrel approx): T_rec = 2*mp*M/(mp+M)^2 * e*(1-mu)
        Trec=2*mp*M/(mp+M)**2*e*(1-mu)
        ep=e-Trec
        Rres=Rng(e);
        cosl=np.cos(thl)
        defic=(Rres-Rng(ep)*np.clip(cosl,0,None))*10  # mm
        dP=nO*sig*dz
        acc['P']+=dP
        acc['tr']+=dP*trapz(wgt*(1-cosl),mu)
        for x in w_mm: accw[x]+=dP*trapz(wgt*(defic>x),mu)
    return acc,accw
for E0 in (150,200):
  for lab,kw in [('table sharp R',{}),('R x1.2',{'Rfac':1.2}),('s=0.5',{'s':0.5}),('s=0.9',{'s':0.9}),('gauss',{'gauss':True})]:
    acc,accw=run(E0,**kw)
    print(E0,f"{lab:14s} P_el(O)={acc['P']:.4f} sum(Sig*<1-cos_lab>)={acc['tr']:.5f}  removal P(deficit>w): "+" ".join(f"{x}mm:{accw[x]:.4f}" for x in accw))
