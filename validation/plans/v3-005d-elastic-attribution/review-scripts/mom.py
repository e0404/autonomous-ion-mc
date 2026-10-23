import numpy as np, math, sys
d=np.load(sys.argv[1])
E=d['grid_e_mev']; Z=d['target_z']; A=d['target_a']
ed=d['edges_mu']
for t in range(len(Z)):
  if Z[t] not in (1,6,8): continue
  for e in (50,100,150,200,250):
    i=np.argmin(abs(E-e))
    edges=ed[t,i]
    m1=[];th2=[];fr=[]
    # uniform in mu within each bin
    a=edges[:-1]; b=edges[1:]
    m=np.mean(1-(a+b)/2)
    # theta^2: integrate numerically
    tot=0; tt=0
    for lo,hi in zip(a,b):
        mu=np.linspace(lo,hi,41); th=np.arccos(np.clip(mu,-1,1))
        tt+=np.trapz(th**2,mu)/(hi-lo)
    th_rms=math.degrees(math.sqrt(tt/len(a)))
    # fractions
    th_edges=np.degrees(np.arccos(edges))
    def frac_gt(deg):
        c=math.cos(math.radians(deg)); return np.interp(c,edges,np.linspace(0,1,len(edges)))
    print(f"Z={Z[t]} A={A[t]} E={E[i]:.1f} sig={d['sigma_barn'][t,i]*1000:.1f}mb R={d['radius_fm'][t,i]:.3f} <1-mu>={m:.4f} th_rms={th_rms:.2f} P(thCM<5)={1-frac_gt(5):.3f} P(<10)={1-frac_gt(10):.3f} P(<20)={1-frac_gt(20):.3f} median={np.degrees(np.arccos(np.interp(0.5,np.linspace(0,1,len(edges)),edges))):.1f}")
