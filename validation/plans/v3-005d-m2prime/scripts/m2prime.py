# Independent optical-limit Glauber check (numpy only). Not repository code.
import numpy as np, math, sys
HBARC=197.3269804; MP=938.27208816; RP=0.8414
def j0(x):  # integral representation, trapezoid on periodic integrand
    t=np.linspace(0,math.pi,161)
    w=np.full(t.size,t[1]-t[0]); w[0]=w[-1]=w[0]/2
    x=np.asarray(x,float)
    out=np.zeros(x.shape)
    for i in range(t.size):
        out+=w[i]*np.cos(x*math.sin(t[i]))
    return out/math.pi
def ho_rho(r,a,al): return (1+al*(r/a)**2)*np.exp(-(r/a)**2)
r=np.linspace(0,20,4001); dr=r[1]-r[0]
def ff_num(rho,q):
    wr=rho*r*r; norm=np.trapezoid(wr,r)
    return np.array([np.trapezoid(rho*r*np.sin(qq*r)/qq,r)/norm if qq>0 else 1.0 for qq in q])
def pcm(T,M):
    s=(MP+M)**2+2*M*T; plab=math.sqrt(T*T+2*T*MP); return plab*M/math.sqrt(s)
def run(rho, Z, N, M, T, nn, unfold=True, label=""):
    q=QG
    F=ff_num(rho,q)
    if unfold: F=F*np.exp(q*q*RP*RP/6)
    b=BG; K=KG
    chi=np.zeros(b.size,complex)
    for (n,sig,al,be) in ((Z,)+nn[0],(N,)+nn[1]):
        env=F*np.exp(-be*q*q/2)
        C=(K*(q*env)).dot(np.full(q.size,q[1]-q[0]))/(2*math.pi)
        chi+=n*sig*(1j+al)/2*C
    k=pcm(T,M)/HBARC
    S=np.exp(1j*chi); wb=np.full(b.size,b[1]-b[0]); wb[0]=wb[-1]=wb[0]/2
    sR=2*math.pi*np.sum(wb*b*(1-abs(S)**2))
    sT=4*math.pi*np.sum(wb*b*(1-S).real)
    th=np.linspace(0,math.pi,1501); qq=2*k*np.sin(th/2)
    Fa=1j*k*(j0(np.outer(qq,b)).dot(wb*b*(1-S)))
    sEl=2*math.pi*np.trapezoid(abs(Fa)**2*np.sin(th),th)
    sTot_ot=4*math.pi/k*Fa[0].imag
    # first min
    d=abs(Fa)**2; i=next((i for i in range(1,d.size-1) if d[i]<d[i-1] and d[i]<=d[i+1]),None)
    qmin=qq[i] if i else float('nan')
    return dict(T=T,lab=label,sR=10*sR,sEl=10*sEl,sTot=10*sT,sTot_ot=10*sTot_ot,qmin=qmin,k=k)
mb=0.1
QG=np.linspace(1e-6,10,801); BG=np.linspace(0,12,481); KG=j0(np.outer(BG,QG))
said={100:((31.64,1.843),(75.38,0.9676)),150:((25.72,1.623),(52.14,0.9435)),200:((24.27,1.304),(43.00,0.7838))}
def nnsaid(T):
    out=[]
    for s,a in said[T]:
        s*=mb; out.append((s,a,(1+a*a)*s/(16*math.pi)))
    return out

RHO0=0.17
BAR={'C-12':{100:253.0,150:125.0,200:84.0},'O-16':{100:352.0,150:159.0,200:102.0}}
sys.path.insert(0, __import__('os').path.dirname(__import__('os').path.abspath(__file__)))
from exfor_tot import measured_sigma_tot
TOT=measured_sigma_tot()  # EXFOR 13569/13753 read from the acquisition cache at run time
TH2={100:13.033388626513915,150:9.108035714119335,200:7.097961439193826}
def unf(a,al):
    c1=al*a*a/(2*(2+3*al)); a2=a*a-2*RP*RP/3; return math.sqrt(a2),4*c1/(a2-6*c1)
def ho_norm(r,a,al):  # int rho d3r = 1
    return (1+al*(r/a)**2)*np.exp(-(r/a)**2)/(math.pi**1.5*a**3*(1+1.5*al))
def ff_raw(rho_eff,q):  # int rho_eff j0(qr) d3r  (not renormalised)
    return np.array([4*math.pi*np.trapezoid(rho_eff*r*np.sin(qq*r)/qq,r) for qq in q])
def bdc_R(rho_tot,E,mult=1.0):
    rt=np.maximum(mult*rho_tot/RHO0,1e-30)**(2/3)
    return np.where(E>46.27*rt,1-37.02*rt/E,E/(231.38*rt))
def chi_of(Fq,Z,N,nn,alscale=1.0,lam=1.0):
    q=QG; chi=np.zeros(BG.size,complex)
    for (n,(sig,al,be)) in ((Z,nn[0]),(N,nn[1])):
        env=Fq*np.exp(-be*q*q/2)
        C=(KG*(q*env)).dot(np.full(q.size,q[1]-q[0]))/(2*math.pi)
        chi+=n*sig*(1j+al*alscale)/2*C
    return lam*chi
WB=np.full(BG.size,BG[1]-BG[0]); WB[0]=WB[-1]=WB[0]/2
def stot(chi): return 10*4*math.pi*np.sum(WB*BG*(1-np.exp(1j*chi)).real)
def full(chi,T,M,th2lab):
    k=pcm(T,M)/HBARC; S=np.exp(1j*chi)
    sR=10*2*math.pi*np.sum(WB*BG*(1-abs(S)**2)); sT=stot(chi)
    th=np.linspace(0,math.pi,1501); qq=2*k*np.sin(th/2)
    Fa=1j*k*(j0(np.outer(qq,BG)).dot(WB*BG*(1-S))); d=abs(Fa)**2
    g=d*np.sin(th); cum=np.concatenate(([0],np.cumsum(0.5*(g[1:]+g[:-1])*np.diff(th))))
    sEl=10*2*math.pi*cum[-1]
    # lab->cm for threshold
    plab=math.sqrt(T*T+2*T*MP); bcm=plab/(T+MP+M); gcm=1/math.sqrt(1-bcm*bcm)
    pst=pcm(T,M); bst=pst/math.sqrt(pst*pst+MP*MP)
    thl=np.arctan2(np.sin(th),gcm*(np.cos(th)+bcm/bst))
    thc=np.interp(math.radians(th2lab),thl,th)
    P=1-np.interp(thc,th,cum)/cum[-1]
    i=next((i for i in range(1,d.size-1) if d[i]<d[i-1] and d[i]<=d[i+1]),None)
    return dict(sR=sR,sEl=sEl,sT=sT,qmin=(qq[i] if i else float('nan')),P2=P)
def solve(f,target,lo,hi):
    flo=f(lo)-target
    for _ in range(50):
        m=(lo+hi)/2; fm=f(m)-target
        if fm*flo>0: lo,flo=m,fm
        else: hi=m
    return (lo+hi)/2
tg={'C-12':(6,6,1.692,1.082,11174.86),'O-16':(8,8,1.833,1.544,14895.08)}
for name,(Z,N,a,al,M) in tg.items():
    ap,alp=unf(a,al); rho=ho_norm(r,ap,alp); A=Z+N
    F0=ff_raw(rho,QG)
    print(f"{name}: point HO a'={ap:.4f} alpha'={alp:.4f} central rho_tot={A*rho[0]:.4f} fm^-3 (max {A*rho.max():.4f})")
    for T in (100,150,200):
        nn=nnsaid(T); rows=[]
        chiA=chi_of(F0,Z,N,nn); rows.append(('(a) M2',chiA,''))
        for mult,lab in ((1.0,'(b) BDC Pauli LDA rho~=rho/rho0'),(2.0,'(b2) BDC rho~=2rho/rho0 [alt reading]')):
            Fb=ff_raw(rho*bdc_R(A*rho,T,mult),QG)
            rows.append((lab,chi_of(Fb,Z,N,nn),''))
        sa=solve(lambda x: stot(chi_of(F0,Z,N,nn,alscale=x)),TOT[name][T],0.0,1.5)
        rows.append(('(c) alpha x %.3f'%sa,chi_of(F0,Z,N,nn,alscale=sa),''))
        lam=solve(lambda x: stot(chi_of(F0,Z,N,nn,lam=x)),TOT[name][T],0.2,1.5)
        rows.append(('(d) chi x %.3f'%lam,chi_of(F0,Z,N,nn,lam=lam),''))
        Fb=ff_raw(rho*bdc_R(A*rho,T,1.0),QG)
        sb=solve(lambda x: stot(chi_of(Fb,Z,N,nn,alscale=x)),TOT[name][T],0.0,1.5)
        rows.append(('(e) BDC + alpha x %.3f'%sb,chi_of(Fb,Z,N,nn,alscale=sb),''))
        for lab,chi,_ in rows:
            x=full(chi,T,M,TH2[T])
            print(f"{name} {T} {lab:38s} sTot {x['sT']:6.1f} ({x['sT']/TOT[name][T]-1:+.1%}) sR {x['sR']:6.1f} sEl {x['sEl']:6.1f} el/Bar {x['sEl']/BAR[name][T]:.3f} qmin {x['qmin']:.3f} P(>th2mm) {x['P2']:.4f}",flush=True)
