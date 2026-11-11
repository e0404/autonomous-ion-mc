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

HB2M=HBARC**2/(2*MP)
def cv_R(rho_tot,E):
    kf=(1.5*math.pi**2*np.maximum(rho_tot,0))**(1/3); EF=HB2M*kf*kf
    x=np.where(EF>0,EF/E,0.0)
    hi=1-1.4*x
    lo=1-1.4*x+0.4*x*np.clip(2-1/np.where(x>0,x,1),0,None)**2.5
    out=np.where(E>=2*EF,hi,lo)
    return np.where(E<=EF,0.0,out)
ASYM_RA=1/(1+1.892)
EF0=HB2M*(1.5*math.pi**2*RHO0)**(2/3)
print(f"E_F(rho0=0.17)={EF0:.2f} MeV, 7E_F/5={1.4*EF0:.2f} MeV; R-A asym factor {ASYM_RA:.4f}")
tg={'C-12':(6,6,1.692,1.082,11174.86),'O-16':(8,8,1.833,1.544,14895.08)}
for name,(Z,N,a,al,M) in tg.items():
    ap,alp=unf(a,al); rho=ho_norm(r,ap,alp); A=Z+N
    for T in (100,150,200):
        nn=nnsaid(T)
        rows=[('M2',ff_raw(rho,QG)),
              ('R-A rho1=0,rho2=rho_tot',ff_raw(rho*ASYM_RA*bdc_R(A*rho,T,1.0),QG)),
              ('R-B rho1=rho_p,rho2=rho_n',ff_raw(rho*bdc_R(A*rho,T,1.0),QG)),
              ('R-C Clementel-Villi LDA',ff_raw(rho*cv_R(A*rho,T),QG))]
        rc=cv_R(np.array([A*rho.max()]),T)[0]; rb=bdc_R(np.array([A*rho.max()]),T)[0]
        for lab,Fq in rows:
            x=full(chi_of(Fq,Z,N,nn),T,M,TH2[T])
            print(f"{name} {T} {lab:28s} sTot {x['sT']:6.1f} ({x['sT']/TOT[name][T]-1:+.1%}) sR {x['sR']:6.1f} sEl {x['sEl']:6.1f} el/Bar {x['sEl']/BAR[name][T]:.3f} qmin {x['qmin']:.3f} P2 {x['P2']:.4f}",flush=True)
        print(f"   factor at peak density {A*rho.max():.3f}: R-B {rb:.3f}  R-C {rc:.3f}  R-A {ASYM_RA*rb:.3f}")
