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
abu={100:((3.16,1.808,0.268),(7.87,0.933,0.293)),200:((2.36,1.221,0.117),(4.54,0.701,0.135))}
hor={100:(5.295,1.435,0.322),150:(3.845,1.245,0.195),200:(3.45,0.953,0.131)}
tg={'C-12':(6,6,1.692,1.082,11174.86),'O-16':(8,8,1.833,1.544,14895.08)}
for name,(Z,N,a,al,M) in tg.items():
    rho=ho_rho(r,a,al)
    for T in (100,150,200):
        res=[run(rho,Z,N,M,T,nnsaid(T),True,'SAID+Eq9 unfolded')]
        res.append(run(rho,Z,N,M,T,nnsaid(T),False,'SAID+Eq9 charge dens'))
        h=hor[T]; res.append(run(rho,Z,N,M,T,(h,h),True,'Horiuchi isospin-mean'))
        if T in abu: res.append(run(rho,Z,N,M,T,abu[T],True,'Abu-Ibrahim TabI'))
        for x in res:
            print(f"{name} {T} {x['lab']:24s} sR {x['sR']:6.1f} sEl {x['sEl']:6.1f} sTot {x['sTot']:6.1f} (OT {x['sTot_ot']:6.1f}) tot-R {x['sTot']-x['sR']:6.1f} qmin {x['qmin']:.3f}",flush=True)
